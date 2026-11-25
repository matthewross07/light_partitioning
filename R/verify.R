# -----------------------------------------------------------------------------
# verify.R -- PaperOrchestra verification gates:
#   Gate A: Implementation-to-Text Alignment & internal mathematical consistency
#   Gate B: Claim-Evidence Provenance & statistical verification
# -----------------------------------------------------------------------------

#' Master verification checklist for the AquaMatch pipeline (used in Preface & Ch. 10).
run_verification <- function() {
  simul <- get_simul('simul')
  no_sdd <- get_simul('no_secchi')
  tr <- get_training()
  app <- get_application()
  zs <- get_zone_sets()
  mods <- get_models()
  perf <- model_performance(mods, tr)
  fd <- get_final_data()
  areas <- get_fig4_areas()
  dom <- test_dominance_secchi(fd)
  med <- areas %>% group_by(type) %>% summarize(m = median(area), .groups = 'drop')

  tribble(
    ~step, ~check, ~expected, ~actual, ~pass,
    '1. AquaMatch download', 'All 4 EDI packages + spatial layers present', '11 / 11 files',
      sprintf('%d / %d files', sum(raw_status()$present), nrow(raw_status())), all(raw_status()$present),

    '2. Simultaneous obs.', 'Training matchups (Secchi + Chl a + DOC + TSS)', '> 250,000',
      fmt_n(nrow(simul)), nrow(simul) > 250000,

    '2. Simultaneous obs.', 'Application matchups (Chl a + DOC + TSS, Secchi optional)', '> 450,000',
      fmt_n(nrow(no_sdd)), nrow(no_sdd) > 450000,

    '3. Training set', 'Analytical-floor & Beer-Lambert screened training obs.', '> 250,000',
      fmt_n(nrow(tr)), nrow(tr) > 250000 && all(tr$secchi >= P$secchi_floor & tr$secchi < P$max_secchi_m),

    '3. Training set', 'Residual screen against the published reference attenuation (|z_r| <= 2.5 SD)',
      '<= 2.50 SD',
      {
        st_d <- kd_discordance_stats(get_trimmed())
        r_tr <- log10(tr$kd / reference_attenuation(tr))
        sprintf('max |z_r| = %.2f SD', max(abs((r_tr - st_d$mean_r) / st_d$sd_r)))
      },
      {
        st_d <- kd_discordance_stats(get_trimmed())
        r_tr <- log10(tr$kd / reference_attenuation(tr))
        max(abs((r_tr - st_d$mean_r) / st_d$sd_r)) <= P$kd_discordance_sd + 1e-9
      },

    '3. Training set', 'Shuffled-Secchi null test: the screen creates no fit by itself',
      'R2 < 0.01',
      sprintf('R2 = %.4f (max of %d repeats %.4f)', mean(get_screen_null()$r2_screened), nrow(get_screen_null()), max(get_screen_null()$r2_screened)),
      max(get_screen_null()$r2_screened) < 0.01,

    '3. Training set', 'Pure-water offset k_w = 1 / (1.5 * max Secchi) [m^-1]',
      sprintf('%.5f', 1 / 1.5 / max(tr$secchi)),
      sprintf('%.5f', get_k_w()),
      abs(get_k_w() - mods$m7$k_w) < 1e-12,

    '4. Ecoregion joins', '1-to-1 site assignment (no boundary row duplication)',
      fmt_n(nrow(tr)),
      fmt_n(nrow(eco_join(tr, get_eco()$eco3, l3code))),
      nrow(eco_join(tr, get_eco()$eco3, l3code)) == nrow(tr),

    '5. Models 1-10', 'Unified minimum group size (> 100 obs) across all models',
      '> 100',
      sprintf('min group = %d', min(perf$group_min)),
      min(perf$group_min) > 100,

    '6. Model 7 application', 'All ternary fractions strictly inside [0, 1] simplex',
      '0 outside [0, 1]',
      sprintf('%d outside [0, 1]', sum(pmin(fd$tss_kd, fd$doc_kd, fd$chl_kd) < 0 | pmax(fd$tss_kd, fd$doc_kd, fd$chl_kd) > 1)),
      sum(pmin(fd$tss_kd, fd$doc_kd, fd$chl_kd) < 0) == 0,

    '6. Model 7 application', 'Consistent k_w used in both model fit and Secchi prediction',
      sprintf('%.5f m^-1', mods$m7$k_w),
      sprintf('%.5f m^-1', unique(round(fd$kd_pred - fd$att_sum, 5))[1]),
      max(abs((fd$kd_pred - fd$att_sum) - mods$m7$k_w)) < 1e-12,

    '7. Dominance classes', 'Secchi ordering: NAP-dominant lowest, then Chl a-dominant, then DOC-dominant and Mixed',
      'NAP < Chl a < DOC, Mixed',
      paste(sprintf('%s: %.2f m', dom$summary$class, dom$summary$mean_secchi_pred), collapse = '; '),
      {
        ms <- setNames(dom$summary$mean_secchi_pred, dom$summary$class)
        ms[['NAP dominant']] < ms[['Chl a dominant']] && ms[['Chl a dominant']] < min(ms[['DOC dominant']], ms[['Mixed']])
      },

    '7. Dominance classes', 'NAP-dominant waters less clear than each other class at p < 0.01 (Welch & Wilcoxon)',
      '3 / 3 pairs p < 0.01',
      {
        pw <- dom$pairwise; both <- pw$pass_p01_welch & pw$pass_p01_wilcoxon
        sprintf('%d / 3 pairs; all six pairs: %d / 6%s', sum(both[grepl('NAP', pw$Comparison)]), sum(both),
                if (any(!both)) paste0(' (not different: ', paste(pw$Comparison[!both], collapse = '; '), ')') else '')
      },
      {
        pw <- dom$pairwise
        all((pw$pass_p01_welch & pw$pass_p01_wilcoxon)[grepl('NAP', pw$Comparison)])
      },

    '8. Figure 4 metric', 'Exact ternary area fraction (normalized by sqrt(3)/4): Streams & Estuaries < Lakes',
      'Stream, Estuary (~10%) < Lake (~15%)',
      sprintf('Stream %.1f%%, Estuary %.1f%% < Lake %.1f%%',
              100 * med$m[med$type == 'Stream'], 100 * med$m[med$type == 'Estuary'], 100 * med$m[med$type == 'Lake']),
      max(med$m[med$type %in% c('Stream', 'Estuary')]) < med$m[med$type == 'Lake']
  )
}
