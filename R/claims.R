# -----------------------------------------------------------------------------
# claims.R -- Formal statistical tests of manuscript claims (Secchi depth by
#             LAC dominance class, p < 0.01) and PaperOrchestra literature
#             verification benchmarks
# -----------------------------------------------------------------------------

#' Classify each record into one of the four LAC dominance classes (>80% threshold).
classify_dominance <- function(fd, threshold = P$dominance_threshold) {
  fd %>%
    mutate(
      class = factor(
        case_when(
          tss_kd > threshold & doc_kd < 1 - threshold ~ 'NAP dominant',
          chl_kd > threshold & tss_kd < 1 - threshold ~ 'Chl a dominant',
          chl_kd < 1 - threshold & doc_kd > threshold ~ 'DOC dominant',
          TRUE ~ 'Mixed'
        ),
        levels = c('NAP dominant', 'Chl a dominant', 'DOC dominant', 'Mixed')
      ),
      secchi_combined = ifelse(is.na(secchi), secchi_pred, secchi)
    )
}

#' Formal statistical test of the Section 3.2 claim that mean Secchi depths
#' across the four LAC dominance classes differ significantly at p < 0.01.
#' Evaluates:
#'   1. Welch heteroscedastic one-way ANOVA + Holm-adjusted pairwise Welch t-tests
#'      (appropriate because SD in Mixed/DOC classes is >2x SD in NAP/Chl-a classes)
#'   2. Non-parametric Kruskal-Wallis + Holm-adjusted pairwise Wilcoxon rank-sum tests
#'   3. Classical equal-variance one-way ANOVA + Tukey HSD
test_dominance_secchi <- function(fd) {
  sz <- classify_dominance(fd)

  summary_tbl <- sz %>%
    group_by(class) %>%
    summarize(
      n_obs = n(),
      n_sites = n_distinct(SiteID),
      mean_secchi_pred = mean(secchi_pred, na.rm = TRUE),
      sd_secchi_pred   = sd(secchi_pred, na.rm = TRUE),
      se_secchi_pred   = sd_secchi_pred / sqrt(n_obs),
      median_secchi_pred = median(secchi_pred, na.rm = TRUE),
      mean_secchi_combined = mean(secchi_combined, na.rm = TRUE),
      mean_secchi_measured = mean(secchi, na.rm = TRUE),
      .groups = 'drop'
    )

  # Omnibus tests on modeled Secchi depth (`secchi_pred`)
  welch_aov <- oneway.test(secchi_pred ~ class, data = sz, var.equal = FALSE)
  kw_test   <- kruskal.test(secchi_pred ~ class, data = sz)
  eq_aov    <- summary(aov(secchi_pred ~ class, data = sz))[[1]]

  omnibus_tbl <- tibble(
    Test = c('Welch heteroscedastic ANOVA (secchi_pred)',
             'Kruskal-Wallis rank-sum test (secchi_pred)',
             'Equal-variance one-way ANOVA (secchi_pred)'),
    Statistic = c(sprintf('F = %.1f', welch_aov$statistic),
                  sprintf('Chi-sq = %.1f', kw_test$statistic),
                  sprintf('F = %.1f', eq_aov$`F value`[1])),
    p_value = c(welch_aov$p.value, kw_test$p.value, eq_aov$`Pr(>F)`[1])
  )

  # All 6 pairwise comparisons
  w_t <- pairwise.t.test(sz$secchi_pred, sz$class, pool.sd = FALSE, p.adjust.method = 'holm')$p.value
  w_x <- pairwise.wilcox.test(sz$secchi_pred, sz$class, p.adjust.method = 'holm')$p.value
  tuk <- as_tibble(TukeyHSD(aov(secchi_pred ~ class, data = sz))$class, rownames = 'pair')

  pairs <- expand_grid(g2 = rownames(w_t), g1 = colnames(w_t)) %>%
    filter(!is.na(map2_dbl(g2, g1, ~ w_t[.x, .y]))) %>%
    transmute(
      pair = paste0(g2, '-', g1),
      Comparison = paste(g2, 'vs.', g1),
      p_welch_holm    = map2_dbl(g2, g1, ~ w_t[.x, .y]),
      p_wilcoxon_holm = map2_dbl(g2, g1, ~ w_x[.x, .y])
    ) %>%
    left_join(tuk %>% select(pair, diff_m = diff, p_tukey_eqvar = `p adj`), by = 'pair') %>%
    select(Comparison, diff_m, p_welch_holm, p_wilcoxon_holm, p_tukey_eqvar) %>%
    mutate(pass_p01_welch    = p_welch_holm < 0.01,
           pass_p01_wilcoxon = p_wilcoxon_holm < 0.01)

  list(summary = summary_tbl, omnibus = omnibus_tbl, pairwise = pairs)
}

#' Compare fitted attenuation coefficients against theoretical and empirical
#' bio-optical benchmarks in the literature (Bannister 1974; Kirk 1994; etc.).
compare_bannister_literature <- function(models = get_models()) {
  c1 <- models$m1$coefs %>% select(term, estimate) %>% deframe()
  c2 <- models$m2$coefs %>% select(type, term, estimate) %>%
    pivot_wider(names_from = type, values_from = estimate) %>%
    column_to_rownames('term')

  tribble(
    ~Constituent, ~Units, ~`AquaMatch Global (Model 1)`, ~`AquaMatch by Type (Model 2)`, ~`Literature Range`, ~`Primary References`,
    'Chl a (alpha_Chl)', 'm^2 (mg Chl a)^-1',
    sprintf('%.4f', c1[['Chl_a']]),
    sprintf('Estuary: %.4f; Lake: %.4f; Stream: %.4f', c2['Chl_a', 'Estuary'], c2['Chl_a', 'Lake'], c2['Chl_a', 'Stream']),
    '0.008 - 0.020 (canonical k_c = 0.016)',
    'Bannister (1974) L&O 19:1-12; Morel (1988) JGR 93:10749; Kirk (1994) Light and Photosynthesis in Aquatic Ecosystems',

    'NAP (alpha_NAP)', 'm^2 g^-1 (L mg^-1 m^-1)',
    sprintf('%.4f', c1[['NAP']]),
    sprintf('Estuary: %.4f; Lake: %.4f; Stream: %.4f', c2['NAP', 'Estuary'], c2['NAP', 'Lake'], c2['NAP', 'Stream']),
    '0.025 - 0.070',
    'Gohin et al. (2005) Int. J. Remote Sens. 26:1639; Gallegos et al. (2011) Estuar. Coast. Shelf Sci. 95:218; Kirk (1994)',

    'DOC (alpha_DOC)', 'm^2 g^-1 (L mg^-1 m^-1)',
    sprintf('%.4f', c1[['DOC']]),
    sprintf('Estuary: %.4f; Lake: %.4f; Stream: %.4f', c2['DOC', 'Estuary'], c2['DOC', 'Lake'], c2['DOC', 'Stream']),
    '0.030 - 0.150',
    'Morris et al. (1995) L&O 40:1381; Spencer et al. (2012) JGR 117:G03001; Griffin et al. (2018) Water Res. 144:719'
  )
}

#' PaperOrchestra Literature Verification Table: documents the peer-reviewed
#' support for every dataset, transformation, and methodological decision.
verified_literature_table <- function() {
  tribble(
    ~Step_or_Decision, ~Equation_or_Parameter, ~Citation, ~Identifier_or_DOI, ~Verification_Notes,
    'In-situ Secchi depth (sdd)', 'AquaMatch harmonized WQP Secchi depth (2.74M records, 1970-2024)',
    'Brousil et al. (2025)', 'EDI package edi.1856.1',
    'Retrieved via AquaMatchr::download_parameters("sdd"); records with sdd_mdl >= 2 (<=0.01 m floor placeholders) set to NA.',

    'In-situ Chlorophyll a (chla)', 'AquaMatch harmonized WQP Chlorophyll a (3.39M records, 1970-2024)',
    'Brousil et al. (2024)', 'EDI package edi.1756.2',
    'Retrieved via AquaMatchr::download_parameters("chla").',

    'In-situ Dissolved Organic Carbon (doc)', 'AquaMatch harmonized WQP DOC (1.84M records, 1970-2024)',
    'Brousil et al. (2024)', 'EDI package edi.1809.1',
    'Retrieved via AquaMatchr::download_parameters("doc").',

    'In-situ Total Suspended Solids & SSC (tss)', 'AquaMatch harmonized WQP TSS & SSC (4.61M records, 1970-2025)',
    'De La Torre et al. (2025)', 'EDI package edi.2048.2',
    'Retrieved via AquaMatchr::download_parameters("tss"); prefers TSS when both TSS and SSC are measured on the same sample.',

    'Diffuse attenuation from Secchi depth', 'K_d = 1 / z_SD (Eq. 1 with a = 1)',
    'Lee et al. (2015, 2018)', '10.1002/2015GL065066; 10.1029/2018JC014018',
    'Contrast-transmittance theory showing K_d(tr) is inversely proportional to Secchi depth without empirical Poole & Atkins (1929) constants.',

    'Algal dry biomass & NAP separation', 'biomass [mg/L] = 0.234 * (Chl a)^0.57; NAP = TSS - biomass',
    'Gohin et al. (2005)', '10.1080/01431160512331338032',
    'Nonlinear power-law conversion accounting for lower C:Chl ratios at high algal concentrations; prevents over-subtraction of algal mass.',

    'Specific Chl a attenuation benchmark', 'k_c = 0.016 m^2 (mg Chl a)^-1 (range 0.008 - 0.020)',
    'Bannister (1974); Kirk (1994)', '10.4319/lo.1974.19.1.0001',
    'Matches AquaMatch global alpha_Chl = 0.0088 m^2/mg and Estuary alpha_Chl = 0.0136 m^2/mg.',

    'A priori WQP Beer-Lambert residual screen', '|z_r| <= 2.5 SD (~4.7-fold [~5x]) on r = log10(K_d / K_d,lit) using published k_w=0.08, k_NAP=0.04, k_DOC=0.07, k_Chl=0.016',
    'Bannister (1974); Kirk (1994); Morris et al. (1995); Gohin et al. (2005)', '10.4319/lo.1974.19.1.0001; 10.4319/lo.1995.40.8.1381',
    'Flags 2.49% discordant records: 1.84% upper-tail shallow bottom hits (Secchi set to NA so valid 3-LAC chemistry is retained in application) and 0.65% lower-tail chemistry unit errors.',

    'DOC as CDOM proxy & regional SUVA variation', 'alpha_DOC varies across ecoregions with DOM aromaticity',
    'Spencer et al. (2012); Griffin et al. (2018)', '10.1029/2011JG001928; 10.1016/j.watres.2018.07.066',
    'Justifies ecoregion-stratified attenuation slopes (Model 7) to capture regional differences in CDOM:DOC stoichiometry.',

    'Hierarchical Ecoregion Stratification', 'EPA / CEC North American Ecoregions Levels I-IV',
    'Omernik & Griffith (2014)', '10.1007/s00267-014-0364-1',
    'Provides continentally nested physiographic, soil, and land-cover boundaries governing sediment and DOM optical properties.',

    'River network longitudinal distances', 'RiverSR modified NHDPlusV2 flowlines (TermnlP = 200005438)',
    'Gardner et al. (2021)', '10.5281/zenodo.4304567',
    'Provides path lengths along the Delaware River mainstem for Figure 3.'
  )
}
