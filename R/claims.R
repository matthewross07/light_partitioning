# -----------------------------------------------------------------------------
# claims.R -- Formal statistical tests of manuscript claims (Secchi depth by
#             LAC dominance class, p < 0.01) and PaperOrchestra literature
#             verification benchmarks
# -----------------------------------------------------------------------------

#' Three ways of saying that one LAC leads, counted per LAC:
#'   dominant            share > P$dominance_threshold (70%)
#'   majority contributor share > P$majority_threshold (50%)
#'   largest contributor  share is the largest of the three
dominance_definitions <- function(fd) {
  f <- cbind(NAP = fd$tss_kd, DOC = fd$doc_kd, `Chl a` = fd$chl_kd)
  top <- colnames(f)[max.col(f, ties.method = 'first')]
  tibble(LAC = colnames(f),
         n_dominant = colSums(f > P$dominance_threshold), n_majority = colSums(f > P$majority_threshold),
         n_largest = as.integer(table(factor(top, colnames(f))))) %>%
    mutate(across(starts_with('n_'), ~ 100 * .x / nrow(fd), .names = 'pct_{.col}')) %>%
    rename_with(~ sub('pct_n_', 'pct_', .x))
}

#' Classify each record into one of the four LAC dominance classes (share > P$dominance_threshold).
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

#' Compare fitted coefficients (per unit of 1/z_SD) with published coefficients (per
#' unit of K_d). Since K_d = a / z_SD, a published k and our alpha are related by
#' k = a * alpha, so each pair implies a conversion factor a = k / alpha, which can be
#' checked against measured values of K_d x z_SD (Koenings and Edmundson 1991).
#' Every published value was read in the source: page or table given.
compare_bannister_literature <- function(models = get_models()) {
  c1 <- models$m1$coefs %>% select(term, estimate) %>% deframe()
  c2 <- models$m2$coefs %>% select(type, term, estimate) %>%
    pivot_wider(names_from = type, values_from = estimate) %>% column_to_rownames('term')
  rows <- tribble(
    ~Constituent, ~term, ~ours, ~published, ~source,
    'Chl a (m² mg⁻¹)', 'Chl_a', c1[['Chl_a']], 0.016,  'Bannister (1974), p. 6; mean of four estimates, 0.013–0.020 (p. 10)',
    'Chl a (m² mg⁻¹)', 'Chl_a', c1[['Chl_a']], 0.014,  'Atlas and Bannister (1980), p. 157; 0.014 ± 0.002 at the surface, 0.005–0.021 overall',
    'Chl a (m² mg⁻¹)', 'Chl_a', c2['Chl_a', 'Estuary'], 0.0154, 'Gallegos (2001), Table 5; estuary coefficient, mean of 0.012–0.019',
    'NAP (m² g⁻¹)',    'NAP',   c1[['NAP']],   0.0625, 'Gohin et al. (2005), Eq. 21; non-living suspended matter',
    'NAP (m² g⁻¹)',    'NAP',   c2['NAP', 'Estuary'], 0.0806, 'Gallegos (2001), Table 5; total suspended solids, estuary, mean of 0.074–0.095',
    'DOC (m² g⁻¹)',    'DOC',   c2['DOC', 'Lake'],    0.22,   'Morris et al. (1995), Table 3; 63 lakes',
    'DOC (m² g⁻¹)',    'DOC',   c2['DOC', 'Estuary'], 0.0471, 'Gallegos (2001), Table 5; estuary, mean of 0.037–0.051')
  which_ours <- c('Model 1, all waters', 'Model 1, all waters', 'Model 2, estuaries', 'Model 1, all waters', 'Model 2, estuaries', 'Model 2, lakes', 'Model 2, estuaries')
  rows %>% mutate(`Our coefficient` = sprintf('%.4f (%s)', ours, which_ours),
                  `Published coefficient` = sprintf('%s', format(published, drop0trailing = TRUE)),
                  `Implied a = published / ours` = sprintf('%.1f', published / ours)) %>%
    select(Constituent, `Our coefficient`, `Published coefficient`, `Implied a = published / ours`, Source = source)
}

#' PaperOrchestra Literature Verification Table: documents the peer-reviewed
#' support for every dataset, transformation, and methodological decision.
verified_literature_table <- function() {
  tribble(
    ~Step_or_Decision, ~Equation_or_Parameter, ~Citation, ~Identifier_or_DOI, ~Verification_Notes,
    'In-situ Secchi depth (sdd)', 'AquaMatch harmonized WQP Secchi depth (2.74M records, 1970-2024)',
    'De La Torre et al. (2025)', 'EDI package edi.1856.1',
    'Retrieved via AquaMatchr::download_parameters("sdd"); records with sdd_mdl >= 2 (<=0.01 m floor placeholders) set to NA.',

    'In-situ Chlorophyll a (chla)', 'AquaMatch harmonized WQP Chlorophyll a (3.39M records, 1970-2024)',
    'Brousil et al. (2024)', 'EDI package edi.1756.2',
    'Retrieved via AquaMatchr::download_parameters("chla").',

    'In-situ Dissolved Organic Carbon (doc)', 'AquaMatch harmonized WQP DOC (1.84M records, 1970-2024)',
    'Brousil et al. (2024)', 'EDI package edi.1809.1',
    'Retrieved via AquaMatchr::download_parameters("doc").',

    'In-situ Total Suspended Solids & SSC (tss)', 'AquaMatch harmonized WQP TSS & SSC (4.61M records, 1970-2025)',
    'Brousil et al. (2025)', 'EDI package edi.2048.2',
    'Retrieved via AquaMatchr::download_parameters("tss"); prefers TSS when both TSS and SSC are measured on the same sample.',

    'Diffuse attenuation from Secchi depth', 'K_d = 1 / z_SD (Eq. 1 with a = 1)',
    'Lee et al. (2015, 2018)', '10.1016/j.rse.2015.08.002; 10.1002/lno.10940',
    'Contrast-transmittance theory showing K_d(tr) is inversely proportional to Secchi depth without empirical Poole & Atkins (1929) constants.',

    'Algal dry biomass & NAP separation', 'biomass [mg/L] = 0.234 * (Chl a)^0.57; NAP = TSS - biomass',
    'Gohin et al. (2005)', '10.1016/j.rse.2004.11.007',
    'Nonlinear power-law conversion accounting for lower C:Chl ratios at high algal concentrations; prevents over-subtraction of algal mass.',

    'Published Chl a attenuation coefficient', 'k_c = 0.016 m^2 (mg Chl a)^-1 for diffuse attenuation of PAR',
    'Bannister (1974); Atlas and Bannister (1980)', '10.4319/lo.1974.19.1.0001; 10.4319/lo.1980.25.1.0157',
    'An estimate, not an upper limit (Bannister p. 6 and p. 10). It is per unit of K_d; ours is per unit of 1/z_SD, so the two differ by the factor a = K_d x z_SD.',

    'Values below reporting limits', 'Chl a, DOC, TSS (and NAP where algal mass >= TSS) below the limit replaced by draws from a fitted lognormal truncated at the limit',
    'Lubin et al. (2004); Helsel (2006)', '10.1289/ehp.7199; 10.1016/j.chemosphere.2006.04.051',
    'Fill-in from a fitted distribution instead of substituting the limit; fewer than 3% of values are affected. Draws are tied to the record so training and application use the same value.',

    'Data screen against a published reference attenuation', 'K_ref = 0.1 + 0.0625 NAP + k_DOC DOC + 0.016 Chl a, k_DOC = 0.22 (lakes, streams) or 0.047 (estuaries); keep |z_r| <= 2.5 on r = log10[(1/z_SD) / K_ref]',
    'Gohin et al. (2005); Morris et al. (1995); Gallegos (2001); Bannister (1974)', '10.1016/j.rse.2004.11.007; 10.4319/lo.1995.40.8.1381; 10.2307/1353240; 10.4319/lo.1974.19.1.0001',
    'Every coefficient is a published value used as printed (Gohin Eq. 21; Morris Table 3; Gallegos Table 5; Bannister p. 6). A constant Secchi conversion factor cancels in z_r. Shuffled-Secchi test and reference sensitivity: Section 4.',

    'Secchi depth and attenuation', 'z_SD = Gamma / (alpha + K), Gamma ~ 8-9; K_d x z_SD = 0.59 to 4.09 across studies',
    'Tyler (1968); Preisendorfer (1986); Koenings and Edmundson (1991)', '10.4319/lo.1968.13.1.0001; 10.4319/lo.1986.31.5.0909; 10.4319/lo.1991.36.1.0091',
    'Theory and measurements behind the conversion factor a in Eq. 1 and the comparison of our coefficients with published ones.',

    'DOC as CDOM proxy & regional SUVA variation', 'alpha_DOC varies across ecoregions with DOM aromaticity',
    'Spencer et al. (2012); Griffin et al. (2018)', '10.1029/2011JG001928; 10.1016/j.watres.2018.08.007',
    'Justifies ecoregion-stratified attenuation slopes (Model 7) to capture regional differences in CDOM:DOC stoichiometry.',

    'Hierarchical Ecoregion Stratification', 'EPA / CEC North American Ecoregions Levels I-IV',
    'Omernik & Griffith (2014)', '10.1007/s00267-014-0364-1',
    'Provides continentally nested physiographic, soil, and land-cover boundaries governing sediment and DOM optical properties.',

    'River network longitudinal distances', 'RiverSR modified NHDPlusV2 flowlines (TermnlP = 200005438)',
    'Gardner et al. (2021)', '10.5281/zenodo.4304567',
    'Provides path lengths along the Delaware River mainstem for Figure 3.'
  )
}
