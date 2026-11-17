# =============================================================================
# 05b_modeles_fonctionnels.R
# Modélisation des indices fonctionnels : saison, profondeur, station
#
# À lancer après 05_diversite_fonctionnelle.R
# Entrée : data/processed/fonctionnel_indices.csv
#
# Démarche (pour chaque indice) :
#   1. Loi adaptée au type de variable (tableau `reponses` ci-dessous)
#   2. Modèle de référence fixé a priori :
#        indice ~ saison + profondeur + (1 | station)
#      -> tests du rapport de vraisemblance (LRT), coefficients,
#         R² de vraisemblance (part expliquée par saison + profondeur)
#   3. Sélection par AICc parmi 5 candidats (nul, saison, profondeur,
#      saison + profondeur, interaction) : quelle structure les données
#      soutiennent-elles ?
#   4. Comparaisons entre saisons (emmeans, ajustement de Tukey)
#   5. Diagnostics DHARMa (loi, dispersion, extrêmes, autocorrélation spatiale)
#   6. Sensibilité : modèle de référence refait sans les traits douteux
#   7. Correction de Holm sur l'ensemble des tests (nombreux indices)
#
# Sorties : data/processed/modeles_fonctionnels/*.csv
#           outputs/figures/modeles_fonctionnels/ (+ diagnostics/)
# =============================================================================

pkgs <- c("here", "dplyr", "tidyr", "ggplot2", "glmmTMB", "DHARMa", "emmeans")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants)) install.packages(manquants)

library(dplyr)
library(ggplot2)

source(here::here("R", "00_config.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "modeles_utils.R"), encoding = "UTF-8")

set.seed(graine_aleatoire)
dossier_out <- file.path(dossiers$processed, "modeles_fonctionnels")
dossier_fig <- here::here("outputs", "figures", "modeles_fonctionnels")
dossier_diag <- file.path(dossier_fig, "diagnostics")
for (d in c(dossier_out, dossier_fig, dossier_diag)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# =============================================================================
# Indices modélisés et loi associée (modifiable)
# =============================================================================
reponses <- tribble(
  ~indice,                                   ~famille,   ~transformation, ~label,
  "S",                                       "genpois",  "aucune", "Richesse en poissons",
  "FDis",                                    "gaussian", "aucune", "Dispersion fonctionnelle (FDis)",
  "RaoQ",                                    "gaussian", "aucune", "Entropie de Rao",
  "SES_FDis",                                "gaussian", "aucune", "SES de FDis",
  "SES_RaoQ",                                "gaussian", "aucune", "SES de Rao",
  "FRic",                                    "beta",     "aucune", "Richesse fonctionnelle (FRic)",
  "FEve",                                    "beta",     "aucune", "Régularité fonctionnelle (FEve)",
  "FDiv",                                    "beta",     "aucune", "Divergence fonctionnelle (FDiv)",
  "CWM_taille_max_cm",                       "gaussian", "log10",  "CWM taille max. (log10 cm)",
  "CWM_niveau_trophique",                    "gaussian", "aucune", "CWM niveau trophique",
  "CWM_guilde_trophique__microbenthivore",   "beta",     "01",     "Part de microbenthivores",
  "CWM_guilde_estuarienne__migrant marin",   "beta",     "01",     "Part de migrants marins",
  "CWM_position_verticale__benthique",       "beta",     "01",     "Part d'espèces benthiques"
)
# transformation : "log10" ; "01" = proportion ramenée dans ]0 ; 1[ pour la
# loi bêta (les 0 et 1 exacts sont impossibles en bêta)

# =============================================================================
# Données
# =============================================================================
indices <- read.csv(file.path(dossiers$processed, "fonctionnel_indices.csv"),
                    encoding = "UTF-8", check.names = FALSE) |>
  mutate(saison = factor(saison, levels = c("printemps", "ete", "automne")),
         station = factor(station))

absents <- setdiff(reponses$indice, names(indices))
if (length(absents)) {
  message("Indices absents (ignorés) : ", paste(absents, collapse = ", "))
  reponses <- filter(reponses, !indice %in% absents)
}

preparer_reponse <- function(d, r) {
  y <- d[[r$indice]]
  y <- switch(r$transformation, log10 = log10(y), `01` = compresser_01(y), y)
  d$y <- y
  d[!is.na(d$y), ]
}

# =============================================================================
# Boucle de modélisation
# =============================================================================
res <- list(selection = list(), tests = list(), coefficients = list(),
            moyennes = list(), contrastes = list(), diagnostics = list(),
            sensibilite = list(), predictions = list(), resume = list())

for (i in seq_len(nrow(reponses))) {
  r <- reponses[i, ]
  message(sprintf("\n--- %s (%s) ---", r$indice, r$famille))
  d <- preparer_reponse(indices, r)
  fam <- famille_modele(r$famille)

  # 1. Candidats et sélection AICc
  cand <- lapply(formules_candidates("y"), ajuster_glmm, data = d, famille = fam)
  sel <- selection_aicc(cand)
  res$selection[[r$indice]] <- cbind(indice = r$indice, sel)

  # 2. Modèle de référence (fixé a priori)
  m <- cand$saison_profondeur
  if (is.null(m)) { warning("Échec du modèle de référence pour ", r$indice); next }
  avert <- attr(m, "avertissements")

  lrt <- tests_lrt(m)
  res$tests[[r$indice]] <- cbind(indice = r$indice, famille = r$famille,
                                 n = nobs(m), lrt)
  res$coefficients[[r$indice]] <- cbind(indice = r$indice, coefficients_modele(m))

  # 3. Comparaisons entre saisons
  cs <- tryCatch(comparaisons_saisons(m), error = function(e) NULL)
  if (!is.null(cs)) {
    res$moyennes[[r$indice]]   <- cbind(indice = r$indice, cs$moyennes)
    res$contrastes[[r$indice]] <- cbind(indice = r$indice, cs$contrastes)
  }

  # 4. Diagnostics
  diag <- diagnostics_dharma(
    m, d, file.path(dossier_diag, paste0(make.names(r$indice), ".png")))
  r2 <- r2_vraisemblance(m, cand$nul)
  res$diagnostics[[r$indice]] <- data.frame(
    indice = r$indice, diag, R2_vraisemblance = r2,
    variance_station = unname(unlist(glmmTMB::VarCorr(m)$cond$station)),
    avertissements = paste(avert, collapse = " | "))

  # 5. Sensibilité aux traits douteux
  d2 <- d[!d$douteux, ]
  m2 <- ajuster_glmm(y ~ saison + profondeur_m + (1 | station), d2, fam)
  if (!is.null(m2)) {
    res$sensibilite[[r$indice]] <- cbind(indice = r$indice, n = nobs(m2),
                                         tests_lrt(m2))
  }

  # 6. Prédictions le long de la profondeur
  pr <- predire_profondeur(m, d)
  res$predictions[[r$indice]] <- cbind(indice = r$indice, label = r$label, pr)

  res$resume[[r$indice]] <- data.frame(
    indice = r$indice, label = r$label, famille = r$famille, n = nobs(m),
    meilleur_modele_AICc = sel$modele[1],
    p_saison = lrt$p[lrt$effet == "saison"],
    p_profondeur = lrt$p[lrt$effet == "profondeur_m"],
    R2 = r2,
    probleme_diagnostic = any(unlist(diag) < 0.05, na.rm = TRUE),
    ajustement_singulier = any(grepl("singular|convergence|Hessian", avert,
                                     ignore.case = TRUE)))
}

lier <- function(x) { x <- do.call(rbind, x); rownames(x) <- NULL; x }
tabs <- lapply(res, lier)

# Correction de Holm sur l'ensemble des tests saison / profondeur
tabs$tests$p_holm <- p.adjust(tabs$tests$p, method = "holm")
tabs$resume <- tabs$resume |>
  left_join(tabs$tests |> filter(effet == "saison") |>
              select(indice, p_saison_holm = p_holm), by = "indice") |>
  left_join(tabs$tests |> filter(effet == "profondeur_m") |>
              select(indice, p_profondeur_holm = p_holm), by = "indice") |>
  relocate(p_saison_holm, .after = p_saison) |>
  relocate(p_profondeur_holm, .after = p_profondeur)

# =============================================================================
# Export
# =============================================================================
for (nom in names(tabs)) {
  write.csv(tabs[[nom]], file.path(dossier_out, paste0(nom, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}

# =============================================================================
# Figures
# =============================================================================
ordre <- function(x) factor(x, levels = reponses$indice)   # ordre du tableau
couleurs_saison <- c(printemps = "#66a61e", ete = "#e6ab02", automne = "#d95f02")
labels <- setNames(reponses$label, reponses$indice)

# a. Moyennes marginales par saison (IC 95 %) sur les données brutes
brut <- lapply(seq_len(nrow(reponses)), function(i) {
  r <- reponses[i, ]
  d <- preparer_reponse(indices, r)
  if (r$transformation == "01") d$y <- d[[r$indice]]   # échelle d'origine
  data.frame(indice = r$indice, saison = d$saison, y = d$y)
}) |> bind_rows() |> mutate(indice = ordre(indice))

moy <- tabs$moyennes |> mutate(indice = ordre(indice))
if (!is.null(moy) && nrow(moy)) {
  # retour à l'échelle d'origine pour les proportions compressées
  p_moy <- ggplot() +
    geom_jitter(data = brut, aes(saison, y, colour = saison),
                width = 0.12, alpha = 0.45, size = 1.4) +
    geom_pointrange(data = moy, aes(saison, moyenne, ymin = IC95_bas,
                                    ymax = IC95_haut, colour = saison),
                    size = 0.5, linewidth = 0.9,
                    position = position_nudge(x = 0.28)) +
    facet_wrap(~ indice, scales = "free_y", ncol = 4,
               labeller = as_labeller(labels)) +
    scale_colour_manual(values = couleurs_saison, guide = "none") +
    labs(title = "Effet de la saison sur les indices fonctionnels",
         subtitle = "Points pâles : traits de chalut ; point + barre : moyenne du modèle (profondeur moyenne) et IC 95 %",
         x = NULL, y = NULL) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(size = 8.5))
  print(p_moy)
  ggsave(file.path(dossier_fig, "effet_saison.png"), p_moy,
         width = 12, height = 9, dpi = 150, bg = "white")
}

# b. Effet de la profondeur (prédictions par saison)
pred <- tabs$predictions |> mutate(indice = ordre(indice))
brut_prof <- lapply(seq_len(nrow(reponses)), function(i) {
  r <- reponses[i, ]
  d <- preparer_reponse(indices, r)
  data.frame(indice = r$indice, saison = d$saison,
             profondeur_m = d$profondeur_m, y = d$y)
}) |> bind_rows() |> mutate(indice = ordre(indice))
p_prof <- ggplot(pred, aes(profondeur_m, prediction, colour = saison,
                           fill = saison)) +
  geom_ribbon(aes(ymin = IC95_bas, ymax = IC95_haut), alpha = 0.12,
              colour = NA) +
  geom_line(linewidth = 0.9) +
  geom_point(data = brut_prof, aes(profondeur_m, y), size = 1.4, alpha = 0.6,
             inherit.aes = FALSE, colour = "grey30") +
  facet_wrap(~ indice, scales = "free_y", ncol = 4,
             labeller = as_labeller(labels)) +
  scale_colour_manual(values = couleurs_saison) +
  scale_fill_manual(values = couleurs_saison) +
  labs(title = "Effet de la profondeur (modèle saison + profondeur, IC 95 %)",
       subtitle = "Échelle du modèle (log10 pour la taille, proportions comprimées pour la loi bêta)",
       x = "Profondeur moyenne du trait (m)", y = NULL, colour = NULL,
       fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8.5))
print(p_prof)
ggsave(file.path(dossier_fig, "effet_profondeur.png"), p_prof,
       width = 12, height = 9, dpi = 150, bg = "white")

# c. Synthèse des p-values (Holm)
syn <- tabs$resume |>
  select(indice, label, p_saison_holm, p_profondeur_holm) |>
  tidyr::pivot_longer(c(p_saison_holm, p_profondeur_holm),
                      names_to = "effet", values_to = "p") |>
  mutate(effet = recode(effet, p_saison_holm = "Saison",
                        p_profondeur_holm = "Profondeur"),
         label = factor(label, levels = rev(reponses$label)),
         signif = cut(p, c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                      labels = c("< 0,001", "< 0,01", "< 0,05", "< 0,1", "n.s.")))
p_syn <- ggplot(syn, aes(effet, label, fill = signif)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = formatC(p, format = "g", digits = 2)), size = 3) +
  scale_fill_manual(values = c("< 0,001" = "#08519c", "< 0,01" = "#3182bd",
                               "< 0,05" = "#6baed6", "< 0,1" = "#c6dbef",
                               "n.s." = "grey95"), drop = FALSE) +
  labs(title = "Synthèse des tests (LRT, p corrigées de Holm)",
       x = NULL, y = NULL, fill = "p") +
  theme_minimal(base_size = 10)
print(p_syn)
ggsave(file.path(dossier_fig, "synthese_tests.png"), p_syn,
       width = 7, height = 5.5, dpi = 150, bg = "white")

# =============================================================================
# Résumé
# =============================================================================
message("\n================ RÉSUMÉ ================")
print(tabs$resume |>
        mutate(across(where(is.numeric), ~ signif(.x, 3))) |>
        select(-label))
message("\nDétails : ", dossier_out)
message("Diagnostics DHARMa (un graphique par indice) : ", dossier_diag)
