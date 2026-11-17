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
source(here::here("R", "fonctions", "figures_modeles.R"), encoding = "UTF-8")

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

# =============================================================================
# Modélisation (fonction commune dans R/fonctions/modeles_utils.R)
# =============================================================================
tabs <- modeliser_indices(indices, reponses, dossier_diag)

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
figures_modeles(indices, reponses, tabs, dossier_fig,
                sujet = "les indices fonctionnels")

# =============================================================================
# Résumé
# =============================================================================
message("\n================ RÉSUMÉ ================")
print(tabs$resume |>
        mutate(across(where(is.numeric), ~ signif(.x, 3))) |>
        select(-label))
message("\nDétails : ", dossier_out)
message("Diagnostics DHARMa (un graphique par indice) : ", dossier_diag)
