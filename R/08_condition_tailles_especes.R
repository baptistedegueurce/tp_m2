# =============================================================================
# 08_condition_tailles_especes.R
# Qu'est-ce qui fait varier la condition et la taille des individus, dans
# l'espace et dans le temps ? Mêmes modèles simples pour chaque espèce.
#
# Entrées : data/processed/preparation.rds          (script 02)
#           data/processed/environnement_traits.csv (script 06)
#
# Deux réponses, à l'échelle de l'INDIVIDU :
#   - condition : log10(poids) ~ log10(taille) + covariables
#       = différence de poids À TAILLE ÉGALE (approche « ANCOVA », préférée
#         à la modélisation directe d'un ratio comme Kn ou K de Fulton)
#   - taille    : log10(taille) ~ covariables
# Effet aléatoire (1 | trait_id) : les poissons d'un même trait de chalut
# partagent le même environnement et ne sont pas indépendants.
#
# Modèles candidats (mêmes pour toutes les espèces, comparés par AICc) :
#   nul        : rien (en plus de la taille pour la condition)
#   spatial    : profondeur (gradient côte-large)
#   saison     : saison + profondeur                   <- modèle de référence
#   env_<var>  : <variable Copernicus> + profondeur    (une seule à la fois :
#                elles sont toutes corrélées entre elles et à la saison)
# Un modèle env_ aussi bon que « saison » (ΔAICc < 2) indique que la
# variable peut expliquer l'effet saison ; il ne prouve pas la causalité.
#
# Covariables continues centrées-réduites sur l'ensemble des traits : les
# effets sont exprimés en % de variation pour +1 écart-type, comparables
# entre espèces et entre variables.
#
# Option : restreindre_juveniles (00_config.R) pour ne comparer que des
# juvéniles entre saisons.
#
# Diagnostics : avec plusieurs centaines d'individus, DHARMa détecte de
# petits écarts ; les tailles mélangent des cohortes (distribution
# multimodale), donc la loi normale n'est qu'approchée pour la réponse
# « taille ». Les effets estimés restent interprétables ; les p-values sont
# approximatives.
#
# Sorties : data/processed/condition_especes/
#           outputs/figures/condition_especes/ (+ diagnostics/)
# =============================================================================

pkgs <- c("here", "dplyr", "tidyr", "ggplot2", "glmmTMB", "DHARMa", "emmeans")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants)) install.packages(manquants)

library(dplyr)
library(ggplot2)

source(here::here("R", "00_config.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "modeles_utils.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "nourricerie_utils.R"), encoding = "UTF-8")

set.seed(graine_aleatoire)
dossier_out  <- file.path(dossiers$processed, "condition_especes")
dossier_fig  <- here::here("outputs", "figures", "condition_especes")
dossier_diag <- file.path(dossier_fig, "diagnostics")
for (d in c(dossier_out, dossier_fig, dossier_diag)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
saisons <- c("printemps", "ete", "automne")
couleurs_saison <- c(printemps = "#66a61e", ete = "#e6ab02", automne = "#d95f02")

prep <- readRDS(file.path(dossiers$processed, "preparation.rds"))
fichier_env <- file.path(dossiers$processed, "environnement_traits.csv")
if (!file.exists(fichier_env)) stop("Lancer d'abord R/06_variables_environnementales.R")
env <- read.csv(fichier_env, encoding = "UTF-8")
noms <- with(read.csv(fichier_traits_fonctionnels, encoding = "UTF-8"),
             setNames(nom_commun, taxon))

# =============================================================================
# 1. Covariables : choix et standardisation
# =============================================================================
env_a <- env |> filter(analysable)

# Covariables spatiales : on retire la seconde si trop corrélée à la première
spatiales <- intersect(covariables_spatiales, names(env))
if (length(spatiales) > 1) {
  r <- cor(env_a[[spatiales[1]]], env_a[[spatiales[2]]], method = "spearman",
           use = "complete.obs")
  if (abs(r) > seuil_correlation) {
    message(sprintf("%s et %s sont corrélées (rho = %.2f) : %s seule est gardée.",
                    spatiales[1], spatiales[2], r, spatiales[1]))
    spatiales <- spatiales[1]
  }
}

# Variables Copernicus disponibles pour tous les traits analysables
cop <- intersect(variables_copernicus_traits, names(env))
cop <- cop[vapply(cop, function(v) all(!is.na(env_a[[v]])) &&
                    sd(env_a[[v]]) > 0, logical(1))]
if (!length(cop)) {
  message("Aucune variable Copernicus : seuls les modèles saison / profondeur ",
          "sont ajustés (lancer 01 puis 06 pour les ajouter).")
}

# Centrage-réduction sur les traits analysables (même échelle pour toutes
# les espèces). Les moyennes / écarts-types sont exportés pour revenir aux
# unités d'origine.
a_reduire <- c(spatiales, cop)
echelles <- data.frame(variable = a_reduire,
                       moyenne = sapply(a_reduire, function(v) mean(env_a[[v]])),
                       ecart_type = sapply(a_reduire, function(v) sd(env_a[[v]])))
for (v in a_reduire) {
  env[[paste0(v, "_z")]] <- (env[[v]] - echelles$moyenne[echelles$variable == v]) /
    echelles$ecart_type[echelles$variable == v]
}

# =============================================================================
# 2. Individus
# =============================================================================
ind <- prep$individus |>
  filter(groupe == "poisson", !aberrant_taille_poids,
         taille_mm > 0, poids_g > 0) |>
  inner_join(env |> filter(analysable) |>
               select(trait_id, station, lat, lon, ends_with("_z"),
                      all_of(a_reduire)),
             by = c("trait_id", "station")) |>
  mutate(saison = factor(saison, levels = saisons),
         trait_id = factor(trait_id),
         log_poids = log10(poids_g),
         log_taille = log10(taille_mm))

if (restreindre_juveniles) {
  lm_mat <- read.csv(fichier_traits_fonctionnels, encoding = "UTF-8") |>
    select(taxon, taille_maturite_cm)
  avant <- nrow(ind)
  ind <- ind |> left_join(lm_mat, by = "taxon") |>
    filter(!is.na(taille_maturite_cm), taille_mm < 10 * taille_maturite_cm)
  message(sprintf("Juvéniles seulement : %d individus sur %d", nrow(ind), avant))
}

especes <- ind |> count(taxon) |> filter(n >= n_min_individus_espece) |> pull(taxon)
message("Espèces analysées (>= ", n_min_individus_espece, " individus) : ",
        paste(especes, collapse = ", "))

# Kn (descriptif, pour les figures)
ind <- condition_le_cren(ind |> mutate(aberrant_taille_poids = FALSE),
                         n_min = 20)$individus

# =============================================================================
# 3. Modèles par espèce
# =============================================================================
terme_spatial <- paste(paste0(spatiales, "_z"), collapse = " + ")

formules <- function(reponse) {
  base <- if (reponse == "condition") "log_poids ~ log_taille_c" else "log_taille ~ 1"
  f <- c(nul     = base,
         spatial = paste(base, "+", terme_spatial),
         saison  = paste(base, "+ saison +", terme_spatial))
  for (v in cop) f[paste0("env_", v)] <- paste(base, "+", paste0(v, "_z"), "+", terme_spatial)
  lapply(f, function(x) as.formula(paste(x, "+ (1 | trait_id)")))
}

# Effet multiplicatif : coefficient b sur log10 -> % de variation
pct <- function(b) 100 * (10^b - 1)

res <- list(selection = list(), tests = list(), effets = list(),
            saisons = list(), diagnostics = list())

for (tx in especes) for (reponse in c("condition", "taille")) {
  cle <- paste(tx, reponse, sep = " | ")
  message("--- ", cle)
  d <- ind |> filter(taxon == tx) |>
    mutate(saison = droplevels(saison),
           log_taille_c = log_taille - mean(log_taille))
  f <- formules(reponse)
  if (nlevels(d$saison) < 2) f$saison <- NULL

  mods <- lapply(f, ajuster_glmm, data = d, famille = gaussian())
  sel <- selection_aicc(mods)
  res$selection[[cle]] <- cbind(taxon = tx, reponse = reponse, n = nrow(d), sel)

  # Modèle de référence : saison + spatial
  m <- mods$saison
  if (is.null(m)) next
  t <- tests_lrt(m)
  res$tests[[cle]] <- cbind(taxon = tx, reponse = reponse, n = nobs(m),
                            n_traits = nlevels(droplevels(d$trait_id)), t)

  # Effets continus (% pour +1 écart-type), modèle saison + chaque modèle env
  for (nm in grep("^saison$|^env_", names(mods), value = TRUE)) {
    mm <- mods[[nm]]
    if (is.null(mm)) next
    co <- coefficients_modele(mm)
    co <- co[grepl("_z$", co$terme), ]
    if (nm != "saison") co <- co[co$terme == paste0(sub("^env_", "", nm), "_z"), ]
    if (!nrow(co)) next
    res$effets[[paste(cle, nm)]] <- data.frame(
      taxon = tx, reponse = reponse, modele = nm,
      variable = sub("_z$", "", co$terme),
      effet_pct_par_ET = pct(co$estimation),
      IC95_bas = pct(co$estimation - 1.96 * co$erreur_type),
      IC95_haut = pct(co$estimation + 1.96 * co$erreur_type),
      p = co$p,
      delta_AICc_modele = sel$delta_AICc[sel$modele == nm])
  }

  # Saisons : moyennes marginales -> % d'écart à la moyenne de l'espèce
  em <- as.data.frame(emmeans::emmeans(m, ~ saison))
  ref <- mean(em$emmean)
  res$saisons[[cle]] <- data.frame(
    taxon = tx, reponse = reponse, saison = em$saison,
    ecart_pct = pct(em$emmean - ref),
    IC95_bas = pct(em[[grep("lower|LCL", names(em))[1]]] - ref),
    IC95_haut = pct(em[[grep("upper|UCL", names(em))[1]]] - ref),
    n = as.vector(table(d$saison)[as.character(em$saison)]))

  diag <- diagnostics_dharma(
    m, d, file.path(dossier_diag, paste0(make.names(tx), "_", reponse, ".png")),
    groupe = "trait_id")
  res$diagnostics[[cle]] <- data.frame(
    taxon = tx, reponse = reponse, diag,
    variance_trait = unname(unlist(glmmTMB::VarCorr(m)$cond$trait_id)),
    avertissements = paste(attr(m, "avertissements"), collapse = " | "))
}

lier <- function(x) { x <- do.call(rbind, x); rownames(x) <- NULL; x }
tabs <- lapply(res, lier)
# Holm sur les tests d'intérêt seulement (saison, covariables) : l'effet de
# la taille sur le poids, évident, n'entre pas dans la correction
interet <- tabs$tests$effet != "log_taille_c"
tabs$tests$p_holm <- NA_real_
tabs$tests$p_holm[interet] <- p.adjust(tabs$tests$p[interet], method = "holm")
if (!is.null(tabs$effets)) tabs$effets$p_holm <- p.adjust(tabs$effets$p, method = "holm")

# Résumé : une ligne par espèce × réponse
tabs$resume <- tabs$selection |>
  group_by(taxon, reponse, n) |>
  summarise(meilleur_modele = modele[which.min(AICc)],
            delta_AICc_saison = delta_AICc[modele == "saison"][1],
            .groups = "drop") |>
  left_join(tabs$tests |> filter(effet != "log_taille_c") |>
              select(taxon, reponse, effet, p_holm) |>
              tidyr::pivot_wider(names_from = effet, values_from = p_holm,
                                 names_prefix = "p_holm_"),
            by = c("taxon", "reponse"))

for (nom in names(tabs)) {
  write.csv(tabs[[nom]], file.path(dossier_out, paste0(nom, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}
write.csv(echelles, file.path(dossier_out, "echelles_covariables.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# 4. Figures
# =============================================================================
lab_esp <- function(tx) ifelse(is.na(noms[tx]), tx, paste0(noms[tx], "\n(", tx, ")"))
lab_rep <- c(condition = "Condition (poids à taille égale)",
             taille = "Taille")
sauver <- function(p, nom, l, h) {
  print(p)
  ggsave(file.path(dossier_fig, nom), p, width = l, height = h, dpi = 150,
         bg = "white")
}

# 4a. Effet de la saison (% d'écart à la moyenne de l'espèce)
p_sais <- tabs$saisons |>
  mutate(taxon = lab_esp(taxon)) |>
  ggplot(aes(saison, ecart_pct, colour = saison)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = IC95_bas, ymax = IC95_haut)) +
  facet_grid(reponse ~ taxon, scales = "free_y",
             labeller = labeller(reponse = lab_rep)) +
  scale_colour_manual(values = couleurs_saison, guide = "none") +
  labs(title = "Effet de la saison par espèce (modèle saison + profondeur)",
       subtitle = "Écart à la moyenne de l'espèce (%), IC 95 % ; condition = poids à taille égale",
       x = NULL, y = "Écart (%)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
sauver(p_sais, "effet_saison_especes.png", 12, 6.5)

# 4b. Effets des variables continues (% pour +1 écart-type)
if (!is.null(tabs$effets) && nrow(tabs$effets)) {
  eff <- tabs$effets |>
    mutate(source = ifelse(modele == "saison", "modèle saison + profondeur",
                           "modèle variable env. + profondeur"),
           taxon = lab_esp(taxon),
           significatif = p_holm < 0.05)
  p_eff <- ggplot(eff, aes(effet_pct_par_ET, taxon, colour = variable,
                           shape = significatif)) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_pointrange(aes(xmin = IC95_bas, xmax = IC95_haut),
                    position = position_dodge(width = 0.6)) +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       labels = c(`TRUE` = "p Holm < 0,05", `FALSE` = "n.s.")) +
    facet_wrap(~ reponse, scales = "free_x", labeller = as_labeller(lab_rep)) +
    labs(title = "Effet de +1 écart-type de chaque variable, par espèce",
         subtitle = "Profondeur : modèle saison + profondeur ; variables Copernicus : chacune seule + profondeur",
         x = "Variation (%) pour +1 écart-type", y = NULL, colour = NULL,
         shape = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", legend.box = "vertical")
  sauver(p_eff, "effets_variables_especes.png", 11, 6.5)
}

# 4c. Comparaison des modèles (ΔAICc)
p_aic <- tabs$selection |>
  mutate(taxon = lab_esp(taxon),
         modele = factor(modele, levels = unique(c("nul", "spatial", "saison",
                                                   paste0("env_", cop))))) |>
  ggplot(aes(modele, taxon, fill = pmin(delta_AICc, 20))) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(delta_AICc, 1)), size = 3) +
  scale_fill_gradient(low = "#08519c", high = "grey95", limits = c(0, 20),
                      name = "Delta AICc") +
  facet_wrap(~ reponse, labeller = as_labeller(lab_rep)) +
  labs(title = "Comparaison des modèles par espèce (Delta AICc : 0 = meilleur ; < 2 = équivalent)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
sauver(p_aic, "comparaison_modeles.png", 11, 5.5)

# 4d. Données brutes : Kn selon la profondeur, par saison
p_brut <- ind |>
  filter(taxon %in% especes, !is.na(Kn)) |>
  mutate(taxon = lab_esp(taxon)) |>
  ggplot(aes(.data[[spatiales[1]]], Kn, colour = saison)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey60") +
  geom_point(alpha = 0.25, size = 0.9) +
  geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  facet_wrap(~ taxon) +
  scale_colour_manual(values = couleurs_saison) +
  coord_cartesian(ylim = c(0.6, 1.4)) +
  labs(title = "Condition relative (Kn) des individus selon la profondeur",
       x = spatiales[1], y = "Kn de Le Cren", colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
sauver(p_brut, "kn_profondeur_brut.png", 10, 6.5)

# =============================================================================
# 5. Résumé
# =============================================================================
message("\n================ RÉSUMÉ ================")
print(tabs$resume |> mutate(across(where(is.numeric), ~ signif(.x, 3))) |>
        as.data.frame())
message("\nDétails : ", dossier_out)
