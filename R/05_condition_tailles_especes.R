# =============================================================================
# 05_condition_tailles_especes.R
# Qu'est-ce qui fait varier la condition et la taille des poissons, dans
# l'espace (profondeur) et dans le temps (saison) ? Même modèle simple pour
# chaque espèce.
#
# Entrées : data/processed/preparation.rds          (script 02)
#           data/processed/environnement_traits.csv (script 03)
#           data/ref/traits_fonctionnels.csv        (taille de 1re maturité)
#
# Espèces : au moins `n_min_individus_espece` individus mesurés et pesés.
# Saisons : pour chaque espèce, seules les saisons où elle a été prise dans
#   au moins `n_min_traits_saison` traits entrent dans les tests (sinon
#   l'effet « saison » se confond avec l'effet d'un ou deux traits).
#
# 1. Description : relation taille-poids W = a L^b, condition relative de
#    Le Cren (Kn), distributions de tailles
# 2. Modèles linéaires mixtes (lme4), à l'échelle de l'individu :
#      condition : log10(W) ~ log10(L) + saison + profondeur + (1 | trait)
#                  -> différence de poids À TAILLE ÉGALE (ANCOVA)
#      taille    : log10(L) ~ saison + profondeur + (1 | trait)
#    (1 | trait) : les poissons d'un même trait de chalut ne sont pas
#    indépendants (même lieu, même jour, même passage du chalut).
#    Profondeur centrée, par tranche de 10 m -> effets en % pour +10 m.
#    Test de chaque effet (rapport de vraisemblance), Tukey entre saisons,
#    résidus vérifiés avec DHARMa.
# 3. Croissance apparente des juvéniles entre l'été et l'automne : décalage
#    de la taille médiane / nombre de jours
#
# Sorties : data/processed/condition_especes/
#           outputs/figures/condition_especes/ (+ diagnostics/)
# =============================================================================

pkgs <- c("here", "dplyr", "tidyr", "ggplot2", "lme4", "DHARMa", "emmeans")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants)) install.packages(manquants)

library(dplyr)
library(ggplot2)

source(here::here("R", "00_config.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "modeles_simples.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "nourricerie_utils.R"), encoding = "UTF-8")

set.seed(graine_aleatoire)
dossier_out  <- file.path(dossiers$processed, "condition_especes")
dossier_fig  <- here::here("outputs", "figures", "condition_especes")
dossier_diag <- file.path(dossier_fig, "diagnostics")
unlink(list.files(dossier_out, full.names = TRUE))   # anciennes sorties
unlink(list.files(dossier_fig, pattern = "\\.png$", full.names = TRUE))
for (d in c(dossier_out, dossier_fig, dossier_diag)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
saisons <- c("printemps", "ete", "automne")
couleurs_saison <- c(printemps = "#66a61e", ete = "#e6ab02", automne = "#d95f02")

prep <- readRDS(file.path(dossiers$processed, "preparation.rds"))
fichier_env <- file.path(dossiers$processed, "environnement_traits.csv")
if (!file.exists(fichier_env)) stop("Lancer d'abord R/03_variables_environnementales.R")
env <- read.csv(fichier_env, encoding = "UTF-8")
traits_ref <- read.csv(fichier_traits_fonctionnels, encoding = "UTF-8")
noms <- setNames(traits_ref$nom_commun, traits_ref$taxon)
lab_esp <- function(tx) ifelse(is.na(noms[tx]), tx, paste0(noms[tx], "\n(", tx, ")"))

# =============================================================================
# 1. Individus
# =============================================================================
prof_moy <- mean(env[[covariable_spatiale]][env$analysable])

ind <- prep$individus |>
  filter(groupe == "poisson", !aberrant_taille_poids,
         taille_mm > 0, poids_g > 0) |>
  inner_join(env |> filter(analysable) |>
               select(trait_id, station, prof = all_of(covariable_spatiale)),
             by = c("trait_id", "station")) |>
  left_join(traits_ref |> select(taxon, taille_maturite_cm), by = "taxon") |>
  mutate(saison = factor(saison, levels = saisons),
         trait_id = factor(trait_id),
         prof10 = (prof - prof_moy) / 10,        # +1 = 10 m plus profond
         log_poids = log10(poids_g),
         log_taille = log10(taille_mm),
         juvenile = taille_mm < 10 * taille_maturite_cm)

especes <- ind |> count(taxon) |> filter(n >= n_min_individus_espece) |> pull(taxon)
message("Espèces analysées (>= ", n_min_individus_espece, " individus) : ",
        paste(especes, collapse = ", "))
ind <- ind |> filter(taxon %in% especes)
ind <- condition_le_cren(ind, n_min = 20)$individus   # Kn (descriptif)

# Effectifs : individus et traits par espèce et saison
effectifs <- ind |>
  group_by(taxon, saison) |>
  summarise(n_individus = n(), n_traits = n_distinct(trait_id),
            part_juveniles = mean(juvenile), L_min = min(taille_mm),
            L_max = max(taille_mm), .groups = "drop") |>
  mutate(saison_testee = n_traits >= n_min_traits_saison)

# Relation taille-poids (toutes saisons) : W = a L^b
rel_tp <- bind_rows(lapply(especes, function(tx) {
  d <- ind[ind$taxon == tx, ]
  fit <- lm(log_poids ~ log_taille, d)
  ic <- confint(fit)[2, ]
  data.frame(taxon = tx, nom_commun = noms[[tx]], n = nrow(d),
             L_min_mm = min(d$taille_mm), L_max_mm = max(d$taille_mm),
             a = 10^coef(fit)[[1]], b = coef(fit)[[2]],
             b_IC95_bas = ic[[1]], b_IC95_haut = ic[[2]],
             r2 = summary(fit)$r.squared)
}))

# =============================================================================
# 2. Modèles
# =============================================================================
pct <- function(b) 100 * (10^b - 1)   # coefficient sur log10 -> % de variation

res <- list(tests = list(), saisons = list(), contrastes = list(),
            profondeur = list(), diagnostics = list())

for (tx in especes) {
  sais_ok <- effectifs$saison[effectifs$taxon == tx & effectifs$saison_testee]
  d <- ind |>
    filter(taxon == tx, saison %in% sais_ok) |>
    mutate(saison = droplevels(saison), trait_id = droplevels(trait_id),
           log_taille_c = log_taille - mean(log_taille))
  avec_saison <- nlevels(d$saison) >= 2

  for (reponse in c("condition", "taille")) {
    cle <- paste(tx, reponse, sep = " | ")
    f <- if (reponse == "condition") "log_poids ~ log_taille_c" else "log_taille ~ 1"
    if (avec_saison) f <- paste(f, "+ saison")
    f <- as.formula(paste(f, "+ prof10 + (1 | trait_id)"))
    m <- ajuster(f, d)
    if (is.null(m)) { warning("Échec du modèle : ", cle); next }

    t <- tester(m) |> filter(effet != "log_taille_c")
    res$tests[[cle]] <- data.frame(
      taxon = tx, reponse = reponse,
      saisons_testees = paste(levels(d$saison), collapse = "/"),
      n_individus = nrow(d), n_traits = nlevels(d$trait_id),
      t |> mutate(effet = ifelse(effet == "prof10", "profondeur", effet)),
      avertissements = paste(attr(m, "avertissements"), collapse = " | "))

    # Profondeur : % pour +10 m (IC 95 % : estimation ± 1,96 erreur-type)
    co <- summary(m)$coefficients["prof10", ]
    res$profondeur[[cle]] <- data.frame(
      taxon = tx, reponse = reponse,
      effet_pct_10m = pct(co[[1]]),
      IC95_bas = pct(co[[1]] - 1.96 * co[[2]]),
      IC95_haut = pct(co[[1]] + 1.96 * co[[2]]),
      p = t$p[t$effet == "prof10"],
      gamme_profondeur_m = paste(round(range(d$prof)), collapse = "-"))

    # Saisons : écart à la moyenne de l'espèce (%), contrastes de Tukey (%)
    if (avec_saison) {
      em <- emmeans::emmeans(m, ~ saison)
      emd <- as.data.frame(em)
      ref <- mean(emd$emmean)
      res$saisons[[cle]] <- data.frame(
        taxon = tx, reponse = reponse, saison = emd$saison,
        ecart_pct = pct(emd$emmean - ref),
        IC95_bas = pct(emd[[grep("lower|LCL", names(emd))[1]]] - ref),
        IC95_haut = pct(emd[[grep("upper|UCL", names(emd))[1]]] - ref),
        valeur_mm = if (reponse == "taille") 10^emd$emmean else NA,
        n_individus = as.vector(table(d$saison)[as.character(emd$saison)]),
        n_traits = as.vector(tapply(d$trait_id, d$saison,
                                    function(x) length(unique(x)))[as.character(emd$saison)]))
      ctr <- as.data.frame(pairs(em, adjust = "tukey"))
      res$contrastes[[cle]] <- data.frame(
        taxon = tx, reponse = reponse, contraste = ctr$contrast,
        difference_pct = pct(ctr$estimate), p_tukey = ctr$p.value)
    }

    res$diagnostics[[cle]] <- data.frame(
      taxon = tx, reponse = reponse,
      diagnostic(m, file.path(dossier_diag, paste0(make.names(tx), "_", reponse, ".png"))),
      variance_trait = as.data.frame(lme4::VarCorr(m))$vcov[1],
      variance_residuelle = sigma(m)^2)
  }
}
tabs <- lapply(res, lier)
tabs$tests$significatif <- tabs$tests$p < 0.05

# =============================================================================
# 3. Croissance apparente des juvéniles (été -> automne)
# =============================================================================
# Hypothèse : on suit la même cohorte, sans arrivée ni départ sélectif en
# taille entre les deux campagnes -> « croissance apparente ».
# Calculée seulement si l'espèce est testée en été ET en automne.
mediane_juv <- function(dd) quantile_pond(dd$taille_mm, dd$facteur_elevation, 0.5)
croissance <- bind_rows(lapply(especes, function(tx) {
  ok <- effectifs |> filter(taxon == tx, saison %in% c("ete", "automne"), saison_testee)
  if (nrow(ok) < 2) return(NULL)
  dj <- ind |> filter(taxon == tx, juvenile)
  e <- dj[dj$saison == "ete", ]; a <- dj[dj$saison == "automne", ]
  jours <- as.numeric(mean(as.Date(a$date)) - mean(as.Date(e$date)))
  data.frame(taxon = tx, nom_commun = noms[[tx]],
             n_juv_ete = nrow(e), n_juv_automne = nrow(a),
             L_med_ete_mm = mediane_juv(e), L_med_automne_mm = mediane_juv(a),
             jours = round(jours),
             gain_mm = mediane_juv(a) - mediane_juv(e),
             croissance_mm_jour = (mediane_juv(a) - mediane_juv(e)) / jours)
}))

tabs$effectifs <- effectifs
tabs$relation_taille_poids <- rel_tp
tabs$croissance <- croissance
for (nom in names(tabs)) {
  write.csv(tabs[[nom]], file.path(dossier_out, paste0(nom, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}

# =============================================================================
# 4. Figures
# =============================================================================
lab_rep <- c(condition = "Condition (poids à taille égale)", taille = "Taille")
sauver <- function(p, nom, l, h) {
  print(p)
  ggsave(file.path(dossier_fig, nom), p, width = l, height = h, dpi = 150,
         bg = "white")
}
etiquette_p <- function(p) ifelse(p < 0.001, "p < 0,001",
                                  paste0("p = ", formatC(p, digits = 2, format = "fg",
                                                         decimal.mark = ",")))

# 4a. Effet de la saison (% d'écart à la moyenne de l'espèce)
d_sais <- tabs$saisons |>
  mutate(saison = factor(as.character(saison), levels = saisons),
         espece = lab_esp(taxon))
etiquettes <- tabs$tests |> filter(effet == "saison") |>
  mutate(espece = lab_esp(taxon), texte = paste("saison :", etiquette_p(p)))
p_sais <- ggplot(d_sais, aes(saison, ecart_pct, colour = saison)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = IC95_bas, ymax = IC95_haut)) +
  geom_text(aes(y = IC95_haut, label = paste0(n_traits, " tr.")), vjust = -0.6,
            size = 2.6, colour = "grey30") +
  geom_text(data = etiquettes, aes(x = 0.6, y = Inf, label = texte),
            inherit.aes = FALSE, hjust = 0, vjust = 1.4, size = 3, colour = "grey20") +
  facet_grid(reponse ~ espece, scales = "free_y",
             labeller = labeller(reponse = lab_rep)) +
  scale_colour_manual(values = couleurs_saison, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.2))) +
  labs(title = "Effet de la saison par espèce (modèle saison + profondeur)",
       subtitle = paste0("Écart à la moyenne de l'espèce (%), IC 95 % ; « tr. » = traits de chalut ; ",
                         "saisons avec < ", n_min_traits_saison, " traits non testées"),
       x = NULL, y = "Écart (%)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
sauver(p_sais, "effet_saison_especes.png", 12, 6.5)

# 4b. Effet de la profondeur (% pour +10 m)
p_prof <- tabs$profondeur |>
  mutate(espece = lab_esp(taxon),
         signif = ifelse(p < 0.05, "p < 0,05", "non significatif")) |>
  ggplot(aes(effet_pct_10m, espece, colour = signif)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = IC95_bas, xmax = IC95_haut)) +
  facet_wrap(~ reponse, scales = "free_x", labeller = as_labeller(lab_rep)) +
  scale_colour_manual(values = c(`p < 0,05` = "#1b7837", `non significatif` = "grey55")) +
  labs(title = "Effet de la profondeur par espèce (modèle saison + profondeur)",
       subtitle = "Variation (%) pour 10 m de profondeur en plus, IC 95 %",
       x = "Variation pour +10 m (%)", y = NULL, colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
sauver(p_prof, "effet_profondeur_especes.png", 10, 5)

# 4c. Condition relative (Kn) par saison, saisons non testées en pâle
p_kn <- ind |>
  left_join(effectifs |> select(taxon, saison, saison_testee), by = c("taxon", "saison")) |>
  mutate(espece = lab_esp(taxon)) |>
  ggplot(aes(saison, Kn, colour = saison, alpha = saison_testee)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.7) +
  facet_wrap(~ espece, nrow = 1) +
  scale_colour_manual(values = couleurs_saison, guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 0.6, `FALSE` = 0.15), guide = "none") +
  coord_cartesian(ylim = c(0.6, 1.4)) +
  labs(title = "Condition relative de Le Cren (Kn) par saison",
       subtitle = paste0("Kn = 1 : poids attendu pour la taille. Pâle : saison non testée (< ",
                         n_min_traits_saison, " traits)"),
       x = NULL, y = "Kn") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
sauver(p_kn, "condition_par_saison.png", 12, 4.5)

# 4d. Distributions de tailles par saison, taille de maturité
lm_lignes <- ind |> distinct(taxon, taille_maturite_cm) |>
  mutate(espece = lab_esp(taxon), L50 = 10 * taille_maturite_cm)
p_dist <- ind |>
  group_by(taxon, saison) |>
  mutate(poids_densite = facteur_elevation / sum(facteur_elevation)) |>
  ungroup() |>
  mutate(espece = lab_esp(taxon)) |>
  ggplot(aes(taille_mm, colour = saison, fill = saison)) +
  geom_density(aes(weight = poids_densite), alpha = 0.15, adjust = 0.8) +
  geom_vline(data = lm_lignes, aes(xintercept = L50), linetype = 2, colour = "grey30") +
  facet_wrap(~ espece, scales = "free", nrow = 1) +
  scale_colour_manual(values = couleurs_saison) +
  scale_fill_manual(values = couleurs_saison) +
  labs(title = "Distribution des tailles par saison (effectifs élevés)",
       subtitle = "Tirets : taille de première maturité (FishBase) ; à gauche = juvéniles",
       x = "Longueur (mm)", y = "Densité", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
sauver(p_dist, "distributions_tailles.png", 13, 4.5)

# =============================================================================
# 5. Résumé
# =============================================================================
message("\n================ RÉSUMÉ ================")
print(effectifs |> mutate(across(where(is.numeric), ~ round(.x, 2))) |> as.data.frame())
print(tabs$tests |> select(taxon, reponse, effet, saisons_testees, n_traits, chi2, ddl, p) |>
        mutate(across(where(is.numeric), ~ signif(.x, 2))) |> as.data.frame())
print(tabs$contrastes |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
print(tabs$profondeur |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
print(tabs$diagnostics |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
print(croissance |> mutate(across(where(is.numeric), ~ round(.x, 2))))
