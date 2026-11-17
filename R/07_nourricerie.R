# =============================================================================
# 07_nourricerie.R
# Structure en taille, part de juvéniles, condition et croissance au cours de
# la saison de nourricerie
#
# Question : comment la structure fonctionnelle et en taille des communautés
# de juvéniles évolue-t-elle au cours de la saison de nourricerie ?
#
# Pourquoi des traits MESURÉS : les poissons capturés sont surtout des
# juvéniles (ex. soles à 11 % de leur taille maximale). Les traits tirés de la
# littérature décrivent des adultes ; la taille et le poids mesurés sur les
# 1 537 individus décrivent les animaux réellement présents.
#
# Entrées : data/processed/preparation.rds (script 02)
#           data/ref/traits_fonctionnels.csv (taille de première maturité)
#
# 1. Stade de chaque individu : juvénile si taille < taille de 1re maturité
# 2. Indicateurs par trait de chalut : part et densité de juvéniles, taille
#    moyenne, taille relative (L / L maturité), pente du spectre de taille,
#    condition moyenne
# 3. Modèles mixtes de ces indicateurs (même démarche que 05b)
# 4. Par espèce : distributions de tailles, condition (Le Cren) par saison,
#    progression des tailles des juvéniles (croissance apparente)
#
# Sorties : data/processed/nourricerie/, outputs/figures/nourricerie/
# =============================================================================

pkgs <- c("here", "dplyr", "tidyr", "ggplot2", "glmmTMB", "DHARMa", "emmeans")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants)) install.packages(manquants)

library(dplyr)
library(ggplot2)

source(here::here("R", "00_config.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "modeles_utils.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "figures_modeles.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "nourricerie_utils.R"), encoding = "UTF-8")

set.seed(graine_aleatoire)
dossier_out  <- file.path(dossiers$processed, "nourricerie")
dossier_fig  <- here::here("outputs", "figures", "nourricerie")
dossier_diag <- file.path(dossier_fig, "diagnostics")
for (d in c(dossier_out, dossier_fig, dossier_diag)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
saisons <- c("printemps", "ete", "automne")
couleurs_saison <- c(printemps = "#66a61e", ete = "#e6ab02", automne = "#d95f02")

prep <- readRDS(file.path(dossiers$processed, "preparation.rds"))
traits_ref <- read.csv(fichier_traits_fonctionnels, encoding = "UTF-8",
                       stringsAsFactors = FALSE)

# =============================================================================
# 1. Individus : stade et condition
# =============================================================================
ind <- prep$individus |>
  filter(groupe == "poisson", !is.na(taille_mm), !is.na(poids_g)) |>
  left_join(traits_ref |> select(taxon, nom_commun, taille_maturite_cm),
            by = "taxon") |>
  left_join(prep$traits |> select(trait_id, distance_km, analysable),
            by = "trait_id") |>
  filter(analysable) |>
  mutate(
    saison = factor(saison, levels = saisons),
    L_maturite_mm = 10 * taille_maturite_cm,
    L_rel  = taille_mm / L_maturite_mm,
    stade  = ifelse(is.na(L_rel), NA, ifelse(L_rel < 1, "juvénile", "adulte"))
  )

sans_lm <- unique(ind$taxon[is.na(ind$L_maturite_mm)])
if (length(sans_lm)) {
  message("Taxons sans taille de maturité (stade non attribué) : ",
          paste(sans_lm, collapse = ", "))
}

cond <- condition_le_cren(ind, n_min = 20)
ind <- cond$individus
ind_ok <- filter(ind, !aberrant_taille_poids)

# =============================================================================
# 2. Indicateurs par trait de chalut
# =============================================================================
par_trait <- ind_ok |>
  group_by(trait_id) |>
  summarise(
    n_mesures      = n(),
    effectif_eleve = sum(facteur_elevation),
    part_juveniles = moy_pond(stade == "juvénile", facteur_elevation),
    dens_juv_km    = sum(facteur_elevation[stade %in% "juvénile"]) / first(distance_km),
    L_moy_mm       = 10^moy_pond(log10(taille_mm), facteur_elevation),
    L_rel_moy      = moy_pond(L_rel, facteur_elevation),
    S_juv          = n_distinct(taxon[stade %in% "juvénile"]),
    Kn_moy         = moy_pond(Kn, facteur_elevation),
    pente_nbss     = pente_nbss(poids_g, facteur_elevation, first(distance_km))[["pente"]],
    .groups = "drop"
  ) |>
  left_join(prep$traits |> select(trait_id, saison, station, date, profondeur_m,
                                  lat, lon, douteux, fune_m),
            by = "trait_id") |>
  mutate(saison = factor(saison, levels = saisons), station = factor(station)) |>
  relocate(trait_id, saison, station, date, profondeur_m)

# =============================================================================
# 3. Modèles mixtes des indicateurs
# =============================================================================
reponses <- tribble(
  ~indice,          ~famille,   ~transformation, ~label,
  "part_juveniles", "beta",     "01",      "Part de juvéniles (abondance)",
  "dens_juv_km",    "gaussian", "log10p1", "Densité de juvéniles (log10 ind/km + 1)",
  "S_juv",          "genpois",  "aucune",  "Espèces présentes en juvéniles",
  "L_moy_mm",       "gaussian", "log10",   "Taille moyenne observée (log10 mm)",
  "L_rel_moy",      "gaussian", "aucune",  "Taille relative (L / L maturité)",
  "pente_nbss",     "gaussian", "aucune",  "Pente du spectre de taille (NBSS)",
  "Kn_moy",         "gaussian", "aucune",  "Condition moyenne (Kn de Le Cren)"
)

tabs <- modeliser_indices(par_trait, reponses, dossier_diag)
figures_modeles(par_trait, reponses, tabs, dossier_fig,
                sujet = "la structure en taille et les juvéniles")

# =============================================================================
# 4. Par espèce
# =============================================================================
n_min_espece <- 50
especes <- ind_ok |> count(taxon) |> filter(n >= n_min_espece) |> pull(taxon)
noms <- setNames(traits_ref$nom_commun, traits_ref$taxon)

# 4a. Tailles par espèce et saison (pondérées par l'élévation)
tailles <- ind_ok |>
  filter(taxon %in% especes) |>
  group_by(taxon, saison) |>
  summarise(
    n_mesures    = n(),
    date_moyenne = mean(as.Date(date)),
    L_q25        = quantile_pond(taille_mm, facteur_elevation, 0.25),
    L_mediane    = quantile_pond(taille_mm, facteur_elevation, 0.50),
    L_q75        = quantile_pond(taille_mm, facteur_elevation, 0.75),
    L_med_juveniles = quantile_pond(taille_mm[stade %in% "juvénile"],
                                    facteur_elevation[stade %in% "juvénile"], 0.5),
    part_juveniles = moy_pond(stade == "juvénile", facteur_elevation),
    Kn_moyen     = mean(Kn, na.rm = TRUE),
    .groups = "drop"
  )

# 4b. Croissance apparente des juvéniles entre saisons successives :
# déplacement de la taille médiane des juvéniles / nombre de jours.
# Hypothèse forte : même cohorte suivie (pas de recrutement, de départ ni de
# mortalité sélective en taille entre les deux dates). À présenter comme
# une "croissance apparente", pas une croissance individuelle.
croissance <- tailles |>
  arrange(taxon, saison) |>
  group_by(taxon) |>
  mutate(saison_precedente = lag(saison),
         jours = as.numeric(date_moyenne - lag(date_moyenne)),
         gain_mm = L_med_juveniles - lag(L_med_juveniles),
         croissance_mm_jour = gain_mm / jours) |>
  ungroup() |>
  filter(!is.na(saison_precedente)) |>
  mutate(interpretation = ifelse(
    gain_mm < 0,
    "Baisse de taille : arrivée d'une nouvelle cohorte ou départ des plus grands, pas une croissance",
    "Croissance apparente (si même cohorte)")) |>
  select(taxon, de = saison_precedente, a = saison, jours, gain_mm,
         croissance_mm_jour, interpretation)

# 4c. Condition par saison : log10(W) ~ log10(L) + saison + (1 | trait)
# Le trait de chalut en effet aléatoire : les poissons d'un même trait ne
# sont pas indépendants. L'effet saison = différence de poids à taille égale.
condition_saison <- list()
condition_contrastes <- list()
for (tx in especes) {
  d <- ind_ok |> filter(taxon == tx) |> mutate(saison = droplevels(saison))
  if (nlevels(d$saison) < 2) next
  m <- ajuster_glmm(log10(poids_g) ~ log10(taille_mm) + saison + (1 | trait_id),
                    d, gaussian())
  if (is.null(m)) next
  t <- tests_lrt(m)
  condition_saison[[tx]] <- data.frame(taxon = tx, n = nobs(m),
                                       saisons = paste(levels(d$saison), collapse = "/"),
                                       chi2 = t$chi2[t$effet == "saison"],
                                       ddl = t$ddl[t$effet == "saison"],
                                       p = t$p[t$effet == "saison"])
  ctr <- as.data.frame(pairs(emmeans::emmeans(m, ~ saison), adjust = "tukey"))
  condition_contrastes[[tx]] <- data.frame(
    taxon = tx, contraste = ctr$contrast,
    # différence de log10(poids) -> % de poids en plus à taille égale
    difference_poids_pct = 100 * (10^ctr$estimate - 1),
    p_tukey = ctr$p.value)
}
lier <- function(x) { x <- do.call(rbind, x); rownames(x) <- NULL; x }
condition_saison <- lier(condition_saison)
if (!is.null(condition_saison)) {
  condition_saison$p_holm <- p.adjust(condition_saison$p, method = "holm")
}
condition_contrastes <- lier(condition_contrastes)

# =============================================================================
# 5. Export
# =============================================================================
ecrire <- function(x, nom) {
  if (is.null(x)) return(invisible())
  write.csv(x, file.path(dossier_out, nom), row.names = FALSE,
            fileEncoding = "UTF-8")
}
ecrire(par_trait, "indicateurs_par_trait.csv")
for (nom in names(tabs)) ecrire(tabs[[nom]], paste0("modeles_", nom, ".csv"))
ecrire(tailles, "tailles_par_espece_saison.csv")
ecrire(croissance, "croissance_apparente.csv")
ecrire(cond$parametres, "relation_taille_poids.csv")
ecrire(condition_saison, "condition_tests.csv")
ecrire(condition_contrastes, "condition_contrastes.csv")

# =============================================================================
# 6. Figures par espèce
# =============================================================================
sauver <- function(p, nom, l, h) {
  print(p)
  ggsave(file.path(dossier_fig, nom), p, width = l, height = h, dpi = 150,
         bg = "white")
}
lab_esp <- function(tx) paste0(noms[tx], " (", tx, ")")

# 6a. Distributions de tailles par saison (pondérées), taille de maturité
dist <- ind_ok |>
  filter(taxon %in% especes) |>
  group_by(taxon, saison) |>
  mutate(poids_densite = facteur_elevation / sum(facteur_elevation)) |>
  ungroup()
lm_lignes <- dist |> distinct(taxon, L_maturite_mm) |> filter(!is.na(L_maturite_mm))
p_dist <- ggplot(dist, aes(taille_mm, colour = saison, fill = saison)) +
  geom_density(aes(weight = poids_densite), alpha = 0.15, adjust = 0.8) +
  geom_vline(data = lm_lignes, aes(xintercept = L_maturite_mm),
             linetype = 2, colour = "grey30") +
  facet_wrap(~ taxon, scales = "free", labeller = as_labeller(lab_esp)) +
  scale_colour_manual(values = couleurs_saison) +
  scale_fill_manual(values = couleurs_saison) +
  labs(title = "Distribution des tailles par saison (abondances élevées)",
       subtitle = "Tirets : taille de première maturité (FishBase) ; à gauche = juvéniles",
       x = "Longueur (mm)", y = "Densité", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
sauver(p_dist, "distributions_tailles.png", 11, 7)

# 6b. Condition (Kn) par espèce et saison
p_kn <- ind_ok |>
  filter(taxon %in% especes, !is.na(Kn)) |>
  ggplot(aes(saison, Kn, colour = saison)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
  facet_wrap(~ taxon, labeller = as_labeller(lab_esp)) +
  scale_colour_manual(values = couleurs_saison, guide = "none") +
  coord_cartesian(ylim = c(0.6, 1.4)) +
  labs(title = "Condition relative de Le Cren (Kn) par saison",
       subtitle = "Kn = 1 : poids attendu pour la taille ; > 1 : individus plus lourds",
       x = NULL, y = "Kn") +
  theme_minimal(base_size = 10)
sauver(p_kn, "condition_par_saison.png", 10, 6.5)

# 6c. Taille médiane (et intervalle interquartile) par saison
p_med <- ggplot(tailles, aes(saison, L_mediane, colour = taxon, group = taxon)) +
  geom_pointrange(aes(ymin = L_q25, ymax = L_q75),
                  position = position_dodge(width = 0.5)) +
  geom_line(position = position_dodge(width = 0.5), alpha = 0.6) +
  scale_colour_discrete(labels = lab_esp) +
  labs(title = "Taille médiane par campagne (barres : 25e - 75e centiles)",
       x = NULL, y = "Longueur (mm)", colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "right")
sauver(p_med, "tailles_medianes.png", 9, 5)

# =============================================================================
# 7. Résumé
# =============================================================================
message("\n================ RÉSUMÉ ================")
message(sprintf("%d poissons mesurés (aberrants exclus) dans %d traits ; %.0f %% de juvéniles (abondance élevée)",
                nrow(ind_ok), nrow(par_trait),
                100 * moy_pond(ind_ok$stade == "juvénile", ind_ok$facteur_elevation)))
message("\nIndicateurs par trait (modèles saison + profondeur + station) :")
print(tabs$resume |> select(indice, famille, n, meilleur_modele_AICc, p_saison_holm,
                            p_profondeur_holm, R2, probleme_diagnostic) |>
        mutate(across(where(is.numeric), ~ signif(.x, 3))))
message("\nTailles par espèce et saison :")
print(tailles |> select(taxon, saison, n_mesures, L_mediane, part_juveniles, Kn_moyen) |>
        mutate(across(where(is.numeric), ~ round(.x, 2))) |> as.data.frame())
message("\nCroissance apparente des juvéniles :")
print(croissance |> select(-interpretation) |>
        mutate(across(where(is.numeric), ~ round(.x, 2))) |> as.data.frame())
message("\nCondition : effet saison (à taille égale) :")
print(condition_saison |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
print(condition_contrastes |> mutate(across(where(is.numeric), ~ round(.x, 3))))
