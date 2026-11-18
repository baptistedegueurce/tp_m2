# =============================================================================
# 07_nourricerie.R
# Arrivée des juvéniles au cours de la saison de nourricerie (échelle du
# trait de chalut). Sert de CONTEXTE aux analyses par espèce du script 08.
#
# Entrées : data/processed/preparation.rds (script 02)
#           data/ref/traits_fonctionnels.csv (taille de première maturité)
#
# 1. Stade de chaque individu : juvénile si taille < taille de 1re maturité
#    (L50, FishBase)
# 2. Par trait : densité de juvéniles (ind/km, effectifs élevés), nombre
#    d'espèces présentes en juvéniles, part de juvéniles
# 3. Modèle mixte de référence : indicateur ~ saison + profondeur + (1 | station)
#    + test de sensibilité sans les traits douteux
#    + test à profondeur égale (traits côtiers seulement, Kruskal-Wallis),
#      car au printemps seuls des traits peu profonds sont analysables
# 4. Tableau descriptif : densité de juvéniles par espèce et par saison
#
# Questions retirées (méthode pas assez fiable avec nos effectifs) :
#   - pente du spectre de taille (NBSS) : souvent < 5 individus par trait,
#     résidus sur-dispersés ;
#   - condition et taille moyennes « toutes espèces » : elles mélangent
#     les espèces et changent dès que la composition change. Condition et
#     taille sont analysées ESPÈCE PAR ESPÈCE dans le script 08.
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
unlink(list.files(dossier_out, full.names = TRUE))   # anciennes sorties
for (d in c(dossier_out, dossier_fig, dossier_diag)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
saisons <- c("printemps", "ete", "automne")
couleurs_saison <- c(printemps = "#66a61e", ete = "#e6ab02", automne = "#d95f02")

prep <- readRDS(file.path(dossiers$processed, "preparation.rds"))
traits_ref <- read.csv(fichier_traits_fonctionnels, encoding = "UTF-8",
                       stringsAsFactors = FALSE)
noms <- setNames(traits_ref$nom_commun, traits_ref$taxon)

# =============================================================================
# 1. Individus : stade
# =============================================================================
# Les individus non pesés comptent aussi (seule la taille sert ici)
ind <- prep$individus |>
  filter(groupe == "poisson", !is.na(taille_mm), !aberrant_taille_poids) |>
  left_join(traits_ref |> select(taxon, taille_maturite_cm), by = "taxon") |>
  inner_join(prep$traits |> filter(analysable) |> select(trait_id, distance_km),
             by = "trait_id") |>
  mutate(saison = factor(saison, levels = saisons),
         stade = ifelse(is.na(taille_maturite_cm), NA,
                        ifelse(taille_mm < 10 * taille_maturite_cm,
                               "juvénile", "adulte")))

sans_lm <- unique(ind$taxon[is.na(ind$stade)])
if (length(sans_lm)) {
  message("Taxons sans taille de maturité (stade non attribué) : ",
          paste(sans_lm, collapse = ", "))
}

# =============================================================================
# 2. Indicateurs par trait de chalut
# =============================================================================
par_trait <- ind |>
  group_by(trait_id) |>
  summarise(
    n_mesures      = n(),
    part_juveniles = moy_pond(stade == "juvénile", facteur_elevation),
    dens_juv_km    = sum(facteur_elevation[stade %in% "juvénile"]) / first(distance_km),
    S_juv          = n_distinct(taxon[stade %in% "juvénile"]),
    .groups = "drop") |>
  left_join(prep$traits |> select(trait_id, saison, station, date, profondeur_m,
                                  lat, lon, douteux, fune_m),
            by = "trait_id") |>
  mutate(saison = factor(saison, levels = saisons), station = factor(station)) |>
  relocate(trait_id, saison, station, date, profondeur_m)

# =============================================================================
# 3. Modèles mixtes (même démarche que 05b)
# =============================================================================
reponses <- tribble(
  ~indice,          ~famille,   ~transformation, ~label,
  "dens_juv_km",    "gaussian", "log10p1", "Densité de juvéniles (log10 ind/km + 1)",
  "S_juv",          "genpois",  "aucune",  "Espèces présentes en juvéniles",
  "part_juveniles", "beta",     "01",      "Part de juvéniles (abondance)"
)

tabs <- modeliser_indices(par_trait, reponses, dossier_diag)
figures_modeles(par_trait, reponses, tabs, dossier_fig,
                sujet = "les juvéniles")

# Sensibilité à la profondeur : au printemps, les traits analysables sont
# presque tous côtiers. On compare les saisons sur les seuls traits peu
# profonds (Kruskal-Wallis : pas d'hypothèse de loi, adapté aux petits n).
cotiers <- par_trait |> filter(profondeur_m <= profondeur_max_cotiers)
sens_prof <- bind_rows(lapply(reponses$indice, function(v) {
  k <- kruskal.test(cotiers[[v]] ~ droplevels(cotiers$saison))
  med <- tapply(cotiers[[v]], cotiers$saison, median)
  n <- table(cotiers$saison)
  data.frame(indice = v, profondeur_max_m = profondeur_max_cotiers,
             n_printemps = n[["printemps"]], n_ete = n[["ete"]],
             n_automne = n[["automne"]],
             mediane_printemps = med[["printemps"]], mediane_ete = med[["ete"]],
             mediane_automne = med[["automne"]],
             chi2_kruskal = unname(k$statistic), p_kruskal = k$p.value)
}))

# =============================================================================
# 4. Densité de juvéniles par espèce et saison (descriptif)
# =============================================================================
traits_a <- prep$traits |> filter(analysable) |> select(trait_id, saison, distance_km)
dens_esp <- ind |>
  filter(stade %in% "juvénile") |>
  group_by(taxon, trait_id) |>
  summarise(n_eleve = sum(facteur_elevation), .groups = "drop") |>
  tidyr::complete(taxon, trait_id = traits_a$trait_id, fill = list(n_eleve = 0)) |>
  left_join(traits_a, by = "trait_id") |>
  mutate(dens = n_eleve / distance_km,
         saison = factor(saison, levels = saisons)) |>
  group_by(taxon, saison) |>
  summarise(n_traits = n(), traits_avec_juveniles = sum(n_eleve > 0),
            densite_moyenne_ind_km = mean(dens), densite_mediane_ind_km = median(dens),
            .groups = "drop") |>
  mutate(nom_commun = noms[taxon]) |>
  arrange(taxon, saison)

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
ecrire(sens_prof, "sensibilite_traits_cotiers.csv")
ecrire(dens_esp, "densite_juveniles_par_espece.csv")

# =============================================================================
# 6. Figure : densité de juvéniles des principales espèces
# =============================================================================
principales <- ind |> count(taxon) |> filter(n >= n_min_individus_espece) |> pull(taxon)
p_esp <- dens_esp |>
  filter(taxon %in% principales) |>
  mutate(espece = paste0(noms[taxon], "\n(", taxon, ")")) |>
  ggplot(aes(saison, densite_moyenne_ind_km, fill = saison)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(traits_avec_juveniles, "/", n_traits)),
            vjust = -0.3, size = 3) +
  facet_wrap(~ espece, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = couleurs_saison, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Densité moyenne de juvéniles par espèce et par campagne",
       subtitle = "Étiquette : traits où l'espèce est présente en juvéniles / traits analysables",
       x = NULL, y = "Juvéniles par km chaluté") +
  theme_minimal(base_size = 10)
print(p_esp)
ggsave(file.path(dossier_fig, "densite_juveniles_especes.png"), p_esp,
       width = 11, height = 4.5, dpi = 150, bg = "white")

# =============================================================================
# 7. Résumé
# =============================================================================
message("\n================ RÉSUMÉ ================")
message(sprintf("%d poissons mesurés dans %d traits ; %.0f %% de juvéniles (abondance élevée)",
                nrow(ind), nrow(par_trait),
                100 * moy_pond(ind$stade == "juvénile", ind$facteur_elevation)))
print(tabs$resume |> select(indice, famille, n, p_saison_holm, p_profondeur_holm,
                            R2, probleme_diagnostic) |>
        mutate(across(where(is.numeric), ~ signif(.x, 3))))
message("\nSans les traits douteux :")
print(tabs$sensibilite |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
message("\nTraits côtiers seulement (<= ", profondeur_max_cotiers, " m) :")
print(sens_prof |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
