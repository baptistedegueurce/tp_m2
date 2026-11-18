# =============================================================================
# 04_communautes.R
# Structure des communautés de poissons par trait de chalut : abondance,
# diversité (richesse, Shannon, Piélou) et arrivée des juvéniles.
#
# Entrées : data/processed/preparation.rds          (script 02)
#           data/processed/environnement_traits.csv (script 03)
#           data/ref/traits_fonctionnels.csv        (taille de 1re maturité)
#
# 1. Indices par trait analysable (poissons) :
#      S        richesse : nombre de taxons
#      H        indice de Shannon, H = -somme(p_i ln p_i), p_i = part du
#               taxon i dans l'abondance (CPUE élevée)
#      J        équitabilité de Piélou, J = H / ln(S) (si S >= 2)
#      N        densité totale (ind/km)
#      N_juv    densité de juvéniles (ind/km) ; juvénile = taille < taille
#               de 1re maturité (FishBase)
# 2. Même modèle pour chaque indice :
#      indice ~ saison + profondeur + (1 | station)
#    (1 | station) : les 14 stations sont visitées jusqu'à 3 fois.
#    Loi normale (lmer) ; loi de Poisson pour la richesse (glmer).
#    Densités en log10(x + 1).
#    Test de chaque effet (rapport de vraisemblance), Tukey entre saisons,
#    résidus vérifiés avec DHARMa.
# 3. Description : composition par saison, juvéniles par espèce
#
# Sorties : data/processed/communautes/, outputs/figures/communautes/
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
dossier_out  <- file.path(dossiers$processed, "communautes")
dossier_fig  <- here::here("outputs", "figures", "communautes")
dossier_diag <- file.path(dossier_fig, "diagnostics")
for (d in c(dossier_out, dossier_fig, dossier_diag)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
saisons <- c("printemps", "ete", "automne")
couleurs_saison <- c(printemps = "#66a61e", ete = "#e6ab02", automne = "#d95f02")

prep <- readRDS(file.path(dossiers$processed, "preparation.rds"))
fichier_env <- file.path(dossiers$processed, "environnement_traits.csv")
if (!file.exists(fichier_env)) stop("Lancer d'abord R/03_variables_environnementales.R")
env <- read.csv(fichier_env, encoding = "UTF-8") |> filter(analysable)
traits_ref <- read.csv(fichier_traits_fonctionnels, encoding = "UTF-8")
noms <- setNames(traits_ref$nom_commun, traits_ref$taxon)

prof_moy <- mean(env[[covariable_spatiale]])
traits <- env |>
  transmute(trait_id, station = factor(station),
            saison = factor(saison, levels = saisons),
            profondeur_m = .data[[covariable_spatiale]],
            prof10 = (profondeur_m - prof_moy) / 10)   # +1 = 10 m plus profond

# =============================================================================
# 1. Indices par trait
# =============================================================================
captures <- prep$captures |>
  filter(groupe %in% groupes_communaute, trait_id %in% traits$trait_id,
         cpue_n_km > 0)

diversite <- captures |>
  group_by(trait_id) |>
  summarise(S = n_distinct(taxon),
            N = sum(cpue_n_km),
            H = -sum((cpue_n_km / N) * log(cpue_n_km / N)),
            .groups = "drop") |>
  mutate(J = ifelse(S >= 2, H / log(S), NA))

# Juvéniles : individus mesurés, élevés à l'effectif du trait
juv <- prep$individus |>
  filter(groupe %in% groupes_communaute, !aberrant_taille_poids,
         !is.na(taille_mm), trait_id %in% traits$trait_id) |>
  left_join(traits_ref |> select(taxon, taille_maturite_cm), by = "taxon") |>
  mutate(juvenile = taille_mm < 10 * taille_maturite_cm)

par_trait <- traits |>
  left_join(diversite, by = "trait_id") |>
  left_join(juv |> group_by(trait_id) |>
              summarise(n_juv = sum(facteur_elevation[juvenile %in% TRUE]),
                        part_juveniles = sum(facteur_elevation[juvenile %in% TRUE]) /
                          sum(facteur_elevation[!is.na(juvenile)]),
                        .groups = "drop"),
            by = "trait_id") |>
  left_join(prep$traits |> select(trait_id, distance_km), by = "trait_id") |>
  mutate(S = coalesce(S, 0L), N = coalesce(N, 0), n_juv = coalesce(n_juv, 0),
         N_juv = n_juv / distance_km) |>
  select(-n_juv)

# =============================================================================
# 2. Modèles
# =============================================================================
indices <- tribble(
  ~indice, ~label,                               ~loi,      ~transfo,
  "S",     "Richesse (nombre de taxons)",        "poisson", "aucune",
  "H",     "Indice de Shannon H'",               "normale", "aucune",
  "J",     "Équitabilité de Piélou J'",          "normale", "aucune",
  "N",     "Densité totale (log10 ind/km + 1)",  "normale", "log10p1",
  "N_juv", "Densité de juvéniles (log10 ind/km + 1)", "normale", "log10p1"
)

res <- list(tests = list(), moyennes = list(), contrastes = list(),
            diagnostics = list(), points = list())
for (i in seq_len(nrow(indices))) {
  r <- indices[i, ]
  d <- par_trait
  d$y <- if (r$transfo == "log10p1") log10(d[[r$indice]] + 1) else d[[r$indice]]
  d <- d[!is.na(d$y), ]
  m <- ajuster(y ~ saison + prof10 + (1 | station), d,
               if (r$loi == "poisson") poisson() else NULL)
  t <- tester(m)
  s <- saisons_tukey(m)
  co <- summary(m)$coefficients["prof10", ]
  res$tests[[r$indice]] <- data.frame(
    indice = r$indice, label = r$label, loi = r$loi, n_traits = nrow(d),
    t, effet_profondeur_10m = co[[1]], ES_profondeur = co[[2]],
    avertissements = paste(attr(m, "avertissements"), collapse = " | "))
  res$moyennes[[r$indice]] <- cbind(indice = r$indice, label = r$label, s$moyennes)
  res$contrastes[[r$indice]] <- cbind(indice = r$indice, s$contrastes)
  res$diagnostics[[r$indice]] <- cbind(
    indice = r$indice,
    diagnostic(m, file.path(dossier_diag, paste0(r$indice, ".png"))))
  res$points[[r$indice]] <- data.frame(indice = r$indice, label = r$label,
                                       saison = d$saison, y = d$y)
}
tabs <- lapply(res, lier)

# =============================================================================
# 3. Description
# =============================================================================
# Composition : CPUE moyenne par taxon et par saison (tous traits analysables,
# zéros compris)
composition <- captures |>
  select(trait_id, taxon, cpue_n_km) |>
  tidyr::complete(trait_id = traits$trait_id, taxon, fill = list(cpue_n_km = 0)) |>
  left_join(traits |> select(trait_id, saison), by = "trait_id") |>
  group_by(saison, taxon) |>
  summarise(cpue_moyenne = mean(cpue_n_km), occurrence = mean(cpue_n_km > 0),
            .groups = "drop") |>
  mutate(nom_commun = noms[taxon])

# Densité de juvéniles par espèce et saison
juv_esp <- juv |>
  filter(juvenile %in% TRUE) |>
  group_by(taxon, trait_id) |>
  summarise(n = sum(facteur_elevation), .groups = "drop") |>
  tidyr::complete(taxon, trait_id = traits$trait_id, fill = list(n = 0)) |>
  left_join(prep$traits |> select(trait_id, distance_km), by = "trait_id") |>
  left_join(traits |> select(trait_id, saison), by = "trait_id") |>
  group_by(taxon, saison) |>
  summarise(n_traits = n(), traits_presents = sum(n > 0),
            densite_moyenne = mean(n / distance_km), .groups = "drop") |>
  mutate(nom_commun = noms[taxon])

# =============================================================================
# Export
# =============================================================================
tabs$indices_par_trait <- par_trait
tabs$composition <- composition
tabs$juveniles_par_espece <- juv_esp
points <- tabs$points; tabs$points <- NULL
for (nom in names(tabs)) {
  write.csv(tabs[[nom]], file.path(dossier_out, paste0(nom, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}

# =============================================================================
# 4. Figures
# =============================================================================
sauver <- function(p, nom, l, h) {
  print(p)
  ggsave(file.path(dossier_fig, nom), p, width = l, height = h, dpi = 150,
         bg = "white")
}

# 4a. Indices par saison : traits (points) + moyenne du modèle (IC 95 %)
p_sig <- tabs$tests |> filter(effet == "saison") |>
  mutate(texte = paste("saison :", ifelse(p < 0.001, "p < 0,001",
         paste0("p = ", formatC(p, digits = 2, format = "fg", decimal.mark = ",")))))
ordre <- function(x) factor(x, levels = indices$label)
points$label <- ordre(points$label)
p_sig$label <- ordre(p_sig$label)
# moyennes du modèle sur la même échelle que les points (log10(x + 1) pour les densités)
moy <- tabs$moyennes |>
  mutate(saison = factor(saison, levels = saisons), label = ordre(label))
p_ind <- ggplot(points,
                aes(saison, y, colour = saison)) +
  geom_jitter(width = 0.12, alpha = 0.45, size = 1.6) +
  geom_pointrange(data = moy, aes(y = moyenne, ymin = IC95_bas, ymax = IC95_haut),
                  position = position_nudge(x = 0.28), size = 0.5) +
  geom_text(data = p_sig, aes(x = 0.5, y = Inf, label = texte), inherit.aes = FALSE,
            hjust = 0, vjust = 1.3, size = 3, colour = "grey25") +
  facet_wrap(~ label, scales = "free_y", nrow = 1,
             labeller = label_wrap_gen(22)) +
  scale_colour_manual(values = couleurs_saison, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(title = "Indices de communauté par saison",
       subtitle = "Points : traits de chalut ; point + barre : moyenne du modèle à profondeur moyenne (IC 95 %)",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
sauver(p_ind, "indices_par_saison.png", 13, 4.5)

# 4b. Composition : CPUE moyenne des principaux taxons par saison
principaux <- composition |> group_by(taxon) |>
  summarise(c = sum(cpue_moyenne)) |> slice_max(c, n = 6) |> pull(taxon)
p_comp <- composition |>
  mutate(groupe = ifelse(taxon %in% principaux,
                         paste0(nom_commun, " (", taxon, ")"), "autres")) |>
  group_by(saison, groupe) |>
  summarise(cpue = sum(cpue_moyenne), .groups = "drop") |>
  ggplot(aes(factor(saison, levels = saisons), cpue, fill = groupe)) +
  geom_col(width = 0.65) +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Composition des captures de poissons par campagne",
       subtitle = "Densité moyenne par trait (ind/km, effectifs élevés)",
       x = NULL, y = "ind/km", fill = NULL) +
  theme_minimal(base_size = 10)
sauver(p_comp, "composition_par_saison.png", 8, 4.5)

# 4c. Densité de juvéniles des principales espèces
especes <- juv |> count(taxon) |> filter(n >= n_min_individus_espece) |> pull(taxon)
p_juv <- juv_esp |>
  filter(taxon %in% especes) |>
  mutate(espece = paste0(noms[taxon], "\n(", taxon, ")"),
         saison = factor(saison, levels = saisons)) |>
  ggplot(aes(saison, densite_moyenne, fill = saison)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0(traits_presents, "/", n_traits)), vjust = -0.3, size = 3) +
  facet_wrap(~ espece, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = couleurs_saison, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Densité moyenne de juvéniles par espèce et par campagne",
       subtitle = "Étiquette : traits où l'espèce est présente en juvéniles / traits analysables",
       x = NULL, y = "Juvéniles par km chaluté") +
  theme_minimal(base_size = 10)
sauver(p_juv, "juveniles_par_espece.png", 11, 4.5)

# =============================================================================
# 5. Résumé
# =============================================================================
message("\n================ RÉSUMÉ ================")
print(tabs$tests |> select(indice, effet, chi2, ddl, p, avertissements) |>
        mutate(across(where(is.numeric), ~ signif(.x, 3))))
print(tabs$moyennes |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
print(tabs$contrastes |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
print(tabs$diagnostics |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
