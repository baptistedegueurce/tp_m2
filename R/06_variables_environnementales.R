# =============================================================================
# 06_variables_environnementales.R
# Une ligne par trait de chalut avec ses variables environnementales
#
# Entrées : data/processed/preparation.rds          (script 02)
#           data/processed/copernicus_mensuel.rds   (script 01, facultatif)
#
# Variables :
#   - mesurées sur le terrain : profondeur moyenne, longueur de fune, jour
#     de l'année
#   - calculée : distance à l'embouchure (km, depuis `point_embouchure`)
#   - Copernicus : valeur du MOIS du trait, au pixel le plus proche du
#     milieu du trait (distance au pixel conservée pour contrôle)
#
# Attention : à 1/12° ou 1/4°, les variables Copernicus varient très peu
# d'une station à l'autre sur une zone de 17 × 10 km. Elles portent surtout
# la variation ENTRE CAMPAGNES (avril, juillet, septembre) : deux variables
# Copernicus seront presque toujours corrélées entre elles et à la saison.
#
# Sorties : data/processed/environnement_traits.csv
#           data/processed/environnement_correlations.csv
#           outputs/figures/environnement/
# =============================================================================

pkgs <- c("here", "dplyr", "tidyr", "ggplot2")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants)) install.packages(manquants)

library(dplyr)
library(ggplot2)

source(here::here("R", "00_config.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "preparation_utils.R"), encoding = "UTF-8")

dossier_fig <- here::here("outputs", "figures", "environnement")
dir.create(dossier_fig, recursive = TRUE, showWarnings = FALSE)

prep <- readRDS(file.path(dossiers$processed, "preparation.rds"))

# =============================================================================
# 1. Variables de terrain et distance à l'embouchure
# =============================================================================
env <- prep$traits |>
  select(trait_id, saison, station, date, lat, lon, profondeur_m, fune_m,
         analysable, douteux) |>
  mutate(
    date = as.Date(date),
    jour_annee = as.integer(format(date, "%j")),
    distance_embouchure_km = haversine_km(lat, lon, point_embouchure[["lat"]],
                                          point_embouchure[["lon"]])
  )

# =============================================================================
# 2. Copernicus : mois du trait, pixel le plus proche
# =============================================================================
fichier_cop <- file.path(dossiers$processed, "copernicus_mensuel.rds")

if (file.exists(fichier_cop)) {
  cop <- readRDS(fichier_cop) |>
    filter(variable %in% variables_copernicus_traits) |>
    mutate(mois_cle = format(date, "%Y-%m"))

  absentes <- setdiff(variables_copernicus_traits, unique(cop$variable))
  if (length(absentes)) {
    message("Variables Copernicus absentes du fichier (relancer le script 01 ?) : ",
            paste(absentes, collapse = ", "))
  }

  for (v in intersect(variables_copernicus_traits, unique(cop$variable))) {
    cv <- cop[cop$variable == v, ]
    # si plusieurs niveaux de profondeur, on garde le plus superficiel
    cv <- cv |> group_by(lon, lat, mois_cle) |>
      slice_min(profondeur, n = 1, with_ties = FALSE) |> ungroup()
    valeurs <- numeric(nrow(env)); dist_px <- numeric(nrow(env))
    for (i in seq_len(nrow(env))) {
      px <- cv[cv$mois_cle == format(env$date[i], "%Y-%m"), ]
      if (!nrow(px)) { valeurs[i] <- NA; dist_px[i] <- NA; next }
      d <- haversine_km(env$lat[i], env$lon[i], px$lat, px$lon)
      j <- which.min(d)
      valeurs[i] <- px$valeur[j]; dist_px[i] <- d[j]
    }
    env[[v]] <- valeurs
    env[[paste0(v, "_dist_pixel_km")]] <- dist_px
    n_na <- sum(is.na(valeurs))
    message(sprintf("%-18s : %d traits renseignés, pixel à %.1f km en moyenne%s",
                    v, sum(!is.na(valeurs)), mean(dist_px, na.rm = TRUE),
                    if (n_na) sprintf(" (%d sans valeur : mois non extrait ?)", n_na) else ""))
  }
} else {
  message("Pas de données Copernicus (", fichier_cop, ").\n",
          "-> Lancer R/01_extraction_copernicus.R (zone et mois du TP déjà réglés",
          " dans 00_config.R). En attendant, seules les variables de terrain sont",
          " disponibles.")
}

# =============================================================================
# 3. Corrélations entre covariables (traits analysables)
# =============================================================================
covs <- intersect(c("profondeur_m", "distance_embouchure_km", "fune_m",
                    "jour_annee", variables_copernicus_traits), names(env))
x <- env |> filter(analysable) |> select(all_of(covs))
x <- x[, vapply(x, function(v) sum(!is.na(v)) > 3 && sd(v, na.rm = TRUE) > 0,
               logical(1)), drop = FALSE]
cor_mat <- cor(x, use = "pairwise.complete.obs", method = "spearman")

cor_long <- as.data.frame(as.table(cor_mat)) |>
  setNames(c("var1", "var2", "rho")) |>
  filter(as.integer(var1) < as.integer(var2)) |>
  arrange(desc(abs(rho)))

# Part de la variance de chaque covariable expliquée par la saison
# (R² de lm(variable ~ saison)) : proche de 1 = la variable ne fait que
# répéter la saison
r2_saison <- sapply(names(x), function(v) {
  d <- env |> filter(analysable)
  summary(lm(d[[v]] ~ d$saison))$r.squared
})

# Conditions de chaque campagne (moyenne, min, max sur les traits
# analysables) : c'est ainsi que les variables Copernicus sont utilisées
# dans le rapport (description), pas comme covariables des modèles
campagnes <- env |>
  filter(analysable) |>
  tidyr::pivot_longer(all_of(names(x)), names_to = "variable", values_to = "valeur") |>
  group_by(saison, variable) |>
  summarise(n_traits = sum(!is.na(valeur)), moyenne = mean(valeur, na.rm = TRUE),
            min = min(valeur, na.rm = TRUE), max = max(valeur, na.rm = TRUE),
            .groups = "drop")
write.csv(campagnes, file.path(dossiers$processed, "environnement_campagnes.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

message("\nCorrélations (Spearman) les plus fortes :")
print(head(cor_long |> mutate(rho = round(rho, 2)), 10))
message("\nPart de chaque variable expliquée par la saison (R²) :")
print(round(r2_saison, 2))

# =============================================================================
# 4. Export et figures
# =============================================================================
write.csv(env, file.path(dossiers$processed, "environnement_traits.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(cor_long,
          file.path(dossiers$processed, "environnement_correlations.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(data.frame(variable = names(r2_saison), r2_saison = r2_saison),
          file.path(dossiers$processed, "environnement_part_saison.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

couleurs_saison <- c(printemps = "#66a61e", ete = "#e6ab02", automne = "#d95f02")

# 4a. Matrice de corrélation
p_cor <- as.data.frame(as.table(cor_mat)) |>
  setNames(c("var1", "var2", "rho")) |>
  ggplot(aes(var1, var2, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 3) +
  scale_fill_gradient2(limits = c(-1, 1), low = "#2166ac", high = "#b2182b") +
  labs(title = "Corrélations entre covariables (Spearman, traits analysables)",
       x = NULL, y = NULL, fill = "rho") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
print(p_cor)
ggsave(file.path(dossier_fig, "correlations_covariables.png"), p_cor,
       width = 7.5, height = 6, dpi = 150, bg = "white")

# 4b. Covariables par saison et selon la position
p_cov <- env |>
  filter(analysable) |>
  select(trait_id, saison, distance_embouchure_km, all_of(names(x))) |>
  tidyr::pivot_longer(all_of(setdiff(names(x), "distance_embouchure_km")),
                      names_to = "variable", values_to = "valeur") |>
  ggplot(aes(distance_embouchure_km, valeur, colour = saison)) +
  geom_point(size = 2) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_colour_manual(values = couleurs_saison) +
  labs(title = "Covariables selon la distance à l'embouchure",
       subtitle = "Une variable qui ne bouge qu'entre saisons forme des bandes horizontales",
       x = "Distance à l'embouchure (km)", y = NULL, colour = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
print(p_cov)
ggsave(file.path(dossier_fig, "covariables_distance.png"), p_cov,
       width = 10, height = 6.5, dpi = 150, bg = "white")
