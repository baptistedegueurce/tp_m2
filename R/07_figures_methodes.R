# =============================================================================
# 07_figures_methodes.R — Figures de la partie Matériel et méthodes
#
# Figure 1 : zone d'étude et stations (A : situation dans le golfe de
#            Gascogne ; B : stations, profondeur moyenne et nombre de traits
#            analysables sur les 3 campagnes)
# Figure 2 : démarche d'analyse, du chalutage aux modèles
#
# Entrée  : data/processed/environnement_traits.csv (script 03)
# Sorties : outputs/figures/article/
# =============================================================================

pkgs <- c("here", "dplyr", "ggplot2", "maps", "mapdata", "patchwork", "ggrepel")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants)) install.packages(manquants)

library(dplyr)
library(ggplot2)
library(patchwork)
library(mapdata)   # base "worldHires" (trait de côte détaillé)

source(here::here("R", "00_config.R"), encoding = "UTF-8")
dossier_fig <- here::here("outputs", "figures", "article")
dir.create(dossier_fig, recursive = TRUE, showWarnings = FALSE)

env <- read.csv(file.path(dossiers$processed, "environnement_traits.csv"),
                encoding = "UTF-8")

# =============================================================================
# Figure 1 : carte
# =============================================================================
stations <- env |>
  group_by(station) |>
  summarise(lat = mean(lat), lon = mean(lon),
            profondeur = mean(profondeur_m, na.rm = TRUE),
            traits_analysables = sum(analysable), .groups = "drop")

trait_cote <- function(xlim, ylim) {
  m <- maps::map("worldHires", xlim = xlim + c(-1, 1), ylim = ylim + c(-1, 1),
                 plot = FALSE, fill = TRUE)
  data.frame(lon = m$x, lat = m$y,
             groupe = cumsum(is.na(m$x)))[!is.na(m$x), ]
}

# A : situation
lim_a <- list(x = c(-6, 0), y = c(43.3, 48.8))
cote_a <- trait_cote(lim_a$x, lim_a$y)
p_a <- ggplot() +
  geom_polygon(data = cote_a, aes(lon, lat, group = groupe),
               fill = "grey85", colour = "grey55", linewidth = 0.2) +
  annotate("rect", xmin = -2.62, xmax = -2.12, ymin = 47.05, ymax = 47.32,
           fill = NA, colour = "#b2182b", linewidth = 0.6) +
  annotate("text", x = -4.6, y = 45.3, label = "Golfe de\nGascogne",
           fontface = "italic", size = 3, colour = "grey30") +
  annotate("text", x = -1.2, y = 47.6, label = "Loire", size = 3, colour = "grey30") +
  coord_quickmap(xlim = lim_a$x, ylim = lim_a$y, expand = FALSE) +
  labs(x = NULL, y = NULL, tag = "A") +
  theme_bw(base_size = 9) +
  theme(panel.background = element_rect(fill = "#eaf2f8"),
        panel.grid = element_blank())

# B : stations
lim_b <- list(x = c(-2.62, -2.12), y = c(47.05, 47.32))
cote_b <- trait_cote(lim_b$x, lim_b$y)
p_b <- ggplot() +
  geom_polygon(data = cote_b, aes(lon, lat, group = groupe),
               fill = "grey85", colour = "grey55", linewidth = 0.2) +
  geom_point(data = stations,
             aes(lon, lat, fill = profondeur, size = factor(traits_analysables)),
             shape = 21, colour = "grey15") +
  ggrepel::geom_text_repel(data = stations, aes(lon, lat, label = station),
                           size = 3, min.segment.length = 0.2, seed = 1,
                           point.padding = 0.3) +
  geom_point(aes(x = point_embouchure[["lon"]], y = point_embouchure[["lat"]]),
             shape = 17, size = 2.5, colour = "#b2182b") +
  annotate("text", x = point_embouchure[["lon"]] - 0.012,
           y = point_embouchure[["lat"]] + 0.008, label = "Saint-Nazaire",
           hjust = 1, size = 3, colour = "#b2182b") +
  scale_fill_gradient(low = "#c6dbef", high = "#08306b", name = "Profondeur\nmoyenne (m)") +
  scale_size_manual(values = c(`0` = 1.5, `1` = 2.5, `2` = 3.5, `3` = 4.5),
                    name = "Traits\nanalysables\n(sur 3)") +
  coord_quickmap(xlim = lim_b$x, ylim = lim_b$y, expand = FALSE) +
  labs(x = "Longitude (°)", y = "Latitude (°)", tag = "B") +
  theme_bw(base_size = 9) +
  theme(panel.background = element_rect(fill = "#eaf2f8"),
        panel.grid = element_blank())

fig1 <- p_a + p_b + plot_layout(widths = c(1, 1.9))
ggsave(file.path(dossier_fig, "figure1_zone_etude.png"), fig1,
       width = 17, height = 8.5, units = "cm", dpi = 300, bg = "white")

# =============================================================================
# Figure 2 : démarche d'analyse
# =============================================================================
boites <- tribble(
  ~id, ~x,  ~y,  ~texte,                                                              ~type,
  1,   2,   9,   "Chalutage de fond\n14 stations × 3 campagnes (avril, juillet, septembre)\n27 traits analysables", "terrain",
  2,   2,   7.2, "Tri, comptage, mesure et pesée\n(≤ 20 individus par taxon et par trait, Gwaleen)\n1 495 poissons mesurés", "terrain",
  3,   2,   5.4, "Élévation à l'effectif du trait\npuis à la distance chalutée (CPUE, ind/km)", "preparation",
  4,   0.2, 3.4, "Par trait de chalut\nrichesse, Shannon, Piélou,\ndensité totale et de juvéniles", "donnees",
  5,   3.8, 3.4, "Par individu (5 taxons ≥ 50 ind.)\nlongueur, masse, stade\n(L < L50 = juvénile)", "donnees",
  6,   0.2, 1.3, "indice ~ saison + profondeur\n+ (1 | station)\nlmer / glmer (Poisson)", "modele",
  7,   3.8, 1.3, "log W ~ log L + saison + profondeur\nlog L ~ saison + profondeur\n+ (1 | trait), saisons ≥ 3 traits", "modele",
  8,   2,   -0.5, "Test du rapport de vraisemblance · Tukey · résidus DHARMa", "test"
)
fleches <- tribble(
  ~de, ~a,
  1, 2, 2, 3, 3, 4, 3, 5, 4, 6, 5, 7, 6, 8, 7, 8
) |>
  left_join(boites |> select(de = id, x0 = x, y0 = y), by = "de") |>
  left_join(boites |> select(a = id, x1 = x, y1 = y), by = "a") |>
  mutate(y0 = y0 - 0.62, y1 = y1 + 0.62)
fleches$y1[fleches$a == 8] <- -0.5 + 0.35
fleches$y0[fleches$de %in% 6:7] <- 1.3 - 0.62

couleurs <- c(terrain = "#e5f4e4", preparation = "#f9eeda", donnees = "#e7f1fb",
              modele = "#efedff", test = "#f0efec")
fig2 <- ggplot() +
  geom_segment(data = fleches, aes(x0, y0, xend = x1, yend = y1),
               arrow = arrow(length = unit(2, "mm"), type = "closed"),
               colour = "grey40", linewidth = 0.4) +
  geom_label(data = boites, aes(x, y, label = texte, fill = type),
             size = 2.7, lineheight = 0.95, label.size = 0.3,
             label.padding = unit(2.2, "mm"), colour = "grey10") +
  annotate("text", x = -1.35, y = c(8.1, 5.4, 3.4, 1.3, -0.5),
           label = c("Terrain", "Préparation", "Données", "Modèles", "Tests"),
           angle = 90, size = 3, fontface = "bold", colour = "grey35") +
  annotate("text", x = c(0.2, 3.8), y = 4.55,
           label = c("Communautés (script 04)", "Espèces (script 05)"),
           size = 3, fontface = "bold", colour = "grey25") +
  scale_fill_manual(values = couleurs, guide = "none") +
  coord_cartesian(xlim = c(-1.6, 5.6), ylim = c(-1.1, 9.7)) +
  theme_void()
ggsave(file.path(dossier_fig, "figure2_demarche.png"), fig2,
       width = 17, height = 13, units = "cm", dpi = 300, bg = "white")

message("Figures : ", dossier_fig)
