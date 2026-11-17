# =============================================================================
# 01_extraction_copernicus.R
# Extraction mensuelle des données environnementales Copernicus Marine
# (température, salinité, bathymétrie, chlorophylle, nutriments, O2, ...)
#
# Paramètres (zone, années, mois, profondeur, variables) : R/00_config.R
#
# Sorties :
#   data/raw/copernicus/*.nc                    fichiers bruts téléchargés
#   data/processed/copernicus_mensuel.rds       tableau long (toutes variables)
#   data/processed/copernicus_mensuel.csv       idem en CSV
#   outputs/figures/copernicus/*.png            une carte par variable + résumé
# =============================================================================

# --- Packages -----------------------------------------------------------------
pkgs <- c("here", "ncdf4", "dplyr", "ggplot2")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants)) install.packages(manquants)

source(here::here("R", "00_config.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "copernicus_utils.R"), encoding = "UTF-8")

invisible(lapply(dossiers, dir.create, recursive = TRUE, showWarnings = FALSE))

# --- 1. Préparation des requêtes ----------------------------------------------
suffixe <- if (source_copernicus == "reanalyse") "reanalyse" else "prevision"
stopifnot(source_copernicus %in% c("reanalyse", "analyse_prevision"))

inconnues <- setdiff(variables_a_extraire, catalogue_variables$nom)
if (length(inconnues)) {
  stop("Variables inconnues dans la config : ", paste(inconnues, collapse = ", "))
}

selection <- catalogue_variables[catalogue_variables$nom %in% variables_a_extraire, ]
selection$var     <- selection[[paste0("var_", suffixe)]]
selection$dataset <- selection[[paste0("dataset_", suffixe)]]

periodes <- periodes_contigues(annees, mois)
exe      <- trouver_copernicusmarine(chemin_copernicusmarine)

message("Zone : lon [", zone$lon_min, " ; ", zone$lon_max, "], lat [",
        zone$lat_min, " ; ", zone$lat_max, "]")
message("Périodes : ", paste(periodes$debut, "→", periodes$fin, collapse = " | "))
message("Source : ", source_copernicus, "\n")

# --- 2. Téléchargement (une requête par dataset et par période) ---------------
# On regroupe les variables d'un même dataset dans une seule requête.
requetes <- split(selection, selection$dataset)
fichiers <- list()   # dataset -> vecteur de fichiers .nc

for (ds in names(requetes)) {
  vars     <- requetes[[ds]]$var
  statique <- all(requetes[[ds]]$statique)

  if (statique) {
    fichiers[[ds]] <- telecharger_subset(
      exe, ds, vars, zone, dossier = dossiers$raw, ecraser = ecraser_fichiers
    )
  } else {
    fichiers[[ds]] <- vapply(seq_len(nrow(periodes)), function(i) {
      telecharger_subset(
        exe, ds, vars, zone,
        debut = periodes$debut[i], fin = periodes$fin[i],
        prof_min = profondeur$min, prof_max = profondeur$max,
        dossier = dossiers$raw, ecraser = ecraser_fichiers
      )
    }, character(1))
  }
}

# --- 3. Lecture et mise en forme ----------------------------------------------
message("\nLecture des NetCDF…")

donnees <- lapply(seq_len(nrow(selection)), function(i) {
  s <- selection[i, ]
  df <- do.call(rbind, lapply(fichiers[[s$dataset]], lire_variable_nc, var = s$var))
  unite <- unite_nc(fichiers[[s$dataset]][1], s$var)

  if (!s$statique) {
    # garde uniquement les mois demandés (sécurité) puis agrège la profondeur
    df <- df[as.integer(format(df$date, "%Y")) %in% annees &
             as.integer(format(df$date, "%m")) %in% mois, ]
    df <- agreger_profondeur(df, profondeur$agregation)
  }

  df$variable <- s$nom
  df$var_nc   <- s$var
  df$unite    <- unite
  df$label    <- s$label
  message(sprintf("  %-20s %8d valeurs", s$nom, nrow(df)))
  df
})
donnees <- do.call(rbind, donnees)

donnees$annee <- as.integer(format(donnees$date, "%Y"))
donnees$mois  <- as.integer(format(donnees$date, "%m"))
donnees <- donnees[, c("variable", "var_nc", "label", "unite", "lon", "lat",
                       "profondeur", "date", "annee", "mois", "valeur")]

saveRDS(donnees, file.path(dossiers$processed, "copernicus_mensuel.rds"))
write.csv(donnees, file.path(dossiers$processed, "copernicus_mensuel.csv"),
          row.names = FALSE)
message("\n✓ Données enregistrées dans ", dossiers$processed)

# --- 4. Graphiques de contrôle ------------------------------------------------
message("\nGraphiques…")

for (i in seq_len(nrow(selection))) {
  s  <- selection[i, ]
  df <- donnees[donnees$variable == s$nom, ]
  if (!nrow(df)) {
    warning("Aucune donnée pour ", s$nom, " (zone à terre ou hors couverture ?)")
    next
  }
  p <- carte_variable(df, s$label, s$palette)
  print(p)

  n_mois <- max(1, length(unique(df$date)))
  ggplot2::ggsave(file.path(dossiers$figures, paste0(s$nom, ".png")), p,
                  width = min(12, 4 + 2.5 * min(n_mois, 4)),
                  height = 3.5 + 2.5 * (ceiling(n_mois / 4) - 1),
                  dpi = 150, bg = "white")
}

if (any(!is.na(donnees$date))) {
  p_resume <- plot_resume_mensuel(donnees)
  print(p_resume)
  ggplot2::ggsave(file.path(dossiers$figures, "resume_mensuel.png"), p_resume,
                  width = 11, height = 7, dpi = 150, bg = "white")
}

message("✓ Figures enregistrées dans ", dossiers$figures)
