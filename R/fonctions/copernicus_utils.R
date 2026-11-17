# =============================================================================
# copernicus_utils.R — Fonctions pour télécharger, lire et tracer les données
# Copernicus Marine. Sourcé par 01_extraction_copernicus.R
# =============================================================================

# -----------------------------------------------------------------------------
# Toolbox Copernicus Marine (CLI Python `copernicusmarine`)
# -----------------------------------------------------------------------------

#' Trouve l'exécutable `copernicusmarine`, ou s'arrête avec un message clair.
trouver_copernicusmarine <- function(chemin = "") {
  if (nzchar(chemin)) {
    if (!file.exists(chemin)) stop("Exécutable introuvable : ", chemin)
    return(chemin)
  }
  exe <- Sys.which("copernicusmarine")
  if (!nzchar(exe)) {
    stop(
      "La toolbox `copernicusmarine` est introuvable.\n",
      "  1) Installer : pip install copernicusmarine\n",
      "  2) Se connecter une fois : copernicusmarine login\n",
      "     (ou renseigner COPERNICUSMARINE_SERVICE_USERNAME / _PASSWORD dans .Renviron)\n",
      "  3) Si R ne la trouve toujours pas, renseigner `chemin_copernicusmarine`\n",
      "     dans R/00_config.R",
      call. = FALSE
    )
  }
  unname(exe)
}

#' Regroupe les couples (année, mois) en périodes de mois consécutifs,
#' pour limiter le nombre de requêtes sans télécharger de mois inutiles.
#' Retourne un data.frame (debut, fin) au format "YYYY-MM-DD".
periodes_contigues <- function(annees, mois) {
  stopifnot(all(mois %in% 1:12))
  grille <- expand.grid(mois = sort(unique(mois)), annee = sort(unique(annees)))
  idx    <- grille$annee * 12 + grille$mois           # index mensuel continu
  groupe <- cumsum(c(1, diff(idx) != 1))
  dates  <- as.Date(sprintf("%d-%02d-01", grille$annee, grille$mois))
  data.frame(
    debut = format(dates[!duplicated(groupe)]),
    fin   = format(dates[!duplicated(groupe, fromLast = TRUE)]),
    stringsAsFactors = FALSE
  )
}

#' Télécharge un sous-ensemble d'un dataset via `copernicusmarine subset`.
#' Retourne le chemin du fichier NetCDF.
telecharger_subset <- function(exe, dataset_id, variables, zone,
                               debut = NULL, fin = NULL,
                               prof_min = NULL, prof_max = NULL,
                               dossier, ecraser = FALSE) {
  periode  <- if (is.null(debut)) "statique" else paste0(debut, "_", fin)
  prof_tag <- if (is.null(prof_min)) "" else sprintf("_z%g-%g", prof_min, prof_max)
  nom      <- sprintf("%s_%s%s.nc", dataset_id, periode, prof_tag)
  fichier  <- file.path(dossier, nom)

  if (file.exists(fichier) && !ecraser) {
    message("  ✓ déjà présent : ", nom)
    return(fichier)
  }

  args <- c(
    "subset",
    "--dataset-id", dataset_id,
    as.vector(rbind("--variable", variables)),
    "--minimum-longitude", zone$lon_min, "--maximum-longitude", zone$lon_max,
    "--minimum-latitude",  zone$lat_min, "--maximum-latitude",  zone$lat_max,
    "--coordinates-selection-method", "outside",   # englobe toute la zone
    "--output-directory", shQuote(dossier),
    "--output-filename", nom,
    "--overwrite",
    "--disable-progress-bar",
    "--log-level", "WARN"
  )
  if (!is.null(debut)) {
    args <- c(args, "--start-datetime", debut, "--end-datetime", fin)
  }
  if (!is.null(prof_min)) {
    args <- c(args, "--minimum-depth", prof_min, "--maximum-depth", prof_max)
  }

  message("  ↓ ", dataset_id, " [", paste(variables, collapse = ", "), "] ", periode)
  statut <- system2(exe, args)
  if (!identical(statut, 0L) || !file.exists(fichier)) {
    stop("Échec du téléchargement de ", dataset_id, " (code ", statut, ").\n",
         "Vérifier identifiants, zone, dates (couverture du dataset) et connexion.",
         call. = FALSE)
  }
  fichier
}

# -----------------------------------------------------------------------------
# Lecture NetCDF -> data.frame long
# -----------------------------------------------------------------------------

#' Convertit l'axe temps d'un NetCDF ("<unité> since <origine>") en Date.
convertir_temps <- function(valeurs, unites) {
  m <- regmatches(unites, regexec("^\\s*(\\w+)\\s+since\\s+(.+)$", unites))[[1]]
  if (length(m) != 3) stop("Unité de temps non reconnue : ", unites)
  facteur <- switch(tolower(m[2]),
                    seconds = 1, second = 1,
                    minutes = 60, minute = 60,
                    hours = 3600, hour = 3600,
                    days = 86400, day = 86400,
                    stop("Unité de temps non gérée : ", m[2]))
  origine <- as.POSIXct(substr(trimws(m[3]), 1, 19), tz = "UTC",
                        tryFormats = c("%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S",
                                       "%Y-%m-%d"))
  as.Date(origine + valeurs * facteur)
}

#' Lit une variable d'un NetCDF Copernicus et renvoie un data.frame long :
#' lon, lat, profondeur, date, valeur.
lire_variable_nc <- function(fichier, var) {
  nc <- ncdf4::nc_open(fichier)
  on.exit(ncdf4::nc_close(nc))

  if (!var %in% names(nc$var)) {
    stop("Variable '", var, "' absente de ", basename(fichier))
  }
  v      <- nc$var[[var]]
  valeur <- ncdf4::ncvar_get(nc, var, collapse_degen = FALSE)

  # Coordonnées dans l'ordre des dimensions de la variable
  coords <- lapply(v$dim, function(d) {
    if (d$name == "time") convertir_temps(d$vals, d$units) else d$vals
  })
  names(coords) <- vapply(v$dim, `[[`, "", "name")

  # expand.grid fait varier la 1re dimension le plus vite = ordre des tableaux R
  df <- expand.grid(coords, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  df$valeur <- as.vector(valeur)

  renommer <- c(longitude = "lon", latitude = "lat", depth = "profondeur",
                time = "date")
  names(df) <- ifelse(names(df) %in% names(renommer),
                      renommer[names(df)], names(df))
  if (!"profondeur" %in% names(df)) df$profondeur <- NA_real_
  if (!"date" %in% names(df))       df$date       <- as.Date(NA)

  df[!is.na(df$valeur), c("lon", "lat", "profondeur", "date", "valeur")]
}

#' Unité d'une variable dans un NetCDF.
unite_nc <- function(fichier, var) {
  nc <- ncdf4::nc_open(fichier)
  on.exit(ncdf4::nc_close(nc))
  nc$var[[var]]$units
}

#' Agrège la dimension verticale selon la méthode choisie dans la config.
agreger_profondeur <- function(df, methode = c("surface", "moyenne", "aucune")) {
  methode <- match.arg(methode)
  if (methode == "aucune" || all(is.na(df$profondeur))) return(df)

  if (methode == "surface") {
    # niveau le moins profond disponible pour chaque pixel / date
    df <- df |>
      dplyr::group_by(lon, lat, date) |>
      dplyr::slice_min(profondeur, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
  } else {
    df <- df |>
      dplyr::group_by(lon, lat, date) |>
      dplyr::summarise(profondeur = mean(profondeur),
                       valeur = mean(valeur), .groups = "drop")
  }
  as.data.frame(df)
}

# -----------------------------------------------------------------------------
# Graphiques
# -----------------------------------------------------------------------------

#' Carte d'une variable : un panneau par mois (ou une carte si statique).
carte_variable <- function(df, label, palette = "viridis", zone = NULL) {
  statique <- all(is.na(df$date))
  if (!statique) {
    df$mois_label <- factor(format(df$date, "%Y-%m"))
  }
  # si plusieurs niveaux par pixel (agregation = "aucune"), on montre le plus
  # superficiel
  if (anyDuplicated(df[, c("lon", "lat", "date")])) {
    df <- df[df$profondeur == min(df$profondeur, na.rm = TRUE), ]
    label <- sprintf("%s — niveau %.1f m", label, min(df$profondeur))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(lon, lat, fill = valeur)) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_viridis_c(option = palette, name = NULL) +
    ggplot2::coord_quickmap(expand = FALSE) +
    ggplot2::labs(title = label, x = "Longitude (°)", y = "Latitude (°)") +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "grey85",
                                                            colour = NA),
                   panel.grid = ggplot2::element_blank())

  if (!statique) p <- p + ggplot2::facet_wrap(~ mois_label)
  p
}

#' Résumé : moyenne spatiale (± écart-type) par mois, toutes variables.
plot_resume_mensuel <- function(donnees) {
  resume <- donnees |>
    dplyr::filter(!is.na(date)) |>
    dplyr::group_by(label, date) |>
    dplyr::summarise(moy = mean(valeur), et = stats::sd(valeur),
                     .groups = "drop")

  ggplot2::ggplot(resume, ggplot2::aes(date, moy)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = moy - et, ymax = moy + et),
                         fill = "steelblue", alpha = 0.25) +
    ggplot2::geom_line(colour = "steelblue4") +
    ggplot2::geom_point(colour = "steelblue4", size = 1.5) +
    ggplot2::facet_wrap(~ label, scales = "free_y") +
    ggplot2::labs(title = "Moyenne spatiale mensuelle sur la zone (± écart-type)",
                  x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 10)
}
