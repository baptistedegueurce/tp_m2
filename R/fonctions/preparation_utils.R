# =============================================================================
# preparation_utils.R — Fonctions utilisées par 02_preparation_donnees.R
# =============================================================================

#' Nombre saisi avec virgule ou point décimal -> numeric ("" -> NA).
lire_nombre <- function(x) {
  x <- trimws(gsub(",", ".", x))
  x[x %in% c("", "NA", "-", "\u2013", "\u2014")] <- NA   # vide, tirets – et —
  suppressWarnings(as.numeric(x))
}

#' Saison saisie ("Printemps (avril)", "Été (juillet)", "p"...) -> code p / e / a.
lire_saison <- function(x) {
  x <- tolower(trimws(x))
  out <- rep(NA_character_, length(x))
  out[grepl("print|avr", x) | x == "p"] <- "p"
  out[grepl("juil", x) | grepl("^.t.$", x) | x %in% c("e", "ete")] <- "e"
  out[grepl("aut|sept", x) | x == "a"] <- "a"
  out
}

#' Heure "7h05", "07:05", "7h5" -> "07:05".
lire_heure <- function(x) {
  m <- regmatches(x, regexec("^\\s*(\\d{1,2})\\s*[hH:]\\s*(\\d{1,2})", x))
  vapply(m, function(v) if (length(v) == 3)
    sprintf("%02d:%02d", as.integer(v[2]), as.integer(v[3])) else NA_character_,
    character(1))
}

#' Longueur de fune (m) extraite d'un commentaire "Fune : 40 m, ...".
lire_fune <- function(x) {
  m <- regmatches(x, regexec("[Ff]une\\s*:?\\s*(\\d+)", x))
  vapply(m, function(v) if (length(v) == 2) as.numeric(v[2]) else NA_real_,
         numeric(1))
}

#' Date jj/mm/aaaa ou jj/mm/aa -> Date.
lire_date <- function(x) {
  x <- trimws(x)
  d <- as.Date(x, format = "%d/%m/%Y")
  courte <- grepl("^\\d{1,2}/\\d{1,2}/\\d{2}$", x)
  d[courte] <- as.Date(x[courte], format = "%d/%m/%y")
  d
}

#' Texte de sous-échantillonnage ("1/4", "1/4 (1 sachet)", "") -> fraction.
lire_fraction <- function(x) {
  x <- trimws(x)
  out <- rep(1, length(x))
  m <- regmatches(x, regexec("(\\d+)\\s*/\\s*(\\d+)", x))
  ok <- lengths(m) == 3
  out[ok] <- vapply(m[ok], function(v) as.numeric(v[2]) / as.numeric(v[3]),
                    numeric(1))
  out
}

#' Distance orthodromique (km) entre deux points en degrés décimaux.
haversine_km <- function(lat1, lon1, lat2, lon2) {
  rad <- pi / 180
  dlat <- (lat2 - lat1) * rad
  dlon <- (lon2 - lon1) * rad
  a <- sin(dlat / 2)^2 + cos(lat1 * rad) * cos(lat2 * rad) * sin(dlon / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}

#' Découpe les lignes d'un couple saison × station × espèce en "blocs".
#'
#' Règle de saisie : l'effectif, s'il est noté, donne le nombre TOTAL de
#' poissons de l'échantillon (mesurés + non mesurés). Il peut être noté sur la
#' 1re ligne seulement, ou recopié sur toutes les lignes. S'il est absent, le
#' nombre de lignes donne l'effectif.
#'
#' Un nouveau bloc commence quand apparaît un effectif DIFFÉRENT du précédent
#' (cas observé de plusieurs sachets saisis à la suite). Les lignes sans
#' effectif appartiennent au bloc en cours.
#'
#' @param effectif vecteur d'effectifs (NA si non noté), dans l'ordre de saisie
#' @return vecteur d'indices de bloc (1, 1, 1, 2, 2, ...)
numeroter_blocs <- function(effectif) {
  bloc <- integer(length(effectif))
  courant <- 1L
  dernier <- NA
  for (i in seq_along(effectif)) {
    e <- effectif[i]
    if (!is.na(e) && !is.na(dernier) && e != dernier) courant <- courant + 1L
    if (!is.na(e)) dernier <- e
    bloc[i] <- courant
  }
  bloc
}

#' Accumulateur d'anomalies : renvoie une fonction `signaler()` et un getter.
creer_registre_anomalies <- function() {
  registre <- list()
  list(
    signaler = function(type, detail, saison = NA, station = NA, taxon = NA,
                        action = NA) {
      registre[[length(registre) + 1]] <<- data.frame(
        type = type, saison = saison, station = station, taxon = taxon,
        detail = detail, action = action, stringsAsFactors = FALSE
      )
      invisible(NULL)
    },
    tableau = function() {
      if (!length(registre)) {
        return(data.frame(type = character(), saison = character(),
                          station = character(), taxon = character(),
                          detail = character(), action = character()))
      }
      do.call(rbind, registre)
    }
  )
}
