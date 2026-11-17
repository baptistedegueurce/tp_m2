# =============================================================================
# nourricerie_utils.R — Structure en taille, stades et condition des individus
# Utilisé par 07_nourricerie.R
# =============================================================================

#' Moyenne pondérée sans NA.
moy_pond <- function(x, w) {
  ok <- !is.na(x) & !is.na(w)
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

#' Quantile pondéré (fonction de répartition empirique pondérée).
quantile_pond <- function(x, w, p) {
  ok <- !is.na(x) & !is.na(w)
  x <- x[ok]; w <- w[ok]
  if (!length(x)) return(rep(NA_real_, length(p)))
  o <- order(x); x <- x[o]; w <- cumsum(w[o]) / sum(w)
  vapply(p, function(q) x[which(w >= q)[1]], numeric(1))
}

#' Relation taille-poids par taxon et condition relative de Le Cren (1951).
#' log10(W) = log10(a) + b log10(L), ajustée sur tous les individus du taxon
#' (aberrants exclus) ; Kn = W / (a L^b). Kn = 1 : poids attendu pour la
#' taille ; Kn > 1 : individu plus lourd que la moyenne de l'espèce.
#' Kn ne dépend pas de la taille (contrairement au K de Fulton si b != 3),
#' ce qui permet de comparer des saisons où les tailles diffèrent.
condition_le_cren <- function(ind, n_min = 20) {
  ind$Kn <- NA_real_
  params <- list()
  for (tx in unique(ind$taxon)) {
    i <- which(ind$taxon == tx & !ind$aberrant_taille_poids &
                 ind$taille_mm > 0 & ind$poids_g > 0)
    if (length(i) < n_min) next
    fit <- lm(log10(poids_g) ~ log10(taille_mm), data = ind[i, ])
    co <- coef(fit)
    ind$Kn[i] <- ind$poids_g[i] / 10^(co[1] + co[2] * log10(ind$taille_mm[i]))
    ic <- confint(fit)[2, ]
    params[[tx]] <- data.frame(taxon = tx, n = length(i), a = 10^co[1],
                               b = co[2], b_IC95_bas = ic[1], b_IC95_haut = ic[2],
                               r2 = summary(fit)$r.squared)
  }
  list(individus = ind, parametres = do.call(rbind, params))
}

#' Pente du spectre de biomasse normalisé (NBSS) d'un trait de chalut.
#' Les individus sont rangés en classes de poids de largeur doublante
#' (log2) ; la biomasse de chaque classe (élevée, par km) est divisée par la
#' largeur de la classe. La pente de log10(biomasse normalisée) en fonction
#' de log10(poids) résume la structure en taille : plus elle est négative,
#' plus la communauté est dominée par les petits individus.
pente_nbss <- function(poids_g, facteur, distance_km, n_classes_min = 4) {
  ok <- poids_g > 0 & !is.na(facteur)
  if (sum(ok) < 5) return(c(pente = NA_real_, n_classes = 0))
  p <- poids_g[ok]; f <- facteur[ok]
  bornes <- 2^(floor(log2(min(p))):ceiling(log2(max(p)) + 1e-9))
  if (length(bornes) < 2) return(c(pente = NA_real_, n_classes = 0))
  classe <- cut(p, bornes, right = FALSE, include.lowest = TRUE)
  biomasse <- tapply(p * f, classe, sum) / distance_km
  largeur <- diff(bornes)
  centre <- sqrt(bornes[-1] * bornes[-length(bornes)])
  nb <- biomasse / largeur
  garde <- !is.na(nb) & nb > 0
  if (sum(garde) < n_classes_min) return(c(pente = NA_real_, n_classes = sum(garde)))
  fit <- lm(log10(nb[garde]) ~ log10(centre[garde]))
  c(pente = unname(coef(fit)[2]), n_classes = sum(garde))
}
