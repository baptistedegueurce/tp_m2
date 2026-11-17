# =============================================================================
# fonctionnel_utils.R — Indices de diversité fonctionnelle
# Utilisé par 05_diversite_fonctionnelle.R
#
# Notations : P = matrice communautés × taxons d'abondances RELATIVES
#             (lignes sommant à 1), D = distance fonctionnelle (Gower, 0-1)
#             X = coordonnées des taxons dans l'espace fonctionnel (PCoA)
# =============================================================================

#' Met les traits au bon format pour la distance de Gower.
#' - taille maximale : log10 (une différence 10 -> 20 cm compte autant que
#'   50 -> 100 cm)
#' - position verticale : facteur ORDONNÉ du fond vers la surface
#' - autres traits texte : facteurs
preparer_traits <- function(df, traits) {
  x <- as.data.frame(df)[, traits, drop = FALSE]
  if ("taille_max_cm" %in% traits) x$taille_max_cm <- log10(x$taille_max_cm)
  if ("position_verticale" %in% traits) {
    x$position_verticale <- factor(
      x$position_verticale,
      levels = c("benthique", "démersal", "bentho-pélagique", "pélagique"),
      ordered = TRUE)
  }
  for (v in names(x)) if (is.character(x[[v]])) x[[v]] <- factor(x[[v]])
  x
}

#' PCoA d'une distance de Gower. On travaille sur sqrt(D), qui est euclidienne
#' (Gower 1966 ; Legendre & Legendre 2012) : pas d'axes à valeur propre
#' négative, donc pas de correction de Cailliez nécessaire.
pcoa_gower <- function(D) {
  n <- attr(D, "Size")
  pc <- cmdscale(sqrt(D), k = n - 1, eig = TRUE)
  pos <- pc$eig > 1e-10
  X <- pc$points[, pos[seq_len(ncol(pc$points))], drop = FALSE]
  colnames(X) <- paste0("PCoA", seq_len(ncol(X)))
  list(coords = X, var_expliquee = pc$eig[pos] / sum(pc$eig[pos]))
}

#' Qualité de l'espace : écart absolu moyen (mAD) entre la distance des traits
#' et la distance euclidienne dans un espace à k axes (Maire et al. 2015,
#' méthode reprise par mFD). Plus c'est bas, plus l'espace est fidèle.
qualite_espace <- function(D, X, kmax = 6) {
  ref <- as.vector(sqrt(D))
  kmax <- min(kmax, ncol(X))
  data.frame(
    n_axes = seq_len(kmax),
    mAD = vapply(seq_len(kmax), function(k)
      mean(abs(as.vector(dist(X[, seq_len(k), drop = FALSE])) - ref)),
      numeric(1))
  )
}

#' Dispersion fonctionnelle FDis (Laliberté & Legendre 2010) :
#' distance moyenne pondérée des taxons au centroïde pondéré de la communauté.
fdis <- function(X, P) {
  apply(P, 1, function(p) {
    centre <- colSums(X * p)
    sum(p * sqrt(rowSums(sweep(X, 2, centre)^2)))
  })
}

#' Entropie quadratique de Rao : Q = somme_ij d_ij p_i p_j
#' (dissimilarité fonctionnelle attendue entre deux individus tirés au hasard).
rao_q <- function(Dm, P) {
  apply(P, 1, function(p) as.numeric(t(p) %*% Dm %*% p))
}

#' Moyennes pondérées des traits (CWM).
#' Trait numérique -> moyenne pondérée ; trait catégoriel -> proportion de
#' l'abondance dans chaque modalité.
cwm <- function(traits, P) {
  out <- list()
  for (v in names(traits)) {
    t <- traits[[v]]
    if (is.numeric(t)) {
      out[[v]] <- as.vector(P %*% t)
    } else {
      for (niv in levels(t)) {
        out[[paste(v, niv, sep = "__")]] <- as.vector(P %*% as.numeric(t == niv))
      }
    }
  }
  res <- as.data.frame(out, check.names = FALSE)
  rownames(res) <- rownames(P)
  res
}

#' Effet standardisé (SES) de FDis et Rao par rapport à un modèle nul :
#' les traits sont permutés entre taxons (les abondances restent fixes).
#' SES = (observé - moyenne nulle) / écart-type nul.
#' SES < -2 : communauté plus homogène fonctionnellement que prévu par hasard
#' (filtrage environnemental) ; SES > 2 : plus dispersée (différenciation).
ses_fonctionnel <- function(X, Dm, P, n_perm = 999) {
  obs <- cbind(FDis = fdis(X, P), RaoQ = rao_q(Dm, P))
  n <- nrow(X)
  nul <- array(NA_real_, dim = c(nrow(P), 2, n_perm))
  for (k in seq_len(n_perm)) {
    o <- sample(n)
    nul[, 1, k] <- fdis(X[o, , drop = FALSE], P)
    nul[, 2, k] <- rao_q(Dm[o, o], P)
  }
  moy <- apply(nul, c(1, 2), mean)
  et  <- apply(nul, c(1, 2), sd)
  ses <- (obs - moy) / et
  ses[et < 1e-12] <- NA                 # 1 seul taxon : SES non défini
  # p unilatérale (rang de l'observé dans la distribution nulle)
  p_inf <- (apply(sweep(nul, c(1, 2), obs, "<="), c(1, 2), sum) + 1) / (n_perm + 1)
  data.frame(SES_FDis = ses[, 1], SES_RaoQ = ses[, 2],
             p_FDis = p_inf[, 1], p_RaoQ = p_inf[, 2],
             row.names = rownames(P))
}
