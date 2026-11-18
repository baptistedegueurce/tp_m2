# =============================================================================
# modeles_simples.R — Modèles mixtes avec lme4 (scripts 04 et 05)
#
# Même démarche partout :
#   1. un modèle de référence fixé à l'avance (saison + profondeur + effet
#      aléatoire) ;
#   2. test de chaque effet par comparaison du modèle avec / sans cet effet
#      (test du rapport de vraisemblance, drop1) ;
#   3. comparaisons des saisons deux à deux (Tukey, emmeans) ;
#   4. vérification graphique des résidus (DHARMa).
# =============================================================================

#' Ajuste un modèle lme4 en gardant les avertissements (ajustement
#' singulier, convergence) au lieu de les afficher en vrac.
#' famille = NULL -> lmer (loi normale, REML = FALSE pour pouvoir comparer
#' des modèles) ; sinon glmer avec cette famille (ex. poisson()).
ajuster <- function(formule, data, famille = NULL) {
  avert <- character()
  m <- withCallingHandlers(
    tryCatch(
      # do.call : formule et données sont stockées dans l'appel du modèle,
      # ce qui permet à drop1() de réajuster le modèle sans l'effet
      if (is.null(famille)) do.call(lme4::lmer, list(formule, data = data, REML = FALSE))
      else do.call(lme4::glmer, list(formule, data = data, family = famille)),
      error = function(e) { avert <<- c(avert, conditionMessage(e)); NULL }),
    warning = function(w) {
      avert <<- c(avert, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(msg) {
      avert <<- c(avert, conditionMessage(msg))
      invokeRestart("muffleMessage")
    })
  if (!is.null(m)) attr(m, "avertissements") <- unique(avert)
  m
}

#' Test du rapport de vraisemblance de chaque effet : on retire l'effet du
#' modèle et on regarde si l'ajustement se dégrade (khi² à ddl = nombre de
#' paramètres retirés).
tester <- function(m) {
  t <- drop1(m, test = "Chisq")
  data.frame(effet = rownames(t)[-1], chi2 = t$LRT[-1], ddl = t$npar[-1],
             p = t$`Pr(Chi)`[-1], row.names = NULL)
}

#' Moyennes par saison (à profondeur moyenne) et comparaisons de Tukey,
#' sur l'échelle de la réponse.
saisons_tukey <- function(m) {
  em <- emmeans::emmeans(m, ~ saison, type = "response")
  moy <- as.data.frame(em)
  moy <- data.frame(saison = moy$saison, moyenne = moy[[2]],
                    IC95_bas = moy[[grep("lower|LCL", names(moy))[1]]],
                    IC95_haut = moy[[grep("upper|UCL", names(moy))[1]]])
  ctr <- as.data.frame(pairs(em, adjust = "tukey"))
  ctr <- data.frame(contraste = ctr$contrast, effet = ctr[[2]],
                    p_tukey = ctr$p.value)
  list(moyennes = moy, contrastes = ctr)
}

#' Résidus simulés (DHARMa) : figure + p-values des tests d'uniformité
#' (la loi choisie convient-elle ?) et de dispersion.
diagnostic <- function(m, fichier_png) {
  sim <- DHARMa::simulateResiduals(m, n = 500)
  png(fichier_png, width = 1600, height = 800, res = 150)
  plot(sim)
  dev.off()
  pv <- function(x) tryCatch(x$p.value, error = function(e) NA_real_)
  data.frame(p_uniformite = pv(DHARMa::testUniformity(sim, plot = FALSE)),
             p_dispersion = pv(DHARMa::testDispersion(sim, plot = FALSE)))
}

lier <- function(x) { x <- do.call(rbind, x); rownames(x) <- NULL; x }
