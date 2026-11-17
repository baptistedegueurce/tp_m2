# =============================================================================
# modeles_utils.R — Modèles mixtes (glmmTMB) : ajustement, sélection,
# tests, comparaisons, diagnostics, prédictions.
# Réutilisable pour les indices taxonomiques (03) et fonctionnels (05b).
# =============================================================================

#' Loi de probabilité glmmTMB à partir d'un nom simple.
#' gaussian : indice continu non borné (SES, CWM de traits continus)
#' poisson / genpois / nbinom2 : comptages (richesse). genpois accepte la
#'   sous-dispersion fréquente pour la richesse ; nbinom2 la sur-dispersion.
#' beta : proportion ou indice borné dans ]0 ; 1[ (FRic, FEve, FDiv, CWM catégoriel)
famille_modele <- function(nom) {
  switch(nom,
         gaussian = gaussian(),
         poisson  = poisson(),
         genpois  = glmmTMB::genpois(),
         nbinom2  = glmmTMB::nbinom2(),
         beta     = glmmTMB::beta_family(),
         stop("Famille inconnue : ", nom))
}

#' Ramène une proportion dans ]0 ; 1[ pour une loi bêta
#' (Smithson & Verkuilen 2006) : y' = (y (n - 1) + 0,5) / n
compresser_01 <- function(y) {
  n <- sum(!is.na(y))
  (y * (n - 1) + 0.5) / n
}

#' AICc d'un modèle (AIC corrigé pour petits effectifs, recommandé quand
#' n / k < 40, ce qui est notre cas).
aicc <- function(m) {
  ll <- logLik(m)
  k <- attr(ll, "df")
  n <- nobs(m)
  -2 * as.numeric(ll) + 2 * k + 2 * k * (k + 1) / max(n - k - 1, 1)
}

#' Ajuste un glmmTMB en capturant les avertissements (convergence, ajustement
#' singulier) au lieu de les laisser s'afficher en vrac.
ajuster_glmm <- function(formule, data, famille) {
  avert <- character()
  m <- withCallingHandlers(
    tryCatch(glmmTMB::glmmTMB(formule, data = data, family = famille),
             error = function(e) { avert <<- c(avert, conditionMessage(e)); NULL }),
    warning = function(w) {
      avert <<- c(avert, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  if (!is.null(m)) attr(m, "avertissements") <- unique(avert)
  m
}

#' Jeu de modèles candidats, du plus simple au plus complet.
#' L'effet aléatoire station tient compte des 3 passages sur une même station.
formules_candidates <- function(reponse) {
  re <- "+ (1 | station)"
  f <- c(
    nul               = "1",
    saison            = "saison",
    profondeur        = "profondeur_m",
    saison_profondeur = "saison + profondeur_m",
    interaction       = "saison * profondeur_m"
  )
  lapply(f, function(x) as.formula(paste(reponse, "~", x, re)))
}

#' Ajuste les candidats et les classe par AICc.
selection_aicc <- function(modeles) {
  ok <- !vapply(modeles, is.null, logical(1))
  tab <- data.frame(
    modele = names(modeles)[ok],
    k      = vapply(modeles[ok], function(m) attr(logLik(m), "df"), numeric(1)),
    AICc   = vapply(modeles[ok], aicc, numeric(1))
  )
  tab$delta_AICc <- tab$AICc - min(tab$AICc)
  tab$poids_akaike <- exp(-tab$delta_AICc / 2) / sum(exp(-tab$delta_AICc / 2))
  tab[order(tab$AICc), ]
}

#' Tests du rapport de vraisemblance (chaque terme retiré du modèle complet).
tests_lrt <- function(m) {
  t <- drop1(m, test = "Chisq")
  data.frame(effet = rownames(t)[-1], chi2 = t$LRT[-1], ddl = t$Df[-1],
             p = t$`Pr(>Chi)`[-1], row.names = NULL)
}

#' Coefficients du modèle (échelle du lien).
coefficients_modele <- function(m) {
  co <- summary(m)$coefficients$cond
  data.frame(terme = rownames(co), estimation = co[, 1], erreur_type = co[, 2],
             z = co[, 3], p = co[, 4], row.names = NULL)
}

#' Moyennes marginales par saison et comparaisons deux à deux (Tukey),
#' sur l'échelle de la réponse.
comparaisons_saisons <- function(m) {
  em <- emmeans::emmeans(m, ~ saison, type = "response")
  moy <- as.data.frame(em)
  # les noms de colonnes d'emmeans dépendent de la loi (emmean, rate,
  # response, prob… ; lower.CL ou asymp.LCL) : on les normalise
  moy <- data.frame(saison = moy$saison, moyenne = moy[[2]], SE = moy$SE,
                    IC95_bas = moy[[grep("lower|LCL", names(moy))[1]]],
                    IC95_haut = moy[[grep("upper|UCL", names(moy))[1]]])
  ctr <- as.data.frame(pairs(em, adjust = "tukey"))
  ctr <- data.frame(contraste = ctr$contrast,
                    type_effet = names(ctr)[2],     # estimate, ratio, odds.ratio
                    effet = ctr[[2]], SE = ctr$SE,
                    statistique = ctr[[grep("^(z|t)\\.ratio$", names(ctr))[1]]],
                    p_tukey = ctr$p.value)
  list(moyennes = moy, contrastes = ctr)
}

#' Diagnostics par résidus simulés (DHARMa, Hartig 2022).
#' Un p < 0,05 signale un problème : mauvaise loi (uniformité), sur- ou
#' sous-dispersion, valeurs extrêmes, autocorrélation spatiale résiduelle.
#' @param groupe nom d'une colonne de `data` (ex. "trait_id") : pour des
#'   données individuelles, l'autocorrélation spatiale est testée sur les
#'   résidus regroupés par trait de chalut (individus d'un même trait =
#'   mêmes coordonnées).
diagnostics_dharma <- function(m, data, fichier_png = NULL, groupe = NULL) {
  # Simulations faites par glmmTMB::simulate puis passées à DHARMa : évite les
  # incompatibilités entre versions de DHARMa et de glmmTMB.
  sims <- as.matrix(simulate(m, nsim = 1000))
  sim <- DHARMa::createDHARMa(
    simulatedResponse = sims,
    observedResponse = model.response(model.frame(m)),
    fittedPredictedResponse = predict(m, type = "response"),
    integerResponse = family(m)$family %in% c("poisson", "genpois", "nbinom1",
                                              "nbinom2"))
  pv <- function(expr) tryCatch(expr$p.value, error = function(e) NA_real_)
  res <- data.frame(
    p_uniformite = pv(DHARMa::testUniformity(sim, plot = FALSE)),
    p_dispersion = pv(DHARMa::testDispersion(sim, plot = FALSE)),
    p_extremes   = pv(DHARMa::testOutliers(sim, type = "bootstrap", plot = FALSE)),
    p_autocorrelation_spatiale = pv(if (is.null(groupe)) {
      DHARMa::testSpatialAutocorrelation(sim, x = data$lon, y = data$lat,
                                         plot = FALSE)
    } else {
      g <- factor(data[[groupe]])
      simg <- DHARMa::recalculateResiduals(sim, group = g)
      DHARMa::testSpatialAutocorrelation(simg, x = tapply(data$lon, g, mean),
                                         y = tapply(data$lat, g, mean),
                                         plot = FALSE)
    })
  )
  if (!is.null(fichier_png)) {
    png(fichier_png, width = 1600, height = 800, res = 150)
    plot(sim)
    dev.off()
  }
  res
}

#' R² de vraisemblance (Magee 1990 ; Nakagawa & Schielzeth 2013, éq. 1) :
#' R² = 1 - exp(-2 (logL_modèle - logL_nul) / n), le modèle nul gardant
#' l'effet station. Valable quelle que soit la loi (gaussienne, bêta, comptage).
#' Part de la variation expliquée par les effets fixes (saison, profondeur).
r2_vraisemblance <- function(m, m_nul) {
  if (is.null(m_nul)) return(NA_real_)
  1 - exp(-2 * (as.numeric(logLik(m)) - as.numeric(logLik(m_nul))) / nobs(m))
}

#' Prédictions (effets fixes seuls, IC 95 %) le long de la profondeur,
#' pour chaque saison.
predire_profondeur <- function(m, data, n = 50) {
  grille <- expand.grid(
    profondeur_m = seq(min(data$profondeur_m), max(data$profondeur_m), length.out = n),
    saison = levels(droplevels(data$saison)),
    station = NA
  )
  grille$saison <- factor(grille$saison, levels = levels(data$saison))
  pr <- predict(m, newdata = grille, type = "link", se.fit = TRUE, re.form = NA,
                allow.new.levels = TRUE)
  inv <- family(m)$linkinv
  grille$prediction <- inv(pr$fit)
  grille$IC95_bas   <- inv(pr$fit - 1.96 * pr$se.fit)
  grille$IC95_haut  <- inv(pr$fit + 1.96 * pr$se.fit)
  grille
}

#' Prépare la variable réponse d'un indice : transformation (log10, ou
#' proportion ramenée dans ]0 ; 1[ pour la loi bêta), lignes NA retirées.
preparer_reponse <- function(d, r) {
  y <- d[[r$indice]]
  y <- switch(r$transformation, log10 = log10(y), log10p1 = log10(y + 1),
              `01` = compresser_01(y), y)
  d$y <- y
  d[!is.na(d$y), ]
}

#' Modélise une série d'indices (une ligne de `reponses` = un indice avec sa
#' loi, sa transformation et son libellé) avec la même démarche pour tous :
#' sélection AICc, modèle de référence saison + profondeur + (1 | station),
#' LRT (Holm), emmeans, DHARMa, R², sensibilité aux traits douteux,
#' prédictions. Les données doivent contenir saison, station, profondeur_m,
#' lon, lat, douteux.
#' @return liste de tableaux (resume, tests, selection, coefficients,
#'   moyennes, contrastes, diagnostics, sensibilite, predictions)
modeliser_indices <- function(indices, reponses, dossier_diag) {
  res <- list(selection = list(), tests = list(), coefficients = list(),
              moyennes = list(), contrastes = list(), diagnostics = list(),
              sensibilite = list(), predictions = list(), resume = list())

  for (i in seq_len(nrow(reponses))) {
    r <- reponses[i, ]
    message(sprintf("\n--- %s (%s) ---", r$indice, r$famille))
    d <- preparer_reponse(indices, r)
    fam <- famille_modele(r$famille)

    # 1. Candidats et sélection AICc
    cand <- lapply(formules_candidates("y"), ajuster_glmm, data = d, famille = fam)
    sel <- selection_aicc(cand)
    res$selection[[r$indice]] <- cbind(indice = r$indice, sel)

    # 2. Modèle de référence (fixé a priori)
    m <- cand$saison_profondeur
    if (is.null(m)) { warning("Échec du modèle de référence pour ", r$indice); next }
    avert <- attr(m, "avertissements")

    lrt <- tests_lrt(m)
    res$tests[[r$indice]] <- cbind(indice = r$indice, famille = r$famille,
                                   n = nobs(m), lrt)
    res$coefficients[[r$indice]] <- cbind(indice = r$indice, coefficients_modele(m))

    # 3. Comparaisons entre saisons
    cs <- tryCatch(comparaisons_saisons(m), error = function(e) NULL)
    if (!is.null(cs)) {
      res$moyennes[[r$indice]]   <- cbind(indice = r$indice, cs$moyennes)
      res$contrastes[[r$indice]] <- cbind(indice = r$indice, cs$contrastes)
    }

    # 4. Diagnostics
    diag <- diagnostics_dharma(
      m, d, file.path(dossier_diag, paste0(make.names(r$indice), ".png")))
    r2 <- r2_vraisemblance(m, cand$nul)
    res$diagnostics[[r$indice]] <- data.frame(
      indice = r$indice, diag, R2_vraisemblance = r2,
      variance_station = unname(unlist(glmmTMB::VarCorr(m)$cond$station)),
      avertissements = paste(avert, collapse = " | "))

    # 5. Sensibilité aux traits douteux
    d2 <- d[!d$douteux, ]
    m2 <- ajuster_glmm(y ~ saison + profondeur_m + (1 | station), d2, fam)
    if (!is.null(m2)) {
      res$sensibilite[[r$indice]] <- cbind(indice = r$indice, n = nobs(m2),
                                           tests_lrt(m2))
    }

    # 6. Prédictions le long de la profondeur
    pr <- predire_profondeur(m, d)
    res$predictions[[r$indice]] <- cbind(indice = r$indice, label = r$label, pr)

    res$resume[[r$indice]] <- data.frame(
      indice = r$indice, label = r$label, famille = r$famille, n = nobs(m),
      meilleur_modele_AICc = sel$modele[1],
      p_saison = lrt$p[lrt$effet == "saison"],
      p_profondeur = lrt$p[lrt$effet == "profondeur_m"],
      R2 = r2,
      probleme_diagnostic = any(unlist(diag) < 0.05, na.rm = TRUE),
      ajustement_singulier = any(grepl("singular|convergence|Hessian", avert,
                                       ignore.case = TRUE)))
  }

  lier <- function(x) { x <- do.call(rbind, x); rownames(x) <- NULL; x }
  tabs <- lapply(res, lier)

  # Correction de Holm sur l'ensemble des tests saison / profondeur
  tabs$tests$p_holm <- p.adjust(tabs$tests$p, method = "holm")
  tabs$resume <- tabs$resume |>
    dplyr::left_join(tabs$tests |> dplyr::filter(effet == "saison") |>
                dplyr::select(indice, p_saison_holm = p_holm), by = "indice") |>
    dplyr::left_join(tabs$tests |> dplyr::filter(effet == "profondeur_m") |>
                dplyr::select(indice, p_profondeur_holm = p_holm), by = "indice") |>
    dplyr::relocate(p_saison_holm, .after = p_saison) |>
    dplyr::relocate(p_profondeur_holm, .after = p_profondeur)
  tabs
}
