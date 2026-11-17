# =============================================================================
# figures_modeles.R — Figures communes des scripts de modélisation
# (effet saison, effet profondeur, synthèse des tests)
# Nécessite ggplot2 chargé et modeles_utils.R sourcé (preparer_reponse).
# =============================================================================

#' @param indices données (une ligne par trait de chalut)
#' @param reponses tableau des indices modélisés (indice, famille,
#'   transformation, label)
#' @param tabs sortie de modeliser_indices()
#' @param dossier_fig dossier où enregistrer les PNG
#' @param sujet fin du titre, ex. "les indices fonctionnels"
figures_modeles <- function(indices, reponses, tabs, dossier_fig, sujet) {
  ordre <- function(x) factor(x, levels = reponses$indice)   # ordre du tableau
  couleurs_saison <- c(printemps = "#66a61e", ete = "#e6ab02", automne = "#d95f02")
  labels <- setNames(reponses$label, reponses$indice)

  # a. Moyennes marginales par saison (IC 95 %) sur les données brutes
  brut <- lapply(seq_len(nrow(reponses)), function(i) {
    r <- reponses[i, ]
    d <- preparer_reponse(indices, r)
    if (r$transformation == "01") d$y <- d[[r$indice]]   # échelle d'origine
    data.frame(indice = r$indice, saison = d$saison, y = d$y)
  }) |> dplyr::bind_rows() |> dplyr::mutate(indice = ordre(indice))

  moy <- tabs$moyennes |> dplyr::mutate(indice = ordre(indice))
  if (!is.null(moy) && nrow(moy)) {
    # retour à l'échelle d'origine pour les proportions compressées
    p_moy <- ggplot() +
      geom_jitter(data = brut, aes(saison, y, colour = saison),
                  width = 0.12, alpha = 0.45, size = 1.4) +
      geom_pointrange(data = moy, aes(saison, moyenne, ymin = IC95_bas,
                                      ymax = IC95_haut, colour = saison),
                      size = 0.5, linewidth = 0.9,
                      position = position_nudge(x = 0.28)) +
      facet_wrap(~ indice, scales = "free_y", ncol = 4,
                 labeller = as_labeller(labels)) +
      scale_colour_manual(values = couleurs_saison, guide = "none") +
      labs(title = paste("Effet de la saison sur", sujet),
           subtitle = "Points pâles : traits de chalut ; point + barre : moyenne du modèle (profondeur moyenne) et IC 95 %",
           x = NULL, y = NULL) +
      theme_minimal(base_size = 10) +
      theme(strip.text = element_text(size = 8.5))
    print(p_moy)
    ggsave(file.path(dossier_fig, "effet_saison.png"), p_moy,
           width = 12, height = 9, dpi = 150, bg = "white")
  }

  # b. Effet de la profondeur (prédictions par saison)
  pred <- tabs$predictions |> dplyr::mutate(indice = ordre(indice))
  brut_prof <- lapply(seq_len(nrow(reponses)), function(i) {
    r <- reponses[i, ]
    d <- preparer_reponse(indices, r)
    data.frame(indice = r$indice, saison = d$saison,
               profondeur_m = d$profondeur_m, y = d$y)
  }) |> dplyr::bind_rows() |> dplyr::mutate(indice = ordre(indice))
  p_prof <- ggplot(pred, aes(profondeur_m, prediction, colour = saison,
                             fill = saison)) +
    geom_ribbon(aes(ymin = IC95_bas, ymax = IC95_haut), alpha = 0.12,
                colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(data = brut_prof, aes(profondeur_m, y), size = 1.4, alpha = 0.6,
               inherit.aes = FALSE, colour = "grey30") +
    facet_wrap(~ indice, scales = "free_y", ncol = 4,
               labeller = as_labeller(labels)) +
    scale_colour_manual(values = couleurs_saison) +
    scale_fill_manual(values = couleurs_saison) +
    labs(title = "Effet de la profondeur (modèle saison + profondeur, IC 95 %)",
         subtitle = "Échelle du modèle (variables transformées : log10, proportions comprimées pour la loi bêta)",
         x = "Profondeur moyenne du trait (m)", y = NULL, colour = NULL,
         fill = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", strip.text = element_text(size = 8.5))
  print(p_prof)
  ggsave(file.path(dossier_fig, "effet_profondeur.png"), p_prof,
         width = 12, height = 9, dpi = 150, bg = "white")

  # c. Synthèse des p-values (Holm)
  syn <- tabs$resume |>
    dplyr::select(indice, label, p_saison_holm, p_profondeur_holm) |>
    tidyr::pivot_longer(c(p_saison_holm, p_profondeur_holm),
                        names_to = "effet", values_to = "p") |>
    dplyr::mutate(effet = dplyr::recode(effet, p_saison_holm = "Saison",
                          p_profondeur_holm = "Profondeur"),
           label = factor(label, levels = rev(reponses$label)),
           signif = cut(p, c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
                        labels = c("< 0,001", "< 0,01", "< 0,05", "< 0,1", "n.s.")))
  p_syn <- ggplot(syn, aes(effet, label, fill = signif)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = formatC(p, format = "g", digits = 2)), size = 3) +
    scale_fill_manual(values = c("< 0,001" = "#08519c", "< 0,01" = "#3182bd",
                                 "< 0,05" = "#6baed6", "< 0,1" = "#c6dbef",
                                 "n.s." = "grey95"), drop = FALSE) +
    labs(title = "Synthèse des tests (LRT, p corrigées de Holm)",
         x = NULL, y = NULL, fill = "p") +
    theme_minimal(base_size = 10)
  print(p_syn)
  ggsave(file.path(dossier_fig, "synthese_tests.png"), p_syn,
         width = 7, height = 5.5, dpi = 150, bg = "white")
  invisible(NULL)
}
