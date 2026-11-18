# =============================================================================
# 06_diversite_fonctionnelle.R
# Structure et diversité fonctionnelles des communautés
#
# Entrées : data/processed/preparation.rds (script 02)
#           data/ref/traits_fonctionnels.csv
# Paramètres : R/00_config.R (traits utilisés, taxons gardés, pondération…)
#
# Étapes :
#   1. Taxons × traits (poissons, traits attribués avec une confiance suffisante)
#   2. Distance de Gower + PCoA, qualité de l'espace fonctionnel
#   3. Indices par trait de chalut : FDis, Rao, FRic, FEve, FDiv, CWM
#   4. Modèle nul : SES de FDis et Rao (indépendants de la richesse)
#   5. Figures (la modélisation est dans 06b_modeles_fonctionnels.R)
#
# Sorties : data/processed/fonctionnel_*.csv, outputs/figures/fonctionnel/
# =============================================================================

pkgs <- c("here", "dplyr", "tidyr", "ggplot2", "FD")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants)) install.packages(manquants)

library(dplyr)
library(ggplot2)

source(here::here("R", "00_config.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "fonctionnel_utils.R"), encoding = "UTF-8")

set.seed(graine_aleatoire)
dossier_fig <- here::here("outputs", "figures", "fonctionnel")
dir.create(dossier_fig, recursive = TRUE, showWarnings = FALSE)

prep <- readRDS(file.path(dossiers$processed, "preparation.rds"))

# =============================================================================
# 1. Taxons et traits
# =============================================================================
traits_ref <- read.csv(fichier_traits_fonctionnels, encoding = "UTF-8",
                       stringsAsFactors = FALSE)

manque <- setdiff(prep$taxons$taxon, traits_ref$taxon)
if (length(manque)) {
  stop("Taxons sans ligne dans traits_fonctionnels.csv : ",
       paste(manque, collapse = ", "))
}

taxons_f <- prep$taxons |>
  select(taxon, groupe, effectif_total) |>
  inner_join(traits_ref, by = "taxon") |>
  mutate(garde = groupe %in% groupes_analyse_fonctionnelle &
           !confiance_espece %in% confiance_exclue)

# Un taxon gardé doit avoir tous ses traits renseignés
incomplets <- taxons_f |>
  filter(garde) |>
  filter(if_any(all_of(traits_utilises), ~ is.na(.x) | .x == ""))
if (nrow(incomplets)) {
  stop("Traits manquants pour : ", paste(incomplets$taxon, collapse = ", "))
}

exclus <- taxons_f |> filter(!garde)
message("Taxons exclus de l'analyse fonctionnelle : ",
        paste0(exclus$taxon, " (", exclus$groupe, ", confiance ",
               exclus$confiance_espece, ")", collapse = "; "))

taxons_f <- taxons_f |> filter(garde)

# =============================================================================
# 2. Matrice d'abondance (traits de chalut × taxons)
# =============================================================================
colonne_poids <- c(abondance = "cpue_n_km", biomasse = "cpue_g_km")[ponderation_fonctionnelle]
if (is.na(colonne_poids)) stop("ponderation_fonctionnelle : 'abondance' ou 'biomasse'")

abond <- prep$captures |>
  filter(taxon %in% taxons_f$taxon) |>
  select(trait_id, taxon, valeur = all_of(unname(colonne_poids)))
if (anyNA(abond$valeur)) stop("Valeurs manquantes de ", colonne_poids)

A <- abond |>
  tidyr::pivot_wider(names_from = taxon, values_from = valeur, values_fill = 0) |>
  as.data.frame()
rownames(A) <- A$trait_id
A <- as.matrix(A[, -1, drop = FALSE])

vides <- rownames(A)[rowSums(A) == 0]
if (length(vides)) {
  message("Traits sans poisson retenu (exclus) : ", paste(vides, collapse = ", "))
}
A <- A[rowSums(A) > 0, colSums(A) > 0, drop = FALSE]
taxons_f <- taxons_f |> filter(taxon %in% colnames(A))
A <- A[, taxons_f$taxon, drop = FALSE]
P <- A / rowSums(A)                     # abondances relatives

message(sprintf("%d traits de chalut × %d taxons", nrow(P), ncol(P)))

# =============================================================================
# 3. Espace fonctionnel
# =============================================================================
tr <- preparer_traits(taxons_f, traits_utilises)
rownames(tr) <- taxons_f$taxon

D  <- FD::gowdis(tr)                    # Gower : traits mixtes, ordinal géré
Dm <- as.matrix(D)
esp <- pcoa_gower(D)
X   <- esp$coords

qualite <- qualite_espace(D, X, kmax = 6)
n_axes <- qualite$n_axes[which.min(qualite$mAD > 0.05)]   # 1er espace avec mAD < 0,05
if (all(qualite$mAD > 0.05)) n_axes <- max(qualite$n_axes)
message(sprintf("Espace fonctionnel : %d axes retenus (mAD = %.3f, %.0f %% de variance)",
                n_axes, qualite$mAD[n_axes],
                100 * sum(esp$var_expliquee[seq_len(n_axes)])))
Xk <- X[, seq_len(n_axes), drop = FALSE]

# =============================================================================
# 4. Indices par trait de chalut
# =============================================================================
indices <- data.frame(
  trait_id = rownames(P),
  S        = rowSums(P > 0),
  FDis     = fdis(Xk, P),
  RaoQ     = rao_q(Dm, P)
)

# FRic, FEve, FDiv : package FD (calculés seulement s'il y a assez de taxons
# dans la communauté ; NA sinon)
fd <- tryCatch(
  suppressWarnings(FD::dbFD(D, A, w.abun = TRUE, corr = "sqrt", m = n_axes,
                            calc.FRic = TRUE, calc.FDiv = TRUE, calc.CWM = FALSE,
                            stand.FRic = TRUE, messages = FALSE)),
  error = function(e) { warning("dbFD : ", conditionMessage(e)); NULL })
if (!is.null(fd)) {
  # FD réduit le nombre d'axes pour FRic quand les communautés ont peu de
  # taxons (il faut plus de taxons que d'axes) : qualité reportée ici
  message(sprintf("FRic calculée dans un espace réduit (qualité R² = %.2f)",
                  fd$qual.FRic))
  indices$FRic <- fd$FRic[indices$trait_id]
  indices$FEve <- fd$FEve[indices$trait_id]
  indices$FDiv <- fd$FDiv[indices$trait_id]
}

# CWM
cw <- cwm(tr, P)
names(cw) <- paste0("CWM_", names(cw))
if ("CWM_taille_max_cm" %in% names(cw)) {
  # moyenne de log10 -> moyenne géométrique en cm, plus lisible
  cw$CWM_taille_max_cm <- 10^cw$CWM_taille_max_cm
}

# Modèle nul
message(sprintf("Modèle nul : %d permutations…", n_permutations_ses))
ses <- ses_fonctionnel(Xk, Dm, P, n_permutations_ses)

indices <- indices |>
  bind_cols(ses, cw) |>
  left_join(prep$traits |> select(trait_id, saison, station, date, profondeur_m,
                                  lat, lon, douteux),
            by = "trait_id") |>
  relocate(trait_id, saison, station, date, profondeur_m, lat, lon, douteux) |>
  mutate(saison = factor(saison, levels = c("printemps", "ete", "automne")))

# =============================================================================
# 5. Modélisation : voir 06b_modeles_fonctionnels.R
# =============================================================================

# =============================================================================
# 6. Export
# =============================================================================
espace <- data.frame(taxon = rownames(X), X[, seq_len(min(4, ncol(X)))],
                     row.names = NULL) |>
  left_join(taxons_f |> select(taxon, nom_commun, espece_reference,
                               all_of(traits_utilises)), by = "taxon")

write.csv(indices, file.path(dossiers$processed, "fonctionnel_indices.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(espace, file.path(dossiers$processed, "fonctionnel_espace.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(qualite, file.path(dossiers$processed, "fonctionnel_qualite_espace.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# =============================================================================
# 7. Figures
# =============================================================================
couleurs_saison <- c(printemps = "#66a61e", ete = "#e6ab02", automne = "#d95f02")
sauver <- function(p, nom, l = 9, h = 6) {
  print(p)
  ggsave(file.path(dossier_fig, nom), p, width = l, height = h, dpi = 150,
         bg = "white")
}

# 7a. Espace fonctionnel (axes 1-2)
ab_tot <- colSums(A)
espace$abondance <- ab_tot[espace$taxon]
pct <- round(100 * esp$var_expliquee[1:2])
p_esp <- ggplot(espace, aes(PCoA1, PCoA2)) +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(size = abondance, colour = guilde_trophique), alpha = 0.8) +
  {if (requireNamespace("ggrepel", quietly = TRUE))
    ggrepel::geom_text_repel(aes(label = nom_commun), size = 3)
   else geom_text(aes(label = nom_commun), size = 3, vjust = -1.1)} +
  scale_size_area(max_size = 12, guide = "none") +
  labs(title = "Espace fonctionnel des taxons",
       subtitle = "Taille des points proportionnelle à l'abondance totale (CPUE)",
       x = sprintf("PCoA 1 (%d %%)", pct[1]), y = sprintf("PCoA 2 (%d %%)", pct[2]),
       colour = "Guilde trophique") +
  theme_minimal(base_size = 11)
sauver(p_esp, "espace_fonctionnel.png", 9, 7)

# 7b. Qualité de l'espace
p_q <- ggplot(qualite, aes(n_axes, mAD)) +
  geom_line() + geom_point(size = 2) +
  geom_point(data = qualite[n_axes, ], colour = "red", size = 3.5) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "grey50") +
  labs(title = "Qualité de l'espace fonctionnel",
       subtitle = "Écart moyen entre distance des traits et distance dans l'espace (rouge : retenu)",
       x = "Nombre d'axes de la PCoA", y = "mAD") +
  theme_minimal(base_size = 11)
sauver(p_q, "qualite_espace.png", 7, 4.5)

# 7c. Indices par saison
ind_long <- indices |>
  select(trait_id, saison, any_of(c("S", "FDis", "RaoQ", "FRic", "FEve",
                                    "FDiv", "SES_FDis", "SES_RaoQ"))) |>
  tidyr::pivot_longer(-c(trait_id, saison), names_to = "indice",
                      values_to = "valeur") |>
  filter(!is.na(valeur)) |>
  mutate(indice = factor(indice, levels = c("S", "FDis", "RaoQ", "FRic", "FEve",
                                            "FDiv", "SES_FDis", "SES_RaoQ")))
ref_ses <- data.frame(indice = factor(c("SES_FDis", "SES_RaoQ"),
                                      levels = levels(ind_long$indice)),
                      y = 0)
p_ind <- ggplot(ind_long, aes(saison, valeur, colour = saison)) +
  geom_hline(data = ref_ses, aes(yintercept = y), linetype = 2,
             colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.7) +
  facet_wrap(~ indice, scales = "free_y", nrow = 2) +
  scale_colour_manual(values = couleurs_saison, guide = "none") +
  labs(title = "Indices fonctionnels par saison (un point = un trait de chalut)",
       subtitle = "SES : 0 = attendu au hasard compte tenu de la richesse",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11)
sauver(p_ind, "indices_par_saison.png", 11, 6.5)

# 7d. CWM numériques
num <- intersect(c("CWM_taille_max_cm", "CWM_niveau_trophique"), names(indices))
if (length(num)) {
  p_num <- indices |>
    select(trait_id, saison, profondeur_m, all_of(num)) |>
    tidyr::pivot_longer(all_of(num), names_to = "trait", values_to = "valeur") |>
    mutate(trait = recode(trait,
                          CWM_taille_max_cm = "Taille max. moyenne (cm, moy. géom.)",
                          CWM_niveau_trophique = "Niveau trophique moyen")) |>
    ggplot(aes(saison, valeur, colour = saison)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 1.6, alpha = 0.7) +
    facet_wrap(~ trait, scales = "free_y") +
    scale_colour_manual(values = couleurs_saison, guide = "none") +
    labs(title = "Moyennes pondérées des traits (CWM)", x = NULL, y = NULL) +
    theme_minimal(base_size = 11)
  sauver(p_num, "cwm_numeriques.png", 9, 4.5)
}

# 7e. CWM catégoriels : composition moyenne par saison
cat_cols <- grep("__", names(indices), value = TRUE)
if (length(cat_cols)) {
  p_cat <- indices |>
    select(saison, all_of(cat_cols)) |>
    tidyr::pivot_longer(-saison, names_to = "var", values_to = "prop") |>
    tidyr::separate(var, into = c("trait", "modalite"), sep = "__") |>
    mutate(trait = sub("^CWM_", "", trait)) |>
    group_by(saison, trait, modalite) |>
    summarise(prop = mean(prop), .groups = "drop") |>
    ggplot(aes(prop, modalite, fill = saison)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    facet_wrap(~ trait, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = couleurs_saison) +
    scale_x_continuous(labels = function(x) paste0(round(100 * x), " %")) +
    labs(title = "Composition fonctionnelle moyenne par saison",
         subtitle = "Part moyenne de l'abondance dans chaque modalité (CWM)",
         x = NULL, y = NULL, fill = NULL) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom")
  sauver(p_cat, "cwm_categoriels.png", 11, 8)
}

# =============================================================================
# 8. Résumé
# =============================================================================
message("\n================ RÉSUMÉ ================")
message(sprintf("%d taxons, %d traits de chalut, pondération : %s",
                ncol(P), nrow(P), ponderation_fonctionnelle))
print(indices |>
        group_by(saison) |>
        summarise(n = n(), S = mean(S), FDis = mean(FDis), RaoQ = mean(RaoQ),
                  SES_FDis = mean(SES_FDis, na.rm = TRUE),
                  SES_RaoQ = mean(SES_RaoQ, na.rm = TRUE)) |>
        mutate(across(where(is.numeric), ~ round(.x, 2))) |>
        as.data.frame())
message("\nModélisation (saison, profondeur, station) : lancer R/06b_modeles_fonctionnels.R")
