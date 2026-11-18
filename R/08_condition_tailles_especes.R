# =============================================================================
# 08_condition_tailles_especes.R
# Qu'est-ce qui fait varier la condition et la taille des poissons, dans
# l'espace (profondeur) et dans le temps (saison) ? Mêmes modèles simples
# pour chaque espèce, puis vérification de la robustesse de chaque effet.
#
# Entrées : data/processed/preparation.rds          (script 02)
#           data/processed/environnement_traits.csv (script 06)
#           data/ref/traits_fonctionnels.csv        (taille de 1re maturité)
#
# Espèces : au moins `n_min_individus_espece` individus mesurés et pesés.
# Saisons : pour chaque espèce, seules les saisons où elle a été prise dans
#   au moins `n_min_traits_saison` traits entrent dans les tests (sinon
#   l'effet « saison » se confond avec l'effet d'un ou deux traits).
#
# 1. Description : relation taille-poids (a, b), effectifs par saison,
#    condition relative de Le Cren (Kn), distributions de tailles
# 2. Modèles de référence (modèles linéaires mixtes, à l'échelle de l'individu)
#      condition : log10(W) ~ log10(L) + saison + profondeur + (1 | trait)
#                  = différence de poids À TAILLE ÉGALE (ANCOVA ; préférée à
#                    l'analyse directe d'un ratio comme Kn, Packard & Boardman 1988)
#      taille    : log10(L) ~ saison + profondeur + (1 | trait)
#    Profondeur centrée et exprimée par tranche de 10 m : effets en % de
#    variation pour +10 m, comparables entre espèces.
#    Tests du rapport de vraisemblance (LRT), correction de Holm sur
#    l'ensemble espèces × réponses × effets, comparaisons de Tukey.
# 3. Vérifications de robustesse de chaque effet :
#    a. au niveau du trait de chalut (unité réellement indépendante) :
#       moyenne par trait des résidus de condition / du log de la taille,
#       puis lm(moyenne ~ saison + profondeur), test F ;
#    b. retrait d'un trait à tour de rôle (p maximal obtenu) ;
#    c. sans les traits douteux.
#    Verdict : « robuste » si l'effet reste significatif dans toutes les
#    versions ; « fragile » s'il est significatif dans le modèle de référence
#    mais pas dans au moins une vérification ; « non significatif » sinon.
# 4. Croissance apparente des juvéniles entre l'été et l'automne (décalage
#    de la taille médiane, IC 95 % par bootstrap des traits)
#
# Variables Copernicus : non utilisées comme covariables (voir 00_config.R).
#
# Sorties : data/processed/condition_especes/
#           outputs/figures/condition_especes/ (+ diagnostics/)
# =============================================================================

pkgs <- c("here", "dplyr", "tidyr", "ggplot2", "glmmTMB", "DHARMa", "emmeans")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants)) install.packages(manquants)

library(dplyr)
library(ggplot2)

source(here::here("R", "00_config.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "modeles_utils.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "nourricerie_utils.R"), encoding = "UTF-8")

set.seed(graine_aleatoire)
dossier_out  <- file.path(dossiers$processed, "condition_especes")
dossier_fig  <- here::here("outputs", "figures", "condition_especes")
dossier_diag <- file.path(dossier_fig, "diagnostics")
unlink(list.files(dossier_out, full.names = TRUE))   # anciennes sorties
unlink(list.files(dossier_fig, pattern = "\\.png$", full.names = TRUE))
for (d in c(dossier_out, dossier_fig, dossier_diag)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
saisons <- c("printemps", "ete", "automne")
couleurs_saison <- c(printemps = "#66a61e", ete = "#e6ab02", automne = "#d95f02")

prep <- readRDS(file.path(dossiers$processed, "preparation.rds"))
fichier_env <- file.path(dossiers$processed, "environnement_traits.csv")
if (!file.exists(fichier_env)) stop("Lancer d'abord R/06_variables_environnementales.R")
env <- read.csv(fichier_env, encoding = "UTF-8")
traits_ref <- read.csv(fichier_traits_fonctionnels, encoding = "UTF-8")
noms <- setNames(traits_ref$nom_commun, traits_ref$taxon)
lab_esp <- function(tx) ifelse(is.na(noms[tx]), tx, paste0(noms[tx], "\n(", tx, ")"))

# =============================================================================
# 1. Individus
# =============================================================================
prof_moy <- mean(env[[covariable_spatiale]][env$analysable])

ind <- prep$individus |>
  filter(groupe == "poisson", !aberrant_taille_poids,
         taille_mm > 0, poids_g > 0) |>
  inner_join(env |> filter(analysable) |>
               select(trait_id, station, lat, lon, douteux,
                      prof = all_of(covariable_spatiale)),
             by = c("trait_id", "station")) |>
  left_join(traits_ref |> select(taxon, taille_maturite_cm), by = "taxon") |>
  mutate(saison = factor(saison, levels = saisons),
         trait_id = factor(trait_id),
         prof10 = (prof - prof_moy) / 10,        # +1 = 10 m plus profond
         log_poids = log10(poids_g),
         log_taille = log10(taille_mm),
         juvenile = taille_mm < 10 * taille_maturite_cm)

especes <- ind |> count(taxon) |> filter(n >= n_min_individus_espece) |> pull(taxon)
message("Espèces analysées (>= ", n_min_individus_espece, " individus) : ",
        paste(especes, collapse = ", "))
ind <- ind |> filter(taxon %in% especes)
ind <- condition_le_cren(ind, n_min = 20)$individus   # Kn (descriptif)

# Effectifs : individus et traits par espèce et saison
effectifs <- ind |>
  group_by(taxon, saison) |>
  summarise(n_individus = n(), n_traits = n_distinct(trait_id),
            part_juveniles = mean(juvenile), L_min = min(taille_mm),
            L_max = max(taille_mm), .groups = "drop") |>
  mutate(saison_testee = n_traits >= n_min_traits_saison)

# Relation taille-poids (toutes saisons) : W = a L^b
rel_tp <- bind_rows(lapply(especes, function(tx) {
  d <- ind[ind$taxon == tx, ]
  fit <- lm(log_poids ~ log_taille, d)
  ic <- confint(fit)[2, ]
  data.frame(taxon = tx, nom_commun = noms[[tx]], n = nrow(d),
             L_min_mm = min(d$taille_mm), L_max_mm = max(d$taille_mm),
             a = 10^coef(fit)[[1]], b = coef(fit)[[2]],
             b_IC95_bas = ic[[1]], b_IC95_haut = ic[[2]],
             r2 = summary(fit)$r.squared)
}))

# =============================================================================
# 2. Modèles de référence + 3. vérifications
# =============================================================================
pct <- function(b) 100 * (10^b - 1)   # coefficient sur log10 -> % de variation

formule <- function(reponse, avec_saison) {
  f <- if (reponse == "condition") "log_poids ~ log_taille_c" else "log_taille ~ 1"
  if (avec_saison) f <- paste(f, "+ saison")
  as.formula(paste(f, "+ prof10 + (1 | trait_id)"))
}
p_effet <- function(m, effet) {
  if (is.null(m)) return(NA_real_)
  t <- tryCatch(tests_lrt(m), error = function(e) NULL)
  if (is.null(t) || !effet %in% t$effet) return(NA_real_)
  t$p[t$effet == effet]
}

res <- list(tests = list(), saisons = list(), contrastes = list(),
            profondeur = list(), diagnostics = list())

for (tx in especes) {
  sais_ok <- effectifs$saison[effectifs$taxon == tx & effectifs$saison_testee]
  d <- ind |>
    filter(taxon == tx, saison %in% sais_ok) |>
    mutate(saison = droplevels(saison), trait_id = droplevels(trait_id),
           log_taille_c = log_taille - mean(log_taille))
  avec_saison <- nlevels(d$saison) >= 2
  effets <- c(if (avec_saison) "saison", "prof10")

  # Résidus de condition (relation taille-poids commune aux saisons testées)
  d$residu_condition <- resid(lm(log_poids ~ log_taille, d))

  for (reponse in c("condition", "taille")) {
    cle <- paste(tx, reponse, sep = " | ")
    message("--- ", cle, " (saisons : ", paste(levels(d$saison), collapse = ", "), ")")
    f <- formule(reponse, avec_saison)
    m <- ajuster_glmm(f, d, gaussian())
    if (is.null(m)) { warning("Échec du modèle : ", cle); next }
    lrt <- tests_lrt(m)

    # 3a. Niveau du trait : une valeur moyenne par trait
    y_ind <- if (reponse == "condition") "residu_condition" else "log_taille"
    pt <- d |>
      group_by(trait_id, saison) |>
      summarise(y = mean(.data[[y_ind]]), n = n(), prof10 = first(prof10),
                .groups = "drop") |>
      filter(n >= n_min_individus_trait) |>
      mutate(saison = droplevels(saison))
    f_trait <- as.formula(paste("y ~", paste(c(if (avec_saison && nlevels(pt$saison) > 1)
      "saison", "prof10"), collapse = " + ")))
    p_trait <- tryCatch({
      dt <- drop1(lm(f_trait, pt), test = "F")
      setNames(dt$`Pr(>F)`[-1], rownames(dt)[-1])
    }, error = function(e) c(saison = NA, prof10 = NA))

    # 3b. Retrait d'un trait à tour de rôle
    loo <- sapply(levels(d$trait_id), function(tr) {
      dd <- d[d$trait_id != tr, ] |> mutate(saison = droplevels(saison))
      if (avec_saison && nlevels(dd$saison) < 2) return(c(saison = NA, prof10 = NA))
      mm <- ajuster_glmm(f, dd, gaussian())
      sapply(effets, function(e) p_effet(mm, e))
    })
    loo <- matrix(loo, nrow = length(effets), dimnames = list(effets, NULL))

    # 3c. Sans les traits douteux
    dd <- d[!d$douteux, ] |> mutate(saison = droplevels(saison))
    m_sd <- if (avec_saison && nlevels(dd$saison) < 2) NULL else
      ajuster_glmm(f, dd, gaussian())

    for (e in effets) {
      res$tests[[paste(cle, e)]] <- data.frame(
        taxon = tx, reponse = reponse,
        effet = ifelse(e == "prof10", "profondeur", e),
        saisons_testees = paste(levels(d$saison), collapse = "/"),
        n_individus = nobs(m), n_traits = nlevels(d$trait_id),
        chi2 = lrt$chi2[lrt$effet == e], ddl = lrt$ddl[lrt$effet == e],
        p = lrt$p[lrt$effet == e],
        n_traits_moyennes = nrow(pt),
        p_niveau_trait = unname(p_trait[e]),
        p_max_retrait_un_trait = max(loo[e, ], na.rm = TRUE),
        trait_le_plus_influent = levels(d$trait_id)[which.max(loo[e, ])],
        p_sans_douteux = p_effet(m_sd, e))
    }

    # Effet de la profondeur : % pour +10 m
    co <- coefficients_modele(m)
    co <- co[co$terme == "prof10", ]
    res$profondeur[[cle]] <- data.frame(
      taxon = tx, reponse = reponse,
      effet_pct_10m = pct(co$estimation),
      IC95_bas = pct(co$estimation - 1.96 * co$erreur_type),
      IC95_haut = pct(co$estimation + 1.96 * co$erreur_type),
      gamme_profondeur_m = paste(round(range(d$prof)), collapse = "-"))

    # Saisons : écart à la moyenne de l'espèce (%), contrastes de Tukey (%)
    if (avec_saison) {
      em <- emmeans::emmeans(m, ~ saison)
      emd <- as.data.frame(em)
      ref <- mean(emd$emmean)
      res$saisons[[cle]] <- data.frame(
        taxon = tx, reponse = reponse, saison = emd$saison,
        ecart_pct = pct(emd$emmean - ref),
        IC95_bas = pct(emd[[grep("lower|LCL", names(emd))[1]]] - ref),
        IC95_haut = pct(emd[[grep("upper|UCL", names(emd))[1]]] - ref),
        valeur_mm = if (reponse == "taille") 10^emd$emmean else NA,
        n_individus = as.vector(table(d$saison)[as.character(emd$saison)]),
        n_traits = as.vector(tapply(d$trait_id, d$saison,
                                    function(x) length(unique(x)))[as.character(emd$saison)]))
      ctr <- as.data.frame(pairs(em, adjust = "tukey"))
      res$contrastes[[cle]] <- data.frame(
        taxon = tx, reponse = reponse, contraste = ctr$contrast,
        difference_pct = pct(ctr$estimate), p_tukey = ctr$p.value)
    }

    # Diagnostics DHARMa (autocorrélation spatiale testée par trait)
    dg <- diagnostics_dharma(
      m, d, file.path(dossier_diag, paste0(make.names(tx), "_", reponse, ".png")),
      groupe = "trait_id")
    res$diagnostics[[cle]] <- data.frame(
      taxon = tx, reponse = reponse, dg,
      variance_trait = unname(unlist(glmmTMB::VarCorr(m)$cond$trait_id)),
      variance_residuelle = sigma(m)^2,
      avertissements = paste(attr(m, "avertissements"), collapse = " | "))
  }
}

lier <- function(x) { x <- do.call(rbind, x); rownames(x) <- NULL; x }
tabs <- lapply(res, lier)

# Correction de Holm sur tous les tests saison / profondeur, puis verdict
tabs$tests <- tabs$tests |>
  mutate(p_holm = p.adjust(p, method = "holm"),
         verdict = case_when(
           p_holm >= 0.05 ~ "non significatif",
           p_niveau_trait < 0.05 & p_max_retrait_un_trait < 0.05 &
             p_sans_douteux < 0.05 ~ "robuste",
           TRUE ~ "fragile"),
         raison_fragile = case_when(
           verdict != "fragile" ~ "",
           TRUE ~ paste0(
             ifelse(is.na(p_niveau_trait) | p_niveau_trait >= 0.05,
                    "non significatif au niveau du trait ; ", ""),
             ifelse(p_max_retrait_un_trait >= 0.05,
                    paste0("dépend du trait ", trait_le_plus_influent, " ; "), ""),
             ifelse(is.na(p_sans_douteux) | p_sans_douteux >= 0.05,
                    "disparaît sans les traits douteux", "")))) |>
  relocate(p_holm, .after = p)

# =============================================================================
# 4. Croissance apparente des juvéniles (été -> automne)
# =============================================================================
# Décalage de la taille médiane des juvéniles (pondérée par l'élévation)
# divisé par le nombre de jours. Hypothèse : même cohorte, sans arrivée ni
# départ sélectif en taille entre les deux campagnes -> « croissance
# apparente ». IC 95 % : bootstrap des traits dans chaque saison.
mediane_juv <- function(dd) quantile_pond(dd$taille_mm, dd$facteur_elevation, 0.5)
croissance <- bind_rows(lapply(especes, function(tx) {
  ok <- effectifs |> filter(taxon == tx, saison %in% c("ete", "automne"),
                             saison_testee)
  if (nrow(ok) < 2) return(NULL)
  dj <- ind |> filter(taxon == tx, juvenile, saison %in% c("ete", "automne"))
  e <- dj[dj$saison == "ete", ]; a <- dj[dj$saison == "automne", ]
  jours <- as.numeric(mean(as.Date(a$date)) - mean(as.Date(e$date)))
  gain <- mediane_juv(a) - mediane_juv(e)
  boot <- replicate(n_bootstrap, {
    tirer <- function(x) {
      tr <- sample(unique(as.character(x$trait_id)), replace = TRUE)
      do.call(rbind, lapply(tr, function(t) x[x$trait_id == t, ]))
    }
    mediane_juv(tirer(a)) - mediane_juv(tirer(e))
  })
  data.frame(taxon = tx, nom_commun = noms[[tx]],
             n_juv_ete = nrow(e), n_juv_automne = nrow(a),
             L_med_ete_mm = mediane_juv(e), L_med_automne_mm = mediane_juv(a),
             jours = jours, gain_mm = gain,
             croissance_mm_jour = gain / jours,
             IC95_bas = quantile(boot, 0.025) / jours,
             IC95_haut = quantile(boot, 0.975) / jours)
}))

# =============================================================================
# Export
# =============================================================================
tabs$effectifs <- effectifs
tabs$relation_taille_poids <- rel_tp
tabs$croissance <- croissance
for (nom in names(tabs)) {
  write.csv(tabs[[nom]], file.path(dossier_out, paste0(nom, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}

# =============================================================================
# 5. Figures
# =============================================================================
lab_rep <- c(condition = "Condition (poids à taille égale)", taille = "Taille")
sauver <- function(p, nom, l, h) {
  print(p)
  ggsave(file.path(dossier_fig, nom), p, width = l, height = h, dpi = 150,
         bg = "white")
}
verdicts <- tabs$tests |> select(taxon, reponse, effet, verdict)

# 5a. Effet de la saison (% d'écart à la moyenne de l'espèce) + verdict
d_sais <- tabs$saisons |>
  mutate(saison = factor(as.character(saison), levels = saisons),
         espece = lab_esp(taxon))
etiquettes <- d_sais |>
  group_by(taxon, espece, reponse) |>
  summarise(y = max(IC95_haut), .groups = "drop") |>
  left_join(verdicts |> filter(effet == "saison"), by = c("taxon", "reponse"))
p_sais <- ggplot(d_sais, aes(saison, ecart_pct, colour = saison)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_pointrange(aes(ymin = IC95_bas, ymax = IC95_haut)) +
  geom_text(aes(y = IC95_haut, label = paste0(n_traits, " tr.")), vjust = -0.6,
            size = 2.6, colour = "grey30") +
  geom_text(data = etiquettes, aes(x = 0.6, y = Inf, label = verdict),
            inherit.aes = FALSE, hjust = 0, vjust = 1.4, size = 3,
            fontface = "italic", colour = "grey20") +
  facet_grid(reponse ~ espece, scales = "free_y",
             labeller = labeller(reponse = lab_rep)) +
  scale_x_discrete(drop = TRUE) +
  scale_colour_manual(values = couleurs_saison, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.15))) +
  labs(title = "Effet de la saison par espèce (modèle saison + profondeur)",
       subtitle = paste0("Écart à la moyenne de l'espèce (%), IC 95 % ; « tr. » = traits de chalut ; ",
                         "saisons avec < ", n_min_traits_saison, " traits non testées"),
       x = NULL, y = "Écart (%)") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
sauver(p_sais, "effet_saison_especes.png", 12, 6.5)

# 5b. Effet de la profondeur (% pour +10 m) + verdict
p_prof <- tabs$profondeur |>
  left_join(verdicts |> filter(effet == "profondeur"), by = c("taxon", "reponse")) |>
  mutate(espece = lab_esp(taxon),
         verdict = factor(verdict, levels = c("robuste", "fragile", "non significatif"))) |>
  ggplot(aes(effet_pct_10m, espece, colour = verdict)) +
  geom_vline(xintercept = 0, colour = "grey60") +
  geom_pointrange(aes(xmin = IC95_bas, xmax = IC95_haut)) +
  facet_wrap(~ reponse, scales = "free_x", labeller = as_labeller(lab_rep)) +
  scale_colour_manual(values = c(robuste = "#1b7837", fragile = "#e08214",
                                 `non significatif` = "grey55"), drop = FALSE) +
  labs(title = "Effet de la profondeur par espèce (modèle saison + profondeur)",
       subtitle = "Variation (%) pour 10 m de profondeur en plus, IC 95 %",
       x = "Variation pour +10 m (%)", y = NULL, colour = "Verdict") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
sauver(p_prof, "effet_profondeur_especes.png", 10, 5)

# 5c. Condition relative (Kn) par saison, saisons testées en plein
p_kn <- ind |>
  left_join(effectifs |> select(taxon, saison, saison_testee),
            by = c("taxon", "saison")) |>
  mutate(espece = lab_esp(taxon)) |>
  ggplot(aes(saison, Kn, colour = saison, alpha = saison_testee)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  geom_boxplot(outlier.shape = NA, width = 0.55) +
  geom_jitter(width = 0.15, size = 0.7) +
  facet_wrap(~ espece, nrow = 1) +
  scale_colour_manual(values = couleurs_saison, guide = "none") +
  scale_alpha_manual(values = c(`TRUE` = 0.6, `FALSE` = 0.15), guide = "none") +
  coord_cartesian(ylim = c(0.6, 1.4)) +
  labs(title = "Condition relative de Le Cren (Kn) par saison",
       subtitle = paste0("Kn = 1 : poids attendu pour la taille. Pâle : saison non testée (< ",
                         n_min_traits_saison, " traits)"),
       x = NULL, y = "Kn") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
sauver(p_kn, "condition_par_saison.png", 12, 4.5)

# 5d. Distributions de tailles par saison, taille de maturité
lm_lignes <- ind |> distinct(taxon, taille_maturite_cm) |>
  mutate(espece = lab_esp(taxon), L50 = 10 * taille_maturite_cm)
p_dist <- ind |>
  group_by(taxon, saison) |>
  mutate(poids_densite = facteur_elevation / sum(facteur_elevation)) |>
  ungroup() |>
  mutate(espece = lab_esp(taxon)) |>
  ggplot(aes(taille_mm, colour = saison, fill = saison)) +
  geom_density(aes(weight = poids_densite), alpha = 0.15, adjust = 0.8) +
  geom_vline(data = lm_lignes, aes(xintercept = L50), linetype = 2, colour = "grey30") +
  facet_wrap(~ espece, scales = "free", nrow = 1) +
  scale_colour_manual(values = couleurs_saison) +
  scale_fill_manual(values = couleurs_saison) +
  labs(title = "Distribution des tailles par saison (effectifs élevés)",
       subtitle = "Tirets : taille de première maturité (FishBase) ; à gauche = juvéniles",
       x = "Longueur (mm)", y = "Densité", colour = NULL, fill = NULL) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
sauver(p_dist, "distributions_tailles.png", 13, 4.5)

# 5e. Moyennes par trait (vérification) : condition et taille selon la profondeur
par_trait <- ind |>
  group_by(taxon, trait_id, saison) |>
  summarise(n = n(), prof = first(prof), Kn = mean(Kn),
            L_moy = 10^mean(log_taille), .groups = "drop") |>
  filter(n >= n_min_individus_trait) |>
  mutate(espece = lab_esp(taxon))
write.csv(par_trait, file.path(dossier_out, "moyennes_par_trait.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
p_pt <- par_trait |>
  tidyr::pivot_longer(c(Kn, L_moy), names_to = "variable", values_to = "valeur") |>
  mutate(variable = c(Kn = "Kn moyen", L_moy = "Taille moyenne (mm)")[variable]) |>
  ggplot(aes(prof, valeur, colour = saison, size = n)) +
  geom_point(alpha = 0.8) +
  facet_grid(variable ~ espece, scales = "free_y") +
  scale_colour_manual(values = couleurs_saison) +
  scale_size_area(max_size = 4) +
  labs(title = "Moyennes par trait de chalut selon la profondeur",
       subtitle = paste0("Traits avec au moins ", n_min_individus_trait,
                         " individus de l'espèce ; taille du point = individus"),
       x = "Profondeur (m)", y = NULL, colour = NULL, size = "n") +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
sauver(p_pt, "moyennes_par_trait_profondeur.png", 13, 6)

# =============================================================================
# 6. Résumé
# =============================================================================
message("\n================ RÉSUMÉ ================")
message("Effectifs (individus / traits) :")
print(effectifs |> mutate(across(where(is.numeric), ~ round(.x, 2))) |> as.data.frame())
message("\nTests et verdicts :")
print(tabs$tests |>
        select(taxon, reponse, effet, saisons_testees, n_traits, p_holm,
               p_niveau_trait, p_max_retrait_un_trait, p_sans_douteux, verdict) |>
        mutate(across(where(is.numeric), ~ signif(.x, 2))) |> as.data.frame())
message("\nCroissance apparente :")
print(croissance |> mutate(across(where(is.numeric), ~ round(.x, 2))))
message("\nDétails : ", dossier_out)
