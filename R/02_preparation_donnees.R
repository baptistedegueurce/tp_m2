# =============================================================================
# 02_preparation_donnees.R
# Préparation et standardisation des données de chalutage
#
# Entrées (voir R/00_config.R) :
#   data/raw/terrain/taille_poids_commun.csv   feuille "commun" (mensurations)
#   data/raw/terrain/metadonnees.csv           fiches de trait
#   data/ref/taxons.csv                        nom saisi -> taxon
#
# Sorties (data/processed/) :
#   traits.csv       une ligne par trait : effort, position, statut
#   captures.csv     trait × taxon : effectifs élevés, biomasse, CPUE (zéros inclus)
#   individus.csv    une ligne par individu mesuré, avec facteur d'élévation
#   taxons.csv       table des taxons présents
#   matrice_cpue.csv matrice traits × taxons (individus / km) pour vegan & co
#   preparation.rds  les tables ci-dessus dans une liste
#   controle_anomalies.csv  TOUT ce qui a été corrigé / supposé : à relire !
# Figures (outputs/figures/preparation/) : carte des traits, taille-poids,
#   effectifs par taxon.
#
# Règles d'effectif (feuille "commun") :
#   - l'effectif noté = nombre TOTAL de l'échantillon (mesurés + non mesurés),
#     ex. 227 gobies = 20 mesurés + 207 non mesurés ;
#   - effectif absent -> nombre de lignes du couple saison × station × espèce ;
#   - effectif recopié sur toutes les lignes -> compté une seule fois ;
#   - plusieurs effectifs différents pour un même couple -> blocs (sachets ?)
#     additionnés, et signalés dans controle_anomalies.csv.
# =============================================================================

pkgs <- c("here", "dplyr", "tidyr", "ggplot2")
manquants <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(manquants)) install.packages(manquants)

library(dplyr)

source(here::here("R", "00_config.R"), encoding = "UTF-8")
source(here::here("R", "fonctions", "preparation_utils.R"), encoding = "UTF-8")

dossier_fig <- here::here("outputs", "figures", "preparation")
dir.create(dossier_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(dossiers$processed, recursive = TRUE, showWarnings = FALSE)

anomalies <- creer_registre_anomalies()
signaler  <- anomalies$signaler

saisons <- c(p = "printemps", e = "ete", a = "automne")

# =============================================================================
# 1. Métadonnées des traits : effort de pêche
# =============================================================================
meta <- read.csv(fichiers_terrain$metadonnees, colClasses = "character",
                 encoding = "UTF-8", fileEncoding = "UTF-8",
                 check.names = FALSE, na.strings = c("", "NA"))

# Les colonnes sont lues dans l'ordre du fichier (les intitulés avec accents
# et unités varient d'un export à l'autre).
colonnes_meta <- c("saison", "station", "date", "vitesse_noeuds",
                   "filage_heure", "filage_prof_m", "filage_lat_N",
                   "filage_lon_O", "virage_heure", "virage_prof_m",
                   "virage_lat_N", "virage_lon_O", "sous_echantillonnage",
                   "nb_sachets", "commentaire")
if (ncol(meta) != length(colonnes_meta)) {
  stop("Le fichier de métadonnées doit avoir ", length(colonnes_meta),
       " colonnes dans cet ordre :\n  ", paste(colonnes_meta, collapse = ", "),
       "\n(", ncol(meta), " trouvées)")
}
names(meta) <- colonnes_meta

meta <- meta |>
  mutate(
    saison  = lire_saison(saison),
    station = toupper(trimws(station)),
    filage_heure = lire_heure(filage_heure),
    virage_heure = lire_heure(virage_heure),
    nb_sachets   = as.character(lire_nombre(nb_sachets)),   # "—" -> NA
    fune_m  = lire_fune(commentaire)
  )
if (anyNA(meta$saison)) stop("Saison non reconnue dans les métadonnées")
if (any(duplicated(paste(meta$saison, meta$station)))) {
  stop("Station en double dans les métadonnées")
}

traits <- meta |>
  mutate(
    trait_id   = paste(saison, station, sep = "_"),
    saison     = factor(saisons[saison], levels = saisons),
    date       = lire_date(date),
    across(c(vitesse_noeuds, filage_prof_m, filage_lat_N, filage_lon_O,
             virage_prof_m, virage_lat_N, virage_lon_O, fune_m), lire_nombre),
    nb_sachets = as.integer(nb_sachets),
    filage_dt  = as.POSIXct(paste(date, filage_heure), tz = "Europe/Paris"),
    virage_dt  = as.POSIXct(paste(date, virage_heure), tz = "Europe/Paris"),
    duree_min  = as.numeric(difftime(virage_dt, filage_dt, units = "mins")),
    # longitudes Ouest -> négatives (convention WGS84)
    filage_lon = -filage_lon_O,
    virage_lon = -virage_lon_O,
    lat        = (filage_lat_N + virage_lat_N) / 2,
    lon        = (filage_lon + virage_lon) / 2,
    profondeur_m = (filage_prof_m + virage_prof_m) / 2,
    distance_gps_km     = haversine_km(filage_lat_N, filage_lon,
                                       virage_lat_N, virage_lon),
    distance_vitesse_km = vitesse_noeuds * 1.852 * duree_min / 60,
    fraction_sous_ech   = lire_fraction(coalesce(sous_echantillonnage, "")),
    douteux        = trait_id %in% names(traits_douteux),
    raison_douteux = unname(traits_douteux[trait_id])
  )

# Contrôles de l'effort
for (i in which(traits$duree_min <= 0 | traits$duree_min > 60)) {
  signaler("duree_anormale", sprintf("Durée de %.0f min", traits$duree_min[i]),
           traits$saison[i], traits$station[i], action = "À vérifier sur la fiche")
}
# Deux traits qui se chevauchent dans le temps = erreur de date ou d'heure
chev <- traits |>
  select(trait_id, saison, station, filage_dt, virage_dt) |>
  merge(select(traits, autre = trait_id, f2 = filage_dt, v2 = virage_dt),
        by = NULL) |>
  filter(trait_id < autre, filage_dt < v2, f2 < virage_dt)
for (i in seq_len(nrow(chev))) {
  signaler("traits_simultanes",
           sprintf("%s et %s se chevauchent (%s)", chev$trait_id[i],
                   chev$autre[i], format(chev$filage_dt[i], "%d/%m %H:%M")),
           chev$saison[i], chev$station[i],
           action = "Erreur de date ou d'heure sur la fiche : à corriger")
}

ratio <- traits$distance_gps_km / traits$distance_vitesse_km
for (i in which(ratio < 0.6 | ratio > 1.4)) {
  signaler("distance_incoherente",
           sprintf("GPS : %.2f km ; vitesse × durée : %.2f km",
                   traits$distance_gps_km[i], traits$distance_vitesse_km[i]),
           traits$saison[i], traits$station[i],
           action = "Distance GPS conservée ; vérifier les coordonnées")
}
for (i in which(!is.na(traits$commentaire) & grepl("corrig", traits$commentaire))) {
  signaler("correction_fiche", traits$commentaire[i],
           traits$saison[i], traits$station[i], action = "Valeur corrigée")
}
for (i in which(traits$fraction_sous_ech < 1)) {
  signaler("sous_echantillonnage",
           sprintf("Fraction %s : effectifs multipliés par %g",
                   traits$sous_echantillonnage[i], 1 / traits$fraction_sous_ech[i]),
           traits$saison[i], traits$station[i],
           action = if (grepl("sachet", traits$sous_echantillonnage[i]))
             "Appliqué à tout le trait alors qu'un seul sachet est concerné : à vérifier"
           else "Hypothèse : les effectifs saisis sont ceux du sous-échantillon")
}

# Distance retenue : GPS (à défaut, vitesse × durée)
traits$distance_km <- coalesce(traits$distance_gps_km, traits$distance_vitesse_km)
traits$surface_m2  <- traits$distance_km * 1000 * ouverture_chalut_m

# =============================================================================
# 2. Base taille-poids (feuille "commun")
# =============================================================================
tp_brut <- read.csv(fichiers_terrain$taille_poids, colClasses = "character",
                    encoding = "UTF-8", check.names = FALSE)

tp <- tp_brut |>
  mutate(
    ligne_saisie = row_number(),                 # l'ordre compte (blocs)
    date      = lire_date(date),
    code_saison = trimws(tolower(saison)),
    station   = toupper(trimws(station)),
    nom_saisi = tolower(trimws(espece)),
    effectif  = as.integer(lire_nombre(effectif)),
    taille_mm = lire_nombre(taille),
    poids_g   = lire_nombre(poids),
    trait_id  = paste(code_saison, station, sep = "_")
  ) |>
  select(-saison, -espece, -taille, -poids)

if (any(!tp$code_saison %in% names(saisons))) {
  stop("Code saison inconnu : ",
       paste(unique(tp$code_saison[!tp$code_saison %in% names(saisons)]),
             collapse = ", "))
}
inconnus <- setdiff(unique(tp$trait_id), traits$trait_id)
if (length(inconnus)) {
  stop("Stations absentes des métadonnées : ", paste(inconnus, collapse = ", "))
}

# Dates : cohérence avec la fiche de trait
dates_tp <- tp |> distinct(trait_id, date) |>
  left_join(select(traits, trait_id, date_fiche = date), by = "trait_id") |>
  filter(date != date_fiche)
for (i in seq_len(nrow(dates_tp))) {
  signaler("date_differente",
           sprintf("%s dans la base, %s sur la fiche",
                   dates_tp$date[i], dates_tp$date_fiche[i]),
           taxon = NA, station = dates_tp$trait_id[i],
           action = "Date de la fiche de trait conservée")
}

# Taxonomie
taxons_ref <- read.csv(fichiers_terrain$taxons, encoding = "UTF-8",
                       stringsAsFactors = FALSE)
taxons_ref$nom_saisi <- tolower(trimws(taxons_ref$nom_saisi))
non_reconnus <- setdiff(unique(tp$nom_saisi), taxons_ref$nom_saisi)
if (length(non_reconnus)) {
  stop("Noms d'espèces absents de data/ref/taxons.csv : ",
       paste(non_reconnus, collapse = ", "),
       "\n-> les ajouter dans ce fichier puis relancer.")
}
tp <- left_join(tp, select(taxons_ref, nom_saisi, taxon, groupe),
                by = "nom_saisi")

# =============================================================================
# 3. Effectifs : blocs et facteurs d'élévation
# =============================================================================
tp <- tp |>
  arrange(ligne_saisie) |>
  group_by(trait_id, nom_saisi) |>
  mutate(bloc = numeroter_blocs(effectif)) |>
  ungroup()

blocs <- tp |>
  group_by(trait_id, code_saison, station, nom_saisi, taxon, groupe, bloc) |>
  summarise(
    n_lignes   = n(),
    n_mesures  = sum(!is.na(taille_mm) & !is.na(poids_g)),
    effectif_note = if (all(is.na(effectif))) NA_integer_
                    else effectif[!is.na(effectif)][1],
    n_valeurs_effectif = sum(!is.na(effectif)),
    poids_mesures_g = sum(poids_g[!is.na(taille_mm) & !is.na(poids_g)]),
    .groups = "drop"
  ) |>
  mutate(effectif_bloc = coalesce(effectif_note, n_lignes))

# Contrôles et corrections
for (i in seq_len(nrow(blocs))) {
  b <- blocs[i, ]
  qui <- list(saison = saisons[b$code_saison], station = b$station,
              taxon = b$nom_saisi)
  if (!is.na(b$effectif_note) && b$effectif_note < b$n_lignes) {
    do.call(signaler, c(list(
      "effectif_inferieur_aux_lignes",
      sprintf("Effectif noté %d mais %d lignes saisies (bloc %d)",
              b$effectif_note, b$n_lignes, b$bloc),
      action = "Nombre de lignes retenu"), qui))
    blocs$effectif_bloc[i] <- b$n_lignes
  }
  if (is.na(b$effectif_note) && b$n_mesures >= 20) {
    do.call(signaler, c(list(
      "effectif_absent_20_mesures",
      sprintf("%d individus mesurés, pas d'effectif total noté", b$n_mesures),
      action = "Effectif = nombre de lignes ; vérifier s'il y avait plus de poissons"),
      qui))
  }
}
multi <- blocs |> count(trait_id, code_saison, station, nom_saisi) |> filter(n > 1)
for (i in seq_len(nrow(multi))) {
  d <- blocs |> filter(trait_id == multi$trait_id[i], nom_saisi == multi$nom_saisi[i])
  signaler("plusieurs_effectifs",
           paste(sprintf("bloc %d : effectif %s, %d lignes", d$bloc,
                         ifelse(is.na(d$effectif_note), "absent",
                                d$effectif_note), d$n_lignes),
                 collapse = " | "),
           saisons[multi$code_saison[i]], multi$station[i], multi$nom_saisi[i],
           action = sprintf("Blocs additionnés (total %d) : sachets différents ?",
                            sum(d$effectif_bloc)))
}

# Facteurs : élévation au bloc × sous-échantillonnage du trait
blocs <- blocs |>
  left_join(select(traits, trait_id, fraction_sous_ech), by = "trait_id") |>
  mutate(
    facteur_bloc   = ifelse(n_mesures > 0, effectif_bloc / n_mesures, NA_real_),
    effectif_eleve = effectif_bloc / fraction_sous_ech,
    biomasse_g     = ifelse(n_mesures > 0,
                            poids_mesures_g * facteur_bloc / fraction_sous_ech,
                            NA_real_)
  )

# =============================================================================
# 4. Table individus (+ détection d'aberrants taille-poids)
# =============================================================================
individus <- tp |>
  filter(!is.na(taille_mm) | !is.na(poids_g)) |>
  left_join(select(blocs, trait_id, nom_saisi, bloc, facteur_bloc,
                   fraction_sous_ech),
            by = c("trait_id", "nom_saisi", "bloc")) |>
  left_join(select(traits, trait_id, saison, date_trait = date), by = "trait_id") |>
  mutate(facteur_elevation = facteur_bloc / fraction_sous_ech,
         aberrant_taille_poids = FALSE,
         residu_taille_poids   = NA_real_)

for (tx in unique(individus$taxon)) {
  idx <- which(individus$taxon == tx & individus$taille_mm > 0 &
                 individus$poids_g > 0)
  if (length(idx) < 10) next
  fit <- lm(log(poids_g) ~ log(taille_mm), data = individus[idx, ])
  r <- rstudent(fit)
  individus$residu_taille_poids[idx] <- r
  out <- idx[abs(r) > seuil_aberrant_taille_poids]
  individus$aberrant_taille_poids[out] <- TRUE
  for (j in out) {
    signaler("aberrant_taille_poids",
             sprintf("%g mm pour %g g (résidu studentisé %.1f)",
                     individus$taille_mm[j], individus$poids_g[j],
                     r[match(j, idx)]),
             individus$saison[j], individus$station[j], tx,
             action = "Conservé mais marqué : vérifier la saisie")
  }
}
for (i in which(individus$taille_mm <= 0 | individus$poids_g <= 0)) {
  signaler("valeur_non_positive", "Taille ou poids <= 0",
           individus$saison[i], individus$station[i], individus$taxon[i],
           action = "Marqué aberrant")
  individus$aberrant_taille_poids[i] <- TRUE
}

individus <- individus |>
  transmute(trait_id, saison, station, date = date_trait, taxon, groupe,
            nom_saisi, taille_mm, poids_g, bloc, facteur_elevation,
            aberrant_taille_poids, residu_taille_poids, ligne_saisie)

# =============================================================================
# 5. Statut des traits
# =============================================================================
traits_avec_donnees <- unique(tp$trait_id)
traits <- traits |>
  mutate(statut = case_when(
    trait_id %in% traits_avec_donnees    ~ "echantillonne",
    trait_id %in% traits_vides_confirmes ~ "vide_confirme",
    is.na(nb_sachets)                    ~ "pas_d_echantillon",
    TRUE                                 ~ "donnees_indisponibles"
  ))
for (i in which(traits$statut %in% c("pas_d_echantillon", "donnees_indisponibles"))) {
  signaler(traits$statut[i],
           if (traits$statut[i] == "donnees_indisponibles")
             sprintf("%d sachet(s) sur la fiche, données non disponibles",
                     traits$nb_sachets[i])
           else "Aucun sachet (« – » sur la fiche) et aucune ligne dans la base",
           traits$saison[i], traits$station[i],
           action = "Trait exclu (NA). Si vide confirmé : traits_vides_confirmes")
}
traits <- traits |>
  mutate(analysable = statut %in% c("echantillonne", "vide_confirme") &
           !(douteux & exclure_traits_douteux))

# =============================================================================
# 6. Captures par trait × taxon (avec zéros) et CPUE
# =============================================================================
captures <- blocs |>
  group_by(trait_id, taxon, groupe) |>
  summarise(
    effectif_saisi = sum(effectif_bloc),
    effectif       = sum(effectif_eleve),
    n_mesures      = sum(n_mesures),
    # biomasse inconnue si un bloc n'a aucun individu pesé
    biomasse_g     = if (any(is.na(biomasse_g))) NA_real_ else sum(biomasse_g),
    .groups = "drop"
  )

traits_analyse <- traits |> filter(analysable)
groupes_taxons <- distinct(captures, taxon, groupe)

captures <- captures |>
  filter(trait_id %in% traits_analyse$trait_id) |>
  tidyr::complete(trait_id = traits_analyse$trait_id,
                  tidyr::nesting(taxon, groupe),
                  fill = list(effectif_saisi = 0, effectif = 0, n_mesures = 0,
                              biomasse_g = 0)) |>
  left_join(select(traits, trait_id, saison, station, date, distance_km,
                   surface_m2, profondeur_m, lat, lon, douteux),
            by = "trait_id") |>
  mutate(
    cpue_n_km    = effectif / distance_km,
    cpue_g_km    = biomasse_g / distance_km,
    dens_n_1000m2 = effectif / surface_m2 * 1000,
    dens_g_1000m2 = biomasse_g / surface_m2 * 1000
  ) |>
  arrange(saison, station, desc(effectif))

if (all(is.na(ouverture_chalut_m))) {
  captures <- select(captures, -dens_n_1000m2, -dens_g_1000m2)
}

# Matrice communautaire (traits en lignes, taxons en colonnes)
matrice_cpue <- captures |>
  select(trait_id, taxon, cpue_n_km) |>
  tidyr::pivot_wider(names_from = taxon, values_from = cpue_n_km,
                     values_fill = 0) |>
  as.data.frame()
rownames(matrice_cpue) <- matrice_cpue$trait_id

# Table des taxons présents
taxons <- taxons_ref |>
  filter(taxon %in% captures$taxon) |>
  group_by(taxon, rang, groupe) |>
  summarise(noms_saisis = paste(sort(unique(nom_saisi)), collapse = " / "),
            nom_commun  = first(nom_commun),
            remarque    = paste(unique(remarque[nzchar(remarque)]), collapse = " ; "),
            .groups = "drop") |>
  left_join(captures |> group_by(taxon) |>
              summarise(effectif_total = sum(effectif),
                        n_traits_presents = sum(effectif > 0),
                        n_mesures = sum(n_mesures)),
            by = "taxon") |>
  arrange(desc(effectif_total))

# =============================================================================
# 7. Export
# =============================================================================
traits_export <- traits |>
  select(trait_id, saison, station, date, filage_heure, virage_heure, duree_min,
         vitesse_noeuds, fune_m, filage_lat_N, filage_lon, virage_lat_N,
         virage_lon, lat, lon, filage_prof_m, virage_prof_m, profondeur_m,
         distance_gps_km, distance_vitesse_km, distance_km, surface_m2,
         nb_sachets, fraction_sous_ech, statut, douteux, raison_douteux,
         analysable, commentaire)

controle <- anomalies$tableau()

ecrire <- function(x, nom) {
  write.csv(x, file.path(dossiers$processed, nom), row.names = FALSE,
            fileEncoding = "UTF-8")
}
ecrire(traits_export, "traits.csv")
ecrire(captures,      "captures.csv")
ecrire(individus,     "individus.csv")
ecrire(taxons,        "taxons.csv")
ecrire(matrice_cpue,  "matrice_cpue.csv")
ecrire(controle,      "controle_anomalies.csv")
saveRDS(list(traits = traits_export, captures = captures, individus = individus,
             taxons = taxons, matrice_cpue = matrice_cpue, controle = controle),
        file.path(dossiers$processed, "preparation.rds"))

# =============================================================================
# 8. Figures de contrôle
# =============================================================================
library(ggplot2)

# 8a. Carte des traits (segment filage -> virage), couleur = statut
p_carte <- ggplot(traits) +
  geom_segment(aes(x = filage_lon, y = filage_lat_N, xend = virage_lon,
                   yend = virage_lat_N, colour = statut),
               linewidth = 1.1,
               arrow = arrow(length = unit(0.12, "cm"), type = "closed")) +
  geom_text(aes(lon, lat, label = station), size = 3, nudge_y = 0.004) +
  facet_wrap(~ saison) +
  scale_colour_manual(values = c(echantillonne = "#1b7837",
                                 vide_confirme = "#5aae61",
                                 pas_d_echantillon = "grey60",
                                 donnees_indisponibles = "#d73027")) +
  coord_quickmap() +
  labs(title = "Traits de chalut par saison", x = "Longitude (°)",
       y = "Latitude (°)", colour = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
print(p_carte)
ggsave(file.path(dossier_fig, "carte_traits.png"), p_carte,
       width = 11, height = 4.5, dpi = 150, bg = "white")

# 8b. Relation taille-poids (échelles log), aberrants en rouge
taxons_mesures <- individus |> count(taxon) |> filter(n >= 10) |> pull(taxon)
p_tp <- individus |>
  filter(taxon %in% taxons_mesures, taille_mm > 0, poids_g > 0) |>
  ggplot(aes(taille_mm, poids_g)) +
  geom_point(aes(colour = saison), alpha = 0.5, size = 1) +
  geom_point(data = ~ filter(.x, aberrant_taille_poids), colour = "red",
             shape = 21, size = 3, stroke = 1) +
  scale_x_log10() + scale_y_log10() +
  facet_wrap(~ taxon, scales = "free") +
  labs(title = "Relation taille-poids (aberrants entourés en rouge)",
       x = "Taille (mm, log)", y = "Poids (g, log)", colour = NULL) +
  theme_minimal(base_size = 10)
print(p_tp)
ggsave(file.path(dossier_fig, "taille_poids.png"), p_tp,
       width = 10, height = 7, dpi = 150, bg = "white")

# 8c. Effectifs élevés par taxon et saison
p_eff <- captures |>
  group_by(saison, taxon) |>
  summarise(effectif = sum(effectif), .groups = "drop") |>
  filter(effectif > 0) |>
  ggplot(aes(effectif, reorder(taxon, effectif, FUN = sum), colour = saison)) +
  geom_point(size = 2.5, alpha = 0.85) +
  scale_x_log10() +
  labs(title = "Effectifs élevés par taxon (somme des traits analysables)",
       x = "Effectif (échelle log)", y = NULL, colour = NULL) +
  theme_minimal(base_size = 10)
print(p_eff)
ggsave(file.path(dossier_fig, "effectifs_taxons.png"), p_eff,
       width = 8, height = 6, dpi = 150, bg = "white")

# =============================================================================
# 9. Résumé
# =============================================================================
message("\n================ RÉSUMÉ ================")
message(sprintf("Traits : %d au total, %d analysables", nrow(traits),
                sum(traits$analysable)))
print(table(traits$saison, traits$statut))
message(sprintf("Taxons : %d | individus mesurés : %d (dont %d aberrants marqués)",
                nrow(taxons), nrow(individus),
                sum(individus$aberrant_taille_poids)))
message(sprintf("Anomalies / hypothèses à relire : %d -> %s", nrow(controle),
                file.path(dossiers$processed, "controle_anomalies.csv")))
print(table(controle$type))
