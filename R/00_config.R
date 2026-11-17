# =============================================================================
# 00_config.R — Paramètres du projet (À MODIFIER ICI, nulle part ailleurs)
# =============================================================================
# Ce fichier est sourcé par tous les autres scripts.
# -----------------------------------------------------------------------------

# --- Zone d'étude (degrés décimaux, WGS84) -----------------------------------
# /!\ À ADAPTER à la zone du TP. Les modèles biogéochimiques sont à 0.25°
#     (~25 km) : prendre une boîte d'au moins ~1° de côté pour avoir des pixels.
zone <- list(
  lon_min = -6.0,
  lon_max = -2.0,
  lat_min = 46.5,
  lat_max = 49.0
)

# --- Période : années × mois à extraire --------------------------------------
# Toutes les combinaisons annee × mois sont extraites.
# Ex. : mois <- 1:12 (année complète), mois <- c(3, 4, 5) (printemps), etc.
annees <- c(2025)
mois   <- c(3, 4, 5, 6)

# --- Profondeur --------------------------------------------------------------
# Plage de profondeur (m) téléchargée pour les variables 3D (température,
# salinité, chlorophylle, nutriments, oxygène).
# agregation :
#   "surface" -> garde uniquement le niveau le moins profond de la plage
#   "moyenne" -> moyenne verticale sur les niveaux de la plage
#   "aucune"  -> garde tous les niveaux (profils)
profondeur <- list(
  min        = 0,
  max        = 10,
  agregation = "surface"
)

# --- Source Copernicus --------------------------------------------------------
# "reanalyse"         : GLOBAL_MULTIYEAR_*  (1993 -> ~2 mois avant aujourd'hui)
#                       -> recommandé, cohérent dans le temps.
# "analyse_prevision" : GLOBAL_ANALYSISFORECAST_* (mi-2022 / fin 2021 -> présent)
#                       -> pour les mois les plus récents.
source_copernicus <- "reanalyse"

# --- Variables à extraire -----------------------------------------------------
# Noms disponibles : voir le catalogue `catalogue_variables` ci-dessous.
# Commenter / supprimer ce dont on n'a pas besoin.
variables_a_extraire <- c(
  "temperature",          # température potentielle (°C)
  "salinite",             # salinité
  "temperature_fond",     # température au fond (°C)
  "couche_melange",       # profondeur de la couche de mélange (m)
  "bathymetrie",          # profondeur du fond (m) — statique
  "chlorophylle",         # chlorophylle a (mg/m3)
  "nitrate",              # NO3 (mmol/m3)
  "phosphate",            # PO4 (mmol/m3)
  "silicate",             # Si (mmol/m3)
  "oxygene",              # O2 dissous (mmol/m3)
  "production_primaire"   # production primaire nette (mg C/m3/j)
)

# --- Options de téléchargement -----------------------------------------------
ecraser_fichiers <- FALSE   # TRUE = re-télécharge même si le .nc existe déjà
# Chemin vers l'exécutable `copernicusmarine` si R ne le trouve pas tout seul
# (ex. "C:/Users/moi/miniconda3/Scripts/copernicusmarine.exe")
chemin_copernicusmarine <- ""

# --- Dossiers ----------------------------------------------------------------
dossiers <- list(
  raw       = here::here("data", "raw", "copernicus"),
  processed = here::here("data", "processed"),
  figures   = here::here("outputs", "figures", "copernicus")
)

# =============================================================================
# Catalogue des variables (normalement pas besoin d'y toucher)
# Identifiants vérifiés sur le catalogue Copernicus Marine (toolbox v2).
# =============================================================================
catalogue_variables <- data.frame(
  nom = c("temperature", "salinite", "temperature_fond", "couche_melange",
          "bathymetrie", "chlorophylle", "nitrate", "phosphate", "silicate",
          "oxygene", "production_primaire"),

  # nom de la variable dans le NetCDF
  var_reanalyse = c("thetao", "so", "bottomT", "mlotst", "deptho",
                    "chl", "no3", "po4", "si", "o2", "nppv"),
  var_prevision = c("thetao", "so", "tob", "mlotst", "deptho",
                    "chl", "no3", "po4", "si", "o2", "nppv"),

  # dataset_id
  dataset_reanalyse = c(
    rep("cmems_mod_glo_phy_my_0.083deg_P1M-m", 4),
    "cmems_mod_glo_phy_my_0.083deg_static",
    rep("cmems_mod_glo_bgc_my_0.25deg_P1M-m", 6)
  ),
  dataset_prevision = c(
    "cmems_mod_glo_phy-thetao_anfc_0.083deg_P1M-m",
    "cmems_mod_glo_phy-so_anfc_0.083deg_P1M-m",
    "cmems_mod_glo_phy_anfc_0.083deg_P1M-m",
    "cmems_mod_glo_phy_anfc_0.083deg_P1M-m",
    "cmems_mod_glo_phy_anfc_0.083deg_static",
    "cmems_mod_glo_bgc-pft_anfc_0.25deg_P1M-m",
    rep("cmems_mod_glo_bgc-nut_anfc_0.25deg_P1M-m", 3),
    rep("cmems_mod_glo_bgc-bio_anfc_0.25deg_P1M-m", 2)
  ),

  statique = c(FALSE, FALSE, FALSE, FALSE, TRUE, rep(FALSE, 6)),

  label = c("Température (°C)", "Salinité", "Température au fond (°C)",
            "Couche de mélange (m)", "Bathymétrie (m)",
            "Chlorophylle a (mg/m³)", "Nitrate (mmol/m³)",
            "Phosphate (mmol/m³)", "Silicate (mmol/m³)",
            "Oxygène dissous (mmol/m³)", "Production primaire (mg C/m³/j)"),

  # palette viridis utilisée pour les cartes
  palette = c("plasma", "viridis", "plasma", "cividis", "mako",
              "viridis", "rocket", "rocket", "rocket", "turbo", "viridis"),

  stringsAsFactors = FALSE
)

# =============================================================================
# Préparation des données de terrain (script 02)
# =============================================================================
fichiers_terrain <- list(
  taille_poids = here::here("data", "raw", "terrain", "taille_poids_commun.csv"),
  metadonnees  = here::here("data", "raw", "terrain", "metadonnees_traits.csv"),
  taxons       = here::here("data", "ref", "taxons.csv")
)

# Ouverture horizontale du chalut (m). NA = inconnue -> CPUE par km chalutés.
# Si renseignée -> CPUE aussi pour 1 000 m² balayés.
ouverture_chalut_m <- NA

# Traits sans aucune ligne dans la base taille-poids que l'on CONFIRME vides
# (rien pêché). Ils entrent alors dans les analyses avec des zéros.
# Format "saison_station", ex. c("p_F", "p_M"). Les autres restent NA.
traits_vides_confirmes <- c()

# Traits signalés comme douteux (engin modifié, trait abrégé…)
traits_douteux <- c(
  p_A = "Chalut alourdi (manilles et chaînes, 5 kg)",
  p_O = "Ligne de sonde trop courte",
  e_A = "Fune allongée en cours de trait",
  e_N = "Chaîne gréée au chalut",
  e_G = "Fune trop longue",
  e_L = "Trait abrégé, chalut rempli de cténaires"
)
exclure_traits_douteux <- FALSE  # TRUE = retirés des tables d'analyse

# Seuil de détection des individus aberrants sur la relation taille-poids
# (résidus studentisés de log(poids) ~ log(taille), par taxon)
seuil_aberrant_taille_poids <- 4
