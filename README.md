# Analyse des résultats de TP — données environnementales

## Structure

```
TP_environnement.Rproj
R/
  00_config.R                  <- TOUS les paramètres (zone, années, mois, profondeur, variables)
  01_extraction_copernicus.R   <- extraction Copernicus Marine + cartes de contrôle
  fonctions/
    copernicus_utils.R         <- fonctions (téléchargement, lecture NetCDF, graphiques)
data/
  raw/copernicus/              <- NetCDF bruts téléchargés
  processed/                   <- tableaux propres (.rds / .csv)
outputs/figures/copernicus/    <- cartes et résumé
```

## Installation (une seule fois)

1. **Compte Copernicus Marine** (gratuit) : https://data.marine.copernicus.eu/register
2. **Toolbox Copernicus Marine** (Python, ≥ 2.0) dans un terminal :
   ```bash
   pip install copernicusmarine
   ```
   ```bash
   copernicusmarine login
   ```
   (ou copier `.Renviron.example` en `.Renviron` et y mettre ses identifiants)
3. Packages R : `here`, `ncdf4`, `dplyr`, `ggplot2` (installés automatiquement par le script si besoin).

Si R ne trouve pas `copernicusmarine` (fréquent sous Windows / conda), renseigner son chemin
dans `chemin_copernicusmarine` (`R/00_config.R`). Pour le trouver : `where copernicusmarine`
(Windows) ou `which copernicusmarine` (Mac/Linux).

## Utilisation

1. Ouvrir `TP_environnement.Rproj` dans RStudio.
2. Modifier `R/00_config.R` : zone, `annees`, `mois`, profondeur, variables.
3. Lancer `R/01_extraction_copernicus.R` (bouton *Source*).

Les fichiers déjà téléchargés ne sont pas re-téléchargés (`ecraser_fichiers = FALSE`).
Les mois non consécutifs sont regroupés en périodes continues : on ne télécharge que les mois demandés.

## Variables disponibles

| nom (config)          | variable | produit Copernicus (réanalyse)        | résolution |
|-----------------------|----------|----------------------------------------|-----------|
| temperature           | thetao   | GLOBAL_MULTIYEAR_PHY_001_030           | 1/12°     |
| salinite              | so       | GLOBAL_MULTIYEAR_PHY_001_030           | 1/12°     |
| temperature_fond      | bottomT  | GLOBAL_MULTIYEAR_PHY_001_030           | 1/12°     |
| couche_melange        | mlotst   | GLOBAL_MULTIYEAR_PHY_001_030           | 1/12°     |
| bathymetrie           | deptho   | GLOBAL_MULTIYEAR_PHY_001_030 (static)  | 1/12°     |
| chlorophylle          | chl      | GLOBAL_MULTIYEAR_BGC_001_029           | 1/4°      |
| nitrate / phosphate / silicate | no3 / po4 / si | GLOBAL_MULTIYEAR_BGC_001_029 | 1/4° |
| oxygene               | o2       | GLOBAL_MULTIYEAR_BGC_001_029           | 1/4°      |
| production_primaire   | nppv     | GLOBAL_MULTIYEAR_BGC_001_029           | 1/4°      |

`source_copernicus = "analyse_prevision"` bascule sur les produits GLOBAL_ANALYSISFORECAST_*
(PHY_001_024 et BGC_001_028), qui couvrent mi-2022 / fin 2021 → aujourd'hui.

## Sortie

`data/processed/copernicus_mensuel.rds` : un tableau long avec les colonnes
`variable, var_nc, label, unite, lon, lat, profondeur, date, annee, mois, valeur`
(une ligne par pixel × mois × variable). C'est ce fichier que liront les scripts suivants.
