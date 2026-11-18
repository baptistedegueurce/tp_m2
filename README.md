# Analyse des résultats de TP — données environnementales

## Structure

```
TP_environnement.Rproj
R/
  00_config.R                  <- TOUS les paramètres (zone, années, mois, profondeur, variables)
  01_extraction_copernicus.R   <- extraction Copernicus Marine + cartes de contrôle
  02_preparation_donnees.R     <- nettoyage, effectifs élevés, effort, CPUE
  05_diversite_fonctionnelle.R <- espace fonctionnel, FDis, Rao, FRic, CWM, SES
  05b_modeles_fonctionnels.R   <- modèles mixtes des indices fonctionnels (glmmTMB)
  06_variables_environnementales.R <- table environnementale par trait (terrain + Copernicus)
  07_nourricerie.R             <- juvéniles, structure en taille, condition, croissance
  08_condition_tailles_especes.R <- condition et taille par espèce ~ saison, profondeur, environnement
  fonctions/
    copernicus_utils.R         <- fonctions (téléchargement, lecture NetCDF, graphiques)
    preparation_utils.R        <- fonctions de préparation des données de terrain
    fonctionnel_utils.R        <- indices fonctionnels (FDis, Rao, CWM, modèle nul)
    modeles_utils.R            <- modèles mixtes : sélection, tests, emmeans, DHARMa
    figures_modeles.R          <- figures communes (effet saison, profondeur, synthèse)
    nourricerie_utils.R        <- condition de Le Cren, spectre de taille, quantiles pondérés
data/
  raw/copernicus/              <- NetCDF bruts téléchargés
  raw/terrain/                 <- feuille "commun" (taille-poids) + fiches de trait saisies
  ref/taxons.csv               <- correspondance nom saisi -> taxon (à compléter)
  ref/traits_fonctionnels.csv  <- traits par taxon (FishBase / SeaLifeBase v24.07)
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

## Script 02 — préparation des données de terrain

Lit `data/raw/terrain/taille_poids_commun.csv` (feuille « commun »), `metadonnees.csv`
(fiches de trait recopiées) et `data/ref/taxons.csv`, puis écrit dans `data/processed/` :

| fichier | contenu |
|---|---|
| `traits.csv` | une ligne par trait : durée, distance, profondeur, position, statut |
| `captures.csv` | trait × taxon : effectif élevé, biomasse, CPUE par km (zéros inclus) |
| `individus.csv` | individus mesurés + facteur d'élévation + aberrants marqués |
| `taxons.csv` | taxons présents |
| `matrice_cpue.csv` | matrice traits × taxons (ind/km), prête pour `vegan` |
| `controle_anomalies.csv` | **tout ce qui a été corrigé ou supposé : à relire** |

Règles d'effectif : l'effectif noté = total de l'échantillon (mesurés + non mesurés) ;
s'il est absent, on compte les lignes ; s'il est recopié sur toutes les lignes, il compte une fois ;
plusieurs effectifs différents pour un même couple saison × station × espèce sont additionnés et signalés.
Un nouveau nom d'espèce dans la base doit d'abord être ajouté à `data/ref/taxons.csv`.

## Script 05 — diversité fonctionnelle

À lancer après le script 02. Lit `data/processed/preparation.rds` et `data/ref/traits_fonctionnels.csv`.
Traits, taxons gardés (poissons par défaut), pondération et nombre de permutations : `R/00_config.R`.

1. Distance de Gower entre taxons, puis PCoA ; nombre d'axes choisi par l'écart moyen (mAD < 0,05).
2. Par trait de chalut : richesse S, FDis, entropie de Rao, FRic / FEve / FDiv (package `FD`), CWM.
3. Modèle nul (traits permutés entre taxons) : SES de FDis et Rao, indépendants de la richesse.

Sorties : `data/processed/fonctionnel_indices.csv`, `fonctionnel_espace.csv`,
`fonctionnel_qualite_espace.csv` ; figures dans `outputs/figures/fonctionnel/`.

## Script 05b — modélisation des indices fonctionnels

À lancer après le script 05. Pour chaque indice (liste et loi en tête du script) :

| étape | méthode |
|---|---|
| loi | gaussienne (FDis, Rao, SES, CWM continus), bêta (FRic, FEve, FDiv, proportions), Poisson généralisée (richesse) |
| modèle de référence | `indice ~ saison + profondeur + (1 \| station)` (`glmmTMB`) |
| tests | rapport de vraisemblance (LRT) par effet, p corrigées de Holm sur tous les indices |
| sélection | AICc parmi 5 modèles (nul, saison, profondeur, saison + profondeur, interaction) |
| saisons | moyennes marginales et comparaisons deux à deux (`emmeans`, Tukey) |
| diagnostics | résidus simulés `DHARMa` : loi, dispersion, extrêmes, autocorrélation spatiale |
| sensibilité | même modèle sans les traits douteux |

Sorties : `data/processed/modeles_fonctionnels/` (`resume.csv` d'abord, puis `tests`, `selection`,
`coefficients`, `moyennes`, `contrastes`, `diagnostics`, `sensibilite`, `predictions`) ;
figures dans `outputs/figures/modeles_fonctionnels/` (+ un graphique DHARMa par indice dans `diagnostics/`).

## Script 07 — arrivée des juvéniles (contexte)

À lancer après le script 02. Un individu est juvénile si sa taille est inférieure à la taille
de première maturité (`taille_maturite_cm` dans `data/ref/traits_fonctionnels.csv`, FishBase).

Par trait de chalut : densité de juvéniles (ind/km), nombre de taxons présents en juvéniles,
part de juvéniles. Modèles `indicateur ~ saison + profondeur + (1 | station)` (même démarche
que 05b), sensibilité sans les traits douteux, et test à profondeur égale (Kruskal-Wallis sur
les traits de moins de `profondeur_max_cotiers` m), car au printemps seuls des traits côtiers
sont analysables. Tableau descriptif de la densité de juvéniles par espèce et par saison.

Retirés (pas assez fiables avec nos effectifs) : pente du spectre de taille (NBSS), condition
et taille moyennes toutes espèces. Condition et taille sont analysées par espèce dans le 08.

Sorties : `data/processed/nourricerie/`, figures dans `outputs/figures/nourricerie/`.

## Script 06 — variables environnementales par trait de chalut

Une ligne par trait : profondeur, fune, jour de l'année, distance à l'embouchure
(depuis `point_embouchure`, Saint-Nazaire par défaut) et, si le script 01 a été lancé,
les variables Copernicus du **mois du trait** au **pixel le plus proche**.
Exporte aussi les corrélations entre covariables et la part de chaque variable expliquée
par la saison. À 1/12° ou 1/4°, les variables Copernicus varient surtout entre campagnes,
très peu entre stations.

Exporte aussi `environnement_campagnes.csv` (moyenne, min, max de chaque variable par
campagne) : c'est ainsi que les variables Copernicus servent dans le rapport.

Ordre : `01` (Copernicus) → `02` → `06` → `07` → `08`.

## Script 08 — condition et tailles par espèce, avec vérification de robustesse

Pour chaque espèce avec au moins `n_min_individus_espece` individus mesurés et pesés :

- condition : `log10(poids) ~ log10(taille) + saison + profondeur + (1 | trait)` (poids à taille égale)
- taille : `log10(taille) ~ saison + profondeur + (1 | trait)`

Une saison n'entre dans les tests d'une espèce que si l'espèce y a été prise dans au moins
`n_min_traits_saison` traits (3). Profondeur par tranche de 10 m : effets en % pour +10 m.
LRT, correction de Holm sur tous les tests, Tukey entre saisons, DHARMa.

Chaque effet est ensuite vérifié de trois façons : test sur les moyennes par trait (traits
d'au moins `n_min_individus_trait` individus), retrait d'un trait à tour de rôle, et sans les
traits douteux. Verdict dans `tests.csv` : « robuste », « fragile » ou « non significatif ».
Croissance apparente été → automne des juvéniles (IC 95 % par bootstrap des traits).

Les variables Copernicus ne sont pas des covariables : à cette échelle, elles sont confondues
avec la saison (voir `00_config.R`).

Sorties : `data/processed/condition_especes/` (`tests`, `saisons`, `contrastes`, `profondeur`,
`effectifs`, `relation_taille_poids`, `croissance`, `diagnostics`, `moyennes_par_trait`) ;
figures dans `outputs/figures/condition_especes/`. Durée : environ 6 minutes.
