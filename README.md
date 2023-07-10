# Plant size, latitude, and phylogeny explain within-population variability in herbivory
---

This dataset or analysis was produced by the members of The Herbivory Variability Network. It examines global and phylogenetic patterns of herbivory via 790 suveys of 503 plant species. 

Our website is [herbvar.org](https://herbvar.org). We are funded by the US National Science Foundation Research Coordination Networks program. The Principal Investigators (PIs) of the project are William Wetzel (william.wetzel@montana.edu), Phil Hahn, Brian Inouye, Nora Underwood, and Susan Whitehead. The steering committee is the PIs and Karen Abbott, Emilio Bruna, N. Ivalú Cacho, and Lee Dyer. Our contact information is [here](https://herbvar.org/leadership.html). For the full list of HerbVar members, please see [this webpage](https://herbvar.org/CollaboratorDirectory.html).

This version of the dataset/analysis is associated with our paper "Plant size, latitude, and phylogeny explain within-population variability in herbivory."


## Description of the data and file structure

There are three csv data files in this dataset. Each represents a different biological scale: plant individuals, plant populations, and plant species.

data_plant_level_prep_2023-01-01 10-15-19.csv -- each row is a plant individual
data_survey_level_prep_2023-01-01 10-15-19.csv -- each row is a plant population
data_species_level_prep_2023-01-01 10-15-19.csv -- each row is a plant species

## Column name definitions

**data_plant_level_prep_2023-01-01 10-15-19.csv**

| Col name | Definition |
| --- | --- |
| surveyID | ID for the population survey |
| Genus_sp | Genus and species of the plant surveyed |
| plantFamily | Taxonomic family of the plant species |
| percHerbPlant | The percent of leaf area damaged by herbivores |
| Lat | The latitude of the survey site |

**data_survey_level_prep_2023-01-01 10-15-19.csv**

| Col name | Definition |
| --- | --- |
| surveyID | ID for the population survey |
| date | The date of the survey |
| Lat | The latitude of the survey site |
| Lat_abs | The absolute value of the latitude of the survey site |
| Lat_abs_scale | The scaled absolute value of the latitude of the survey site |
| Long | The longitude of the survey site |
| hemi | The hemisphere of the survey site |
| Genus_sp | Genus and species of the plant surveyed |
| plantFamily | Taxonomic family of the plant species |
| Biome | The biome of the survey site |
| growthForm | The growth form of the plant species |
| growthForm_simp | A simplified version of the growth form of the plant species |
| sizeDiameterMean | The mean height for most plant species or mean diameter for prostrate species |
| sizeDiameterMeanLog | Log of sizeDiameterMean |
| sizeDiameterMeanLogScale | Scaled log of sizeDiameterMean |
| plantMean | The mean percent herbivory in the survey |
| plantMeanProp | The mean proportion herbivory in the survey |
| plantMeanPropLogit | The logit mean proportion herbivory in the survey |
| plantGini | The Gini coefficient of herbivory in the survey |
| plantGiniAdj | The adjusted Gini coefficient of herbivory in the survey (see methods) |
| plantGiniAdjLogit | The logit adjusted Gini coefficient of herbivory in the survey (see methods) |
| propPlantsHerb50 | The minimum proportion of plant individuals with >= 50% of the population's proportion herbivory |
| focalPlantCoverMean | The average cover of the focal plant species within the survey |

**data_species_level_prep_2023-01-01 10-15-19.csv**

| Col name | Definition |
| --- | --- |
| Genus_sp | Genus and species of the plant surveyed |
| plantFamily | Taxonomic family of the plant species |
| plantMean | The mean percent herbivory averaged across the surveys of each species |
| plantMeanProp | The mean proportion herbivory averaged across the surveys of each species |
| plantMeanPropLogit | The logit mean proportion herbivory averaged across the surveys of each species |
| plantGini | The Gini coefficient of herbivory averaged across the surveys of each species |
| plantGiniAdj | The adjusted Gini coefficient of herbivory averaged across the surveys of each species (see methods) |
| plantGiniAdjLogit | The logit adjusted Gini coefficient of herbivory averaged across the surveys of each species (see methods) |


## Code/Software

Our analysis scripts are available at 
[https://github.com/HerbVar-Network/HV-Large-Patterns-MS-public](https://github.com/HerbVar-Network/HV-Large-Patterns-MS-public)
