# DynaChicks_LargeGroup_Movement-Dominance
### Analysis of social behaviour and movements of laying hens within an aviary

Data and scripts for the paper: Exploring implications of laying hens’ tendencies for aggression in aviaries


## Data files:  
Only columns are listed which were used in the analysis.  

1. **filtered_30s_full_combined.csv**  
   Output from raw tracking after validated filtering technique (see Candelotto et al., 2022).  

   **Columns:**  
   |----------|-------------|
   | **Time**  | Time and date stamp of transition between zones of the aviary  |
   | **Zone**  | New zone after transition (Wintergarten, Litter, Tier_2, Ramp_Nestbox, Tier_4)  |
   | **Hen**  | Unique identifier of animal ("Hen_00")  |
   | **PackID**  | Pen ID in addition to backpack colour (non-unique identifier)  |

2. **socialData.csv**  
   Observed social interactions of focal individuals.  

   **Columns:**  
   |----------|-------------|
   | **Observer**  | Who performed the observation (3 distinct observers)  |
   | **Reliability**  | Boolean to indicate whether to use the data for the reliability analysis  |
   | **Pen**  | Animal group out of 3, 4, 5, 10, 11, 12  |
   | **ID**  | Pen ID in addition to backpack colour (non-unique identifier)  |
   | **Count**  | Indicator to signal repeats for the reliability analysis  |
   | **Sum_Sub**  | Sum of interactions where an animal showed submissive behaviours  |
   | **Sum_Dom**  | Sum of interactions where an animal showed dominant behaviours  |
   | **Sum_Actions**  | Sum of all interactions an animal engaged in during observations  |
   | **Exclusion**  | Whether the observation did not fulfill criteria to classify for analysis  |

3. **focalSelection.csv**  
   Contains comb sizes of focal individuals, among other data.  

   **Columns:**  
   |----------|-------------|
   | **Pen**  | Animal group out of 3, 4, 5, 10, 11, 12  |
   | **ID**  | Backpack colour (non-unique identifier)  |
   | **Date**  | Date a hen was outfitted with a backpack  |
   | **HenID**  | Unique numeric identifier of a focal animal  |
   | **TagID**  | Numeric identifier of tracking tag  |
   | **Mean5_7**  | Comb size (cm²)  |

4. **ControlsComb.csv**  
   Data on control cohort, to compare to focal individuals.  

   **Columns:**  
   |----------|-------------|
   | **Pen**  | Animal group out of 3, 4, 5, 10, 11, 12  |
   | **Backpack**  | Unique numeric identifier of controls  |
   | **Area**  | Comb size (cm²)  |

5. **HA_all.csv**  
   Data of all physical assessments performed.  

   **Columns:**  
   |----------|-------------|
   | **Date**  | Date of physical assessment  |
   | **Pen**  | Animal group out of 3, 4, 5, 10, 11, 12  |
   | **Backpack**  | Backpack colour (non-unique identifier)  |
   | **Weight**  | Body mass of an animal (g)  |
   | **Comb**  | Severity degree of comb wounds (0-100, 100 = highest severity)  |
   | **Neck/Wings/Tail/Cloaca/Breast**  | Severity degree of feather loss on different areas (0-100, 100 = highest severity)  |
   | **Wounds**  | Severity degree of wounds (0-100, 100 = highest severity)  |
   | **\*_podo / \*_bumble / \*_injure**  | Severity degree of foot problems (0-100, 100 = highest severity)  |

6. **KBF_scores.csv**  
   Data of radiograph assessments of keel bone fracture severity.  

   **Columns:**  
   |----------|-------------|
   | **Date**  | Date of radiograph  |
   | **HenID**  | Pen ID in addition to backpack colour (non-unique identifier)  |
   | **Severity**  | Keel fracture severity degree (0-10, 10 = highest severity)  |

7. **trackingData.csv**  
   Output of `prepareTracking()` function. The tracking data used for the analysis after preprocessing.  

   **Columns:**  
   |----------|-------------|
   | **Time**  | Time and date stamp of transition between zones of the aviary  |
   | **Zone**  | New zone after transition (Wintergarten, Litter, Tier_2, Ramp_Nestbox, Tier_4)  |
   | **Hen**  | Unique identifier of animal ("Hen_00")  |
   | **PackID**  | Pen ID in addition to backpack colour (non-unique identifier)  |
   | **HenID**  | Unique numeric identifier of a focal animal  |
   | **Date**  | Date stamp  |
   | **Pen**  | Animal group out of 3, 4, 5, 10, 11, 12  |
   | **Light**  | Boolean to indicate whether lights are on (TRUE) or off (FALSE)  |
   | **TrueTransition**  | Boolean to indicate whether the entry is a transition into a new zone  |
   | **DayIndic**  | Boolean to indicate whether the entry represents the first or last entry of a new day  |
   | **LightIndic**  | Boolean to indicate whether the entry represents the moment the lights go on/off  |
   | **Duration**  | Duration within the zone from the current to the next transition of the animal  |
   | **NightCycle**  | Date stamp shifted to fit a whole dark cycle within one date  |
   | **WoA**  | Weeks of age of the animals  |

---

## R scripts:  

1. **analysis_Movement-Dominance.R**  
   Calculation of the Elo ratings for the chickens and analysis of hierarchy characteristics, interactions, and  
   other factors to compare small and large group sizes.  
   Creates data file `IndividualsOutput.csv` used in `Fear_Recognition_Badges`.  

2. **helper_functions.R**  
   Contains functions to ease working with the tracking data.  

3. **prepareTracking.R**  
   Function to clean and prepare tracking data for analysis.  