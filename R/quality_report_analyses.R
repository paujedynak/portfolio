## Analyses for the Quality Report for the compounds study
## Author: Paulina Jedynak
## Date 14/04/2023

options(warn = -1)

#Load packages
if (!"pacman" %in% installed.packages()[,"Package"]) {
  install.packages("pacman")
}

# Install and load packages
pacman::p_load(car, corrplot, factoextra, ggpubr, here, Hmisc, janitor, 
               magrittr, mycor, naniar, rio, tidyverse, stats, tableone)

# Source functions
source(here("R/functions.R"))

### Data preparation
## Load data
# Original data on samples
sample_orig <- import(here("data/data_sample/dataset_1.RDS"))

# Data on samples with corrected IDs
sample <- import(here("data/data_sample/dataset_2.xlsx"))

# Data on covariates (par, BMI, active smoking)
covariates <- import(here("data/data_sample/dataset_3.RDS"))

# Information on conc and del dates
sample_dates <- import(here("data/data_sample/dataset_4.RDS"))

# Data on compounds returned from the lab
results_compounds_cohort_1_cohort_2 <- import(here("data/data_compounds/dataset_5.xlsx"))

# Data on LODs and LOQs
LOD_LOQ <- import(here("data/data_compounds/LOD_LOQ_compounds.xlsx"))


## Create compounds names and clean names list
compounds_names <- paste0("comp", seq(1:9))

compounds_names_long <- paste0("Compound ", seq(1:9))

#=====================

## Merge info on sample and compounds data
# Preprocess sample variables
sample <- sample %>% 
  mutate(ID = as.character(subjectID),
         shipment_orig = ifelse(shipment == "pilot", 0, shipment),
         shipment_orig = factor(shipment_orig),
         
         shipment = ifelse(shipment == "pilot", 1, shipment),
         shipment = factor(shipment),
         
         # Rename variable where comments are introduced to distinguish from comments on the results
         comment_sample = comment,
         coll_time = as.Date(date_collection_corrected)) %>%
  
  # Separate Cohorts IDs
  separate(sampleID, c("cohort_1_id", "cohort_2_id"), sep = " ") %>% # NAs are produced 
  # where cohort_2 ID is missing
  
  # Select variables of interest
  select(ID, cohort_1_id, cohort_2_id, shipment, comment_sample, coll_time) %>%
  
  # Add conc, sample collection, bi, and sample analysis dates
  left_join(sample_dates, by = c("ID", "cohort_1_id")) # n = 1009

# Merge data on hospital and sample treatments
sample_orig <- sample_orig %>% 
  select(idparticipant, hospital, col_dye, col_treat) %>% # n = 852
  mutate(ID = as.character(idparticipant), 
         col_dye = factor(col_dye),
         col_treat = factor(col_treat))

sample <- sample %>% 
  left_join(sample_orig, by = "ID") %>% 
  
  # Fill missing hospital based on the sample ID
  mutate(hospital = ifelse(is.na(hospital), str_sub(ID, start = 1L, end = 2L), hospital),
         hospital = factor(hospital)) # n = 1009

#=====================



# Make the final file with all results and sample info
results_compounds_cohort_1_cohort_2 <- results_compounds_cohort_1_cohort_2  %>% 
  as.data.frame() %>% 
  mutate(ID = as.character(ID),
         batch_orig = ifelse(batch_results == "pilot", 0, batch_results),
         batch_orig = factor(batch_orig),
         batch = ifelse(batch_results == "pilot", "1", batch_results),
         batch = factor(batch)) %>% 
  
  # Merge with sample analysis results
  left_join(sample, by = c("ID", "cohort_1_id", "cohort_2_id")) %>%
  
  # Caluclate GW at sample collection
  mutate(GW_coll = round(difftime(coll_time, conc_estimate, units = "weeks"), 1),
         GW_coll = as.numeric(GW_coll),
         
         # Calucate gest age
         gest_age = round(difftime(p_part, conc_estimate, units = "weeks"), 1),
         gest_age = as.numeric(gest_age),
         
         # Change analysis date format
         analysis_date = as.Date(analysis_date, origin = "1970-01-01"), 
         
         # Calculate time elapsed from sample collection to analysis (in days)
         time_col_analysis = difftime(analysis_date, coll_time, units = "days"),
         time_col_analysis = as.numeric(time_col_analysis),
         
         # Calculate season of sample collection
         coll_season = case_when(
           lubridate::month(coll_time) %in% c(1:3) ~ "Winter",
           lubridate::month(coll_time) %in% c(4:6) ~ "Spring",
           lubridate::month(coll_time) %in% c(7:9) ~ "Summer",
           TRUE ~ "Autumn"),
         coll_season = factor(coll_season))

results_compounds_cohort_1_cohort_2 <- results_compounds_cohort_1_cohort_2 %>% 
  left_join(covariates, by = "ID") %>% 
  
  # Select variables of interest
  select(ID:cohort_2_id,
         batch_orig,
         batch,
         comment_results, 
         hospital,
         shipment,
         comment_sample,
         conc_estimate,
         coll_time,
         GW_coll, 
         p_part, 
         gest_age, 
         analysis_date,
         time_col_analysis,
         sample_weight,
         col_dye, 
         col_treat,
         coll_season,
         par,
         BMI_cat,
         active_smoking,
         comp1:comp9) %>% 
  
  # Remove samples that did not pass quality control
  filter(!str_detect(comment_sample, "exclude") | is.na(comment_sample))


#=====================


## Create the final file with compounds results
results_compounds <- results_compounds_cohort_1_cohort_2 %>% 
  
  # Add the cohort_2 sample collected at period 1 that could be added to Cohort 1
  filter(is.na(cohort_2_id) | cohort_2_id == "(C_52)") %>% 
  select(-cohort_2_id) 

length(results_compounds$ID) # n = 830

# Rename and select variables of interest
results_compounds <- results_compounds %>% 
  rename_at(vars(comp1:comp9), function(x) paste0(x, "_samplew"))

results_compounds_LOD <- results_compounds %>% 
  select(ID, batch, batch_orig, shipment, hospital, conc_estimate, coll_time, GW_coll, p_part, gest_age, time_col_analysis:comp9_samplew) %>% 
  
  # Remove adjustment for sample weight
  mutate(across(contains("samplew"), list(nadj = function(x) ifelse(x == 0.00001, x, x * sample_weight)))) %>% 
  rename_at(vars(contains("_nadj")), function(x) str_replace(x, "_nadj", "nadj"))



# For plots on data in wide format, replace values <LOD with NA
results_compounds <- results_compounds_LOD %>% 
  replace_with_na_at(data = ., .vars = compounds_names, condition = ~.x == 0.00001)


#=====================

# Create a dataset in a long format with present <LOD values
results_compounds_LOD_long <- results_compounds_LOD %>% 
  pivot_longer(cols = -c(ID:active_smoking),
               names_to = c("compound", "samplew_adj"),
               values_to = "value",
               names_pattern = "(.*)_(.*)") %>% 
  mutate(samplew_adj = factor(samplew_adj, 
                              levels = c("samplewnadj", "samplew"), 
                              labels = c("nadjusted", "adjusted")))


# Remove adjustment for sample weight for LODs/LOQs caluclations
results_compounds_LOD_long <- results_compounds_LOD_long %>% 
  left_join(LOD_LOQ, by = "compound") %>% 
  mutate(LOD = ifelse(samplew_adj == "nadjusted", LOD * 50, LOD),
         LOQ_low = ifelse(samplew_adj == "nadjusted", LOQ_low * 50, LOQ_low),
         LOQ_high = ifelse(samplew_adj == "nadjusted", LOQ_high * 50, LOQ_high),
         compound = factor(compound, levels = compounds_names),
         compound_name = factor(compound_name, levels = compounds_names_long)) %>% 
  select(ID, compound_name, compound, samplew_adj, value, LOD:LOQ_high, everything())


#=====================

# Create a wide dataset for plots etc., with values <LOD replaced with NA
results_compounds_long <- results_compounds_LOD_long %>% 
  mutate(value = ifelse(value < LOD, NA, value)) %>%
  mutate_if(is.difftime, as.numeric) %>% 
  arrange(compound)


#=====================


### Imputation of missing compounds concentrations values 
## Use a fill-in protocol applied in sepages
ln_i_comp <- results_compounds_LOD_long %>%
  rowwise() %>% 
  
  # Simulate values <LOD and <LOQ
  mutate(value_orig = ifelse(value == 0.00001, NA, value),
         value_orig = ifelse(value_orig > 0.00001 & value_orig < LOQ_low, NA, value),
         value_LOD = ifelse(value == 0.00001, runif(1, 0.00001, LOD), value),
         
         # ln-transform compounds concentrations
         ln_value = log(value_orig),
         ln_value_LOD = log(value_LOD)) %>% 
  select(ID, compound, compound_name, LOD:LOQ_high, value, value_orig:ln_value_LOD, everything()) %>% 
  
  pivot_wider(id_cols = c("ID", "samplew_adj"),
              names_from = compound,
              values_from = c("ln_value", "ln_value_LOD")) %>% 
  rename_all(function(x) str_replace(x, "value_", ""))


# Impute ln-transformed values
seed <- 652765
lim_LOD <- 0.00001

ln_i_comp <- ln_i_comp %>% 
  mutate(across(contains("ln_LOD_comp1"), 
                .fns = list(ln_i = ~fun_imp(., lim_LOD = log(lim_LOD), lim_infer = log(0.06), lim_super = log(0.063), seed = seed)),
                .names = "{fn}_{col}"))

ln_i_comp <- ln_i_comp %>% 
  mutate(across(contains("ln_LOD_comp2"), 
                .fns = list(ln_i = ~fun_imp(., lim_LOD = log(lim_LOD), lim_infer = log(0.07), lim_super = log(0.07), seed = seed)),
                .names = "{fn}_{col}"))

ln_i_comp <- ln_i_comp %>% 
  mutate(across(contains("ln_LOD_comp3"), 
                .fns = list(ln_i = ~fun_imp(., lim_LOD = log(lim_LOD), lim_infer = log(0.01), lim_super = log(0.016), seed = seed)),
                .names = "{fn}_{col}"))

ln_i_comp <- ln_i_comp %>% 
  mutate(across(matches("ln_LOD_comp4|ln_LOD_comp7"), 
                .fns = list(ln_i = ~fun_imp(., lim_LOD = log(lim_LOD), lim_infer = log(0.1), lim_super = log(0.4), seed = seed)),
                .names = "{fn}_{col}"))

ln_i_comp <- ln_i_comp %>% 
  mutate(across(matches("ln_LOD_comp5|ln_LOD_comp6"), 
                .fns = list(ln_i = ~fun_imp(., lim_LOD = log(lim_LOD), lim_infer = log(0.02), lim_super = log(0.08), seed = seed)),
                .names = "{fn}_{col}"))

ln_i_comp <- ln_i_comp %>%
  mutate(across(matches("ln_LOD_comp8|ln_LOD_comp9"),
                .fns = list(ln_i = ~fun_imp(., lim_LOD = log(lim_LOD), lim_infer = log(0.1), lim_super = log(0.12), seed = seed)),
                .names = "{fn}_{col}")) %>% 
  rename_at(vars(contains("ln_i_ln_LOD_")), function(x) str_replace(x, "ln_i_ln_LOD_", "ln_i_"))


ln_comp_long <- ln_i_comp %>% 
  select(ID, samplew_adj, starts_with("ln_"), -contains("LOD")) %>% 
  pivot_longer(cols = starts_with("ln_"),
               names_to = "compound",
               values_to = "ln_value",
               names_pattern = "ln_(.*)")

ln_i_comp_long <- ln_i_comp %>% 
  select(ID, samplew_adj, starts_with("ln_i_")) %>% 
  pivot_longer(cols = starts_with("ln_i_"),
               names_to = "compound",
               values_to = "ln_i_value",
               names_pattern = "ln_i_(.*)")

#=====================

# Merge imputed data with original dataset

results_compounds_i_long <- results_compounds_long %>% 
  left_join(ln_comp_long, by = c("ID", "compound", "samplew_adj")) %>%  
  left_join(ln_i_comp_long, by = c("ID", "compound", "samplew_adj")) %>% 
  
  mutate(compound = factor(compound, levels = compounds_names),
         
         # Add imputed but not log-transformed (exponentiated ln-transformed) variables
         i_value_from_ln = exp(ln_i_value)) %>% 
  select(ID:samplew_adj, LOD:LOQ_high, value, ln_value:i_value_from_ln, everything())


#=====================

# Format and clean the long dataset for plotting
plot_imp <- results_compounds_i_long %>% 
  select(compound_name, samplew_adj, contains("ln_")) %>% 
  pivot_longer(cols = contains("ln_"),
               names_to = "imputation",
               values_to = "value",
               names_pattern = "(.*)_value") %>% 
  mutate(samplew_adj = factor(samplew_adj, 
                              levels = c("nadjusted", "adjusted"), 
                              labels = c("Non-adjusted", "Adjusted")))


#=====================

# Change imputed dataset to wide format
results_compounds_i <- results_compounds_i_long %>% 
  pivot_wider(id_cols = c(ID, batch:active_smoking), 
              names_from = c("compound", "samplew_adj"),
              values_from = contains("value")) %>% 
  rename_all(function(x) str_replace(x, "value_", "")) %>% 
  select(ID:active_smoking, contains("_adjusted"), contains("_nadjusted"))


#=====================


### Cohort 1 compounds data characteristics
## compound concentrations values percentiles
# Detection rates for compound conc. 
compounds_detection_rates <- results_compounds_LOD_long %>% 
  select(compound_name:LOQ_high)

comp_det_rates_table <- compounds_detection_rates %>% 
  group_by(samplew_adj, compound_name) %>% 
  summarise(n_cohort_1 = sum(!is.na(value)),
            pct_below_LOD = round((sum(value < LOD, na.rm = TRUE) * 100 / n_cohort_1), 1),
            n_below_LOD = sum(value < LOD, na.rm = TRUE),
            pct_below_LOQ_low = round((sum(value < LOQ_low, na.rm = TRUE) * 100 / n_cohort_1), 1),
            n_below_LOQ_low = sum(value < LOQ_low, na.rm = TRUE),
            pct_above_LOQ_high = round((sum(value > LOQ_high, na.rm = TRUE) * 100) / n_cohort_1, 1),
            n_above_LOQ_high = sum(value > LOQ_high, na.rm = TRUE), 
            .groups = 'drop')


# Calculate statistics for the original dataset (where real <LOQ values were provided), to determine which values (min, p5, p25... are <LOQ)
results_compounds_LOD_orig_long <- results_compounds_LOD_long %>% 
  group_by(compound, samplew_adj) %>% 
  summarise(min_orig = min(value, na.rm = TRUE), 
            .groups = 'drop') %>% 
  mutate(compound = factor(compound))

# Calculate statistics for the dataset
compounds_descr_perc <- results_compounds_long %>% 
  group_by(compound, samplew_adj) %>% 
  summarise(min = min(value, na.rm = TRUE),              
            q5 = quantile(value, 0.05, na.rm = TRUE),
            q25 = quantile(value, 0.25, na.rm = TRUE),
            median = median(value, na.rm = TRUE),
            mean = mean(value, na.rm = TRUE),
            q75 = quantile(value, 0.75, na.rm = TRUE),
            q95 = quantile(value, 0.95, na.rm = TRUE),
            max = max(value, na.rm = TRUE), 
            .groups = 'drop') %>% 
  left_join(results_compounds_LOD_orig_long, by = c("compound", "samplew_adj")) %>% 
  left_join(LOD_LOQ, by = "compound") %>% 
  
  # Replace values <LOD with a "<LOD" string
  mutate(min = ifelse(min_orig < LOD, "<LOD", min)) %>% 
  pivot_wider(id_cols = c("compound", "compound_name"),
              names_from = "samplew_adj",
              values_from = min:max) %>% 
  mutate(compound = factor(compound, levels = compounds_names)) %>% 
  arrange(compound) %>% 
  select(compound_name, contains("_adjusted"), contains("nadjusted")) %>% 
  mutate_if(is.numeric, round, 1)


#=====================


## PCA
# Perform PCA and compare individual contribution
data_PCA_comp <- results_compounds_i %>% 
  column_to_rownames(var = "ID") %>% 
  select(starts_with("ln_i_")) %>% 
  select(-contains("nadjusted")) %>% 
  na.omit() %>% 
  rename_all(~ str_replace(., "ln_i_", "")) %>% 
  rename_all(~ str_replace(., "_adjusted", "")) 

data_PCA_cov <- results_compounds_i %>% 
  filter(ID %in% rownames(data_PCA_comp)) %>% 
  select(batch:active_smoking)

res_pca <- prcomp(data_PCA_comp, scale = TRUE)


#=====================


## Correlations between compounds
# Effect on compound concentrations - correlations
data_cor_plot <- results_compounds_i_long %>% 
  mutate(compound = compound_name,
         value = ln_i_value)

batch_corr_effects <- data.frame()

for (b_f in c("GW_coll", "time_col_analysis", "sample_weight")) {
  
  batch_corr <- compute_corr(model_data_long = data_cor_plot,
                             b_factor = b_f) %>% 
    mutate(term = b_f,
           corr = round(corr, 2))
  
  batch_corr_effects <- rbind(batch_corr_effects, batch_corr)
}

# Running the defined correlation function for each compound pair
colnames(data_PCA_comp) <- compounds_names_long
corr_compounds <- mycor(data_PCA_comp, method = "spearman", na.action = na.omit)


#=====================


## Effect on compound concentrations - univariate linear regressions
lm_cont_batch_effects <- data.frame()

for (b_f in c("sample_weight")) {
  
  batch_lm <- compute_lm(model_data_long = data_cor_plot,
                         b_factor = b_f)
  
  lm_cont_batch_effects <- rbind(lm_cont_batch_effects, batch_lm)
}



