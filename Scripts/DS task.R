#####################
## Holmusk RWE SDS Task
## Author: Rose Sisk
## Date: 23/01/2024
#####################

### load packages
library(here)
library(tidyverse)
library(gtsummary)
library(broom)
library(broom.helpers)
library(treemapify)

## load data
nm_fld <- 'SDS RWE Data Challenge Dataset'
clin_data <- read.csv(here(nm_fld, 'clinical_data.csv'), header = TRUE)
demo_data <- read.csv(here(nm_fld, 'demographics.csv'), header = TRUE)

#### clean demographic data 
table(demo_data$gender)
table(demo_data$race)
table(demo_data$resident_status)

demo_data <- demo_data %>%
  mutate(gender = case_when(gender == 'f' ~ 'Female',
                            gender == 'm' ~ 'Male',
                            TRUE ~ gender),
         race = case_when(race == 'chinese' ~ 'Chinese',
                          race == 'India' ~ 'Indian',
                          TRUE ~ race),
         resident_status = case_when(resident_status == 'Singaporean' ~ 'Singapore citizen',
                                     resident_status == 'PR' ~ 'Permanent resident',
                                     TRUE ~ resident_status))

### clean clinical data (limited to columns necessary for RQ)
## how many unique IDs?
length(unique(clin_data$id))
## 3000 unique IDs - restrict to first admission for each patient.

### create a minimum dataset with essential cols only for this analysis
min_df_clin <- clin_data %>% 
  select(
  patient_id = id, date_of_admission:trt_oth,
  cgis_adm) %>%
  left_join(demo_data, by = 'patient_id') %>%
  mutate(medical_history_hbp = as.numeric(case_when(medical_history_hbp == 'No' ~ '0',
                                         medical_history_hbp == 'Yes' ~ '1',
                                         TRUE ~ medical_history_hbp)),
         trt_antidep = ifelse(trt_adt == 1 | trt_ssr == 1, 1, 0),
         trt_std = case_when(trt_antidep == 1 & trt_the == 1 ~ 1,
                             TRUE ~ 0),
         cgis_sev_cat = ifelse(cgis_adm >= 4, 'Mod/severe', 'Mild'),
         do_adm = as.Date(date_of_admission, format = '%d/%m/%y'),
         do_disch = as.Date(date_of_discharge, format = '%d/%m/%y'),
         dob = as.Date(date_of_birth),
         age_adm = floor(interval(dob, do_adm) / years(1)),
         medical_history_sud = ifelse(is.na(medical_history_sud), 0, medical_history_sud)) %>%
  group_by(patient_id) %>%
  arrange(patient_id, do_adm) %>%
  slice(1) %>%
  ungroup()

#### create a treemap of composite outcome.
treemapdat <- min_df_clin %>%
  group_by(trt_antidep, trt_the, trt_std) %>%
  summarise(n = n(),
            prop = n/3000) 

treemapdat$label <- c("No treatment", "Therapy only", "AD only", "Treatment standard met")
         
  
 treemap_trt <- ggplot(treemapdat, aes(area = n,  
                                        fill = as.factor(trt_std),
                                        subgroup = trt_std,
                        label = paste0(label, '\n', ' (n = ', n, ', ', sprintf('%.1f%%',prop*100), ')'))) +
  geom_treemap() +
   geom_treemap_text(colour = "grey12", place = 'centre',  reflow = T) +
  geom_treemap_subgroup_border(colour = 'white', size = 3)  +
   scale_fill_manual(values = c('dodgerblue1', 'palegreen3')) +
   theme(legend.position = 'none')
 
 ggsave(plot = treemap_trt, filename = here('Output', 'treemap.png'), dpi = 600)

##########################################
### summarise cohort by treatment status.
### Table
##########################################
gtsummtab <- min_df_clin %>%
  select(age_adm, gender, race, resident_status, 
         cgis_sev_cat, medical_history_sud, medical_history_anx, medical_history_mood, 
         trt_antidep, trt_the, trt_std,
         trt_anx, trt_con, trt_oth) %>%
  mutate(trt_std = ifelse(trt_std == 0, 'No', 'Yes')) %>%
  tbl_summary(by = trt_std,
    statistic = list(all_continuous() ~ "{mean} ({sd})",
                                all_categorical() ~ "{n} ({p}%)"),
    label = list(age_adm ~ "Age at admission",
                 gender ~ "Gender",
                 race ~ "Race",
                 resident_status ~ "Residential status",
                 cgis_sev_cat ~ "CGIS Severity Category",
                 medical_history_sud ~ "Substance abuse disorder",
                 medical_history_anx ~ "Anxiety disorder",
                 medical_history_mood ~ "Other mood disorder",
                 trt_antidep ~ "Treated with antidepressants",
                 trt_the ~ "Treated with therapy",
                 trt_anx ~ "Treated with anxiolytics",
                 trt_con ~ "Treated with anticonvulsants",
                 trt_oth ~ 'Treated with other pysch. meds'))%>%
  modify_header(all_stat_cols() ~ "**{level}**, N = {n} ({style_percent(p)}%)") %>%
  add_overall()  %>%
  modify_spanning_header(c(stat_1, stat_2) ~ "**Treatment Standard Met**")

#gt::gtsave(as_gt(gtsummtab), "tab_1.png")

####### LR for associated factors
mod_trt_std <- glm(trt_std ~ I(age_adm/10) + gender + race + resident_status + cgis_sev_cat  + medical_history_sud  +medical_history_anx + medical_history_mood, data = min_df_clin, family = "binomial")
  
tidy_mod <- mod_trt_std %>%
  broom.helpers::tidy_and_attach(exponentiate = TRUE, conf.int = TRUE) %>%
  # adding in the reference row for categorical variables
  broom.helpers::tidy_add_reference_rows() %>%
  # adding a reference value to appear in plot
  broom.helpers::tidy_add_estimate_to_reference_rows() %>%
  # adding the variable labels
  broom.helpers::tidy_add_term_labels() %>%
  # removing intercept estimate from model
  broom.helpers::tidy_remove_intercept() %>%
  select(variable, reference_row, label, estimate:conf.high) %>%
  mutate(label2 = case_when(label == 'I(age_adm/10)' ~ 'Age (10yr increase)',
                            label == 'medical_history_sud' ~ 'Substance abuse disorder',
                            label == 'medical_history_anx' ~ 'Anxiety disorder',
                            label == 'medical_history_mood' ~ "Other mood disorder",
                            variable == 'cgis_sev_cat' ~ paste('CGIS category:', label),
                            TRUE ~ label),
         label2 = ifelse(reference_row == FALSE | is.na(reference_row), label2, paste(label2, '(ref)'))) 

#### create forest plot of logistic regression results.
forest_or <- tidy_mod %>%
  mutate(label2 = factor(label2, levels = rev(label2))) %>%
  ggplot(aes(y = label2, x = estimate)) +
  geom_vline(xintercept = 1, linetype="dashed", colour = 'grey') +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.4) +
  geom_point(shape = 18, size = 4, colour = "darkblue") +
  scale_x_continuous(limits = c(0.4, 1.6), breaks = seq(0.4, 1.8, by = 0.2)) +
  theme_bw() +
  xlab('Odds Ratio (95% CI)') +
  ylab('')+
  theme(text = element_text(size = 17)) 

ggsave(plot = forest_or, filename = here('Output', 'fig1.png'), dpi = 600)



