-   [Aim](#aim)
-   [Packages](#packages)
    -   [Install packages](#install-packages)
    -   [Load Packages](#load-packages)
-   [Data preparation](#data-preparation)
    -   [Load the data](#load-the-data)
    -   [Calculate the time difference between the MRI-scan and the
        trauma](#calculate-the-time-difference-between-the-mri-scan-and-the-trauma)
    -   [Structure](#structure)
    -   [total TAI](#total-tai)
    -   [Unfavourable %](#unfavourable)
    -   [Calculating the precentage of unfavourable outcomes per
        variable](#calculating-the-precentage-of-unfavourable-outcomes-per-variable)
    -   [Calculating the precentage of mortalities per
        variable](#calculating-the-precentage-of-mortalities-per-variable)
    -   [continue pre-processing the
        data](#continue-pre-processing-the-data)
-   [Stockholm DAI](#stockholm-dai)
-   [Table depicting the frequency and severity of brain lesions seen on
    MRI](#table-depicting-the-frequency-and-severity-of-brain-lesions-seen-on-mri)
    -   [function used to combine tables with different row
        lengths](#function-used-to-combine-tables-with-different-row-lengths)
    -   [combining the results for FLAIR and the hemorrhagic sequences
        with DWI and the
        infarctions](#combining-the-results-for-flair-and-the-hemorrhagic-sequences-with-dwi-and-the-infarctions)
    -   [Table - unfavourable vs
        favourable](#table---unfavourable-vs-favourable)
    -   [Table - mortalities](#table---mortalities)
    -   [Totals](#totals)
    -   [Table - Stockholm MRI](#table---stockholm-mri)
-   [Demography](#demography)
    -   [Loading the full dataset, including the patients with no
        MRI](#loading-the-full-dataset-including-the-patients-with-no-mri)
    -   [Data processing](#data-processing)
    -   [Checking for normality](#checking-for-normality)
-   [*χ*<sup>2</sup> - test](#chi2---test)
-   [Demographics](#demographics)

Aim
===

The aim of this script is to analyse and display demographic data

Packages
========

Install packages
----------------

This step only has to be done once!

``` r
install.packages("tidymodels")
install.packages("tidyverse")
install.packages("knitr")
install.packages("ggfortify")
install.packages("kableExtra")
install.packages("summarytools")
install.packages("readxl")
install.packages("corrr", dependencies = T)
install.packages("RColorBrewer")
install.packages("corrplot")
install.packages('MASS')
install.packages('descr')
install.packages("CaTools")
install.packages("table1", dependencies = T)
install.packages("mice")
install.packages("pscl")
install.packages("naniar")
install.packages("cowplot")
install.packages("glmnet")
install.packages("pROC")
install.packages("rcompanion")
install.packages("DescTools")
install.packages("arsenal")
```

Load Packages
-------------

``` r
library(MASS)
library(tidymodels)
library(tidyverse)
library(knitr)
library(ggfortify)
library(kableExtra)
library(summarytools)
library(corrplot)
library(readxl)
library(descr) 
library(caTools)
library(table1)
library(corrr)
library(lubridate)
library(mice)
library(pscl)
library(naniar)
library(cowplot)
library(glmnet)
library(pROC)
library(DescTools)
library(rcompanion)
library(table1)
library(arsenal)

options(knitr.kable.NA = '')
```

Data preparation
================

Load the data
-------------

``` r
Data <- read_excel('../raw_data/dai_results_mri_hem.xls', na = "NA") 
```

Calculate the time difference between the MRI-scan and the trauma
-----------------------------------------------------------------

calculate the difference, in days, between the date of the trauma and
the date of the MRI-scan

``` r
Data <- Data %>% 
   mutate(mri_date = as_date(mri_date)) %>% 
   mutate(trauma_date = as_date(trauma_date)) %>% 
   select(-GCSOpl) 

Data <- Data %>% 
   mutate(mri_time = difftime(time1 = mri_date, time2 = trauma_date, 
                                     units = "days")) %>% 
  mutate(mri_time = as.numeric(mri_time)) %>% 
  mutate(dich_marshall = as.factor(ifelse(marshall %in% c("1","2", "3", "4"),
                                "Diffuse", 
                                "Focal"))) %>% 
  mutate(MT = as.factor(MT))
```

Structure
---------

keep only the first 406 rows and only patients who have had their
MRI-scans performed within the first 28 days after their trauma. Exclude
patients below the age of 15 (we have 1 14 year old in our cohort)

``` r
Data <- Data %>%
  filter(id < 1593) %>% 
  drop_na(mri_date) %>%
   select(-contains("_date")) %>% 
   filter(mri_time <= 28) %>% 
  filter(Alder >= 15) %>% 
  filter(!id %in% c(1041)) 
```

Examine the structure of the data

``` r
str(Data, list.len=ncol(Data))
```

change the class of certain variables, based on that which was seen
using the str() function above.

``` r
Data <- Data %>% 
   mutate(FinalGOS = factor(FinalGOS, ordered = TRUE)) %>% 
   mutate(Pupreak = as.factor(Pupreak))
```

Blood pressure and oxygen saturation

**Not used in the current study**

``` r
# replace "0" with "NA"

Data <- Data %>% 
  replace_with_na(replace = list(SaO2Opl = 0,
                             BltrOpl = 0, FinalGOS = 0))

# Dichotomize the saturation into hypoxia (1) or no hypoxia (0) 
# Also dichotomize the bloodpressure into hypotension (1) 
# or no hypotension (0)

Data <- Data%>% 
   mutate(hypoxia = as_factor(ifelse(SaO2Opl < 90, "1", "0"))) %>% 
   mutate(hypotension = as_factor(ifelse(BltrOpl < 90, "1", "0")))

# removing the variables for prehospital oxygen saturation and prehospital blood pressure

Data <- Data %>% 
   select(-c(SaO2Opl, BltrOpl))
```

Create a dichotomized GOS variable and a mortality variable

``` r
Data <- Data %>% 
   mutate(dich_gos = as_factor(ifelse(FinalGOS == 4 | FinalGOS == 5, 
                                      "favourable", "unfavourable"))) %>% 
   mutate(mortality = as_factor(ifelse(FinalGOS == 1, "dead", 
                                      "alive"))) 
```

total TAI
---------

``` r
Data %>% 
  group_by(total_total) %>% 
  summarise(Frequency = n()) %>% 
  mutate(Percentage = round(100 * (Frequency/sum(Frequency))))  %>% 
  kable(caption = "Percentage of patients with TAI, including NA's")

Data %>% 
  group_by(total_total) %>% 
  summarise(Frequency = n()) %>% 
  drop_na() %>% 
  mutate(Percentage = round(100 * (Frequency/sum(Frequency))))  %>% 
  kable(caption = "Percentage of patients with TAI, excluding NA's")
```

Unfavourable %
--------------

Identify all the MRI variables

``` r
imaging_vars <- Data %>% 
  select(-total_total) %>% 
   select(hem_subcort_unilat:infarct_cerebellum, 
          contains("total")) %>% 
   names() %>% 
   as_vector()
```

remove variables with fewer than 10 obvervations, as these have
essentially too few observations to be of any value in the analysis
(currently not run, but the code itself is kept for possible future
uses)

``` r
cols_to_remove <- imaging_vars[colSums(Data[imaging_vars], na.rm = TRUE) <= 12]

cols_to_keep <- imaging_vars[colSums(Data[imaging_vars], na.rm = TRUE) >= 12]
Data[!names(Data) %in% imaging_vars || names(Data) %in% cols_to_keep]
```

Calculating the precentage of unfavourable outcomes per variable
----------------------------------------------------------------

``` r
unfavourable <- Data %>% 
   select(dich_gos,imaging_vars) %>% 
  select(-contains("infarct_")) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Unfavourable outcome (%)" = round(unfavourable/Total
          *100, digits = 0)) %>% 
  rename(Frequency = Total) %>% 
  separate(Variable, c("Sequence", "Region", "lesion_type"), "_") %>% 
  mutate(Sequence = str_replace(Sequence, "flair", "FLAIR")) %>% 
  mutate(Sequence = str_replace(Sequence, "diff", "DWI")) %>%
  mutate(Sequence = str_replace(Sequence, "hem", "Hemorrhagic")) %>%
  mutate(Sequence = str_replace(Sequence, "total", "Total")) %>%
  mutate(Region = str_replace(Region, "cc", "Corpus Callosum")) %>% 
  mutate(Region = str_replace(Region, "no", 
                              "No TAI")) %>%
  mutate(Region = str_replace(Region, "mes", "Midbrain")) %>% 
  mutate(Region = str_replace(Region, "pons", "Pons")) %>% 
  mutate(Region = str_replace(Region, "tha", "Thalamus")) %>% 
  mutate(Region = str_replace(Region, "bg", "Basal ganglia")) %>% 
  mutate(Region = str_replace(Region, "ic", "Internal capsule")) %>% 
   mutate(Region = str_replace(Region, "subcort", "Subcortical")) %>% 
  mutate(lesion_type = str_replace(lesion_type, "unilat", 
                                   "Unilateral")) %>% 
    mutate(lesion_type = str_replace(lesion_type,                    "tai", "Total")) %>%
   mutate(lesion_type = str_replace(lesion_type, "bilateral", 
                                   "Bilateral")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "splenium", 
                                   "Splenium")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "genu", 
                                   "Genu and Rostrum")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "corpus", 
                                   "Trunk")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "tegmentum", 
                                   "Tegmentum")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "crus", 
                                   "Cerebral peduncles")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "tectum", 
                                   "Tectum")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "ventral", 
                                   "Ventral")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "dorsal", 
                                   "Dorsal")) %>% 
  select(Sequence, Region,"Lesion type" = lesion_type, 
         Frequency, "Unfavourable outcome (%)") %>% 
    group_by(Sequence, Region, `Lesion type` ) %>% 
  arrange(desc(`Lesion type`))

diff <- unfavourable %>% 
  filter(Sequence == "DWI") %>% 
  ungroup() %>% 
  select(-Sequence) %>% 
  group_by(Region, `Lesion type` ) %>% 
  arrange(Region)

flair <- unfavourable %>% 
  filter(Sequence == "FLAIR") %>% 
  ungroup() %>% 
  select(-Sequence) %>% 
  group_by(Region, `Lesion type` ) %>% 
  arrange(Region)

hem <- unfavourable %>% 
  filter(Sequence == "Hemorrhagic")%>% 
  ungroup() %>% 
  select(-Sequence) %>% 
  group_by(Region, `Lesion type` ) %>% 
  arrange(Region)

total <- unfavourable %>% 
  filter(Sequence == "Total")%>% 
  ungroup() %>% 
  select(-Sequence) %>% 
  group_by(Region, `Lesion type` ) %>% 
  arrange(Region)

# Data regarding infarcts was excluded from the
# current study
  
infarct <- Data %>% 
   select(dich_gos, contains("infarct_")) %>% 
  select(-infarct_bg) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Unfavourable outcome (%)" = round(unfavourable/Total
          *100, digits = 0)) %>% 
  rename(Frequency = Total) %>% 
  arrange(desc(Frequency)) %>% 
  ungroup() %>% 
  mutate(Variable = str_replace(Variable, "infarct_", "")) %>% 
  rename("Region" = Variable) %>% 
  mutate(Region = str_replace(Region, "posterior", 
                              "PCA infarct")) %>% 
  mutate(Region = str_replace(Region, "anterior", 
                              "ACA infarct")) %>% 
  mutate(Region = str_replace(Region, "media", 
                              "MCA infarct")) %>% 
  mutate(Region = str_replace(Region, "cerebellum", 
                              "Cerebellar infarct")) %>% 
  select("Lesion type" = Region, Frequency, "Unfavourable outcome (%)")

diff <- diff %>% 
  ungroup() %>% 
  select(-Region) %>% 
  bind_rows(infarct) 

freq <- bind_cols(flair, hem) %>% 
  ungroup() %>% 
  select(-contains("Region"))
```

Calculating the precentage of mortalities per variable
------------------------------------------------------

``` r
mortalities <- Data %>% 
   select(mortality,imaging_vars) %>% 
  select(-contains("infarct_")) %>% 
   group_by(mortality) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -mortality) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(mortality, Frequency) %>% 
   mutate(Total = dead+alive) %>% 
   mutate("Mortalities (%)" = round(dead/Total
          *100, digits = 0)) %>% 
  rename(Frequency = Total) %>% 
  separate(Variable, c("Sequence", "Region", "lesion_type"), "_") %>% 
  mutate(Sequence = str_replace(Sequence, "flair", "FLAIR")) %>% 
  mutate(Sequence = str_replace(Sequence, "diff", "DWI")) %>%
  mutate(Sequence = str_replace(Sequence, "hem", "Hemorrhagic")) %>%
  mutate(Sequence = str_replace(Sequence, "total", "Total")) %>%
  mutate(Region = str_replace(Region, "cc", "Corpus Callosum")) %>% 
  mutate(Region = str_replace(Region, "no", 
                              "No TAI")) %>% 
  mutate(Region = str_replace(Region, "mes", "Midbrain")) %>% 
  mutate(Region = str_replace(Region, "pons", "Pons")) %>% 
  mutate(Region = str_replace(Region, "tha", "Thalamus")) %>% 
  mutate(Region = str_replace(Region, "bg", "Basal ganglia")) %>% 
  mutate(Region = str_replace(Region, "ic", "Internal capsule")) %>% 
   mutate(Region = str_replace(Region, "subcort", "Subcortical")) %>% 
  mutate(lesion_type = str_replace(lesion_type, "unilat", 
                                   "Unilateral")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "bilateral", 
                                   "Bilateral")) %>% 
  mutate(lesion_type = str_replace(lesion_type, "tai", 
                                   "Total")) %>%
   mutate(lesion_type = str_replace(lesion_type, "splenium", 
                                   "Splenium")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "genu", 
                                   "Genu and Rostrum")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "corpus", 
                                   "Trunk")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "tegmentum", 
                                   "Tegmentum")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "crus", 
                                   "Cerebral peduncles")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "tectum", 
                                   "Tectum")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "ventral", 
                                   "Ventral")) %>% 
   mutate(lesion_type = str_replace(lesion_type, "dorsal", 
                                   "Dorsal")) %>% 
  select(Sequence, Region,"Lesion type" = lesion_type, 
         Frequency, "Mortalities (%)") %>% 
    group_by(Sequence, Region, `Lesion type` ) %>% 
  arrange(desc(`Lesion type`))

diff_mort <- mortalities %>% 
  filter(Sequence == "DWI") %>% 
  ungroup() %>% 
  select(-Sequence) %>% 
  group_by(Region, `Lesion type` ) %>% 
  arrange(Region)

flair_mort <- mortalities %>% 
  filter(Sequence == "FLAIR") %>% 
  ungroup() %>% 
  select(-Sequence) %>% 
  group_by(Region, `Lesion type` ) %>% 
  arrange(Region)

hem_mort <- mortalities %>% 
  filter(Sequence == "Hemorrhagic")%>% 
  ungroup() %>% 
  select(-Sequence) %>% 
  group_by(Region, `Lesion type` ) %>% 
  arrange(Region)

total_mort <- mortalities %>% 
  filter(Sequence == "Total")%>% 
  ungroup() %>% 
  select(-Sequence) %>% 
  group_by(Region, `Lesion type` ) %>% 
  arrange(Region)  

# Data regarding infarcts was excluded from the
# current study

infarct_mort <- Data %>% 
   select(mortality, contains("infarct_")) %>% 
  select(-infarct_bg) %>% 
   group_by(mortality) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -mortality) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(mortality, Frequency) %>% 
   mutate(Total = dead+alive) %>% 
   mutate("Mortalities (%)" = round(dead/Total
          *100, digits = 0)) %>% 
  rename(Frequency = Total) %>% 
  arrange(desc(Frequency)) %>% 
  ungroup() %>% 
  mutate(Variable = str_replace(Variable, "infarct_", "")) %>% 
  rename("Region" = Variable) %>% 
  mutate(Region = str_replace(Region, "posterior", 
                              "PCA infarct")) %>% 
  mutate(Region = str_replace(Region, "anterior", 
                              "ACA infarct")) %>% 
  mutate(Region = str_replace(Region, "media", 
                              "MCA infarct")) %>% 
  mutate(Region = str_replace(Region, "cerebellum", 
                              "Cerebellar infarct")) %>% 
  select("Lesion type" = Region, Frequency, "Mortalities (%)")

diff_mort <- diff_mort %>% 
  ungroup() %>% 
  select(-Region) %>% 
  bind_rows(infarct_mort) 

freq_mort <- bind_cols(flair_mort, hem_mort) %>% 
  ungroup() %>% 
  select(-contains("Region"))
```

continue pre-processing the data
--------------------------------

Coerce all binary imaging variables into factors instead of numeric
variables

``` r
Data <- Data %>% 
   mutate_at(imaging_vars, funs(factor(.)))
```

Drop variables with no identified lesions

``` r
Data <- Data[sapply(Data, function(x) !is.factor(x) | nlevels(x) > 1)]

Data <- Data %>% select_if(~sum(!is.na(.)) > 0)
```

Stockholm DAI
=============

``` r
ct_tai <- Data %>% 
  mutate(Stockholm_DAI = as.factor(Stockholm_DAI)) %>% 
  group_by(Stockholm_DAI) %>% 
  summarise(n = n()) %>% 
  mutate(percent = (n/351)*100) %>% 
  mutate(Stockholm_DAI = ifelse(is.na(Stockholm_DAI), "Missing", Stockholm_DAI))

ct_tai %>% kable(digits = 2)
```

Table depicting the frequency and severity of brain lesions seen on MRI
=======================================================================

function used to combine tables with different row lengths
----------------------------------------------------------

This function works by filling the table which has fewer rows with
additional rows that all say “NA”, such that the two tables are of equal
length and can thus be combined into a single table. As the tidyverse
functions, such as bind\_rows, do not allow columns with identical
names, I am forced to resort to “Base R” functions such as cbind and
rbind….

``` r
combine_tables <- function(x, y) {
    rows.x <- nrow(x)
    rows.y <- nrow(y)
    if (rows.x > rows.y) {
        diff <- rows.x - rows.y
        df.na <- matrix(NA, diff, ncol(y))
        colnames(df.na) <- colnames(y)
        cbind(x, rbind(y, df.na))
    } else {
        diff <- rows.y - rows.x
        df.na <- matrix(NA, diff, ncol(x))
        colnames(df.na) <- colnames(x)
        cbind(rbind(x, df.na), y)
    } 
}
```

combining the results for FLAIR and the hemorrhagic sequences with DWI and the infarctions
------------------------------------------------------------------------------------------

``` r
combined <- combine_tables(diff, freq)

combined_mort <- combine_tables(diff_mort, freq_mort)
```

Table - unfavourable vs favourable
----------------------------------

``` r
kable(combined, format = "html", booktabs = T,  digits = 2, 
             caption = "The frequency and severity of brain 
      lesions seen on MRI",
      col.names = c("Lesion type", "Frequency", 
                    "Unfavourable outcome (%)",  
                    "Lesion type", "Frequency", 
                    "Unfavourable outcome (%)", 
                    "Lesion type", "Frequency", 
                    "Unfavourable outcome (%)"), align= 'c') %>%
pack_rows("Basal ganglia", 1, 2) %>% 
  pack_rows("Corpus Callosum", 3, 5) %>% 
  pack_rows("Internal capsule", 6, 7) %>% 
  pack_rows("Midbrain", 8, 12) %>% 
  pack_rows("No TAI", 13, 13) %>%
  pack_rows("Pons", 14, 17) %>% 
  pack_rows("Subcortical", 18, 19) %>% 
  pack_rows("Thalamus", 20, 21) %>% 
  pack_rows("Infarcts", 22, 25) %>%
  row_spec(25, extra_css = "border-bottom: 1px solid") %>% 
add_header_above(c("DWI" = 3, "FLAIR" = 3, "Haemorrhagic" = 3))
```

Table - mortalities
-------------------

``` r
kable(combined_mort, format = "html", booktabs = T,  digits = 2, 
             caption = "The frequency and severity of brain 
      lesions seen on MRI",
      col.names = c("Lesion type", "Frequency", 
                    "Mortalities (%)",  
                    "Lesion type", "Frequency", 
                    "Mortalities (%)", 
                    "Lesion type", "Frequency", 
                    "Mortalities (%)"), align= 'c') %>%
pack_rows("Basal ganglia", 1, 2) %>% 
  pack_rows("Corpus Callosum", 3, 5) %>% 
  pack_rows("Internal capsule", 6, 7) %>% 
  pack_rows("Midbrain", 8, 12) %>% 
  pack_rows("No TAI", 13, 13) %>% 
  pack_rows("Pons", 14, 17) %>% 
  pack_rows("Subcortical", 18, 19) %>% 
  pack_rows("Thalamus", 20, 21) %>% 
  pack_rows("Infarcts", 22, 25) %>%
  row_spec(25, extra_css = "border-bottom: 1px solid") %>% 
add_header_above(c("DWI" = 3, "FLAIR" = 3, "Haemorrhagic" = 3))
```

Totals
------

``` r
total_both <- bind_cols(total, total_mort) %>% 
  ungroup() %>% 
  select(-contains("Region"))

kable(total_both, format = "html", booktabs = T,  digits = 2, 
             caption = "The frequency and severity of brain 
      lesions seen on MRI",
      col.names = c("Lesion type", "Frequency", 
                    "Unfavourable outcome (%)", 
                    "Lesion type", "Frequency", 
                    "Mortalities (%)"), align= 'c') %>%
pack_rows("Basal ganglia", 1, 2) %>% 
  pack_rows("Corpus Callosum", 3, 5) %>% 
  pack_rows("Internal capsule", 6, 7) %>% 
  pack_rows("Midbrain", 8, 12) %>% 
  pack_rows("No TAI", 13, 13) %>%
  pack_rows("Pons", 14, 17) %>% 
  pack_rows("Subcortical", 18, 19) %>% 
  pack_rows("Thalamus", 20, 21) %>% 
  row_spec(21, 
    extra_css = "border-bottom: 1px solid") %>% 
add_header_above(c("Dichotomised Glasgow outcome scale" = 3, "Mortality" = 3))
```

Table - Stockholm MRI
---------------------

``` r
stockholm_mri_grade <- Data %>% 
  select(dich_gos, mortality, stockholm_mri) %>% 
  drop_na() %>% 
   mutate("Grade 4" = ifelse(stockholm_mri == 4, 1, 0)) %>% 
   mutate("Grade 3" = ifelse(stockholm_mri == 3, 1, 0)) %>% 
  mutate("Grade 2" = ifelse(stockholm_mri == 2, 1, 0)) %>% 
  mutate("Grade 1" = ifelse(stockholm_mri == 1, 1, 0)) %>% 
  select(-stockholm_mri)

stockholm_mri_grade %>% 
  select(-mortality) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Unfavourable outcome (%)" = round(unfavourable/Total
          *100, digits = 0)) %>% 
  rename(Frequency = Total, "Stockholm MRI" = Variable) %>% 
  select("Stockholm MRI", Frequency, "Unfavourable outcome (%)") %>% 
  kable(format = "html", caption = "Unfavourable outcomes among each grade of the stockholm MRI grading system", booktabs = T)

stockholm_mri_grade %>% 
  select(-dich_gos) %>% 
group_by(mortality) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -mortality) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(mortality, Frequency) %>% 
   mutate(Total = dead+alive) %>% 
   mutate("Mortalities (%)" = round(dead/Total
          *100, digits = 0)) %>% 
  rename(Frequency = Total, "Stockholm MRI" = Variable) %>% 
  select("Stockholm MRI", Frequency, "Mortalities (%)") %>% 
  kable(format = "html", caption = "Mortalities among each grade of the stockholm MRI grading system", booktabs = T)
```

Demography
==========

Loading the full dataset, including the patients with no MRI
------------------------------------------------------------

``` r
full_data <- read_excel('../raw_data/AV-id dataset till Jonathan 191020.xlsx', na = "NA") %>% 
  select("id" = Patientnummer, "MRI" = MRT, GOSA, Alder, FinalGOS, GCSKS, VardDNIVA, Kon, Pupreak, Traumad, GOSB, GOSC, MT, penetrating)
```

Data processing
---------------

Keep only the trauma dates which are between 2005-01-01 and 2018-12-31

``` r
full_data <- full_data %>% 
   mutate(Traumad = as_date(Traumad)) %>% 
   filter(Traumad >= as.Date("2005-01-01") & Traumad <= as.Date("2018-12-31")) %>% 
     select(-Traumad)
```

fix NA in the age variable

id 1041 had a septal defect and an iatrogenic air embolus

``` r
full_data <- full_data %>% 
  filter(id != 1041) %>% 
  replace_with_na(replace = list(Alder = 0, FinalGOS = 0)) %>% 
  filter(Alder >= 15) %>% 
  filter(penetrating == "0")
```

change the class of certain variables, based on that which was seen
using the str() function above.

``` r
full_data <- full_data %>% 
   mutate(FinalGOS = factor(FinalGOS, ordered = TRUE)) %>% 
   mutate(Pupreak = as.factor(Pupreak)) %>% 
   mutate(VardDNIVA = round(VardDNIVA, digits = 1)) 
```

Create a dichotomized GOS variable

``` r
full_data <- full_data %>% 
   mutate(dich_gos = as_factor(ifelse(FinalGOS == 4 | FinalGOS == 5, 
                                      "favourable", "unfavourable"))) %>% 
  mutate(MT = as.factor(MT))
```

In-hospital mortality

``` r
full_data <- full_data %>% 
   mutate("In-hospital mortality" = as_factor(ifelse(GOSA == 1, 
                                      "Deceased prior to discharge", "Discharged"))) %>% 
  select(-GOSA)
```

Pupils

``` r
# Fix pupils

# levels(Pupreak)
 

# Princip -> om en sida saknas imupteras det frön den andra
# om en sida saknas pga extern skada räknas som normal

attach(full_data)

Pupiller<-rep(4,length(full_data$Pupreak))

levels(full_data$Pupreak)

Pupreak<-full_data$Pupreak

Pupiller[which(Pupreak=='Norm')]<-2
Pupiller[which(Pupreak=='TrögNorm')]<-2

Pupiller[which(Pupreak=='Trög')]<-2
Pupiller[which(Pupreak=='NormTrög')]<-2
Pupiller[which(Pupreak=='TrögNA')]<-2
Pupiller[which(Pupreak=='NormNorm')]<-2
Pupiller[which(Pupreak=='NATrög')]<-2

Pupiller[which(Pupreak=='TRÖG')]<-2
Pupiller[which(Pupreak=='TRÖ')]<-2
Pupiller[which(Pupreak=='TRÖ')]<-2
Pupiller[which(Pupreak=='TÖG')]<-2
Pupiller[which(Pupreak=='TRÖD')]<-2
Pupiller[which(Pupreak=='norm/stel (skada vä)')]<-2
Pupiller[which(Pupreak=='norm/stel (defekt vä)')]<-2


Pupiller[which(Pupreak=='TrögStel')]<-1
Pupiller[which(Pupreak=='NAStel')]<-1
Pupiller[which(Pupreak=='NormStel')]<-1
Pupiller[which(Pupreak=='StelTrög')]<-1
Pupiller[which(Pupreak=='StelNA')]<-1
Pupiller[which(Pupreak=='StelNorm')]<-1

Pupiller[which(Pupreak=='Stel')]<-0
Pupiller[which(Pupreak=='StelStel')]<-0
Pupiller[which(Pupreak=='Norm')]<-2

Pupiller[which(is.na(Pupreak))]<-NA
   
Pupiller[which(Pupreak=="NA")]<-NA
Pupiller[which(Pupreak=="u/norm(emaljöga hö)")]<-2
Pupiller[which(Pupreak== "" )]<-NA
Pupiller[which(Pupreak== "u" )]<-NA
Pupiller[which(Pupreak== "TöG"  )]<-2
Pupiller[which(Pupreak==  "Trög/Stel" )]<-1
Pupiller[which(Pupreak== "Trög/Norm"  )]<-2
Pupiller[which(Pupreak==  "Trög"  )]<-2
Pupiller[which(Pupreak=="TRöG")]<-2
Pupiller[which(Pupreak=="TRöD")]<-2
Pupiller[which(Pupreak== "StelNK" )]<-0
Pupiller[which(Pupreak=="Stel/Trög")]<-1
Pupiller[which(Pupreak=="stel/trög")]<-1
Pupiller[which(Pupreak=="stel/normal")]<-1
Pupiller[which(Pupreak=="stel/norm.\n" )]<-1
Pupiller[which(Pupreak=="Stel/Norm (skada hö)" )]<-2
Pupiller[which(Pupreak=="Stel/Norm")]<-1
Pupiller[which(Pupreak=="Stel/norm")]<-1
Pupiller[which(Pupreak=="stel/ej bedömbar")]<-0
Pupiller[which(Pupreak== "Stel" )]<-0
Pupiller[which(Pupreak== "stel")]<-0
Pupiller[which(Pupreak== "stel")]<-0
Pupiller[which(Pupreak== "saknar pupiller" )]<-NA
Pupiller[which(Pupreak=="op. bil.")]<-NA
Pupiller[which(Pupreak=="ONOR")]<-NA
Pupiller[which(Pupreak=="NormTrög")]<-2
Pupiller[which(Pupreak=="NormTrög")]<-2

Pupiller[which(Pupreak=="NATr\xf6g")]<-2

Pupiller[which(Pupreak=="StelTr\xf6g")]<-1
Pupiller[which(Pupreak=="Tr\xf6g")]<-2
Pupiller[which(Pupreak=="Tr\xf6gNA")]<-2
Pupiller[which(Pupreak=="Tr\xf6gNorm")]<-2
Pupiller[which(Pupreak=="Tr\xf6gStel")]<-1

Pupiller[which(Pupreak=="NormNK")]<-2
Pupiller[which(Pupreak=="NormNA")]<-2
Pupiller[which(Pupreak=="Norm/Trög")]<-2
Pupiller[which(Pupreak== "norm/trög")]<-2
Pupiller[which(Pupreak=="norm/stel(skada)")]<-2
Pupiller[which(Pupreak== "Norm/Stel(starr)")]<-2
Pupiller[which(Pupreak== "norm/stel (skada vö)" )]<-2
Pupiller[which(Pupreak=="norm/stel (defekt vö)")]<-2
Pupiller[which(Pupreak=="Norm/Stel")]<-1
Pupiller[which(Pupreak=="norm/stel" )]<-1
Pupiller[which(Pupreak=="Norm/NA" )]<-2
Pupiller[which(Pupreak=="NANorm" )]<-2
Pupiller[which(Pupreak=="Norm" )]<-2
Pupiller[which(Pupreak=="NORM" )]<-2
Pupiller[which(Pupreak=="Nej" )]<-0
Pupiller[which(Pupreak=="NK")]<-NA
Pupiller[which(Pupreak=="NKNorm")]<-2
Pupiller[which(Pupreak=="NATrög")]<-2
Pupiller[which(Pupreak=="NAStel")]<-0
Pupiller[which(Pupreak=="NANorm")]<-2
Pupiller[which(Pupreak==  "linsop. bil.  stela"  )]<-NA
Pupiller[which(Pupreak== "NANorm")]<-2
Pupiller[which(Pupreak==  "glaucomop. bil." )]<-NA
Pupiller[which(Pupreak==  "" )]<-NA
Pupiller[which(Pupreak=="0" )]<-0
Pupiller[which(Pupreak=="0,TR" )]<-1
Pupiller[which(Pupreak=="-/NORM" )]<-1
Pupiller[which(Pupreak=="stel/norm" )]<-1
Pupiller[which(Pupreak=="TRö" )]<-2
Pupiller[which(Pupreak=="Norm/Stel (starr)" )]<-2
Pupiller[which(Pupreak=="Normtrög)" )]<-2
Pupiller[which(Pupreak=="Normtrög" )]<-2
Pupiller[which(Pupreak=="NormTr\xf6g" )]<-2
Pupiller[which(Pupreak=="Normtrög" )]<-2
Pupiller[which(Pupreak=="Normtrög" )]<-2
Pupiller[which(Pupreak=="Steltrög" )]<-1
Pupiller[which(Pupreak=="stel/norm.\r\n")]<-1



Pupiller<-as.factor(Pupiller)
full_data$Pupiller_old<-full_data$Pupiller
full_data$Pupiller<-Pupiller


ifelse(any(Pupiller==4),print('Some Pupil names not caught in program'), return(full_data$Pupiller<-Pupiller))

remove (Pupiller)

# To check
# levels(Pupiller)
# which(Pupiller==4)
# Pupreak[which(Pupiller==4)]
# levels(Pupreak[which(Pupiller==4)])
# levels(Data$Pupiller[which(Pupiller==4)])
# 
detach(full_data)

full_data <- full_data %>% 
   select(-Pupreak)
```

Change the pupils from numerically coded to characters

``` r
full_data <- full_data %>% 
mutate(Pupiller = as_factor(case_when(Pupiller == "2" ~ "Responsive",
                        Pupiller == "1"  ~ "Unilaterally unresponsive",
                        Pupiller == "0"  ~ 
                           "Bilaterally unresponsive"))) #%>% 
  # mutate(GCSKS = as_factor(case_when(GCSKS <= 8 ~ "3-8",
  #                       GCSKS > 8 | GCSKS <= 13   ~ "9-13",
  #                       GCSKS > 13  ~ "14-15")))
```

Set MRI = NEJ on all patients in the full data table and create an if
else statement to change MRI to JA only if that patient is present in
the study dataset

``` r
mri_data <- Data %>% 
  select(id) %>% 
  mutate(MRI_study = "Ja")

full_data <- full_data %>% 
  mutate(MRI = "Nej") %>% 
  left_join(mri_data) %>% 
  mutate(MRI = ifelse(!is.na(MRI_study), MRI_study, MRI)) %>% 
  select(-MRI_study) 
```

Add the column which shows the time between the trauma and the MRI-scan,
as well as the dichotomized Marshall CT score

``` r
full_data <- Data %>% 
  select(id,mri_time, dich_marshall) %>% 
  right_join(full_data) %>% 
  mutate(mri_time = as.numeric(mri_time)) %>% 
  mutate(MRI = as_factor(MRI)) %>% 
  select(-penetrating)
```

Checking for normality
----------------------

``` r
full_data %>% 
  select_if(., is.numeric) %>% 
  select(- c(id, mri_time)) %>% 
    gather(key = "variable_name", value = "value", 
           Alder:VardDNIVA) %>% 
    group_by(variable_name) %>% 
   na.omit() %>% 
    do(tidy(shapiro.test(.$value))) %>% 
    ungroup() %>% 
    select(-method) %>% 
    kable(digits = 12, 
          caption = "Shapiro-wilk test for normality", booktabs = TRUE)
```

*χ*<sup>2</sup> - test
======================

$$\\chi^2 = \\sum \\frac {(O - E)^2}{E}$$

Used for *categorical* data

``` r
chi_sq_p_value <- full_data %>%
  select_if(is.factor) %>%
  select(c(-dich_marshall)) %>% 
  drop_na() %>% 
  select(dich_gos, everything()) %>% 
  select(- c(`In-hospital mortality`, MRI)) %>% 
    summarise_at(2:4,funs(chisq.test(., 
           dich_gos, simulate.p.value = TRUE)$p.value)) %>% 
  mutate(Estimate = "p.value")
 

chi_sq_statistic <- full_data %>%
  select_if(is.factor) %>% 
    select(c(-dich_marshall)) %>% 
  drop_na() %>% 
  select(dich_gos, everything()) %>% 
  select(- c(`In-hospital mortality`, MRI)) %>% 
    summarise_at(2:4,funs(chisq.test(., 
           dich_gos)$statistic)) %>% 
  mutate(Estimate = "Statistic")

chi_sq_results <- bind_rows(chi_sq_statistic, chi_sq_p_value) %>% 
  select(Estimate, everything())

chi_sq_results %>% 
  kable(digits = 3, booktabs = TRUE, caption = "$\\chi^2$ - test")
```

mann-whitney

``` r
full_data %>%
  select_if(is.numeric) %>%
  select(-mri_time) %>% 
  mutate(mri = full_data$MRI) %>% 
    summarise_at(2:4,funs(wilcox.test(. ~ mri)$p.value)) %>% 
  mutate(Estimate = "p.value")  %>% 
  kable(digits = 3, booktabs = TRUE, caption = "$\\chi^2$ - test")
```

Demographics
============

``` r
full_data <- full_data %>% 
  select(MRI, Alder, Kon, GCSKS, Pupiller, FinalGOS,VardDNIVA, mri_time) 


full_data <- full_data %>% 
  mutate(MRI = factor(ifelse(MRI == "Ja", "0", "1"), levels =c(0,1),
         labels = c("MRI scan performed", 
                    "MRI scan not performed"))) %>% 
  mutate(Kon = as_factor(ifelse(Kon == "Man", "Male", "Female")))

my_controls <- tableby.control(
  test = T,
  total = F,
  numeric.test = "kwt", cat.test = "chisq",
  numeric.stats = c("meansd", "medianq1q3", "range"),
  cat.stats = c("countpct"),
  ordered.stats = c("countpct"),
  stats.labels = list(
    meansd = "Mean (SD)",
    medianq1q3 = "Median (Q1-Q3)",
    iqr = "IQR",
    range = "Range (Min - Max)"
  )
)


my_labels <- list(
  Alder = "Age (years)",
  Kon = "Gender",
  GCSKS = "GCS",
  Pupiller = "Pupillary response",
  FinalGOS = "Glasgow outcome scale",
  VardDNIVA = "NCCU stay duration (days)",
  mri_time = "Time until the MRI examination (days)"
)
 
tab <- tableby(MRI ~ . ,data = full_data,
  control = my_controls)

summary(tab, labelTranslations = my_labels,
  title = "Patient demographics", pfootnote = TRUE, digits = 1
)
```
