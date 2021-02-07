-   [Aim](#aim)
-   [Packages](#packages)
    -   [Install packages](#install-packages)
    -   [Load Packages](#load-packages)
-   [Data preparation](#data-preparation)
    -   [Load the data](#load-the-data)
    -   [Calculate the time difference between the MRI-scan and the
        trauma](#calculate-the-time-difference-between-the-mri-scan-and-the-trauma)
    -   [Structure](#structure)
    -   [Studying the data more
        thoroughly](#studying-the-data-more-thoroughly)
    -   [continue pre-processing the
        data](#continue-pre-processing-the-data)
    -   [David Nelson’s pupil fix](#david-nelsons-pupil-fix)
    -   [mutate MRI scores as ordered
        “factors”](#mutate-mri-scores-as-ordered-factors)
    -   [Data summary](#data-summary)
-   [Missing](#missing)
-   [Missing - tables and figures for the
    article](#missing---tables-and-figures-for-the-article)
-   [Multiple imputation](#multiple-imputation)
-   [Session information](#session-information)

Aim
===

The aim of this script is to pre-process data for use in the analysis
comparing MRI-based TAI grading systems.

Packages
========

Install packages
----------------

This step only has to be done once!

``` r
install.packages("tidymodels")
install.packages("tidyverse")
install.packages("knitr")
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
```

Load Packages
-------------

``` r
library(MASS)
library(tidymodels)
library(tidyverse)
library(knitr)
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
   ungroup() %>% 
   mutate(mri_date = as_date(mri_date)) %>% 
   mutate(trauma_date = as_date(trauma_date)) %>% 
   select(-c(GCSOpl, contains("_no_tai"), contains("total")))

Data <- Data %>% 
   mutate(mri_time = difftime(time1 = mri_date, time2 = trauma_date, 
                                     units = "days"))
```

Structure
---------

keep only the first 406 rows, corresponding to all patients admitted
prior to December 31 2018, and only patients who have had their
MRI-scans performed within the first 28 days after their trauma and were
above the age 15

``` r
Data <- Data %>%
  filter(id < 1593) %>% 
  drop_na(mri_date) %>%
   select(-contains("_date")) %>% 
   filter(mri_time <= 28) %>% 
  filter(Alder >= 15)
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
   mutate(Pupreak = as.factor(Pupreak)) %>% 
    mutate(MT = as.factor(MT)) %>%
   mutate(rotterdam = factor(rotterdam, ordered = TRUE)) %>%
   mutate(marshall = as.factor(marshall)) %>% 
   mutate(adams = factor(adams, ordered = TRUE)) %>% 
   select(-c(GOSB,GOSC, infarct_bg))
```

Blood pressure and oxygen saturation

**Not used during the current study**

``` r
# replace "0" with "NA"

Data <- Data %>% 
  replace_with_na(replace = list(SaO2Opl = 0,
                             BltrOpl = 0, FinalGOS = 0,
                             BltrKS = 0, BltrPS = 0, GCSPS = 0,
                             GCSKS = 0))

# Dichotomize the saturation into hypoxia (1) or no hypoxia (0) 
# Also dichotomize the bloodpressure into hypotension (1) 
# or no hypotension (0)

Data <- Data%>% 
   group_by(id) %>% 
   gather(key, Bltr, contains("Bltr")) %>% 
   mutate(bp = min(Bltr, na.rm = TRUE)) %>% 
   spread(key, Bltr) %>% 
   mutate(hypoxia = as_factor(ifelse(SaO2Opl < 90, "1", "0"))) %>% 
   mutate(hypotension = as_factor(ifelse(bp < 90, "1", "0"))) %>% 
   ungroup()

# removing the variables for prehospital oxygen saturation and prehospital blood pressure

Data <- Data %>% 
   select(-c(SaO2Opl,contains("Bltr"), bp, hypotension, hypoxia))
```

Create a dichotomized GOS variable and a mortality variable

``` r
Data <- Data %>% 
   mutate(dich_gos = as_factor(ifelse(FinalGOS == 4 | FinalGOS == 5, 
                                      "favourable", 
                                      "unfavourable"))) %>% 
   mutate(mortality = as_factor(ifelse(FinalGOS == 1, "Dead", 
                                      "Alive"))) 
```

Studying the data more thoroughly
---------------------------------

dich\_gos

``` r
Data %>% 
   select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(-contains("pons")) %>% 
   select(-contains("infarct")) %>%
   select(-contains("mes")) %>%
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Percent (unfavourable)" = round(unfavourable/Total
          *100, digits = 0)) %>% 
   arrange(desc(unfavourable)) %>% 
   kable(caption = "Percent of unfavourable outcomes following lesions outside of the brainstem")

Data %>% 
   select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("bg") | contains("_ic_"), dich_gos) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Percent (unfavourable)" = round(unfavourable/Total
          *100, digits = 0)) %>% 
   arrange(desc(unfavourable)) %>% 
   kable(caption = "Percent of unfavourable outcomes following lesions in the basal ganglia and the internal capsule")

Data %>% 
   select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("pons"), dich_gos) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Percent (unfavourable)" = round(unfavourable/Total
          *100, digits = 0)) %>% 
   arrange(desc(unfavourable)) %>% 
   kable(caption = "Percent of unfavourable outcomes following lesions in the pons")

Data %>% 
   select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("mes"), dich_gos) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Percent (unfavourable)" = round(unfavourable/Total
          *100, digits = 0)) %>%
   arrange(desc(unfavourable)) %>% 
   kable(caption = "Percent of unfavourable outcomes following lesions in the midbrain")

Data %>% 
   select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("infarct"), dich_gos) %>% 
   group_by(dich_gos) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(dich_gos, Frequency) %>% 
   mutate(Total = unfavourable+favourable) %>% 
   mutate("Percent (unfavourable)" = round(unfavourable/Total
          *100, digits = 0)) %>% 
   arrange(desc(unfavourable)) %>% 
   kable(caption = "Percent of unfavourable outcomes following infarcts")
```

mortality

``` r
Data %>% 
   select(mortality,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(-contains("pons")) %>% 
   select(-contains("infarct")) %>%
   select(-contains("mes")) %>% 
   group_by(mortality) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -mortality) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(mortality, Frequency) %>% 
   mutate(Total = Dead+Alive) %>% 
   mutate("Percent (Dead)" = round(Dead/Total
          *100, digits = 0)) %>% 
   arrange(desc(Dead)) %>% 
   kable(caption = "Percent mortalities following lesions outside of the brainstem")

Data %>% 
   select(mortality,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("bg") | contains("_ic_"), mortality) %>% 
   group_by(mortality) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -mortality) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(mortality, Frequency) %>% 
   mutate(Total = Dead+Alive) %>% 
   mutate("Percent (Dead)" = round(Dead/Total
          *100, digits = 0)) %>% 
   arrange(desc(Dead)) %>% 
   kable(caption = "Percent mortalities following lesions in the basal ganglia and internal capsule")

Data %>% 
   select(mortality,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("pons"), mortality) %>% 
   group_by(mortality) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -mortality) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(mortality, Frequency) %>% 
   mutate(Total = Dead+Alive) %>% 
   mutate("Percent (Dead)" = round(Dead/Total
          *100, digits = 0)) %>% 
   arrange(desc(Dead)) %>% 
   kable(caption = "Percent mortalities following lesions in the pons")

Data %>% 
   select(mortality,hem_subcort_unilat:infarct_cerebellum) %>% 
   select(contains("mes"), mortality) %>% 
   group_by(mortality) %>% 
   summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
   gather(key = "Variable", value = "Frequency", -mortality) %>% 
   na.omit() %>%
   ungroup() %>% 
   group_by(Variable) %>% 
   spread(mortality, Frequency) %>% 
   mutate(Total = Dead+Alive) %>% 
   mutate("Percent (Dead)" = round(Dead/Total
          *100, digits = 0)) %>%
   arrange(desc(Dead)) %>% 
   kable(caption = "Percent of unfavourable outcomes following lesions in the midbrain")

# Data %>% 
#    select(dich_gos,hem_subcort_unilat:infarct_cerebellum) %>% 
#    select(contains("infarct"), dich_gos) %>% 
#    group_by(dich_gos) %>% 
#    summarize_if(is.numeric, funs(sum(., na.rm = TRUE))) %>% 
#    gather(key = "Variable", value = "Frequency", -dich_gos) %>% 
#    na.omit() %>%
#    ungroup() %>% 
#    group_by(Variable) %>% 
#    spread(dich_gos, Frequency) %>% 
#    mutate(Total = unfavourable+favourable) %>% 
#    mutate("Percent (unfavourable)" = round(unfavourable/Total
#           *100, digits = 0)) %>% 
#    arrange(desc(unfavourable)) %>% 
#    kable(caption = "Percent of unfavourable outcomes following infarcts")
```

continue pre-processing the data
--------------------------------

Identify all the MRI variables

``` r
imaging_vars <- Data %>% 
   select(hem_subcort_unilat:infarct_cerebellum) %>% 
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

David Nelson’s pupil fix
------------------------

``` r
# Fix pupils

# levels(Pupreak)
 

# Princip -> om en sida saknas imupteras det frön den andra
# om en sida saknas pga extern skada räknas som normal

attach(Data)

Pupiller<-rep(4,length(Data$Pupreak))

levels(Data$Pupreak)

Pupreak<-Data$Pupreak

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
Data$Pupiller_old<-Data$Pupiller
Data$Pupiller<-Pupiller


ifelse(any(Pupiller==4),print('Some Pupil names not caught in program'), return(Data$Pupiller<-Pupiller))

remove (Pupiller)

# To check
# levels(Pupiller)
# which(Pupiller==4)
# Pupreak[which(Pupiller==4)]
# levels(Pupreak[which(Pupiller==4)])
# levels(Data$Pupiller[which(Pupiller==4)])
# 
detach(Data)

Data <- Data %>% 
   select(-Pupreak) %>% 
   mutate(Pupiller = factor(Pupiller, ordered = TRUE))
```

mutate MRI scores as ordered “factors”
--------------------------------------

``` r
Data <- Data %>% 
   mutate(adams = factor(adams, ordered = TRUE)) %>% 
   mutate(firsching = factor(firsching, ordered = TRUE)) %>% 
   mutate(hamdeh = factor(hamdeh, ordered = TRUE)) %>%  
   mutate(stockholm_mri = factor(stockholm_mri, ordered = TRUE)) %>% 
   select(-c(mri_time, contains("infarct"), contains("total")))
```

Data summary
------------

``` r
 Data %>% map_df(~(data.frame(class = class(.x),
                              unique_values = n_distinct(.x, na.rm = TRUE))),
                     .id = "variable") %>% 
   kable(caption = 
            "The variables, their class and the number of unique values",
         format = "latex", booktabs = T,linesep = "", 
         longtable = T) %>% 
   kable_styling(full_width = F, latex_options = "repeat_header")
```

Missing
=======

``` r
Data %>%
    gather(key = "key", value = "val") %>%
    mutate(is_missing = is.na(val)) %>%
    group_by(key, is_missing) %>%
    summarise(Missing = n()) %>%
    filter(is_missing==T) %>%
    select(-is_missing) %>%
    arrange(desc(Missing)) %>% 
   kable(caption = "Missing data",
         format = "latex", booktabs = T,linesep = "", 
         longtable = T) %>% 
   kable_styling(full_width = F, latex_options = "repeat_header")
```

``` r
Data %>%
   select(-c(stockholm_mri, hamdeh, firsching, adams)) %>% 
  gather(key = "key", value = "val") %>% 
  mutate(is_missing = is.na(val)) %>%
  group_by(key, is_missing) %>%
  summarise(Missing = n()) %>%
  group_by(is_missing) %>%
  summarise(Missing = sum(Missing)) %>%  
   mutate(sum_missing = sum(Missing)) %>% 
mutate(per=paste0(round(100*Missing/sum_missing,2),'%')) %>% 
   ungroup() %>% 
   filter(is_missing == "TRUE") %>% 
   select(Missing = per) %>% 
   kable()
```

Row plot

``` r
missing_values <- Data %>%
  gather(key = "key", value = "val") %>%
  mutate(is_na = is.na(val)) %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, is_na) %>%
  summarise(num.isna = n()) %>%
  mutate(pct = num.isna / total * 100)


levels <-
    (missing_values  %>% filter(is_na == T) %>% arrange(desc(pct)))$key

print(Data %>%
  mutate(id = row_number()) %>%
  gather(-id, key = "key", value = "val") %>%
  mutate(isna = is.na(val)) %>%
  ggplot(aes(key, id, fill = isna)) +
    geom_raster(alpha=0.8) +
    scale_fill_manual(name = "",
        values = c('steelblue', 'tomato3'),
        labels = c("Present", "Missing")) +
    scale_x_discrete(limits = levels) +
    labs(x = "Variable",
           y = "Row Number", title = "Missing values in rows") +
    coord_flip())
```

Missing - tables and figures for the article
============================================

``` r
missing_datatable <- Data %>%
    select("MRI: Haemorrhagic sequences" = hem_tha_unilat,
           "MRI: DWI sequence" = diff_tha_unilat,
           "MRI: FLAIR sequence" = flair_tha_unilat,
           "GCS" = GCSKS,
           "GOS" = dich_gos,
           "Age" = Alder,
           "Pupillary reactivity" = Pupiller,
           "Computed tomography" = rotterdam,
           "Multitrauma" = MT) %>% 
    gather(key = "key", value = "val") %>% 
    mutate(is_missing = is.na(val))


missing_table_saved <- missing_datatable %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, is_missing) %>%
  summarise(Missing = n()) %>%
  mutate("Proportion" = Missing / total * 100) %>% 
   filter(is_missing==T) %>%
   ungroup() %>% 
    select(-c(is_missing, total)) %>%
    arrange(desc(Missing)) %>% 
   mutate(Proportion = round(Proportion, digits = 0)) %>% 
   mutate(Proportion = paste(Proportion, " %")) %>% 
   mutate(Proportion = str_replace(Proportion, "0", "< 1")) %>% 
   rename(Variable = key) %>% 
     kable(format = "html", booktabs = T, linesep = "",
        caption = "Missing data",
        align = 'c') %>%
   kable_styling(full_width = F)

# save_kable(missing_table_saved, "../Tables/missing_table.html")
```

Row plot

``` r
missing_values_table <- missing_datatable %>%
  group_by(key) %>%
  mutate(total = n()) %>%
  group_by(key, total, is_missing) %>%
  summarise(missing_n = n()) %>%
  mutate(pct = missing_n / total * 100) %>% 
   ungroup()


levels <- (missing_values_table  %>% filter(is_missing == T) %>% arrange(desc(pct)))$key

missing_row_plot <- Data %>%
    select("MRI: Haemorrhagic sequences" = hem_tha_unilat,
           "MRI: DWI sequence" = diff_tha_unilat,
           "MRI: FLAIR sequence" = flair_tha_unilat,
           "GCS" = GCSKS,
           "GOS" = dich_gos,
           "Age" = Alder,
           "Pupillary reactivity" = Pupiller,
           "Computed tomography" = rotterdam,
           "Multitrauma" = MT) %>% 
  mutate(id = row_number()) %>%
  gather(-id, key = "key", value = "val") %>%
  mutate(is_na = is.na(val)) %>%
  ggplot(aes(key, id, fill = is_na)) +
    geom_raster(alpha=0.9) +
    scale_fill_manual(name = "",
        values = c('steelblue', 'tomato3'),
        labels = c("Present", "Missing")) +
    scale_x_discrete(limits = levels) +
    labs(x = "Variable",
           y = "Individual study participants", 
         title = "Missing values") +
    coord_flip() + theme_bw() + theme(axis.title.y = element_blank(), 
                         axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      plot.title = element_text(hjust = 0.5))  
   
# ggsave(filename = "../Figures/row_plot.tiff", plot = missing_row_plot, width = 170, height = 170, units = "mm", dpi = 300)
```

Multiple imputation
===================

``` r
Data <- Data %>% 
   select(-contains("hypo"))

imputed <- mice(data = Data, seed = 123, m = 10, maxit = 35) 

# in order to save the data so that it may be used for e.g. SVM.
saveRDS(imputed,'../derived_data/imputed_data.rds')
```

Session information
===================

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19041)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_Sweden.1252  LC_CTYPE=English_Sweden.1252   
    ## [3] LC_MONETARY=English_Sweden.1252 LC_NUMERIC=C                   
    ## [5] LC_TIME=English_Sweden.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] rcompanion_2.3.25  DescTools_0.99.38  pROC_1.16.1        glmnet_3.0-2      
    ##  [5] Matrix_1.2-18      cowplot_1.1.0      naniar_0.6.0       pscl_1.5.5        
    ##  [9] mice_3.8.0         lubridate_1.7.8    corrr_0.4.2        table1_1.2        
    ## [13] caTools_1.18.0     descr_1.1.4        readxl_1.3.1       corrplot_0.84     
    ## [17] summarytools_0.9.6 kableExtra_1.1.0   knitr_1.29         forcats_0.5.0     
    ## [21] stringr_1.4.0      readr_1.3.1        tidyverse_1.3.0    yardstick_0.0.7   
    ## [25] workflows_0.2.0    tune_0.1.1         tidyr_1.1.2        tibble_3.0.3      
    ## [29] rsample_0.0.8      recipes_0.1.13     purrr_0.3.4        parsnip_0.1.3     
    ## [33] modeldata_0.0.2    infer_0.5.3        ggplot2_3.3.2      dplyr_1.0.2       
    ## [37] dials_0.0.9        scales_1.1.0       broom_0.7.0        tidymodels_0.1.1  
    ## [41] MASS_7.3-51.5     
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.1.5    plyr_1.8.6         splines_3.6.2     
    ##   [4] listenv_0.8.0      TH.data_1.0-10     pryr_0.1.4        
    ##   [7] digest_0.6.25      foreach_1.5.0      htmltools_0.4.0   
    ##  [10] magick_2.4.0       fansi_0.4.1        magrittr_1.5      
    ##  [13] checkmate_2.0.0    globals_0.13.0     modelr_0.1.8      
    ##  [16] gower_0.2.1        matrixStats_0.55.0 sandwich_2.5-1    
    ##  [19] colorspace_1.4-1   blob_1.2.1         rvest_0.3.6       
    ##  [22] haven_2.2.0        xfun_0.16          tcltk_3.6.2       
    ##  [25] crayon_1.3.4       jsonlite_1.6.1     libcoin_1.0-6     
    ##  [28] Exact_2.0          zoo_1.8-7          survival_3.2-3    
    ##  [31] iterators_1.0.12   glue_1.4.2         gtable_0.3.0      
    ##  [34] ipred_0.9-9        webshot_0.5.2      shape_1.4.5       
    ##  [37] rapportools_1.0    mvtnorm_1.1-1      DBI_1.1.0         
    ##  [40] Rcpp_1.0.3         viridisLite_0.3.0  xtable_1.8-4      
    ##  [43] GPfit_1.0-8        Formula_1.2-3      stats4_3.6.2      
    ##  [46] lava_1.6.8         prodlim_2019.11.13 httr_1.4.2        
    ##  [49] ellipsis_0.3.0     modeltools_0.2-23  pkgconfig_2.0.3   
    ##  [52] multcompView_0.1-8 nnet_7.3-12        dbplyr_1.4.4      
    ##  [55] tidyselect_1.1.0   rlang_0.4.7        DiceDesign_1.8-1  
    ##  [58] munsell_0.5.0      cellranger_1.1.0   tools_3.6.2       
    ##  [61] cli_2.0.2          generics_0.0.2     evaluate_0.14     
    ##  [64] EMT_1.1            yaml_2.2.1         fs_1.3.1          
    ##  [67] pander_0.6.3       coin_1.3-1         visdat_0.5.3      
    ##  [70] future_1.19.1      xml2_1.2.2         compiler_3.6.2    
    ##  [73] rstudioapi_0.11    e1071_1.7-3        reprex_0.3.0      
    ##  [76] lhs_1.0.1          stringi_1.4.6      lattice_0.20-38   
    ##  [79] vctrs_0.3.4        pillar_1.4.6       lifecycle_0.2.0   
    ##  [82] furrr_0.1.0        lmtest_0.9-38      bitops_1.0-6      
    ##  [85] lmom_2.8           R6_2.4.1           gld_2.6.2         
    ##  [88] codetools_0.2-16   boot_1.3-23        assertthat_0.2.1  
    ##  [91] nortest_1.0-4      withr_2.3.0        multcomp_1.4-14   
    ##  [94] expm_0.999-5       parallel_3.6.2     hms_0.5.3         
    ##  [97] grid_3.6.2         rpart_4.1-15       timeDate_3043.102 
    ## [100] class_7.3-15       rmarkdown_2.3      base64enc_0.1-3
