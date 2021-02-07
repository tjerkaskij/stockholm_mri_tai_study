-   [Aim](#aim)
-   [Packages](#packages)
    -   [Install packages](#install-packages)
    -   [Load Packages](#load-packages)
-   [Data preparation](#data-preparation)
    -   [Load the data](#load-the-data)
    -   [Calculate the time difference between the MRI-scan and the
        trauma](#calculate-the-time-difference-between-the-mri-scan-and-the-trauma)
    -   [Structure](#structure)
    -   [continue pre-processing the
        data](#continue-pre-processing-the-data)
    -   [David Nelson’s pupil fix](#david-nelsons-pupil-fix)
-   [Select only variables that are to be used in the final
    analysis](#select-only-variables-that-are-to-be-used-in-the-final-analysis)
-   [Multiple imputation](#multiple-imputation)
-   [Session info](#session-info)

Aim
===

The aim of this script is to pre-process data for use in conjunction
with feature selection using a genetic algorithm and validation using
random forest.

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
install.packages("readxl")
#This may be required for the Support Vector Machines analysis
install.packages("mice")
install.packages("pscl")
install.packages("naniar")
```

Load Packages
-------------

``` r
library(data.table)
library(tidymodels)
library(tidyverse)
library(knitr)
library(kableExtra)
library(readxl)
library(lubridate)
library(mice)
library(pscl)
library(naniar)

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
   select(-c(GCSOpl, contains("_no_tai"), contains("total")))

Data <- Data %>% 
   mutate(mri_time = difftime(time1 = mri_date, time2 = trauma_date, 
                                     units = "days"))
```

Structure
---------

keep only the first 406 rows and only patients who have had their
MRI-scans performed within the first 28 days after their trauma

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

Create a dichotomized GOS variable

``` r
Data <- Data %>% 
   mutate(dich_gos = as_factor(ifelse(FinalGOS == 4 | FinalGOS == 5, 
                                      "favourable", 
                                      "unfavourable"))) 
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
```

    ##  [1] "NANorm"   "NAStel"   "NATrög"   "NK"       "Norm"     "NormStel"
    ##  [7] "NormTrög" "Stel"     "StelNorm" "StelTrög" "Trög"     "TrögNA"  
    ## [13] "TrögNorm" "TrögStel"

``` r
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
```

    ## [1] NA

``` r
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

Select only variables that are to be used in the final analysis
===============================================================

``` r
Data <- Data %>% 
   select(hem_subcort_unilat:flair_pons_dorsal, Alder, GCSKS,
         Pupiller, dich_gos)
```

Multiple imputation
===================

``` r
imputed <- mice(data = Data, seed = 123, m = 1, maxit = 50) 

saveRDS(imputed,'../derived_data/imputed_data_ml.rds')
```

Session info
============

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
    ##  [1] naniar_0.6.0      pscl_1.5.5        mice_3.8.0        lubridate_1.7.8  
    ##  [5] readxl_1.3.1      kableExtra_1.1.0  knitr_1.29        forcats_0.5.0    
    ##  [9] stringr_1.4.0     readr_1.3.1       tidyverse_1.3.0   yardstick_0.0.7  
    ## [13] workflows_0.2.0   tune_0.1.1        tidyr_1.1.2       tibble_3.0.3     
    ## [17] rsample_0.0.8     recipes_0.1.13    purrr_0.3.4       parsnip_0.1.3    
    ## [21] modeldata_0.0.2   infer_0.5.3       ggplot2_3.3.2     dplyr_1.0.2      
    ## [25] dials_0.0.9       scales_1.1.0      broom_0.7.0       tidymodels_0.1.1 
    ## [29] data.table_1.13.0
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] fs_1.3.1           webshot_0.5.2      DiceDesign_1.8-1   httr_1.4.2        
    ##  [5] tools_3.6.2        backports_1.1.5    R6_2.4.1           rpart_4.1-15      
    ##  [9] DBI_1.1.0          colorspace_1.4-1   nnet_7.3-12        withr_2.3.0       
    ## [13] tidyselect_1.1.0   compiler_3.6.2     cli_2.0.2          rvest_0.3.6       
    ## [17] xml2_1.2.2         digest_0.6.25      rmarkdown_2.3      pkgconfig_2.0.3   
    ## [21] htmltools_0.4.0    lhs_1.0.1          dbplyr_1.4.4       rlang_0.4.7       
    ## [25] rstudioapi_0.11    generics_0.0.2     jsonlite_1.6.1     magrittr_1.5      
    ## [29] Matrix_1.2-18      Rcpp_1.0.3         munsell_0.5.0      fansi_0.4.1       
    ## [33] GPfit_1.0-8        visdat_0.5.3       lifecycle_0.2.0    furrr_0.1.0       
    ## [37] stringi_1.4.6      pROC_1.16.1        yaml_2.2.1         MASS_7.3-51.5     
    ## [41] plyr_1.8.6         grid_3.6.2         blob_1.2.1         parallel_3.6.2    
    ## [45] listenv_0.8.0      crayon_1.3.4       lattice_0.20-38    haven_2.2.0       
    ## [49] splines_3.6.2      hms_0.5.3          pillar_1.4.6       codetools_0.2-16  
    ## [53] reprex_0.3.0       glue_1.4.2         evaluate_0.14      modelr_0.1.8      
    ## [57] vctrs_0.3.4        foreach_1.5.0      cellranger_1.1.0   gtable_0.3.0      
    ## [61] future_1.19.1      assertthat_0.2.1   xfun_0.16          gower_0.2.1       
    ## [65] prodlim_2019.11.13 viridisLite_0.3.0  class_7.3-15       survival_3.2-3    
    ## [69] timeDate_3043.102  iterators_1.0.12   lava_1.6.8         globals_0.13.0    
    ## [73] ellipsis_0.3.0     ipred_0.9-9
