-   [Aim](#aim)
-   [Packages](#packages)
    -   [Install packages](#install-packages)
    -   [Load Packages](#load-packages)
-   [Data preparation and splitting](#data-preparation-and-splitting)
    -   [Standardizing numeric
        variables](#standardizing-numeric-variables)
    -   [split](#split)
-   [Genetic algorithm](#genetic-algorithm)
    -   [Preparing functions that will be used later
        on](#preparing-functions-that-will-be-used-later-on)
    -   [Running the GA](#running-the-ga)
    -   [results GA](#results-ga)
-   [RF](#rf)
    -   [tune control](#tune-control)
    -   [custom random forest function](#custom-random-forest-function)
-   [Model performance](#model-performance)
    -   [GA RF](#ga-rf)
    -   [Core rf](#core-rf)
-   [ggplot graph](#ggplot-graph)
-   [Session info](#session-info)

Aim
===

The aim of this script is to run the genetic algorithm on the MRI data
and to validate the results using RF.

Packages
========

Install packages
----------------

This step only has to be done once!

``` r
install.packages("tidymodels")
install.packages("tidyverse")
install.packages("knitr")
install.packages("AppliedPredictiveModeling")
install.packages("caret")
install.packages("MLmetrics")
install.packages("ggfortify")
install.packages("kableExtra")
install.packages("summarytools")
install.packages("readxl")
install.packages("corrr", dependencies = T)
install.packages("RColorBrewer")
install.packages("corrplot")
#This may be required for the Support Vector Machines analysis
install.packages("kernlab")
install.packages('modelgrid')
install.packages('e1071')
install.packages('MASS')
install.packages('descr')
install.packages("DescTools")
install.packages("CaTools")
install.packages("table1", dependencies = T)
install.packages("mice")
install.packages("pscl")
install.packages("naniar")
install.packages("glmnet")
install.packages("pROC")
install.packages("stringi")
install.packages("randomForest")
install.packages("mltools")
install.packages("data.table")
install.packages("cowplot")
install.packages("GA")
install.packages("funModeling")
install.packages("desirability")
install.packages("doParallel")
install.packages(c("FactoMineR", "factoextra"))
```

Load Packages
-------------

``` r
library(MASS)
library(data.table)
library(tidymodels)
library(tidyverse)
library(knitr)
library(AppliedPredictiveModeling)
library(caret)
library(MLmetrics)
library(ggfortify)
library(kableExtra)
library(summarytools)
library(corrplot)
library(readxl)
library(modelgrid)
library(e1071)
library(descr)
library(DescTools)
library(caTools)
library(table1)
library(corrr)
library(lubridate)
library(mice)
library(pscl)
library(naniar)
library(stringi)
library(glmnet)
library(pROC)
library(mltools)
library(cowplot)
library(GA)
library(funModeling)
library(doParallel)
library(FactoMineR)
library(factoextra)


options(knitr.kable.NA = '')
```

Data preparation and splitting
==============================

Load the data, which has already undergone pre-processing and multiple
imputation (please examine the pdf document which contains the stepwise
selection procedure).

``` r
# Selecting a single imputed dataset

set.seed(99)

Data <- readRDS('../derived_data/imputed_data_ml.rds') %>% 
  mice::complete("long") %>%
  group_by(.imp) %>%
  nest() %>% 
  ungroup() %>% 
  sample_n(size =  1, replace = F) %>% 
  unnest() %>% 
 select(- c(.imp, .id)) 
```

check how many there are of each category

``` r
Data %>% 
  group_by(dich_gos) %>% 
  summarise(Frequency = n()) %>% 
  mutate(Percentage = round(100 * (Frequency/sum(Frequency))))  %>% 
  kable()
```

Examine the structure of the data

``` r
str(Data, list.len=ncol(Data))
```

Standardizing numeric variables
-------------------------------

``` r
pre_processing_recipe <- recipe(dich_gos ~ ., data = Data) %>% 
  step_normalize(all_numeric()) %>% 
  step_dummy(all_nominal(), - all_outcomes()) 

Data <- pre_processing_recipe %>% 
  prep(training = Data) %>% 
  juice(all_predictors()) %>% 
  mutate(dich_gos = Data$dich_gos)
```

split
-----

``` r
set.seed(99)

train_test_split <- initial_split(Data, strata = "dich_gos", prop = .66)

train_data <- training(train_test_split)
test_data <- testing(train_test_split)
```

Genetic algorithm
=================

Preparing functions that will be used later on
----------------------------------------------

Define paramateres to be used when evaluating population fitness when
running the Genetic Algorithm. Here, a desirability function was used
using the desirability package in R to balance discrimination (AUC) with
complexity (the number of included variables)

``` r
cvIndex <- createMultiFolds(train_data$dich_gos, k = 10, times = 5)

ctrl <- trainControl(method = "repeatedcv", 
                     number = 10,
                     repeats = 5, 
                     classProbs = TRUE, 
                     summaryFunction = twoClassSummary, 
                     allowParallel = FALSE, 
                     index = cvIndex)

Dcv <- function(ind, x, y, cntrl) 
{
    library(caret)
    library(MASS)
    library(desirability)
    ind <- which(ind == 1)
    if (length(ind) == 0) return(0)
    
     mtry = sqrt(ncol(x[, ind]))
tunegrid = expand.grid(.mtry=round(mtry))
    
    out <- train(x[, ind], y, 
                 method = "rf", 
                 metric = "ROC", 
                 trControl = cntrl,
                 tuneGrid = tunegrid)
    
    
    rocVal <- caret:::getTrainPerf(out)[, "TrainROC"]
    dROC <- dMax(0.5, 1)
    dPreds <- dMin(10, ncol(x))
    ## Comnined the two with a geometirc mean
    allD <- dOverall(dROC, dPreds)
    ## Calculate the overall desirability value
    predict(allD, data.frame(ROC = rocVal, NumPred = length(ind)))
}

initialSmall <- function(object, ...)
{
    population <- sample(0:1,
                         replace = TRUE,
                         size = object@nBits * object@popSize,
                         prob = c(0.7, 0.3))
    population <- matrix(population,
                         nrow = object@popSize,
                         ncol = object@nBits)
    return(population)
}
```

preparing for the GA

``` r
data_x <- train_data %>% 
  select(-dich_gos)

param_nBits=ncol(data_x)
col_names=colnames(data_x)

y <- train_data %>% 
  pull(dich_gos) 
```

Running the GA
--------------

``` r
GA <- ga(type = "binary", 
           fitness = Dcv, 
           lower = 0, upper = 1, 
           run = 200,
           maxiter = 2000, 
           nBits = param_nBits, 
           names = col_names, 
           x = data_x, 
           y = y, 
           cntrl = ctrl, 
           population = initialSmall,
           keepBest = TRUE, 
           parallel = FALSE,
           monitor = FALSE,
           seed=99)

saveRDS(GA, "../derived_data/ga_results.rds")
```

results GA
----------

``` r
GA <- readRDS("../derived_data/ga_results.rds")

plot(GA)
```

``` r
core_model <- c("Alder", "GCSKS", "Pupiller_1", "Pupiller_2")

only_letters <- function(x){      

  
  y <- strsplit(unlist(x), "[^a-zA-Z]+") 
  z <- y %>% map(~paste(., collapse = "_")) %>% simplify()
  return(z)}

chosen_vars_ga=col_names[GA@solution[1,]==1]

chosen_vars_ga <- setdiff(chosen_vars_ga, core_model)

chosen_vars_ga %>% 
  tibble() %>% 
  rename("Variable" = ".") %>% 
  mutate(Variable = only_letters(Variable)) %>% 
  mutate(Variable = str_replace(Variable, "_X", "")) %>% 
  kable(caption = "Variables selected by the GA")
```

RF
==

``` r
rec_ga <-recipe(chosen_vars_formula_ga, data = train_data)

core_formula <- paste("dich_gos","~", paste(core_model, 
                                                       collapse = " + "), 
                                 sep = " " ) %>% 
  as.formula()

rec_core <- recipe(core_formula, data = train_data) 
```

tune control
------------

In preparation for the repeated 10-fold cross validation

``` r
control_random <- trainControl(
 method = "repeatedcv",
 number = 10,
 repeats = 10,
 classProbs = TRUE, 
  summaryFunction = twoClassSummary, search = "random")

control_grid <- trainControl(
 method = "repeatedcv",
 number = 10,
 repeats = 10,
 classProbs = TRUE, 
  summaryFunction = twoClassSummary, search = "grid")
```

custom random forest function
-----------------------------

randomForest implements Breiman’s random forest algorithm for
classification and regression. The hyperparameters are:

ntree: the number of trees to grow mtry: Number of variables randomly
sampled as candidates at each split.

A custom fucntion (see below) was created to enable grid search across
both mtry and ntree using the caret package, as the standard function in
caret simply searches across mtry by default.

Other parameters, which were not set, but rather maintain their default
values, include the nodesize, i.e. the minimum size of terminal nodes,
and “maxnodes”, maximum number of terminal nodes trees in the forest can
have.

``` r
customRF <- list(type = "Classification", library = "randomForest", loop = NULL)
customRF$parameters <- data.frame(parameter = c("mtry", "ntree"), class = rep("numeric", 2), label = c("mtry", "ntree"))
customRF$grid <- function(x, y, len = NULL, search = "grid") {}
customRF$fit <- function(x, y, wts, param, lev, last, weights, classProbs, ...) {
  randomForest(x, y, mtry = param$mtry, ntree=param$ntree, ...)
}
customRF$predict <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata)
customRF$prob <- function(modelFit, newdata, preProc = NULL, submodels = NULL)
   predict(modelFit, newdata, type = "prob")
customRF$sort <- function(x) x[order(x[,1]),]
customRF$levels <- function(x) x$classes
```

``` r
grid_rf <- expand.grid(.mtry=c(1:10), .ntree=c(1000, 1500, 2000, 2500, 3000))


set.seed(99)

rf_core <- train(rec_core, data = train_data,
method = customRF,
tuneGrid = grid_rf,
trControl = control_grid,
  metric = "ROC")

# saveRDS(rf_core, file = "../derived_data/rf_core.rds")

set.seed(99)

rf_ga <- train(rec_ga, data = train_data,
method = customRF,
tuneGrid = grid_rf,
trControl = control_grid,
  metric = "ROC")

# saveRDS(rf_ga, file = "../derived_data/rf_ga.rds")
```

Model performance
=================

GA RF
-----

``` r
testing_results <- test_data %>%
        mutate("ga" = predict(rf_ga, test_data),
               "rf_core" = predict(rf_core, test_data))
  
ga <- conf_mat(testing_results, truth = dich_gos, estimate = ga)

autoplot(ga, type = "heatmap") +
  ggtitle("confusion matrix, rf with GA") +
  theme(plot.title = element_text(hjust = 0.5))

ga <- (summary(ga)) %>% 
  select("Metric" = .metric, "Estimate" = .estimate) %>% 
  mutate(Metric = str_replace(string=Metric, pattern='f_meas', replacement='F-score')) %>% 
 mutate(Metric = str_replace(string=Metric, pattern='j_index', replacement="Youden's index")) %>% 
  mutate(Metric = str_replace(string=Metric, pattern='kap', 
                              replacement="kappa")) %>% 
  filter(Metric %in% c("specificity", "sensitivity", "F-score", 
                       "Youden's index", "accuracy", "kappa", "precision",                              "recall"))

kable(ga,digits = 3, caption = "rf with GA")
```

ROC

``` r
prob_ga <- predict(rf_ga, newdata = test_data %>% select(-dich_gos), 
                      type = "prob")
roc_ga <- roc(response = test_data$dich_gos, predictor = prob_ga[, "unfavourable"])

roc_ga_graph <- ggroc(roc_ga) + theme_minimal() + ggtitle("ROC: rf with GA") + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(roc_ga$auc, 2))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)") + ggtitle("ROC: rf with GA")
```

Core rf
-------

``` r
core_rf <- conf_mat(testing_results, truth = dich_gos, estimate = rf_core)

autoplot(core_rf, type = "heatmap") +
  ggtitle("confusion matrix, rf with the core variables") +
  theme(plot.title = element_text(hjust = 0.5))

core_rf <- (summary(core_rf)) %>% 
  select("Metric" = .metric, "Estimate" = .estimate) %>% 
  mutate(Metric = str_replace(string=Metric, pattern='f_meas', replacement='F-score')) %>% 
 mutate(Metric = str_replace(string=Metric, pattern='j_index', replacement="Youden's index")) %>% 
  mutate(Metric = str_replace(string=Metric, pattern='kap', 
                              replacement="kappa")) %>% 
  filter(Metric %in% c("specificity", "sensitivity", "F-score", 
                       "Youden's index", "accuracy", "kappa", "precision",                              "recall"))

kable(core_rf,digits = 3, caption = "rf with the core variables")
```

ROC

``` r
prob_core_rf <- predict(rf_core, newdata = test_data %>% select(-dich_gos), 
                      type = "prob")
roc_core_rf <- roc(response = test_data$dich_gos, predictor = prob_core_rf[, "unfavourable"])

roc_core_rf_graph <- ggroc(roc_core_rf) + theme_minimal() + ggtitle("ROC: rf with the core variables") + 
    geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") +  annotate("text", x = .9, y = .2, 
       label = paste("AUC =", round(roc_core_rf$auc, 2))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")
```

``` r
print(roc_plots <- plot_grid(roc_core_rf_graph,
                       roc_ga_graph, 
                       rows = 1, cols = 2))
```

ggplot graph
============

``` r
ggroc(roc_core_rf, alpha = 0.7, colour = "Blue",
linetype = 1, size = 2) +
theme_minimal() +
geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed") + annotate("text", x= 0.68, y = 0.5, label = paste("AUC =", round(auc(roc_core_rf), digits = 2))) + xlab("False Positive Rate (1-Specificity)") + ylab("True Positive Rate (Sensitivity)")
```

``` r
# roc_list <- roc(test_data$dich_gos ~ prob_ga[, "unfavourable"] + prob_core_rf[, "unfavourable"])

results <- test_data %>%
        mutate("ga" = predict(rf_ga, test_data, type = "prob")[,"unfavourable"],
    "rf_core" = predict(rf_core, test_data, type = "prob")[,"unfavourable"]) %>% 
  select(dich_gos, "Clinical variables + MRI characteristics" = ga, "Clinical variables" = rf_core) %>% 
  gather(key = "Model", value = "fitted", - dich_gos) 


results %>% 
  group_by(Model) %>% 
  roc_curve(truth = dich_gos, fitted) %>% 
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = Model
    )
  ) + # plot with 2 ROC curves for each model
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, size = 0.4) +
  scale_color_manual(values = c("#48466D", "#3D84A8")) +
  coord_fixed() +
  theme_cowplot() 

# ggsave("auc_plot_initial.png")
```

``` r
results %>% 
  group_by(Model) %>%
roc_auc(truth = dich_gos, fitted) %>% 
  select(Model, "AUC" = .estimate) %>% 
  kable(digits = 2, format = "latex", booktabs = T, linesep = "", align = 'c') %>% 
  row_spec(0,  bold = T, color = "Black") %>%
  row_spec(2, bold = T, color = "#3D84A8") %>% 
  row_spec(1, bold = T, color = "#48466D") %>% 
  save_kable(., "../Figures/auc_table_rf.png")
```

``` r
p <- results %>% 
  group_by(Model) %>% 
  roc_curve(truth = dich_gos, fitted) %>% 
  ggplot(
    aes(
      x = 1 - Specificity, 
      y = Sensitivity, 
      color = Model
    )
  ) + # plot with 2 ROC curves for each model
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, size = 0.4) +
  scale_color_manual(values = c("#48466D", "#3D84A8")) +
  coord_fixed() +
  theme_cowplot() 
  
library(png)

auc_table <- readPNG("../Figures/auc_table_rf.png")

roc_plot <- ggdraw(p) + 
  draw_image(auc_table, scale = 0.45, x = 0.1, y = - 0.1)  

ggsave2("../Figures/roc_plot_rf.png", plot = roc_plot,
        dpi = 300)
```

Without a legend

``` r
l <- results %>% 
  group_by(Model) %>% 
  roc_curve(truth = dich_gos, fitted) %>% 
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = Model
    )
  ) + # plot with 2 ROC curves for each model
  geom_line(size = 1.1) +
  geom_abline(slope = 1, intercept = 0, size = 0.4) +
  scale_color_manual(values = c("#48466D", "#3D84A8")) +
  coord_fixed() +
  theme_cowplot() +
  theme(legend.position = "none")

ggsave2("../Figures/roc_plot_no_legend.png", plot = l,
        dpi = 300)
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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] factoextra_1.0.7                FactoMineR_2.3                 
    ##  [3] doParallel_1.0.15               funModeling_1.9.4              
    ##  [5] Hmisc_4.4-1                     Formula_1.2-3                  
    ##  [7] survival_3.2-3                  GA_3.2                         
    ##  [9] iterators_1.0.12                foreach_1.5.0                  
    ## [11] cowplot_1.1.0                   mltools_0.3.5                  
    ## [13] pROC_1.16.1                     glmnet_3.0-2                   
    ## [15] Matrix_1.2-18                   stringi_1.4.6                  
    ## [17] naniar_0.6.0                    pscl_1.5.5                     
    ## [19] mice_3.8.0                      lubridate_1.7.8                
    ## [21] corrr_0.4.2                     table1_1.2                     
    ## [23] caTools_1.18.0                  DescTools_0.99.38              
    ## [25] descr_1.1.4                     e1071_1.7-3                    
    ## [27] modelgrid_1.1.1.0               readxl_1.3.1                   
    ## [29] corrplot_0.84                   summarytools_0.9.6             
    ## [31] kableExtra_1.1.0                ggfortify_0.4.10               
    ## [33] MLmetrics_1.1.1                 caret_6.0-86                   
    ## [35] lattice_0.20-38                 AppliedPredictiveModeling_1.1-7
    ## [37] knitr_1.29                      forcats_0.5.0                  
    ## [39] stringr_1.4.0                   readr_1.3.1                    
    ## [41] tidyverse_1.3.0                 yardstick_0.0.7                
    ## [43] workflows_0.2.0                 tune_0.1.1                     
    ## [45] tidyr_1.1.2                     tibble_3.0.3                   
    ## [47] rsample_0.0.8                   recipes_0.1.13                 
    ## [49] purrr_0.3.4                     parsnip_0.1.3                  
    ## [51] modeldata_0.0.2                 infer_0.5.3                    
    ## [53] ggplot2_3.3.2                   dplyr_1.0.2                    
    ## [55] dials_0.0.9                     scales_1.1.0                   
    ## [57] broom_0.7.0                     tidymodels_0.1.1               
    ## [59] data.table_1.13.0               MASS_7.3-51.5                  
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] backports_1.1.5      plyr_1.8.6           lazyeval_0.2.2      
    ##   [4] splines_3.6.2        entropy_1.2.1        listenv_0.8.0       
    ##   [7] pryr_0.1.4           digest_0.6.25        htmltools_0.4.0     
    ##  [10] magick_2.4.0         fansi_0.4.1          magrittr_1.5        
    ##  [13] checkmate_2.0.0      cluster_2.1.0        ROCR_1.0-11         
    ##  [16] globals_0.13.0       modelr_0.1.8         gower_0.2.1         
    ##  [19] matrixStats_0.55.0   jpeg_0.1-8.1         colorspace_1.4-1    
    ##  [22] ggrepel_0.8.2        blob_1.2.1           rvest_0.3.6         
    ##  [25] haven_2.2.0          xfun_0.16            tcltk_3.6.2         
    ##  [28] crayon_1.3.4         jsonlite_1.6.1       Exact_2.0           
    ##  [31] glue_1.4.2           gtable_0.3.0         ipred_0.9-9         
    ##  [34] webshot_0.5.2        shape_1.4.5          rapportools_1.0     
    ##  [37] mvtnorm_1.1-1        DBI_1.1.0            Rcpp_1.0.3          
    ##  [40] CORElearn_1.54.2     plotrix_3.7-8        htmlTable_2.1.0     
    ##  [43] viridisLite_0.3.0    xtable_1.8-4         flashClust_1.01-2   
    ##  [46] foreign_0.8-72       GPfit_1.0-8          stats4_3.6.2        
    ##  [49] lava_1.6.8           prodlim_2019.11.13   htmlwidgets_1.5.1   
    ##  [52] httr_1.4.2           RColorBrewer_1.1-2   ellipsis_0.3.0      
    ##  [55] pkgconfig_2.0.3      nnet_7.3-12          dbplyr_1.4.4        
    ##  [58] tidyselect_1.1.0     rlang_0.4.7          DiceDesign_1.8-1    
    ##  [61] reshape2_1.4.3       munsell_0.5.0        cellranger_1.1.0    
    ##  [64] tools_3.6.2          cli_2.0.2            moments_0.14        
    ##  [67] generics_0.0.2       evaluate_0.14        yaml_2.2.1          
    ##  [70] ModelMetrics_1.2.2.2 fs_1.3.1             pander_0.6.3        
    ##  [73] visdat_0.5.3         future_1.19.1        nlme_3.1-142        
    ##  [76] leaps_3.1            xml2_1.2.2           compiler_3.6.2      
    ##  [79] rstudioapi_0.11      png_0.1-7            reprex_0.3.0        
    ##  [82] lhs_1.0.1            vctrs_0.3.4          pillar_1.4.6        
    ##  [85] lifecycle_0.2.0      furrr_0.1.0          bitops_1.0-6        
    ##  [88] lmom_2.8             latticeExtra_0.6-29  R6_2.4.1            
    ##  [91] gridExtra_2.3        rpart.plot_3.0.9     gld_2.6.2           
    ##  [94] codetools_0.2-16     boot_1.3-23          assertthat_0.2.1    
    ##  [97] withr_2.3.0          expm_0.999-5         hms_0.5.3           
    ## [100] grid_3.6.2           rpart_4.1-15         timeDate_3043.102   
    ## [103] class_7.3-15         rmarkdown_2.3        scatterplot3d_0.3-41
    ## [106] base64enc_0.1-3      ellipse_0.4.2
