-   [Aim](#aim)
-   [Packages](#packages)
    -   [Install packages](#install-packages)
    -   [Load Packages](#load-packages)
-   [Preparation to fit logistic regression
    models](#preparation-to-fit-logistic-regression-models)
    -   [A function to calculate the pseudo
        R-squared](#a-function-to-calculate-the-pseudo-r-squared)
    -   [Model formulas using the Rotterdam CT-score in the core
        model](#model-formulas-using-the-rotterdam-ct-score-in-the-core-model)
-   [Fitting the models](#fitting-the-models)
    -   [ROC curves](#roc-curves)
    -   [Likelihood ratio test](#likelihood-ratio-test)
    -   [De long’s test](#de-longs-test)
    -   [summary](#summary)
-   [Analysis using the Stockholm
    CT-score](#analysis-using-the-stockholm-ct-score)
-   [Analysis using the Helsinki
    CT-score](#analysis-using-the-helsinki-ct-score)
-   [Analysis using the Marshall
    CT-score](#analysis-using-the-marshall-ct-score)
-   [Mortality](#mortality)
-   [Session info](#session-info)

Aim
===

The aim of this script is to compare TAI based grading systems to each
other and to known outcome predictors in TBI

Packages
========

Install packages
----------------

This step only has to be done once!

``` r
install.packages("tidymodels")
install.packages("tidyverse")
install.packages("knitr")
install.packages("MLmetrics")
install.packages("kableExtra")
install.packages("summarytools")
install.packages("readxl")
install.packages("corrr", dependencies = T)
install.packages("RColorBrewer")
install.packages("corrplot")
install.packages('MASS')
install.packages('descr')
install.packages("DescTools")
install.packages("CaTools")
install.packages("mice")
install.packages("pscl")
install.packages("naniar")
install.packages("glmnet")
install.packages("pROC")
install.packages("stringi")
install.packages("data.table")
install.packages("cowplot")
```

Load Packages
-------------

``` r
library(MASS)
library(data.table)
library(tidymodels)
library(tidyverse)
library(knitr)
library(kableExtra)
library(summarytools)
library(corrplot)
library(readxl)
library(modelgrid)
library(descr)
library(DescTools)
library(caTools)
library(corrr)
library(lubridate)
library(mice)
library(pscl)
library(naniar)
library(stringi)
library(glmnet)
library(pROC)
library(cowplot)

options(knitr.kable.NA = '')
```

Preparation to fit logistic regression models
=============================================

A function to calculate the pseudo R-squared
--------------------------------------------

This function will used to calculate Nagelkerke’s Pseudo-*R*<sup>2</sup>
in the multivariate logistic regression models, which will be fitted
later on in this report.

``` r
# function used to calculate R squared

r2_calc <- function(data, model, nullmodel){
  CoxSnell <- 1 - exp(-(2/nrow(data) * (logLik(model) -  logLik(nullmodel))))
  r2_full <- CoxSnell/(1 - exp((2 * nrow(data)^(-1)) * logLik(nullmodel)))
  return(r2_full)
}
```

Model formulas using the Rotterdam CT-score in the core model
-------------------------------------------------------------

``` r
core_model <- c("Alder", "GCSKS", "Pupiller")

imputed_log <- readRDS('../derived_data/imputed_data.rds') %>% 
  mice::complete("long") %>%
    group_by(.imp) %>%
  select("Imputation" =.imp, id, dich_gos,
         adams, rotterdam, firsching, hamdeh,
          stockholm_mri, marshall, helsinki,
         all_of(core_model), mortality, MT,
         stockholm) %>%
   group_by(Imputation) %>%  nest()

core_vars_formula <- paste("dich_gos ~", paste(core_model, 
                                               collapse = " + "),
                            sep = " " ) %>% 
   as.formula()

core_rotterdam <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "rotterdam") %>% 
   as.formula()

core_stockholm <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "stockholm") %>% 
   as.formula()

core_helsinki <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "helsinki") %>% 
   as.formula()

core_marshall <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "marshall") %>% 
   as.formula()

core_adams <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "adams") %>% 
   as.formula()

core_karolinska <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "stockholm_mri") %>% 
   as.formula()

core_firsching <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "firsching") %>% 
   as.formula()

core_hamdeh <-  paste("dich_gos ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "hamdeh") %>% 
   as.formula()

var_names <- readRDS('../derived_data/imputed_data.rds') %>% 
  mice::complete("long") %>% 
   select(adams, stockholm_mri, hamdeh, firsching) %>%
   names() 

formulas_rotterdam <- paste('dich_gos ~ ', paste(core_model, 
                                                 collapse = " + "),
                            sep = " " , "+", "rotterdam", "+",var_names)

formulas_stockholm <- paste('dich_gos ~ ', paste(core_model, 
                                                 collapse = " + "),
                            sep = " " , "+", "stockholm", "+",var_names)

formulas_helsinki <- paste('dich_gos ~ ', paste(core_model, 
                                                 collapse = " + "),
                            sep = " " , "+", "helsinki", "+",var_names)

formulas_marshall <- paste('dich_gos ~ ', paste(core_model, 
                                                 collapse = " + "),
                            sep = " " , "+", "marshall", "+",var_names)

#formulas to be used in the the glm function

formulas_univariate <- paste('dich_gos ~ ', var_names)
```

Fitting the models
==================

``` r
estimates <- imputed_log %>%
mutate(adams = map(data, ~ glm(formulas_rotterdam[1], 
                               family = binomial(link="logit"), 
                               data = .))) %>%
  mutate(hamdeh = map(data, ~ glm(formulas_rotterdam[3], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
mutate(firsching = map(data, ~ glm(formulas_rotterdam[4], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
  mutate(karolinska = map(data, ~ glm(formulas_rotterdam[2], 
                                      family = binomial(link="logit"), data = .))) %>%
  mutate(null = map(data, ~ glm(dich_gos ~ 1, 
                              family = binomial(link="logit"), 
                              data = .))) %>%
    mutate(rotterdam = map(data, ~ glm(dich_gos ~ rotterdam, 
                              family = binomial(link="logit"), 
                              data = .))) %>% 
   mutate(core = map(data, ~ glm(core_vars_formula, family = binomial(link="logit"), data = .))) %>%
   mutate(core_rotterdam = map(data, ~ glm(core_rotterdam, family = binomial(link="logit"), data = .))) %>%
   mutate(core_adams = map(data, ~ glm(core_adams, family = binomial(link="logit"), data = .))) %>%
   mutate(core_karolinska = map(data, ~ glm(core_karolinska, family = binomial(link="logit"), data = .))) %>%
   mutate(core_firsching = map(data, ~ glm(core_firsching, family = binomial(link="logit"), data = .))) %>%
   mutate(core_hamdeh = map(data, ~ glm(core_hamdeh, family = binomial(link="logit"), data = .))) %>%
   mutate(r_sq_karolinska = pmap_dbl(list(data, karolinska, null), r2_calc)) %>%
  mutate(r_sq_hamdeh = pmap_dbl(list(data, hamdeh, null), 
                                r2_calc)) %>%   
  mutate(r_sq_firsching = pmap_dbl(list(data, firsching, null), 
                                 r2_calc)) %>% 
   mutate(r_sq_core = pmap_dbl(list(data, core, null), r2_calc)) %>%
  mutate(r_sq_rotterdam = pmap_dbl(list(data, rotterdam, null), r2_calc)) %>%
  mutate(r_sq_core_rotterdam = pmap_dbl(list(data, core_rotterdam, null), r2_calc)) %>%
  mutate(r_sq_core_adams = pmap_dbl(list(data, core_adams, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_firsching = pmap_dbl(list(data, core_firsching, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_hamdeh = pmap_dbl(list(data, core_hamdeh, null), 
                                    r2_calc)) %>%
  mutate(r_sq_core_karolinska = pmap_dbl(list(data, core_karolinska, null), r2_calc)) %>%
  mutate(r_sq_adams = pmap_dbl(list(data, adams, null), r2_calc)) %>%
mutate(glance = map(karolinska, ~ glance(.))) %>%
mutate(aug = map2(karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(tid = map(karolinska, ~ tidy(.))) %>%
mutate(glance_core = map(core, ~ glance(.))) %>% 
mutate(glance_core_rotterdam = map(core_rotterdam,
                                   ~ glance(.))) %>% 
  mutate(glance_rotterdam = map(rotterdam,
                                   ~ glance(.))) %>% 
mutate(glance_core_adams = map(core_adams, ~ glance(.))) %>% 
mutate(glance_core_hamdeh = map(core_hamdeh, ~ glance(.))) %>% 
mutate(glance_core_firsching = map(core_firsching, ~ glance(.))) %>% 
mutate(glance_core_karolinska = map(core_karolinska, ~ glance(.))) %>%
mutate(glance_adams = map(adams, ~ glance(.))) %>%
  mutate(glance_hamdeh = map(hamdeh, ~ glance(.))) %>%
   mutate(glance_firsching = map(firsching, ~ glance(.))) %>%
   mutate(aug_core = map2(core, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_rotterdam = map2(core_rotterdam, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
  mutate(aug_rotterdam = map2(rotterdam, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_adams = map2(core_adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_hamdeh = map2(core_hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_firsching = map2(core_firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_karolinska = map2(core_karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
mutate(aug_adams = map2(adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh = map2(hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_firsching = map2(firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(- c(adams, core_adams, core_rotterdam, 
             core_karolinska,
       core, hamdeh, firsching, karolinska, core_firsching, core_hamdeh, rotterdam)) %>% 
    mutate(karolinska_only = map(data, ~ glm(formulas_univariate[2], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(adams_only = map(data, ~ glm(formulas_univariate[1], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(hamdeh_only = map(data, ~ glm(formulas_univariate[3], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(firsching_only = map(data, ~ glm(formulas_univariate[4], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(r_sq_karolinska_only = pmap_dbl(list(data, karolinska_only, null), r2_calc)) %>%
  mutate(r_sq_adams_only = pmap_dbl(list(data, adams_only, null), r2_calc)) %>%
  mutate(r_sq_hamdeh_only = pmap_dbl(list(data, hamdeh_only, null), r2_calc)) %>%
  mutate(r_sq_firsching_only = pmap_dbl(list(data, firsching_only, null), r2_calc)) %>%
  mutate(glance_hamdeh_only = map(hamdeh_only, ~ glance(.))) %>% 
  mutate(glance_firsching_only = map(firsching_only, ~ glance(.))) %>% 
   mutate(glance_karolinska_only = map(karolinska_only,
                                       ~ glance(.))) %>%
  mutate(glance_adams_only = map(adams_only, ~ glance(.))) %>% 
  mutate(aug_firsching_only = map2(firsching_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh_only = map2(hamdeh_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_adams_only = map2(adams_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(-c(null, karolinska_only, hamdeh_only, 
        adams_only, firsching_only)) %>% 
     mutate_at(vars(contains("aug")), funs(map(., ~ roc(.x,dich_gos,
                                              .fitted)))) %>% 
  mutate_at(vars(contains("aug")), funs(map_dbl(., ~ auc(.x))))

names(estimates) <- gsub("aug", "auc", names(estimates))

model_vars <- estimates %>%  ungroup() %>% 
  sample_n(size =  1, replace = F) %>% 
   unnest(tid, .drop = T) %>% 
      pull(term) %>% 
     as_vector()
```

ROC curves
----------

``` r
core_formula_list <- tibble(models = c("core",
                    "core_adams", 
                    "core_stockholm",
                    "core_firsching",
                    "core_hamdeh"),
          formulas = c(core_vars_formula,                 core_adams, core_karolinska, core_firsching,
       core_hamdeh))

set.seed(123)

auroc_data <- imputed_log %>%
  ungroup() %>% 
  sample_n(1) %>% 
  select(data) %>% 
  unnest()

auroc <- core_formula_list %>% 
  mutate(fitted_models = map(formulas, ~
                              glm(., 
                               family = binomial(link="logit"), data = auroc_data))) %>% 
  mutate(predictions = map(fitted_models,
                      ~ as_tibble(predict(., 
        type = "response")))) %>% 
  mutate(predictions = map(predictions, 
                ~ rename(., "fitted" = "value"))) %>%
  mutate(predictions = map(predictions, 
                ~ mutate(., truth = auroc_data$dich_gos))) %>%
  select(-c(fitted_models, formulas))

roc_curve_adams <- auroc %>% 
  filter(models == "core" | models == "core_adams") %>%
  unnest(predictions) %>% 
  mutate(truth = as.factor(ifelse(truth == "unfavourable", 1, 0))) %>% 
  group_by(models) %>% 
  roc_curve(truth, fitted) %>% 
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = models
    )
  ) + # plot with 2 ROC curves for each model
  geom_path(size = 1) +
  geom_abline(lty = 3) +
  coord_equal() + 
  scale_color_manual(values = c("#48466D", 
                                "#3D84A8"),
                     labels = c('Core','Core + TAI grading system')) +
  theme(plot.title = element_text(size = 10, face = "bold")) +
  theme_cowplot() 
  

roc_curve_stockholm <- auroc %>% 
  filter(models == "core" | 
           models == "core_stockholm") %>%
  unnest(predictions) %>% 
  mutate(truth = as.factor(ifelse(truth == "unfavourable", 1, 0))) %>% 
  group_by(models) %>% 
  roc_curve(truth, fitted) %>% 
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = models
    )
  ) + # plot with 2 ROC curves for each model
  geom_path(size = 1) +
  geom_abline(lty = 3) +
  coord_equal() + 
  scale_color_manual(values = c("#48466D", 
                                "#3D84A8"),
                     labels = c('Core','Core + TAI grading system')) +
  theme(plot.title = element_text(size = 10, face = "bold")) +
  theme_cowplot() 

roc_curve_firsching <- auroc %>% 
  filter(models == "core" | 
           models == "core_firsching") %>%
  unnest(predictions) %>% 
  mutate(truth = as.factor(ifelse(truth == "unfavourable", 1, 0))) %>% 
  group_by(models) %>% 
  roc_curve(truth, fitted) %>% 
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = models
    )
  ) + # plot with 2 ROC curves for each model
  geom_path(size = 1) +
  geom_abline(lty = 3) +
  coord_equal() + 
  scale_color_manual(values = c("#48466D", 
                                "#3D84A8"),
                     labels = c('Core','Core + TAI grading system')) +
  theme(plot.title = element_text(size = 10, face = "bold")) +
  theme_cowplot() 

roc_curve_hamdeh <- auroc %>% 
  filter(models == "core" | 
           models == "core_hamdeh") %>%
  unnest(predictions) %>% 
  mutate(truth = as.factor(ifelse(truth == "unfavourable", 1, 0))) %>% 
  group_by(models) %>% 
  roc_curve(truth, fitted) %>% 
  ggplot(
    aes(
      x = 1 - specificity, 
      y = sensitivity, 
      color = models
    )
  ) + # plot with 2 ROC curves for each model
  geom_path(size = 1) +
  geom_abline(lty = 3) +
  coord_equal() + 
  scale_color_manual(values = c("#48466D", 
                                "#3D84A8"),
                     labels = c('Core','Core + TAI grading system')) +
   theme(plot.title = element_text(size = 10, face = "bold")) +
  theme_cowplot() 

roc_plots_no_legend <- plot_grid(roc_curve_adams + theme(legend.position = "none"),
                       roc_curve_firsching+ theme(legend.position = "none"),
                       roc_curve_hamdeh+ theme(legend.position = "none"),
                       roc_curve_stockholm+ theme(legend.position = "none"),
                       rows = 2, cols = 2,
  align = 'vh',
  labels = c("A", "B", "C", "D"),
  hjust = -1)

# extract the legend from one of the plots
# legend <- get_legend(
#   # create some space to the left of the legend
#   roc_curve_stockholm + theme(legend.box.margin = margin(0, 0, 0, 12))
# )

# add the legend to the right of the roc curve plot which was made previously. 

# 
# roc_plots <- plot_grid(roc_plots_no_legend, legend, rel_widths = c(3, .2), cols = 2, rows = 1)

legend <- get_legend(
  roc_curve_stockholm +
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

# add the legend underneath the plot with the roc curves.

roc_plots <- plot_grid(roc_plots_no_legend, legend, ncol = 1, rel_heights = c(1, .1))

# ggsave2("roc_plot.png", roc_plots)

# ggsave2("roc_plot_stockholm.png", roc_curve_stockholm)
```

Likelihood ratio test
---------------------

to compare the addition of MRI to core variables and the rotterdam
CT-score

``` r
data_table <- imputed_log %>% 
  filter(Imputation == 3) %>% 
  select(data) %>% 
  unnest()

core_model <- glm(core_vars_formula, 
                               family = binomial(link="logit"), 
                               data = data_table)

core_stockholm_mri <- glm(core_karolinska, 
                               family = binomial(link="logit"), 
                               data = data_table)

adams <- glm(formulas_rotterdam[1], 
                               family = binomial(link="logit"), 
                               data = data_table)


hamdeh <- glm(formulas_rotterdam[3], 
                               family = binomial(link="logit"), 
                               data = data_table)

firsching <- glm(formulas_rotterdam[4], 
                               family = binomial(link="logit"), 
                               data = data_table) 

karolinska <- glm(formulas_rotterdam[2], 
                                      family = binomial(link="logit"), data = data_table)

core_rotterdam <- glm(core_rotterdam, family = binomial(link="logit"), data = data_table) 


karolinska_only <- glm(formulas_univariate[2], 
                                      family = binomial(link="logit"), data = data_table) 

adams_only <- glm(formulas_univariate[1], 
                                      family = binomial(link="logit"),
                    data = data_table)

hamdeh_only <- glm(formulas_univariate[3], 
                                      family = binomial(link="logit"),
                    data = data_table)

firsching_only <- glm(formulas_univariate[4], 
                                      family = binomial(link="logit"),
                    data = data_table)


core_karolinska <- anova(core_rotterdam, karolinska , test="LRT")

tidy(core_karolinska) %>% 
  kable(digits = 3, caption = "comparison between CT score and Stockholm mri")

core_hamdeh <- anova(core_rotterdam, hamdeh , test="LRT")

tidy(core_hamdeh) %>% 
 kable(digits = 3, caption = "comparison between CT score and hamdeh")

core_firsching <- anova(core_rotterdam, firsching , test="LRT")

tidy(core_firsching) %>% 
  kable(digits = 3, caption = "comparison between CT score and firsching")

core_adams <- anova(core_rotterdam, adams , test="LRT")

tidy(core_adams) %>% 
 kable(digits = 3, caption = "comparison between CT score and adams")

stockholm_mri_core_comparison <- anova(core_model, core_stockholm_mri, test="LRT")

tidy(stockholm_mri_core_comparison) %>% 
  kable(digits = 3, caption = "comparison between Core and Core + Stockholm mri")


Comparison <- c("Core + Rotterdam CT score vs. Core + Rotterdam CT score + MRI grading system of Adams et al.", "Core + Rotterdam CT score vs. Core + Rotterdam CT score + MRI grading system of Firsching et al.",
"Core + Rotterdam CT score vs. Core + Rotterdam CT score + MRI grading system of Abu Hamdeh  et al.", 
"Core + Rotterdam CT score vs. Core + Rotterdam CT score + Stockholm MRI grading system")

log_likelihood_test <- bind_rows(slice(tidy(core_adams),-1),
               slice(tidy(core_firsching),-1),
               slice(tidy(core_hamdeh),-1),
               slice(tidy(core_karolinska),-1)) %>% 
  add_column(Comparison = Comparison) %>% 
  select(Comparison, p_value = p.value) %>% 
  mutate(p_value = round(p_value, 7)) %>% 
  mutate_at(vars(p_value), funs(ifelse(p_value < 0.05, "< 0.001", p_value))) %>% 
  kable(format = "html", booktabs = T, linesep = "",
        caption = "Log-likelihood test",
         align = 'c', col.names = c("Comparison", "p value")) %>%
   kable_styling(full_width = F)



# save_kable(log_likelihood_test, "../Tables/log_likelihood_test.html")
```

De long’s test
--------------

To see if there are any statistically significant differences in the AUC
between different MRI-based TAI grading systems

``` r
roc_karo <- augment(karolinska_only) %>% 
  roc(., dich_gos,.fitted)


roc_adams <- augment(adams_only) %>% 
  roc(., dich_gos,.fitted)

roc_hamdeh <- augment(hamdeh_only) %>% 
  roc(., dich_gos,.fitted)

roc_firsching <- augment(firsching_only) %>% 
  roc(., dich_gos,.fitted)


auc_p_values <- c(roc.test(roc_adams, roc_firsching, method=c("delong"), boot.n=2000, boot.stratified=TRUE)[["p.value"]],
roc.test(roc_adams, roc_hamdeh, method=c("delong"), boot.n=2000, boot.stratified=TRUE )[["p.value"]],

roc.test(roc_adams, roc_karo, method=c("delong"), boot.n=2000, boot.stratified=TRUE )[["p.value"]],

roc.test(roc_firsching, roc_hamdeh, method=c("delong"), boot.n=2000, boot.stratified=TRUE )[["p.value"]],

roc.test(roc_firsching, roc_karo, method=c("delong"), boot.n=2000, boot.stratified=TRUE )[["p.value"]],

roc.test(roc_hamdeh, roc_karo, method=c("delong"), boot.n=2000, boot.stratified=TRUE )[["p.value"]])

compare_auc <- c("MRI grading system of Adams et al. vs. MRI grading system of Firsching et al.",
"MRI grading system of Adams et al. vs. MRI grading system of Abu Hamdeh et al.",
"MRI grading system of Adams et al. vs. Stockholm MRI grading system",
"MRI grading system of Firsching et al. vs. MRI grading system of Abu Hamdeh et al.",
"MRI grading system of Firsching et al. vs. Stockholm MRI grading system",
"MRI grading system of Abu Hamdeh et al. vs. Stockholm MRI grading system")

de_long <- tibble(Comparison = compare_auc, p_value = auc_p_values) %>% 
  mutate(p_value = round(p_value,3)) %>% 
  kable(format = "html", booktabs = T, linesep = "",
        caption = "De Long's test",
         align = 'c', col.names = c("Comparison", "p value")) %>%
   kable_styling(full_width = F)

# save_kable(de_long, "../Tables/de_long_test.html")
```

summary
-------

``` r
core_comparison <- estimates %>%
  unnest(glance) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska, AIC, 
                "AUC" = auc) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + rotterdam CT-score + Stockholm MRI grading system") 

core_comparison <- estimates %>%
  unnest(glance_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh,AIC, 
                "AUC" = auc_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + rotterdam CT-score + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching,AIC, 
                "AUC" = auc_firsching) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + rotterdam CT-score + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams, AIC, "AUC" = auc_adams) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + rotterdam CT-score + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_karolinska) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_karolinska,AIC, 
                "AUC" = auc_core_karolinska) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm MRI grading system") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_hamdeh,AIC, 
                "AUC" = auc_core_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_firsching,AIC, 
                "AUC" = auc_core_firsching) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_adams, AIC, 
                "AUC" = auc_core_adams) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core_rotterdam) %>%
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_rotterdam, AIC, 
                "AUC" = auc_core_rotterdam) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + rotterdam CT-score") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_karolinska_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska_only, AIC, 
                "AUC" = auc_karolinska_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Stockholm MRI grading system") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_hamdeh_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh_only, AIC, 
                "AUC" = auc_hamdeh_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_firsching_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching_only, AIC, 
                "AUC" = auc_firsching_only) %>% 
summarize_all(funs(mean)) %>%
   mutate(Model = "MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_adams_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams_only, AIC, 
                "AUC" = auc_adams_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_rotterdam) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_rotterdam, AIC, 
                "AUC" = auc_rotterdam) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "rotterdam CT-score") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())

core_comparison <- estimates %>%
  unnest(glance_core) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core, AIC, 
                "AUC" = auc_core) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core") %>% 
   bind_rows(core_comparison) %>% 
   select(Model, everything())
  

core_comparison %>%
  kable(format = "latex", escape = FALSE, booktabs = T, linesep = "",
        caption = "Average result for all imputed datasets",
        digits = 2, align = 'c', col.names = c("Model","Pseudo-R<sup>2</sup>",
                                               "AIC",
                                                "AUC")) %>%
   kable_styling(full_width = F)

# results_table_rotterdam_orig <- core_comparison %>%
#   kable(format = "html",escape = FALSE, booktabs = T, linesep = "",
#         caption = "Average result for all imputed datasets",
#          digits = c(0,2,0,2), align = 'c', col.names = c("Model","Pseudo-R<sup>2</sup>",
#         "AIC",
#                                                 "AUC")) %>%
#    kable_styling(full_width = F)
# 
# 
# save_kable(results_table_rotterdam_orig, "../Tables/rotterdam_results.html")
```

Analysis using the Stockholm CT-score
=====================================

fitting the models

``` r
core_model <- c("Alder", "GCSKS", "Pupiller")

estimates_stockholm <- imputed_log %>%
mutate(adams = map(data, ~ glm(formulas_stockholm[1], 
                               family = binomial(link="logit"), 
                               data = .))) %>%
  mutate(hamdeh = map(data, ~ glm(formulas_stockholm[3], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
mutate(firsching = map(data, ~ glm(formulas_stockholm[4], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
  mutate(karolinska = map(data, ~ glm(formulas_stockholm[2], 
                                      family = binomial(link="logit"), data = .))) %>%
  mutate(null = map(data, ~ glm(dich_gos ~ 1, 
                              family = binomial(link="logit"), 
                              data = .))) %>%
    mutate(stockholm = map(data, ~ glm(dich_gos ~ stockholm, 
                              family = binomial(link="logit"), 
                              data = .))) %>% 
   mutate(core = map(data, ~ glm(core_vars_formula, family = binomial(link="logit"), data = .))) %>%
   mutate(core_stockholm = map(data, ~ glm(core_stockholm, family = binomial(link="logit"), data = .))) %>%
   mutate(core_adams = map(data, ~ glm(core_adams, family = binomial(link="logit"), data = .))) %>%
   mutate(core_karolinska = map(data, ~ glm(core_karolinska, family = binomial(link="logit"), data = .))) %>%
   mutate(core_firsching = map(data, ~ glm(core_firsching, family = binomial(link="logit"), data = .))) %>%
   mutate(core_hamdeh = map(data, ~ glm(core_hamdeh, family = binomial(link="logit"), data = .))) %>%
   mutate(r_sq_karolinska = pmap_dbl(list(data, karolinska, null), r2_calc)) %>%
  mutate(r_sq_hamdeh = pmap_dbl(list(data, hamdeh, null), 
                                r2_calc)) %>%   
  mutate(r_sq_firsching = pmap_dbl(list(data, firsching, null), 
                                 r2_calc)) %>% 
   mutate(r_sq_core = pmap_dbl(list(data, core, null), r2_calc)) %>%
  mutate(r_sq_stockholm = pmap_dbl(list(data, stockholm, null), r2_calc)) %>%
  mutate(r_sq_core_stockholm = pmap_dbl(list(data, core_stockholm, null), r2_calc)) %>%
  mutate(r_sq_core_adams = pmap_dbl(list(data, core_adams, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_firsching = pmap_dbl(list(data, core_firsching, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_hamdeh = pmap_dbl(list(data, core_hamdeh, null), 
                                    r2_calc)) %>%
  mutate(r_sq_core_karolinska = pmap_dbl(list(data, core_karolinska, null), r2_calc)) %>%
  mutate(r_sq_adams = pmap_dbl(list(data, adams, null), r2_calc)) %>%
mutate(glance = map(karolinska, ~ glance(.))) %>%
mutate(aug = map2(karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(tid = map(karolinska, ~ tidy(.))) %>%
mutate(glance_core = map(core, ~ glance(.))) %>% 
mutate(glance_core_stockholm = map(core_stockholm,
                                   ~ glance(.))) %>% 
  mutate(glance_stockholm = map(stockholm,
                                   ~ glance(.))) %>% 
mutate(glance_core_adams = map(core_adams, ~ glance(.))) %>% 
mutate(glance_core_hamdeh = map(core_hamdeh, ~ glance(.))) %>% 
mutate(glance_core_firsching = map(core_firsching, ~ glance(.))) %>% 
mutate(glance_core_karolinska = map(core_karolinska, ~ glance(.))) %>%
mutate(glance_adams = map(adams, ~ glance(.))) %>%
  mutate(glance_hamdeh = map(hamdeh, ~ glance(.))) %>%
   mutate(glance_firsching = map(firsching, ~ glance(.))) %>%
   mutate(aug_core = map2(core, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_stockholm = map2(core_stockholm, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
  mutate(aug_stockholm = map2(stockholm, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_adams = map2(core_adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_hamdeh = map2(core_hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_firsching = map2(core_firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_karolinska = map2(core_karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
mutate(aug_adams = map2(adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh = map2(hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_firsching = map2(firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(- c(adams, core_adams, core_stockholm, 
             core_karolinska,
       core, hamdeh, firsching, karolinska, core_firsching, core_hamdeh, stockholm)) %>% 
    mutate(karolinska_only = map(data, ~ glm(formulas_univariate[2], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(adams_only = map(data, ~ glm(formulas_univariate[1], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(hamdeh_only = map(data, ~ glm(formulas_univariate[3], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(firsching_only = map(data, ~ glm(formulas_univariate[4], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(r_sq_karolinska_only = pmap_dbl(list(data, karolinska_only, null), r2_calc)) %>%
  mutate(r_sq_adams_only = pmap_dbl(list(data, adams_only, null), r2_calc)) %>%
  mutate(r_sq_hamdeh_only = pmap_dbl(list(data, hamdeh_only, null), r2_calc)) %>%
  mutate(r_sq_firsching_only = pmap_dbl(list(data, firsching_only, null), r2_calc)) %>%
  mutate(glance_hamdeh_only = map(hamdeh_only, ~ glance(.))) %>% 
  mutate(glance_firsching_only = map(firsching_only, ~ glance(.))) %>% 
   mutate(glance_karolinska_only = map(karolinska_only,
                                       ~ glance(.))) %>%
  mutate(glance_adams_only = map(adams_only, ~ glance(.))) %>% 
  mutate(aug_firsching_only = map2(firsching_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh_only = map2(hamdeh_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_adams_only = map2(adams_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(-c(null, karolinska_only, hamdeh_only, 
        adams_only, firsching_only)) %>% 
     mutate_at(vars(contains("aug")), funs(map(., ~ roc(.x,dich_gos,
                                              .fitted)))) %>% 
  mutate_at(vars(contains("aug")), funs(map_dbl(., ~ auc(.x))))

names(estimates_stockholm) <- gsub("aug", "auc", names(estimates_stockholm))

model_vars_stockholm <- estimates_stockholm %>%  ungroup() %>% 
  sample_n(size =  1, replace = F) %>% 
   unnest(tid, .drop = T) %>% 
      pull(term) %>% 
     as_vector()
```

``` r
core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska, AIC, 
                "AUC" = auc) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + Stockholm CT-score + Stockholm MRI grading system") 

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh,AIC, 
                "AUC" = auc_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm CT-score + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching,AIC, 
                "AUC" = auc_firsching) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm CT-score + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams, AIC, "AUC" = auc_adams) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + Stockholm CT-score + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_core_karolinska) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_karolinska,AIC, 
                "AUC" = auc_core_karolinska) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm MRI grading system") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_core_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_hamdeh,AIC, 
                "AUC" = auc_core_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_core_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_firsching,AIC, 
                "AUC" = auc_core_firsching) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_core_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_adams, AIC, 
                "AUC" = auc_core_adams) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_core_stockholm) %>%
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_stockholm, AIC, 
                "AUC" = auc_core_stockholm) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm CT-score") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_karolinska_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska_only, AIC, 
                "AUC" = auc_karolinska_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Stockholm MRI grading system") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_hamdeh_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh_only, AIC, 
                "AUC" = auc_hamdeh_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_firsching_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching_only, AIC, 
                "AUC" = auc_firsching_only) %>% 
summarize_all(funs(mean)) %>%
   mutate(Model = "MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_adams_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams_only, AIC, 
                "AUC" = auc_adams_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_stockholm) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_stockholm, AIC, 
                "AUC" = auc_stockholm) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Stockholm CT-score") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())

core_comparison_stockholm <- estimates_stockholm %>%
  unnest(glance_core) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core, AIC, 
                "AUC" = auc_core) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core") %>% 
   bind_rows(core_comparison_stockholm) %>% 
   select(Model, everything())
  

core_comparison_stockholm %>% 
  kable(format = "latex", escape = FALSE, booktabs = T, linesep = "",
        caption = "Average result for all imputed datasets",
        digits = 2, align = 'c', col.names = c("Model","Pseudo-R<sup>2</sup>", 
                                               "AIC",
                                                "AUC")) %>%
   kable_styling(full_width = F)

results_table_Stockholm_orig <- core_comparison_stockholm %>%
  kable(format = "html", booktabs = T, linesep = "",
  escape = FALSE,
        caption = "Average result for all imputed datasets",
        digits = c(0,2,0,2), align = 'c', col.names = c("Model","Pseudo-R<sup>2</sup>",
        "AIC",
                                                "AUC")) %>%
   kable_styling(full_width = F)


save_kable(results_table_Stockholm_orig, "../Tables/Stockholm_results.html")
```

Analysis using the Helsinki CT-score
====================================

fitting the models

``` r
core_model <- c("Alder", "GCSKS", "Pupiller")

estimates_helsinki <- imputed_log %>%
mutate(adams = map(data, ~ glm(formulas_helsinki[1], 
                               family = binomial(link="logit"), 
                               data = .))) %>%
  mutate(hamdeh = map(data, ~ glm(formulas_helsinki[3], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
mutate(firsching = map(data, ~ glm(formulas_helsinki[4], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
  mutate(karolinska = map(data, ~ glm(formulas_helsinki[2], 
                                      family = binomial(link="logit"), data = .))) %>%
  mutate(null = map(data, ~ glm(dich_gos ~ 1, 
                              family = binomial(link="logit"), 
                              data = .))) %>%
    mutate(helsinki = map(data, ~ glm(dich_gos ~ helsinki, 
                              family = binomial(link="logit"), 
                              data = .))) %>% 
   mutate(core = map(data, ~ glm(core_vars_formula, family = binomial(link="logit"), data = .))) %>%
   mutate(core_helsinki = map(data, ~ glm(core_helsinki, family = binomial(link="logit"), data = .))) %>%
   mutate(core_adams = map(data, ~ glm(core_adams, family = binomial(link="logit"), data = .))) %>%
   mutate(core_karolinska = map(data, ~ glm(core_karolinska, family = binomial(link="logit"), data = .))) %>%
   mutate(core_firsching = map(data, ~ glm(core_firsching, family = binomial(link="logit"), data = .))) %>%
   mutate(core_hamdeh = map(data, ~ glm(core_hamdeh, family = binomial(link="logit"), data = .))) %>%
   mutate(r_sq_karolinska = pmap_dbl(list(data, karolinska, null), r2_calc)) %>%
  mutate(r_sq_hamdeh = pmap_dbl(list(data, hamdeh, null), 
                                r2_calc)) %>%   
  mutate(r_sq_firsching = pmap_dbl(list(data, firsching, null), 
                                 r2_calc)) %>% 
   mutate(r_sq_core = pmap_dbl(list(data, core, null), r2_calc)) %>%
  mutate(r_sq_helsinki = pmap_dbl(list(data, helsinki, null), r2_calc)) %>%
  mutate(r_sq_core_helsinki = pmap_dbl(list(data, core_helsinki, null), r2_calc)) %>%
  mutate(r_sq_core_adams = pmap_dbl(list(data, core_adams, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_firsching = pmap_dbl(list(data, core_firsching, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_hamdeh = pmap_dbl(list(data, core_hamdeh, null), 
                                    r2_calc)) %>%
  mutate(r_sq_core_karolinska = pmap_dbl(list(data, core_karolinska, null), r2_calc)) %>%
  mutate(r_sq_adams = pmap_dbl(list(data, adams, null), r2_calc)) %>%
mutate(glance = map(karolinska, ~ glance(.))) %>%
mutate(aug = map2(karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(tid = map(karolinska, ~ tidy(.))) %>%
mutate(glance_core = map(core, ~ glance(.))) %>% 
mutate(glance_core_helsinki = map(core_helsinki,
                                   ~ glance(.))) %>% 
  mutate(glance_helsinki = map(helsinki,
                                   ~ glance(.))) %>% 
mutate(glance_core_adams = map(core_adams, ~ glance(.))) %>% 
mutate(glance_core_hamdeh = map(core_hamdeh, ~ glance(.))) %>% 
mutate(glance_core_firsching = map(core_firsching, ~ glance(.))) %>% 
mutate(glance_core_karolinska = map(core_karolinska, ~ glance(.))) %>%
mutate(glance_adams = map(adams, ~ glance(.))) %>%
  mutate(glance_hamdeh = map(hamdeh, ~ glance(.))) %>%
   mutate(glance_firsching = map(firsching, ~ glance(.))) %>%
   mutate(aug_core = map2(core, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_helsinki = map2(core_helsinki, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
  mutate(aug_helsinki = map2(helsinki, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_adams = map2(core_adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_hamdeh = map2(core_hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_firsching = map2(core_firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_karolinska = map2(core_karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
mutate(aug_adams = map2(adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh = map2(hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_firsching = map2(firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(- c(adams, core_adams, core_helsinki, 
             core_karolinska,
       core, hamdeh, firsching, karolinska, core_firsching, core_hamdeh, helsinki)) %>% 
    mutate(karolinska_only = map(data, ~ glm(formulas_univariate[2], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(adams_only = map(data, ~ glm(formulas_univariate[1], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(hamdeh_only = map(data, ~ glm(formulas_univariate[3], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(firsching_only = map(data, ~ glm(formulas_univariate[4], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(r_sq_karolinska_only = pmap_dbl(list(data, karolinska_only, null), r2_calc)) %>%
  mutate(r_sq_adams_only = pmap_dbl(list(data, adams_only, null), r2_calc)) %>%
  mutate(r_sq_hamdeh_only = pmap_dbl(list(data, hamdeh_only, null), r2_calc)) %>%
  mutate(r_sq_firsching_only = pmap_dbl(list(data, firsching_only, null), r2_calc)) %>%
  mutate(glance_hamdeh_only = map(hamdeh_only, ~ glance(.))) %>% 
  mutate(glance_firsching_only = map(firsching_only, ~ glance(.))) %>% 
   mutate(glance_karolinska_only = map(karolinska_only,
                                       ~ glance(.))) %>%
  mutate(glance_adams_only = map(adams_only, ~ glance(.))) %>% 
  mutate(aug_firsching_only = map2(firsching_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh_only = map2(hamdeh_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_adams_only = map2(adams_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(-c(null, karolinska_only, hamdeh_only, 
        adams_only, firsching_only)) %>% 
     mutate_at(vars(contains("aug")), funs(map(., ~ roc(.x,dich_gos,
                                              .fitted)))) %>% 
  mutate_at(vars(contains("aug")), funs(map_dbl(., ~ auc(.x))))

names(estimates_helsinki) <- gsub("aug", "auc", names(estimates_helsinki))

model_vars_helsinki <- estimates_helsinki %>%  ungroup() %>% 
  sample_n(size =  1, replace = F) %>% 
   unnest(tid, .drop = T) %>% 
      pull(term) %>% 
     as_vector()
```

``` r
core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska, AIC, 
                "AUC" = auc) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + Helsinki CT-score + Stockholm MRI grading system") 

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh,AIC, 
                "AUC" = auc_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Helsinki CT-score + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching,AIC, 
                "AUC" = auc_firsching) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Helsinki CT-score + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams, AIC, "AUC" = auc_adams) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + Helsinki CT-score + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_core_karolinska) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_karolinska,AIC, 
                "AUC" = auc_core_karolinska) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm MRI grading system") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_core_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_hamdeh,AIC, 
                "AUC" = auc_core_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_core_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_firsching,AIC, 
                "AUC" = auc_core_firsching) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_core_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_adams, AIC, 
                "AUC" = auc_core_adams) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_core_helsinki) %>%
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_helsinki, AIC, 
                "AUC" = auc_core_helsinki) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Helsinki CT-score") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_karolinska_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska_only, AIC, 
                "AUC" = auc_karolinska_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Stockholm MRI grading system") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_hamdeh_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh_only, AIC, 
                "AUC" = auc_hamdeh_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_firsching_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching_only, AIC, 
                "AUC" = auc_firsching_only) %>% 
summarize_all(funs(mean)) %>%
   mutate(Model = "MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_adams_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams_only, AIC, 
                "AUC" = auc_adams_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_helsinki) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_helsinki, AIC, 
                "AUC" = auc_helsinki) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Helsinki CT-score") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())

core_comparison_helsinki <- estimates_helsinki %>%
  unnest(glance_core) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core, AIC, 
                "AUC" = auc_core) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core") %>% 
   bind_rows(core_comparison_helsinki) %>% 
   select(Model, everything())
  

core_comparison_helsinki %>% 
  kable(format = "latex", escape = FALSE, booktabs = T, linesep = "",
        caption = "Average result for all imputed datasets",
        digits = 2, align = 'c', col.names = c("Model","Pseudo-R<sup>2</sup>", 
                                               "AIC",
                                                "AUC")) %>%
   kable_styling(full_width = F)

results_table_Stockholm_orig <- core_comparison_helsinki %>%
  kable(format = "html", booktabs = T, linesep = "",
  escape = FALSE,
        caption = "Average result for all imputed datasets",
        digits = c(0,2,0,2), align = 'c', col.names = c("Model","Pseudo-R<sup>2</sup>",
        "AIC",
                                                "AUC")) %>%
   kable_styling(full_width = F)


# save_kable(results_table_Stockholm_orig, "../Tables/Helsinki_results.html")
```

Analysis using the Marshall CT-score
====================================

fitting the models

``` r
core_model <- c("Alder", "GCSKS", "Pupiller")

estimates_marshall <- imputed_log %>%
mutate(adams = map(data, ~ glm(formulas_marshall[1], 
                               family = binomial(link="logit"), 
                               data = .))) %>%
  mutate(hamdeh = map(data, ~ glm(formulas_marshall[3], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
mutate(firsching = map(data, ~ glm(formulas_marshall[4], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
  mutate(karolinska = map(data, ~ glm(formulas_marshall[2], 
                                      family = binomial(link="logit"), data = .))) %>%
  mutate(null = map(data, ~ glm(dich_gos ~ 1, 
                              family = binomial(link="logit"), 
                              data = .))) %>%
    mutate(marshall = map(data, ~ glm(dich_gos ~ marshall, 
                              family = binomial(link="logit"), 
                              data = .))) %>% 
   mutate(core = map(data, ~ glm(core_vars_formula, family = binomial(link="logit"), data = .))) %>%
   mutate(core_marshall = map(data, ~ glm(core_marshall, family = binomial(link="logit"), data = .))) %>%
   mutate(core_adams = map(data, ~ glm(core_adams, family = binomial(link="logit"), data = .))) %>%
   mutate(core_karolinska = map(data, ~ glm(core_karolinska, family = binomial(link="logit"), data = .))) %>%
   mutate(core_firsching = map(data, ~ glm(core_firsching, family = binomial(link="logit"), data = .))) %>%
   mutate(core_hamdeh = map(data, ~ glm(core_hamdeh, family = binomial(link="logit"), data = .))) %>%
   mutate(r_sq_karolinska = pmap_dbl(list(data, karolinska, null), r2_calc)) %>%
  mutate(r_sq_hamdeh = pmap_dbl(list(data, hamdeh, null), 
                                r2_calc)) %>%   
  mutate(r_sq_firsching = pmap_dbl(list(data, firsching, null), 
                                 r2_calc)) %>% 
   mutate(r_sq_core = pmap_dbl(list(data, core, null), r2_calc)) %>%
  mutate(r_sq_marshall = pmap_dbl(list(data, marshall, null), r2_calc)) %>%
  mutate(r_sq_core_marshall = pmap_dbl(list(data, core_marshall, null), r2_calc)) %>%
  mutate(r_sq_core_adams = pmap_dbl(list(data, core_adams, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_firsching = pmap_dbl(list(data, core_firsching, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_hamdeh = pmap_dbl(list(data, core_hamdeh, null), 
                                    r2_calc)) %>%
  mutate(r_sq_core_karolinska = pmap_dbl(list(data, core_karolinska, null), r2_calc)) %>%
  mutate(r_sq_adams = pmap_dbl(list(data, adams, null), r2_calc)) %>%
mutate(glance = map(karolinska, ~ glance(.))) %>%
mutate(aug = map2(karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(tid = map(karolinska, ~ tidy(.))) %>%
mutate(glance_core = map(core, ~ glance(.))) %>% 
mutate(glance_core_marshall = map(core_marshall,
                                   ~ glance(.))) %>% 
  mutate(glance_marshall = map(marshall,
                                   ~ glance(.))) %>% 
mutate(glance_core_adams = map(core_adams, ~ glance(.))) %>% 
mutate(glance_core_hamdeh = map(core_hamdeh, ~ glance(.))) %>% 
mutate(glance_core_firsching = map(core_firsching, ~ glance(.))) %>% 
mutate(glance_core_karolinska = map(core_karolinska, ~ glance(.))) %>%
mutate(glance_adams = map(adams, ~ glance(.))) %>%
  mutate(glance_hamdeh = map(hamdeh, ~ glance(.))) %>%
   mutate(glance_firsching = map(firsching, ~ glance(.))) %>%
   mutate(aug_core = map2(core, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_marshall = map2(core_marshall, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
  mutate(aug_marshall = map2(marshall, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_adams = map2(core_adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_hamdeh = map2(core_hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_firsching = map2(core_firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_karolinska = map2(core_karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
mutate(aug_adams = map2(adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh = map2(hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_firsching = map2(firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(- c(adams, core_adams, core_marshall, 
             core_karolinska,
       core, hamdeh, firsching, karolinska, core_firsching, core_hamdeh, marshall)) %>% 
    mutate(karolinska_only = map(data, ~ glm(formulas_univariate[2], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(adams_only = map(data, ~ glm(formulas_univariate[1], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(hamdeh_only = map(data, ~ glm(formulas_univariate[3], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(firsching_only = map(data, ~ glm(formulas_univariate[4], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(r_sq_karolinska_only = pmap_dbl(list(data, karolinska_only, null), r2_calc)) %>%
  mutate(r_sq_adams_only = pmap_dbl(list(data, adams_only, null), r2_calc)) %>%
  mutate(r_sq_hamdeh_only = pmap_dbl(list(data, hamdeh_only, null), r2_calc)) %>%
  mutate(r_sq_firsching_only = pmap_dbl(list(data, firsching_only, null), r2_calc)) %>%
  mutate(glance_hamdeh_only = map(hamdeh_only, ~ glance(.))) %>% 
  mutate(glance_firsching_only = map(firsching_only, ~ glance(.))) %>% 
   mutate(glance_karolinska_only = map(karolinska_only,
                                       ~ glance(.))) %>%
  mutate(glance_adams_only = map(adams_only, ~ glance(.))) %>% 
  mutate(aug_firsching_only = map2(firsching_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh_only = map2(hamdeh_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_adams_only = map2(adams_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(-c(null, karolinska_only, hamdeh_only, 
        adams_only, firsching_only)) %>% 
     mutate_at(vars(contains("aug")), funs(map(., ~ roc(.x,dich_gos,
                                              .fitted)))) %>% 
  mutate_at(vars(contains("aug")), funs(map_dbl(., ~ auc(.x))))

names(estimates_marshall) <- gsub("aug", "auc", names(estimates_marshall))

model_vars_marshall <- estimates_marshall %>%  ungroup() %>% 
  sample_n(size =  1, replace = F) %>% 
   unnest(tid, .drop = T) %>% 
      pull(term) %>% 
     as_vector()
```

``` r
core_comparison_marshall <- estimates_marshall %>%
  unnest(glance) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska, AIC, 
                "AUC" = auc) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + Marshall CT classification + Stockholm MRI grading system") 

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh,AIC, 
                "AUC" = auc_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Marshall CT classification + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching,AIC, 
                "AUC" = auc_firsching) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Marshall CT classification + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams, AIC, "AUC" = auc_adams) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + Marshall CT classification + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_core_karolinska) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_karolinska,AIC, 
                "AUC" = auc_core_karolinska) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm MRI grading system") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_core_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_hamdeh,AIC, 
                "AUC" = auc_core_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_core_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_firsching,AIC, 
                "AUC" = auc_core_firsching) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_core_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_adams, AIC, 
                "AUC" = auc_core_adams) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_core_marshall) %>%
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_marshall, AIC, 
                "AUC" = auc_core_marshall) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Marshall CT classification") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_karolinska_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska_only, AIC, 
                "AUC" = auc_karolinska_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Stockholm MRI grading system") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_hamdeh_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh_only, AIC, 
                "AUC" = auc_hamdeh_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_firsching_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching_only, AIC, 
                "AUC" = auc_firsching_only) %>% 
summarize_all(funs(mean)) %>%
   mutate(Model = "MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_adams_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams_only, AIC, 
                "AUC" = auc_adams_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_marshall) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_marshall, AIC, 
                "AUC" = auc_marshall) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Marshall CT classification") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())

core_comparison_marshall <- estimates_marshall %>%
  unnest(glance_core) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core, AIC, 
                "AUC" = auc_core) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core") %>% 
   bind_rows(core_comparison_marshall) %>% 
   select(Model, everything())
  

core_comparison_marshall %>% 
  kable(format = "latex", escape = FALSE, booktabs = T, linesep = "",
        caption = "Average result for all imputed datasets",
        digits = 2, align = 'c', col.names = c("Model","Pseudo-R<sup>2</sup>", 
                                               "AIC",
                                                "AUC")) %>%
   kable_styling(full_width = F)

results_table_Stockholm_orig <- core_comparison_marshall %>%
  kable(format = "html", booktabs = T, linesep = "",
  escape = FALSE,
        caption = "Average result for all imputed datasets",
        digits = c(0,2,0,2), align = 'c', col.names = c("Model","Pseudo-R<sup>2</sup>",
        "AIC",
                                                "AUC")) %>%
   kable_styling(full_width = F)


save_kable(results_table_Stockholm_orig, "../Tables/Marshall_results.html")
```

Mortality
=========

``` r
core_model <- c("Alder", "GCSKS", "Pupiller")

core_vars_formula_mortality <- paste("mortality ~", paste(core_model, 
                                               collapse = " + "),
                            sep = " " ) %>% 
   as.formula()

core_rotterdam_mortality <-  paste("mortality ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "rotterdam") %>% 
   as.formula()


core_adams_mortality <-  paste("mortality ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "adams") %>% 
   as.formula()

core_karolinska_mortality <-  paste("mortality ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "stockholm_mri") %>% 
   as.formula()

core_firsching_mortality <-  paste("mortality ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "firsching") %>% 
   as.formula()

core_hamdeh_mortality <-  paste("mortality ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "hamdeh") %>% 
   as.formula()

formulas_mortality <- paste('mortality ~ ', paste(core_model, 
                                                 collapse = " + "),
                            sep = " " , "+", "rotterdam", "+",var_names)


formulas_univariate_mortality <- paste('mortality ~ ', var_names)
```

``` r
core_model <- c("Alder", "GCSKS", "Pupiller")

core_vars_formula_mortality <- paste("mortality ~", paste(core_model, 
                                               collapse = " + "),
                            sep = " " ) %>% 
   as.formula()

core_rotterdam_mortality <-  paste("mortality ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "stockholm") %>% 
   as.formula()


core_adams_mortality <-  paste("mortality ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "adams") %>% 
   as.formula()

core_karolinska_mortality <-  paste("mortality ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "stockholm_mri") %>% 
   as.formula()

core_firsching_mortality <-  paste("mortality ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "firsching") %>% 
   as.formula()

core_hamdeh_mortality <-  paste("mortality ~", paste(core_model, 
                                             collapse = " + "),
                            sep = " " , "+", "hamdeh") %>% 
   as.formula()

formulas_mortality <- paste('mortality ~ ', paste(core_model, 
                                                 collapse = " + "),
                            sep = " " , "+", "stockholm", "+",var_names)


formulas_univariate_mortality <- paste('mortality ~ ', var_names)
```

fitting the models

``` r
estimates_mortality <- imputed_log %>%
mutate(adams = map(data, ~ glm(formulas_mortality[1], 
                               family = binomial(link="logit"), 
                               data = .))) %>%
  mutate(hamdeh = map(data, ~ glm(formulas_mortality[3], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
mutate(firsching = map(data, ~ glm(formulas_mortality[4], 
                               family = binomial(link="logit"), 
                               data = .))) %>% 
  mutate(karolinska = map(data, ~ glm(formulas_mortality[2], 
                                      family = binomial(link="logit"), data = .))) %>%
  mutate(null = map(data, ~ glm(mortality ~ 1, 
                              family = binomial(link="logit"), 
                              data = .))) %>%
    mutate(rotterdam = map(data, ~ glm(mortality ~ rotterdam, 
                              family = binomial(link="logit"), 
                              data = .))) %>% 
   mutate(core = map(data, ~ glm(core_vars_formula_mortality, family = binomial(link="logit"), data = .))) %>%
   mutate(core_rotterdam = map(data, ~ glm(core_rotterdam_mortality, family = binomial(link="logit"), data = .))) %>%
   mutate(core_adams = map(data, ~ glm(core_adams_mortality, family = binomial(link="logit"), data = .))) %>%
   mutate(core_karolinska = map(data, ~ glm(core_karolinska_mortality, family = binomial(link="logit"), data = .))) %>%
   mutate(core_firsching = map(data, ~ glm(core_firsching_mortality, family = binomial(link="logit"), data = .))) %>%
   mutate(core_hamdeh = map(data, ~ glm(core_hamdeh_mortality, family = binomial(link="logit"), data = .))) %>%
   mutate(r_sq_karolinska = pmap_dbl(list(data, karolinska, null), r2_calc)) %>%
  mutate(r_sq_hamdeh = pmap_dbl(list(data, hamdeh, null), 
                                r2_calc)) %>%   
  mutate(r_sq_firsching = pmap_dbl(list(data, firsching, null), 
                                 r2_calc)) %>% 
   mutate(r_sq_core = pmap_dbl(list(data, core, null), r2_calc)) %>%
  mutate(r_sq_rotterdam = pmap_dbl(list(data, rotterdam, null), r2_calc)) %>%
  mutate(r_sq_core_rotterdam = pmap_dbl(list(data, core_rotterdam, null), r2_calc)) %>%
  mutate(r_sq_core_adams = pmap_dbl(list(data, core_adams, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_firsching = pmap_dbl(list(data, core_firsching, null), 
                                    r2_calc)) %>%
   mutate(r_sq_core_hamdeh = pmap_dbl(list(data, core_hamdeh, null), 
                                    r2_calc)) %>%
  mutate(r_sq_core_karolinska = pmap_dbl(list(data, core_karolinska, null), r2_calc)) %>%
  mutate(r_sq_adams = pmap_dbl(list(data, adams, null), r2_calc)) %>%
mutate(glance = map(karolinska, ~ glance(.))) %>%
mutate(aug = map2(karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(tid = map(karolinska, ~ tidy(.))) %>%
mutate(glance_core = map(core, ~ glance(.))) %>% 
mutate(glance_core_rotterdam = map(core_rotterdam,
                                   ~ glance(.))) %>% 
  mutate(glance_rotterdam = map(rotterdam,
                                   ~ glance(.))) %>% 
mutate(glance_core_adams = map(core_adams, ~ glance(.))) %>% 
mutate(glance_core_hamdeh = map(core_hamdeh, ~ glance(.))) %>% 
mutate(glance_core_firsching = map(core_firsching, ~ glance(.))) %>% 
mutate(glance_core_karolinska = map(core_karolinska, ~ glance(.))) %>%
mutate(glance_adams = map(adams, ~ glance(.))) %>%
  mutate(glance_hamdeh = map(hamdeh, ~ glance(.))) %>%
   mutate(glance_firsching = map(firsching, ~ glance(.))) %>%
   mutate(aug_core = map2(core, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_rotterdam = map2(core_rotterdam, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
  mutate(aug_rotterdam = map2(rotterdam, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_adams = map2(core_adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_hamdeh = map2(core_hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_firsching = map2(core_firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>% 
mutate(aug_core_karolinska = map2(core_karolinska, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
mutate(aug_adams = map2(adams, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh = map2(hamdeh, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_firsching = map2(firsching, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(- c(adams, core_adams, core_rotterdam, 
             core_karolinska,
       core, hamdeh, firsching, karolinska, core_firsching, core_hamdeh, rotterdam)) %>% 
    mutate(karolinska_only = map(data, ~ glm(formulas_univariate_mortality[2], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(adams_only = map(data, ~ glm(formulas_univariate_mortality[1], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(hamdeh_only = map(data, ~ glm(formulas_univariate_mortality[3], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(firsching_only = map(data, ~ glm(formulas_univariate_mortality[4], 
                                      family = binomial(link="logit"), data = .))) %>% 
  mutate(r_sq_karolinska_only = pmap_dbl(list(data, karolinska_only, null), r2_calc)) %>%
  mutate(r_sq_adams_only = pmap_dbl(list(data, adams_only, null), r2_calc)) %>%
  mutate(r_sq_hamdeh_only = pmap_dbl(list(data, hamdeh_only, null), r2_calc)) %>%
  mutate(r_sq_firsching_only = pmap_dbl(list(data, firsching_only, null), r2_calc)) %>%
  mutate(glance_hamdeh_only = map(hamdeh_only, ~ glance(.))) %>% 
  mutate(glance_firsching_only = map(firsching_only, ~ glance(.))) %>% 
   mutate(glance_karolinska_only = map(karolinska_only,
                                       ~ glance(.))) %>%
  mutate(glance_adams_only = map(adams_only, ~ glance(.))) %>% 
  mutate(aug_firsching_only = map2(firsching_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_hamdeh_only = map2(hamdeh_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_adams_only = map2(adams_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  mutate(aug_karolinska_only = map2(karolinska_only, data, 
                  ~ augment(.x, newdata = .y,
                            type.predict = "response"))) %>%
  select(-c(null, karolinska_only, hamdeh_only, 
        adams_only, firsching_only)) %>% 
     mutate_at(vars(contains("aug")), funs(map(., ~ roc(.x,mortality,
                                              .fitted)))) %>% 
  mutate_at(vars(contains("aug")), funs(map_dbl(., ~ auc(.x))))

names(estimates_mortality) <- gsub("aug", "auc", names(estimates_mortality))

model_vars_mortality <- estimates_mortality %>%  ungroup() %>% 
  sample_n(size =  1, replace = F) %>% 
   unnest(tid, .drop = T) %>% 
      pull(term) %>% 
     as_vector()
```

``` r
core_comparison_mortality <- estimates_mortality %>%
  unnest(glance) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska, AIC, 
                "AUC" = auc) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + Rotterdam CT-score + Stockholm MRI grading system") 

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh,AIC, 
                "AUC" = auc_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Rotterdam CT-score + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching,AIC, 
                "AUC" = auc_firsching) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Rotterdam CT-score + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams, AIC, "AUC" = auc_adams) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + Rotterdam CT-score + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_core_karolinska) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_karolinska,AIC, 
                "AUC" = auc_core_karolinska) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Stockholm MRI grading system") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_core_hamdeh) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_hamdeh,AIC, 
                "AUC" = auc_core_hamdeh) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_core_firsching) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_firsching,AIC, 
                "AUC" = auc_core_firsching) %>% 
summarize_all(funs(mean)) %>% 
mutate(Model = "Core + MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_core_adams) %>% 
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_adams, AIC, 
                "AUC" = auc_core_adams) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_core_rotterdam) %>%
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core_rotterdam, AIC, 
                "AUC" = auc_core_rotterdam) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core + Rotterdam CT-score") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_karolinska_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_karolinska_only, AIC, 
                "AUC" = auc_karolinska_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Stockholm MRI grading system") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_hamdeh_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_hamdeh_only, AIC, 
                "AUC" = auc_hamdeh_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Abu Hamdeh et al.") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_firsching_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_firsching_only, AIC, 
                "AUC" = auc_firsching_only) %>% 
summarize_all(funs(mean)) %>%
   mutate(Model = "MRI grading system of Firsching et al.") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_adams_only) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_adams_only, AIC, 
                "AUC" = auc_adams_only) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "MRI grading system of Adams et al.") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_rotterdam) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_rotterdam, AIC, 
                "AUC" = auc_rotterdam) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Rotterdam CT-score") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())

core_comparison_mortality <- estimates_mortality %>%
  unnest(glance_core) %>%  
  ungroup() %>% 
  dplyr::select("r_sq" = r_sq_core, AIC, 
                "AUC" = auc_core) %>% 
summarize_all(funs(mean)) %>% 
   mutate(Model = "Core") %>% 
   bind_rows(core_comparison_mortality) %>% 
   select(Model, everything())
  

core_comparison_mortality %>%
  kable(format = "latex", escape = FALSE, booktabs = T, linesep = "",
        caption = "Average result for all imputed datasets",
        digits = c(0,2,0,2), align = 'c', col.names = c("Model","Pseudo-R<sup>2</sup>",
                                               "AIC",
                                                "AUC")) %>%
   kable_styling(full_width = F)

# results_table_mortality <- core_comparison_mortality %>%
#   kable(format = "html", booktabs = T, linesep = "",
#   escape = FALSE,
#         caption = "Average result for all imputed datasets",
#         digits = c(0,2,0,2), align = 'c', col.names = c("Model","Pseudo-R<sup>2</sup>",
#         "AIC",
#                                                 "AUC")) %>%
#    kable_styling(full_width = F)
# 
#  save_kable(results_table_mortality, "../Tables/mortality.html")
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
    ##  [1] cowplot_1.1.0      pROC_1.16.1        glmnet_3.0-2       Matrix_1.2-18     
    ##  [5] stringi_1.4.6      naniar_0.6.0       pscl_1.5.5         mice_3.8.0        
    ##  [9] lubridate_1.7.8    corrr_0.4.2        caTools_1.18.0     DescTools_0.99.38 
    ## [13] descr_1.1.4        modelgrid_1.1.1.0  readxl_1.3.1       corrplot_0.84     
    ## [17] summarytools_0.9.6 kableExtra_1.1.0   knitr_1.29         forcats_0.5.0     
    ## [21] stringr_1.4.0      readr_1.3.1        tidyverse_1.3.0    yardstick_0.0.7   
    ## [25] workflows_0.2.0    tune_0.1.1         tidyr_1.1.2        tibble_3.0.3      
    ## [29] rsample_0.0.8      recipes_0.1.13     purrr_0.3.4        parsnip_0.1.3     
    ## [33] modeldata_0.0.2    infer_0.5.3        ggplot2_3.3.2      dplyr_1.0.2       
    ## [37] dials_0.0.9        scales_1.1.0       broom_0.7.0        tidymodels_0.1.1  
    ## [41] data.table_1.13.0  MASS_7.3-51.5     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] backports_1.1.5      plyr_1.8.6           splines_3.6.2       
    ##  [4] listenv_0.8.0        pryr_0.1.4           digest_0.6.25       
    ##  [7] foreach_1.5.0        htmltools_0.4.0      magick_2.4.0        
    ## [10] fansi_0.4.1          magrittr_1.5         checkmate_2.0.0     
    ## [13] globals_0.13.0       modelr_0.1.8         gower_0.2.1         
    ## [16] matrixStats_0.55.0   colorspace_1.4-1     blob_1.2.1          
    ## [19] rvest_0.3.6          haven_2.2.0          xfun_0.16           
    ## [22] tcltk_3.6.2          crayon_1.3.4         jsonlite_1.6.1      
    ## [25] Exact_2.0            survival_3.2-3       iterators_1.0.12    
    ## [28] glue_1.4.2           gtable_0.3.0         ipred_0.9-9         
    ## [31] webshot_0.5.2        shape_1.4.5          rapportools_1.0     
    ## [34] mvtnorm_1.1-1        DBI_1.1.0            Rcpp_1.0.3          
    ## [37] viridisLite_0.3.0    xtable_1.8-4         GPfit_1.0-8         
    ## [40] stats4_3.6.2         lava_1.6.8           prodlim_2019.11.13  
    ## [43] httr_1.4.2           ellipsis_0.3.0       pkgconfig_2.0.3     
    ## [46] nnet_7.3-12          dbplyr_1.4.4         caret_6.0-86        
    ## [49] tidyselect_1.1.0     rlang_0.4.7          DiceDesign_1.8-1    
    ## [52] reshape2_1.4.3       munsell_0.5.0        cellranger_1.1.0    
    ## [55] tools_3.6.2          cli_2.0.2            generics_0.0.2      
    ## [58] evaluate_0.14        yaml_2.2.1           ModelMetrics_1.2.2.2
    ## [61] fs_1.3.1             pander_0.6.3         visdat_0.5.3        
    ## [64] future_1.19.1        nlme_3.1-142         xml2_1.2.2          
    ## [67] compiler_3.6.2       rstudioapi_0.11      e1071_1.7-3         
    ## [70] reprex_0.3.0         lhs_1.0.1            lattice_0.20-38     
    ## [73] vctrs_0.3.4          pillar_1.4.6         lifecycle_0.2.0     
    ## [76] furrr_0.1.0          bitops_1.0-6         lmom_2.8            
    ## [79] R6_2.4.1             gld_2.6.2            codetools_0.2-16    
    ## [82] boot_1.3-23          assertthat_0.2.1     withr_2.3.0         
    ## [85] expm_0.999-5         parallel_3.6.2       hms_0.5.3           
    ## [88] grid_3.6.2           rpart_4.1-15         timeDate_3043.102   
    ## [91] class_7.3-15         rmarkdown_2.3        base64enc_0.1-3
