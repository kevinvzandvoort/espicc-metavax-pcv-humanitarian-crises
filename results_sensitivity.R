#' load/update all libraries
pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve, magrittr, units, ggplot2, BayesianTools, ggplot2, ggh4x, lvplot, binom, tictoc, patchwork)

# Can setwd by clicking on Session > Set Working Directory > To Source File Location in the toolbar in RStudio

#' set arguments if R is running interactively, otherwise arguments should be passed to the Rscript command
#'  1: working directory of pcvmr (make sure to setwd manually if running interactively)
#'  2: working directory of pcvm
#'  3: name of output folder to store objects (e.g. date, or specific scenario name)
.args = if(interactive()) c(getwd(), "../pcvm") else commandArgs(trailingOnly = TRUE)
.args = setNames(.args, c("wd", "metavax_dir"))
setwd(.args["wd"])

#' Load metavax model and setup initial model_params skeleton
METAVAX_FOLDER = .args["metavax_dir"]
source(sprintf("%s/model/metavax_streppneumo_diamond.R", METAVAX_FOLDER))  

source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")
model_params = modelSetup()

vaccination_strategies_plot_settings = list(
  list(model_name = "<0", name = "No vaccination", target = T, dosing = F, colour = lshtm_colours$black),
  list(model_name = "<1", name = "<1y", target = T, dosing = T, colour = lshtm_colours$red),
  list(model_name = "<2", name = "<2y", target = T, dosing = T, colour = lshtm_colours$yellow),
  list(model_name = "<5", name = "<5y", target = T, dosing = T, colour = lshtm_colours$blue),
  list(model_name = "<10", name = "<10y", target = T, dosing = T, colour = lshtm_colours$green),
  list(model_name = "<15", name = "<15y", target = T, dosing = T, colour = lshtm_colours$lightgreen)) %>%
  lapply(as.data.table) %>% rbindlist()
vaccination_strategies_plot_settings[, name := factor(name, vaccination_strategies_plot_settings$name)]

population_size_all = getParameter(model_params, key = "N", level = "population", population = c("idp_non_malnourished", "idp_malnourished")) %>%
  merge(model_params$global_settings$age_groups_model, by = "age") %>%
  .[, .(N = sum(N)), by = c("age", "name", "from", "to")] %>%
  combineAgeBreaks(x = setAgeBreaks(0), y = ., method = "sum", value.var = "N")

population_size_infants = getParameter(model_params, key = "N", level = "population", population = c("idp_non_malnourished", "idp_malnourished")) %>%
  merge(model_params$global_settings$age_groups_model, by = "age") %>%
  .[, .(N = sum(N)), by = c("age", "name", "from", "to")] %>%
  combineAgeBreaks(x = setAgeBreaks(1), y = ., method = "sum", value.var = "N") %>%
  .[to <= set_units(1, "year")]

population_size = rbind(population_size_all, population_size_infants)

upper_age_limits = c(1, 2, 5, 10, 15)
vaccine_doses = lapply(upper_age_limits, function(upper_age_limit){
  getParameter(model_params, key = "N", level = "population", population = c("idp_non_malnourished", "idp_malnourished")) %>%
    merge(model_params$global_settings$age_groups_model, by = "age") %>%
    .[, .(N = sum(N)), by = c("age", "name", "from", "to")] %>%
    .[, coverage := getVaccineCoverage(model_params$global_settings$age_groups_model,
                                       c(set_units(6, "weeks"), set_units(upper_age_limit, "year")), 0.85)] %>%
    .[, .(vaccine_strategy = sprintf("<%s", upper_age_limit), doses = sum(N * coverage))]
}) %>% rbindlist()

sensitivity_analyses = list(mal_trans = list(values = c(2, 1.8, 1.6, 1.4, 1.2, 1),
                                             labels = c("2", "1.8", "1.6", "1.4", "* 1.2", "1"),
                                             main = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)),
                            migration_rate = list(values = c(0, 10, 8, 6, 4, 2),
                                                  labels = c("None", "10", "* 8", "6", "4", "2"),
                                                  main = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE)),
                            contacts_host = list(values = c("0.40", "0.20", "0.10", "0.05", "0.03", "0"),
                                                 labels = c("40%", "20%", "10%", "5%", "* 3%", "None"),
                                                 main = c(FALSE, FALSE, FALSE, FALSE, TRUE, FALSE)))
sensitivity_analyses_labels = c(mal_trans = "Increased acquisition\nmalnourished",
                                migration_rate = "Migration: duration\n of residence (y)",
                                contacts_host = "Prop extra-household\ncontacts with host pop")

getKeyStatistics = function(posterior_runs,
                            dt_population_size = population_size,
                            dt_vaccine_doses = vaccine_doses){
  prevalence = posterior_runs %>%
    .[outcome == "prevalence"] %>%
    dcast(...~vaccine_strategy, value.var = "value") %>%
    melt(measure.vars = vaccination_strategies_plot_settings[dosing == TRUE, model_name],
         variable.name = "vaccine_strategy", value.name = "value") %>%
    .[, .(relative_impact = 1 - value/`<0`), by = c("age_group", "time", "vaccine_strategy", "compartment", "run")]
  
  prevalence %>%
    .[, .SD[which(relative_impact == max(relative_impact)), .(max_impact_prevalence = relative_impact,
                                                              max_impact_prevalence_time = time)],
      by = c("age_group", "vaccine_strategy", "run")] %>%
    melt(measure.vars = c("max_impact_prevalence", "max_impact_prevalence_time")) %>%
    rbind(posterior_runs %>%
            .[outcome == "incidence"] %>%
            .[time <= 365 * 3] %>%
            .[, .(value = sum(value)), by = c("age_group", "vaccine_strategy", "run")] %>%
            dcast(...~vaccine_strategy, value.var = "value") %>%
            melt(measure.vars = vaccination_strategies_plot_settings[dosing == TRUE, model_name],
                 variable.name = "vaccine_strategy", value.name = "value") %>%
            .[, .(absolute_impact_3y = (`<0` - value) * 5.98,
                  relative_impact_3y = 1 - value/`<0`), by = c("age_group", "vaccine_strategy", "run")] %>%
            merge(dt_population_size[, c("name", "N")], by.y = "name", by.x = "age_group") %>%
            merge(dt_vaccine_doses, by = "vaccine_strategy") %>%
            .[, NNV_3y := doses/absolute_impact_3y] %>%
            .[, absolute_impact_p10000_3y := absolute_impact_3y/N * 10000] %>%
            .[, -c("N", "doses")] %>%
            melt(measure.vars = c("absolute_impact_3y", "relative_impact_3y", "NNV_3y", "absolute_impact_p10000_3y")), fill = TRUE)
}

#' include values used in main model run
posterior_runs_main = readRDS(sprintf("./output/digaale/simulations/posterior_runs_main_v2.RDS"))
posterior_runs_main = posterior_runs %>% getKeyStatistics() 

results = lapply(names(sensitivity_analyses), function(s){
  lapply(seq_along(sensitivity_analyses[[s]][["values"]]), function(i){
    value = sensitivity_analyses[[s]][["values"]][i]
    label = sensitivity_analyses[[s]][["labels"]][i]
    
    message(sprintf("%s - %s (%s)", s, value, label))
    
    if(sensitivity_analyses[[s]][["main"]][i]){
      posterior_runs = posterior_runs_main
    } else {
      posterior_runs = readRDS(sprintf("./output/sensitivity/%s/%s/output/simulations/posterior_runs_main_v2.RDS", s, value))
      
      if(s == "vac_cov"){
        upper_age_limits = c(1, 2, 5, 10, 15)
        vaccine_doses = lapply(upper_age_limits, function(upper_age_limit){
          getParameter(model_params, key = "N", level = "population", population = c("idp_non_malnourished", "idp_malnourished")) %>%
            merge(model_params$global_settings$age_groups_model, by = "age") %>%
            .[, .(N = sum(N)), by = c("age", "name", "from", "to")] %>%
            .[, coverage := getVaccineCoverage(model_params$global_settings$age_groups_model,
                                               c(set_units(6, "weeks"), set_units(upper_age_limit, "year")), as.numeric(v))] %>%
            .[, .(vaccine_strategy = sprintf("<%s", upper_age_limit), doses = sum(N * coverage))]
        }) %>% rbindlist()
      }
      posterior_runs = posterior_runs %>% getKeyStatistics(dt_vaccine_doses = vaccine_doses) 
    }
    
    posterior_runs %>%
      .[, c("sens_var", "sens_val", "sens_val_f") :=
          .(s, value, sprintf("%s-%s", s, label))] %>%
      .[]
  }) %>% rbindlist()
}) %>% rbindlist()

results = results[age_group == "[0y, 1y)" & variable == "max_impact_prevalence"] %>%
  .[, -c("age_group", "variable")] %>%
  .[, variable := c("max_impact_prevalence_infants")] %>%
  rbind(results[age_group == "[0y, 1y)" & variable == "max_impact_prevalence_time"] %>%
          .[, -c("age_group", "variable")] %>%
          .[, variable := c("max_impact_prevalence_time_infants")]) %>%
  rbind(results[age_group == "[0y, 1y)" & variable == "relative_impact_3y"] %>%
          .[, -c("age_group", "variable")] %>%
          .[, variable := c("max_impact_incidence_infants")]) %>%
  rbind(results[age_group == "[0y, 120y)" & variable == "relative_impact_3y"] %>%
          .[, -c("age_group", "variable")] %>%
          .[, variable := c("max_impact_incidence_all")]) %>%
  rbind(results[age_group == "[0y, 1y)" & variable == "NNV_3y"] %>%
          .[, -c("age_group", "variable")] %>%
          .[, variable := c("nnv_infants")]) %>%
  rbind(results[age_group == "[0y, 120y)" & variable == "NNV_3y"] %>%
          .[, -c("age_group", "variable")] %>%
          .[, variable := c("nnv_all")])

scales_y = list(
  scale_y_continuous(labels=scales::percent, lim=c(0, 0.6), n.breaks=4),
  scale_y_continuous(lim=c(0, 600), n.breaks=4),
  scale_y_continuous(labels=scales::percent, lim=c(0, NA), n.breaks=4),
  scale_y_continuous(labels=scales::percent, lim=c(0, 0.6), breaks = c(0, 0.2, 0.4, 0.6)),
  scale_y_log10(lim = c(200, 12000), breaks = c(300, 1000, 3000, 10000)),
  scale_y_log10(lim = c(75, 400), breaks = c(75, 150, 300)))

results %>%
  ggplot(aes(x = factor(sens_val_f, as.vector(sapply(names(sensitivity_analyses), function(s) sprintf("%s-%s", s, sensitivity_analyses[[s]][["labels"]])))),
             y = value,
             colour=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             fill=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name)))+
  facet_grid(factor(sens_var, rev(names(sensitivity_analyses_labels)), rev(sensitivity_analyses_labels))~factor(variable,
                             c("max_impact_prevalence_infants",
                               "max_impact_prevalence_time_infants",
                               "max_impact_incidence_infants",
                               "max_impact_incidence_all",
                               "nnv_infants",
                               "nnv_all"),
                             c("Max impact\nprevalence (infants)",
                               "Peak impact (days)\n prevalence infants",
                               "3y disease\nimpact (infants)",
                               "3y disease\nimpact (all ages)",
                               "3y NNV\n(infants)",
                               "3y NNV\n(all ages)")),
             scales = "free")+
  geom_lv(varwidth = TRUE,
          width.method = "height",
          aes(alpha = after_stat(LV)),
          outlier.colour = adjustColBrightness(lshtm_colours$grey, 0.25),
          outlier.shape = 19,
          outlier.size = 1,
          outlier.stroke = 0.25,
          percent=5)+
  scale_colour_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
  scale_fill_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
  scale_alpha_manual(values = c(1, seq(0.6, 0.1, length.out=5)))+
  guides(alpha = "none")+
  theme_minimal()+
  labs(x="Parameter value",
       y="Outcome value",
       colour="Strategy",
       fill="Strategy")+
  theme(legend.position="bottom",
        strip.text = element_text(face="bold", size=13),
        axis.title = element_text(size=12), axis.text = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size = 11),
        panel.border = element_rect(size=1, fill="#FFFFFF00"),
        panel.spacing.x = unit(c(3), "mm"), panel.spacing.y = unit(c(1), "mm"),
        panel.grid.minor = element_blank())+
  scale_x_discrete(labels = function(x) sapply(strsplit(x, "-"), "[[", 2))+
  coord_flip()+
  ggh4x::facetted_pos_scales(y = scales_y)
ggsave(sprintf("./output/sensitivity/parametric_sensitivity_crisis.png"),
       width = PLOT_WIDTH_ONE*1.1, height = PLOT_HEIGHT_ONE*1.375, units = "in", bg = "#FFFFFF")

sensitivity_analyses = list(vac_eff_dur = list(values = c(2, 4, 6, 8, 10),
                                               labels = c("2", "4", "* 6", "8", "10"),
                                               main = c(FALSE, FALSE, FALSE, FALSE, FALSE)),
                            vac_eff = list(values = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
                                           labels = c("0%", "10%", "20%", "30%", "40%", "* 50%", "60%"),
                                           main = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)),
                            vac_cov = list(values = c("0.45", "0.55", "0.65", "0.75", "0.85", "0.95"),
                                           labels = c("45%", "55%", "65%", "75%", "* 85%", "95%"),
                                           main = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)))
sensitivity_analyses_labels = c(vac_eff_dur = "Duration vaccine\nprotection (y)",
                                vac_eff = "VE (transmission)",
                                vac_cov = "Vaccine coverage")

results = lapply(names(sensitivity_analyses), function(s){
  lapply(seq_along(sensitivity_analyses[[s]][["values"]]), function(i){
    value = sensitivity_analyses[[s]][["values"]][i]
    label = sensitivity_analyses[[s]][["labels"]][i]
    
    message(sprintf("%s - %s (%s)", s, value, label))
    
    if(sensitivity_analyses[[s]][["main"]][i]){
      posterior_runs = posterior_runs_main
    } else {
      posterior_runs = readRDS(sprintf("./output/sensitivity/%s/%s/output/simulations/posterior_runs_main_v2.RDS", s, value))
      
      if(s == "vac_cov"){
        upper_age_limits = c(1, 2, 5, 10, 15)
        vaccine_doses = lapply(upper_age_limits, function(upper_age_limit){
          getParameter(model_params, key = "N", level = "population", population = c("idp_non_malnourished", "idp_malnourished")) %>%
            merge(model_params$global_settings$age_groups_model, by = "age") %>%
            .[, .(N = sum(N)), by = c("age", "name", "from", "to")] %>%
            .[, coverage := getVaccineCoverage(model_params$global_settings$age_groups_model,
                                               c(set_units(6, "weeks"), set_units(upper_age_limit, "year")), as.numeric(value))] %>%
            .[, .(vaccine_strategy = sprintf("<%s", upper_age_limit), doses = sum(N * coverage))]
        }) %>% rbindlist()
      }
      posterior_runs = posterior_runs %>% getKeyStatistics(dt_vaccine_doses = vaccine_doses) 
    }
    
    posterior_runs %>%
      .[, c("sens_var", "sens_val", "sens_val_f") :=
          .(s, value, sprintf("%s-%s", s, label))] %>%
      .[]
  }) %>% rbindlist()
}) %>% rbindlist()

results = results[age_group == "[0y, 1y)" & variable == "max_impact_prevalence"] %>%
  .[, -c("age_group", "variable")] %>%
  .[, variable := c("max_impact_prevalence_infants")] %>%
  rbind(results[age_group == "[0y, 1y)" & variable == "max_impact_prevalence_time"] %>%
          .[, -c("age_group", "variable")] %>%
          .[, variable := c("max_impact_prevalence_time_infants")]) %>%
  rbind(results[age_group == "[0y, 1y)" & variable == "relative_impact_3y"] %>%
          .[, -c("age_group", "variable")] %>%
          .[, variable := c("max_impact_incidence_infants")]) %>%
  rbind(results[age_group == "[0y, 120y)" & variable == "relative_impact_3y"] %>%
          .[, -c("age_group", "variable")] %>%
          .[, variable := c("max_impact_incidence_all")]) %>%
  rbind(results[age_group == "[0y, 1y)" & variable == "NNV_3y"] %>%
          .[, -c("age_group", "variable")] %>%
          .[, variable := c("nnv_infants")]) %>%
  rbind(results[age_group == "[0y, 120y)" & variable == "NNV_3y"] %>%
          .[, -c("age_group", "variable")] %>%
          .[, variable := c("nnv_all")])

results %>%
  ggplot(aes(x = factor(sens_val_f, unlist(sapply(names(sensitivity_analyses), function(s) sprintf("%s-%s", s, sensitivity_analyses[[s]][["labels"]])))),
             y = value,
             colour=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             fill=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name)))+
  facet_grid(factor(sens_var, rev(names(sensitivity_analyses_labels)), rev(sensitivity_analyses_labels))~factor(variable,
                                                                                                                c("max_impact_prevalence_infants",
                                                                                                                  "max_impact_prevalence_time_infants",
                                                                                                                  "max_impact_incidence_infants",
                                                                                                                  "max_impact_incidence_all",
                                                                                                                  "nnv_infants",
                                                                                                                  "nnv_all"),
                                                                                                                c("Max impact\nprevalence (infants)",
                                                                                                                  "Peak impact (days)\n prevalence infants",
                                                                                                                  "3y disease\nimpact (infants)",
                                                                                                                  "3y disease\nimpact (all ages)",
                                                                                                                  "3y NNV\n(infants)",
                                                                                                                  "3y NNV\n(all ages)")),
             scales = "free")+
  geom_lv(varwidth = TRUE,
          width.method = "height",
          aes(alpha = after_stat(LV)),
          outlier.colour = adjustColBrightness(lshtm_colours$grey, 0.25),
          outlier.shape = 19,
          outlier.size = 1,
          outlier.stroke = 0.25,
          percent=5)+
  scale_colour_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
  scale_fill_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
  scale_alpha_manual(values = c(1, seq(0.6, 0.1, length.out=5)))+
  guides(alpha = "none")+
  theme_minimal()+
  labs(x="Parameter value",
       y="Outcome value",
       colour="Strategy",
       fill="Strategy")+
  theme(legend.position="bottom",
        strip.text = element_text(face="bold", size=13),
        axis.title = element_text(size=12), axis.text = element_text(size = 11),
        legend.title = element_text(size = 12), legend.text = element_text(size = 11),
        panel.border = element_rect(size=1, fill="#FFFFFF00"),
        panel.spacing.x = unit(c(3), "mm"), panel.spacing.y = unit(c(1), "mm"),
        panel.grid.minor = element_blank())+
  scale_x_discrete(labels = function(x) sapply(strsplit(x, "-"), "[[", 2))+
  coord_flip()+
  ggh4x::facetted_pos_scales(y = scales_y)
ggsave(sprintf("./output/sensitivity/parametric_sensitivity_vaccine.png"),
       width = PLOT_WIDTH_ONE*1.1, height = PLOT_HEIGHT_ONE*1.375, units = "in", bg = "#FFFFFF")
