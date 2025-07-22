#' load/update all libraries
pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve, magrittr, units, ggplot2, BayesianTools, ggplot2, ggh4x, lvplot, binom, tictoc, patchwork)

# Can setwd by clicking on Session > Set Working Directory > To Source File Location in the toolbar in RStudio

#' set arguments if R is running interactively, otherwise arguments should be passed to the Rscript command
#'  1: working directory of pcvmr (make sure to setwd manually if running interactively)
#'  2: working directory of pcvm
#'  3: name of output folder to store objects (e.g. date, or specific scenario name)
.args = if(interactive()) c(getwd(), "../pcvm", "digaale", "simulations") else commandArgs(trailingOnly = TRUE)
.args = setNames(.args, c("wd", "metavax_dir", "output_subdir", "output_simdir"))
setwd(.args["wd"])

#' Load metavax model and setup initial model_params skeleton
METAVAX_FOLDER = .args["metavax_dir"]
source(sprintf("%s/model/metavax_streppneumo_diamond.R", METAVAX_FOLDER))

#' Create folder to store output
OUTPUT_FOLDER = setOutputFolder(.args["wd"], .args["output_subdir"])
OUTPUT_SUBFOLDER = setOutputFolder(.args["wd"], sprintf("%s/%s", .args["output_subdir"], .args["output_simdir"]))
OUTPUT_RESULTSFOLDER = setOutputFolder(.args["wd"], sprintf("%s/%s/%s", .args["output_subdir"], .args["output_simdir"], "results"))

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

#' create model params without vaccination
if(.args["output_subdir"] == "bambari"){
  model_params = modelSetup(contact_data = list(data$demography[[.args["output_subdir"]]]$displaced$contact_data$data),
                            contacts_extra_host = data$demography[[.args["output_subdir"]]]$displaced$contacts_extra_host$data$value,
                            contacts_prop_extra = data$demography[[.args["output_subdir"]]]$displaced$contacts_prop_extra$data$value,
                            population_data = list(data$demography[[.args["output_subdir"]]]$displaced$population_data$data),
                            prop_malnourished = c(data$demography[[.args["output_subdir"]]]$displaced$malnourished$data$value),
                            migration_rate = data$demography[[.args["output_subdir"]]]$displaced$migration_rate$data$value,
                            vaccination_groups = "unvaccinated")
} else {
  model_params = modelSetup(contact_data = list(data$demography[[.args["output_subdir"]]]$displaced$contact_data$data,
                                                data$demography[[.args["output_subdir"]]]$host$contact_data$data),
                            contacts_extra_host = data$demography[[.args["output_subdir"]]]$displaced$contacts_extra_host$data$value,
                            contacts_prop_extra = data$demography[[.args["output_subdir"]]]$displaced$contacts_prop_extra$data$value,
                            population_data = list(data$demography[[.args["output_subdir"]]]$displaced$population_data$data,
                                                   data$demography[[.args["output_subdir"]]]$host$population_data$data),
                            prop_malnourished = c(data$demography[[.args["output_subdir"]]]$displaced$malnourished$data$value,
                                                  data$demography[[.args["output_subdir"]]]$host$malnourished$data$value),
                            migration_rate = data$demography[[.args["output_subdir"]]]$displaced$migration_rate$data$value,
                            vaccination_groups = "unvaccinated")
}

vaccination_strategies_plot_settings = list(
  list(model_name = "<0", name = "No vaccination", target = T, dosing = F, colour = lshtm_colours$black),
  list(model_name = "<1", name = "<1y", target = T, dosing = T, colour = lshtm_colours$red),
  list(model_name = "<2", name = "<2y", target = T, dosing = T, colour = lshtm_colours$yellow),
  list(model_name = "<5", name = "<5y", target = T, dosing = T, colour = lshtm_colours$blue),
  list(model_name = "<10", name = "<10y", target = T, dosing = T, colour = lshtm_colours$green),
  list(model_name = "<15", name = "<15y", target = T, dosing = T, colour = lshtm_colours$lightgreen)) %>%
  lapply(as.data.table) %>% rbindlist()
vaccination_strategies_plot_settings[, name := factor(name, vaccination_strategies_plot_settings$name)]

posterior_runs_1dose = readRDS(sprintf("%s/posterior_runs_main_v2.RDS", OUTPUT_SUBFOLDER))
posterior_runs_2dose = readRDS(sprintf("%s/posterior_runs_2dose_infants.RDS", OUTPUT_SUBFOLDER))

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
vaccine_doses_1dose = lapply(upper_age_limits, function(upper_age_limit){
  getParameter(model_params, key = "N", level = "population", population = c("idp_non_malnourished", "idp_malnourished")) %>%
    merge(model_params$global_settings$age_groups_model, by = "age") %>%
    .[, .(N = sum(N)), by = c("age", "name", "from", "to")] %>%
    .[, coverage := getVaccineCoverage(model_params$global_settings$age_groups_model,
                                       c(set_units(6, "weeks"), set_units(upper_age_limit, "year")), 0.85)] %>%
    .[, .(vaccine_strategy = sprintf("<%s", upper_age_limit), doses = sum(N * coverage))]
}) %>% rbindlist()
vaccine_doses_2dose = lapply(upper_age_limits, function(upper_age_limit){
  getParameter(model_params, key = "N", level = "population", population = c("idp_non_malnourished", "idp_malnourished")) %>%
    merge(model_params$global_settings$age_groups_model, by = "age") %>%
    .[, .(N = sum(N)), by = c("age", "name", "from", "to")] %>%
    .[, coverage := getVaccineCoverage(model_params$global_settings$age_groups_model,
                                       c(set_units(6, "weeks"), set_units(upper_age_limit, "year")), 0.85)] %>%
    .[, coverage := ifelse(to <= set_units(1, "year"), coverage * 2, coverage)] %>%
    .[, .(vaccine_strategy = sprintf("<%s", upper_age_limit), doses = sum(N * coverage))]
}) %>% rbindlist()

upper_time_limits = c(0.5, 1, 2, 3)
incidence_impact_1dose = lapply(upper_time_limits, function(upper_time_limit){
  posterior_runs_1dose %>%
    .[outcome == "incidence"] %>%
    .[time <= 365 * upper_time_limit] %>%
    .[, .(value = sum(value)), by = c("age_group", "vaccine_strategy", "run")] %>%
    dcast(...~vaccine_strategy, value.var = "value") %>%
    melt(measure.vars = vaccination_strategies_plot_settings[dosing == TRUE, model_name],
         variable.name = "vaccine_strategy", value.name = "value") %>%
    .[, .(absolute_impact = (`<0` - value) * 5.98,
          relative_impact = 1 - value/`<0`), by = c("age_group", "vaccine_strategy", "run")] %>%
    merge(population_size[, c("name", "N")], by.y = "name", by.x = "age_group") %>%
    merge(vaccine_doses_1dose, by = "vaccine_strategy") %>%
    .[, NNV := doses/absolute_impact] %>%
    .[, absolute_impact_p10000 := absolute_impact/N * 10000] %>%
    .[, -c("N")] %>%
    melt(measure.vars = c("absolute_impact", "relative_impact", "NNV", "absolute_impact_p10000")) %>%
    .[, time := upper_time_limit] %>%
    return()
}) %>% rbindlist()

incidence_impact_2dose = lapply(upper_time_limits, function(upper_time_limit){
  posterior_runs_2dose %>%
    .[outcome == "incidence"] %>%
    .[time <= 365 * upper_time_limit] %>%
    .[, .(value = sum(value)), by = c("age_group", "vaccine_strategy", "run")] %>%
    dcast(...~vaccine_strategy, value.var = "value") %>%
    melt(measure.vars = vaccination_strategies_plot_settings[dosing == TRUE, model_name],
         variable.name = "vaccine_strategy", value.name = "value") %>%
    .[, .(absolute_impact = (`<0` - value) * 5.98,
          relative_impact = 1 - value/`<0`), by = c("age_group", "vaccine_strategy", "run")] %>%
    merge(population_size[, c("name", "N")], by.y = "name", by.x = "age_group") %>%
    merge(vaccine_doses_2dose, by = "vaccine_strategy") %>%
    .[, NNV := doses/absolute_impact] %>%
    .[, absolute_impact_p10000 := absolute_impact/N * 10000] %>%
    .[, -c("N")] %>%
    melt(measure.vars = c("absolute_impact", "relative_impact", "NNV", "absolute_impact_p10000")) %>%
    .[, time := upper_time_limit] %>%
    return()
}) %>% rbindlist()

incidence_impact = incidence_impact_1dose[, scenario_doses := 1] %>%
  rbind(incidence_impact_2dose[, scenario_doses := 2])

plot_incidence_impact = incidence_impact[variable == "relative_impact"] %>%
  .[time >= 1] %>%
  ggplot(aes(x = factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             y = value,
             colour = factor(as.character(scenario_doses)),
             fill = factor(as.character(scenario_doses))))+
  facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~factor(time, c(0.5, 1, 2, 3), c("0-6m", "0-1y", "0-2y", "0-3y")), scales="free")+
  geom_lv(varwidth = TRUE,
          width.method = "height",
          aes(alpha = after_stat(LV)),
          outlier.colour = adjustColBrightness(lshtm_colours$grey, 0.25),
          outlier.shape = 19,
          outlier.size = 1,
          outlier.stroke = 0.25,
          percent=5)+
  theme_minimal()+
  scale_y_continuous(labels = scales::percent)+
  scale_colour_manual(values = c("1" = lshtm_colours$green, "2" = lshtm_colours$lightgreen))+
  scale_fill_manual(values = c("1" = lshtm_colours$green, "2" = lshtm_colours$lightgreen))+
  labs(colour = "Number of doses in infants",
       fill = "Number of doses in infants",
       x = "Vaccination strategy",
       y= "Impact on cumulative severe pneumococcal\ndisease cases compared to no vaccination")+
  guides(alpha = "none")+
  theme(legend.position="bottom",
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(face="bold", size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_rect(size=1, fill="#FFFFFF00"))

plot_incidence_NNV = incidence_impact[variable == "NNV"] %>%
  .[time >= 1] %>%
  ggplot(aes(x = factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             y = value,
             colour = factor(as.character(scenario_doses)),
             fill = factor(as.character(scenario_doses))))+
  facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~factor(time, c(0.5, 1, 2, 3), c("0-6m", "0-1y", "0-2y", "0-3y")), scales="free")+
  geom_lv(varwidth = TRUE,
          width.method = "height",
          aes(alpha = after_stat(LV)),
          outlier.colour = adjustColBrightness(lshtm_colours$grey, 0.25),
          outlier.shape = 19,
          outlier.size = 1,
          outlier.stroke = 0.25,
          percent=5)+
  theme_minimal()+
  scale_colour_manual(values = c("1" = lshtm_colours$green, "2" = lshtm_colours$lightgreen))+
  scale_fill_manual(values = c("1" = lshtm_colours$green, "2" = lshtm_colours$lightgreen))+
  labs(colour = "Number of doses in infants",
       fill = "Number of doses in infants",
       x = "Vaccination strategy",
       y= "Number of doses needed to prevent one\ncase of severe pneumococcal disease")+
  guides(alpha = "none")+
  theme(legend.position="bottom",
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(face="bold", size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_rect(size=1, fill="#FFFFFF00"))+
  scale_y_continuous(labels = function(x) format(x, scientific=FALSE, big.mark=","), limits = c(0, NA))

(plot_doses_infants = plot_incidence_impact+
    plot_incidence_NNV+
    plot_annotation(tag_levels = "A")+
    plot_layout(guides = "collect")&
    theme(legend.position="bottom"))
ggsave(plot_doses_infants, filename = sprintf("%s/figure_doses_infants.png", OUTPUT_RESULTSFOLDER),
       width = PLOT_WIDTH_ONE, height = PLOT_HEIGHT_ONE, units = "in", bg = "#FFFFFF")

#' save output to combine in Figure
saveRDS(incidence_impact[variable == "relative_impact"],
        sprintf("%s/plotdata_incidence_impact.RDS", OUTPUT_RESULTSFOLDER))

plot_incidence_impact_3y = incidence_impact[variable == "relative_impact"] %>%
  .[time == 3] %>%
  ggplot(aes(x = factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             y = value,
             colour = factor(as.character(scenario_doses)),
             fill = factor(as.character(scenario_doses))))+
  facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~., scales="free")+
  geom_lv(varwidth = TRUE,
          width.method = "height",
          aes(alpha = after_stat(LV)),
          outlier.colour = adjustColBrightness(lshtm_colours$grey, 0.25),
          outlier.shape = 19,
          outlier.size = 1,
          outlier.stroke = 0.25,
          percent=5)+
  theme_minimal()+
  scale_y_continuous(labels = scales::percent)+
  scale_colour_manual(values = c("1" = lshtm_colours$green, "2" = lshtm_colours$lightgreen))+
  scale_fill_manual(values = c("1" = lshtm_colours$green, "2" = lshtm_colours$lightgreen))+
  labs(colour = "Number of doses in infants",
       fill = "Number of doses in infants",
       x = "Vaccination strategy",
       y= "Impact on cumulative severe pneumococcal\ndisease cases compared to no vaccination")+
  guides(alpha = "none")+
  theme(legend.position="bottom",
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(face="bold", size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_rect(size=1, fill="#FFFFFF00"))

plot_incidence_NNV_3y = incidence_impact[variable == "NNV"] %>%
  .[time == 3] %>%
  ggplot(aes(x = factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             y = value,
             colour = factor(as.character(scenario_doses)),
             fill = factor(as.character(scenario_doses))))+
  facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~., scales="free")+
  geom_lv(varwidth = TRUE,
          width.method = "height",
          aes(alpha = after_stat(LV)),
          outlier.colour = adjustColBrightness(lshtm_colours$grey, 0.25),
          outlier.shape = 19,
          outlier.size = 1,
          outlier.stroke = 0.25,
          percent=5)+
  theme_minimal()+
  scale_colour_manual(values = c("1" = lshtm_colours$green, "2" = lshtm_colours$lightgreen))+
  scale_fill_manual(values = c("1" = lshtm_colours$green, "2" = lshtm_colours$lightgreen))+
  labs(colour = "Number of doses in infants",
       fill = "Number of doses in infants",
       x = "Vaccination strategy",
       y= "Number of doses needed to prevent one\ncase of severe pneumococcal disease")+
  guides(alpha = "none")+
  theme(legend.position="bottom",
        panel.grid.minor.y = element_blank(),
        strip.text = element_text(face="bold", size=14),
        axis.title = element_text(size=14),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        panel.border = element_rect(size=1, fill="#FFFFFF00"))+
  scale_y_continuous(labels = function(x) format(x, scientific=FALSE, big.mark=","), limits = c(0, NA))

(plot_doses_infants_3y = plot_incidence_impact_3y+
    plot_incidence_NNV_3y+
    plot_annotation(tag_levels = "A")+
    plot_layout(guides = "collect")&
    theme(legend.position="bottom"))

ggsave(plot_doses_infants_3y, filename = sprintf("%s/figure_doses_infants_3y.png", OUTPUT_RESULTSFOLDER),
       width = PLOT_WIDTH_TWOTHIRD, height = PLOT_WIDTH_ONETWO, units = "in", bg = "#FFFFFF")
