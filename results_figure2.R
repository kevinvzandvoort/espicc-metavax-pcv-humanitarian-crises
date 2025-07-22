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
if(grepl("neutral", .args["output_subdir"])){
  source(sprintf("%s/model/metavax_streppneumo_neutral.R", METAVAX_FOLDER))
} else {
  source(sprintf("%s/model/metavax_streppneumo_diamond.R", METAVAX_FOLDER))  
}

#' Create folder to store output
OUTPUT_FOLDER = setOutputFolder(.args["wd"], .args["output_subdir"])
OUTPUT_SUBFOLDER = setOutputFolder(.args["wd"], sprintf("%s/%s", .args["output_subdir"], .args["output_simdir"]))
OUTPUT_RESULTSFOLDER = setOutputFolder(.args["wd"], sprintf("%s/%s/%s", .args["output_subdir"], .args["output_simdir"], "results"))

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

#' create model params without vaccination
setting = strsplit(.args["output_subdir"], "_", TRUE)[[1]][1]
if(setting == "maiduguri") setting = "maiduguri_acute"
if(setting == "bambari"){
  model_params = modelSetup(contact_data = list(data$demography[[setting]]$displaced$contact_data$data),
                            contacts_extra_host = data$demography[[setting]]$displaced$contacts_extra_host$data$value,
                            contacts_prop_extra = data$demography[[setting]]$displaced$contacts_prop_extra$data$value,
                            population_data = list(data$demography[[setting]]$displaced$population_data$data),
                            prop_malnourished = c(data$demography[[setting]]$displaced$malnourished$data$value),
                            migration_rate = data$demography[[setting]]$displaced$migration_rate$data$value,
                            vaccination_groups = "unvaccinated")
} else {
  model_params = modelSetup(contact_data = list(data$demography[[setting]]$displaced$contact_data$data,
                                                data$demography[[setting]]$host$contact_data$data),
                            contacts_extra_host = data$demography[[setting]]$displaced$contacts_extra_host$data$value,
                            contacts_prop_extra = data$demography[[setting]]$displaced$contacts_prop_extra$data$value,
                            population_data = list(data$demography[[setting]]$displaced$population_data$data,
                                                   data$demography[[setting]]$host$population_data$data),
                            prop_malnourished = c(data$demography[[setting]]$displaced$malnourished$data$value,
                                                  data$demography[[setting]]$host$malnourished$data$value),
                            migration_rate = data$demography[[setting]]$displaced$migration_rate$data$value,
                            vaccination_groups = "unvaccinated")
}

if(grepl("neutral", .args["output_subdir"])){
  #' set initial states for neutral model
  vaccination_groups = "unvaccinated"
  model_params_prevac = model_params_prevac %>%
    setParameter(level = "global", key = "initial_states",
                 value = c(0.8, 0.1, 0.1, 0, 0, 0) %>%
                   rep(each = model_params$global_settings$age_groups_model[, .N]) %>%
                   c(c(0, 0, 0, 0, 0, 0) %>% rep(each = model_params$global_settings$age_groups_model[, .N]) %>%
                       rep(length(vaccination_groups) - 1)) %>%
                   rep(length(model_params$populations)))
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

#' Read plotdata
posterior_runs_summarized_posterior = readRDS(sprintf("%s/plotdata_posterior_runs_summarized_prevalence.RDS", OUTPUT_RESULTSFOLDER))
posterior_runs_summarized_impact = readRDS(sprintf("%s/plotdata_posterior_runs_summarized_impact.RDS", OUTPUT_RESULTSFOLDER))
incidence_impact_relative_impact = readRDS(sprintf("%s/plotdata_incidence_impact.RDS", OUTPUT_RESULTSFOLDER))

(plot_prevalence_fig = posterior_runs_summarized_posterior %>%
    ggplot(aes(x=time, y=med,
               colour=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
               fill=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
               group = paste0(vaccine_strategy, age_group)))+
    facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~., scales="free")+
    geom_ribbon(alpha=0.2, aes(ymin=low95, ymax=high95, colour = NULL))+
    geom_line(linewidth=2)+
    scale_colour_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
    scale_fill_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
    theme_minimal()+
    labs(x="Years since PCV campaign", y="VT carriage prevalence", colour="Vaccine strategy", fill="Vaccine strategy")+
    theme(legend.position="bottom", panel.grid.minor.y = element_blank(), strip.text = element_text(face="bold", size=14),
          axis.title = element_text(size=14), axis.text = element_text(size = 12), legend.title = element_text(size = 14),
          legend.text = element_text(size = 12), panel.spacing.y=unit(5, "mm"))+
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5) * 365, labels = paste0(c("0", 1, 2, 3, 4, 5), "y"))+
    scale_y_continuous(labels = scales::percent, lim=c(0, NA))+
    geom_vline(data=data.table(age_group = factor("[0y, 1y)", c("[0y, 120y)", "[0y, 1y)")), time = 1*365),
               aes(xintercept=time, x=NULL, colour=NULL, fill=NULL, group=NULL),
               colour="#000000", linewidth=1, linetype=2))

(plot_incidence_impact = posterior_runs_summarized_impact %>%
    rbind(data.table(age_group=factor("[0y, 1y)", c("[0y, 120y)", "[0y, 1y)")), vaccine_strategy="<0"), fill=TRUE) %>%
    ggplot(aes(x=time, y=med,
               colour=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
               fill=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
               group = paste0(vaccine_strategy, age_group)))+
    facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~., scales="free")+
    geom_ribbon(alpha=0.2, aes(ymin=low95, ymax=high95, colour = NULL))+
    geom_line(size=2)+
    scale_colour_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
    scale_fill_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
    theme_minimal()+labs(x="Years since PCV campaign", y="Daily impact on severe\npneumococcal disease cases", colour="Vaccine strategy", fill="Vaccine strategy")+
    theme(legend.position="bottom", panel.grid.minor.y = element_blank(), strip.text = element_text(face="bold", size=14),
          axis.title = element_text(size=14), axis.text = element_text(size = 12), legend.title = element_text(size = 14),
          legend.text = element_text(size = 12), panel.spacing.y=unit(5, "mm"))+
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5) * 365, labels = paste0(c("0", 1, 2, 3, 4, 5), "y"))+
    scale_y_continuous(labels = scales::percent, lim=c(0, NA))+
    geom_vline(data=data.table(age_group = factor("[0y, 1y)", c("[0y, 120y)", "[0y, 1y)")), time = 1*365),
               aes(xintercept=time, x=NULL, colour=NULL, fill=NULL, group=NULL),
               colour="#000000", linewidth=1, linetype=2))

(plot_incidence_impact_3y = incidence_impact_relative_impact %>%
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
    scale_colour_manual(values = c("1" = lshtm_colours$purple, "2" = lshtm_colours$pink))+
    scale_fill_manual(values = c("1" = lshtm_colours$purple, "2" = lshtm_colours$pink))+
    labs(colour = "Number of doses in infants",
         fill = "Number of doses in infants",
         x = "Vaccination strategy",
         y= "Cumulative impact on severe\npneumococcal disease cases over 3 years")+
    guides(alpha = "none")+
    theme(legend.position="bottom",
          panel.grid.minor.y = element_blank(),
          strip.text = element_text(face="bold", size=14),
          axis.title = element_text(size=14),
          axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 12),
          panel.border = element_rect(size=1, fill="#FFFFFF00")))

(plot_fig2 = plot_prevalence_fig+
    (plot_incidence_impact)+
    (plot_incidence_impact_3y)+
    plot_annotation(tag_levels = "A")+
    plot_layout(guides = "collect")&
    theme(legend.position="bottom"))
ggsave(plot_fig2, filename = sprintf("%s/figure_prevalence_incidence_doses.png", OUTPUT_RESULTSFOLDER),
       width = PLOT_WIDTH_ONE, height = PLOT_HEIGHT_ONE, units = "in", bg = "#FFFFFF")

