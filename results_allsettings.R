#' load/update all libraries
pacman::p_load(data.table, Rcpp, RcppArmadillo, inline, deSolve, rootSolve, magrittr, units, ggplot2, BayesianTools, ggplot2, ggh4x, lvplot, binom, tictoc, patchwork)

# Can setwd by clicking on Session > Set Working Directory > To Source File Location in the toolbar in RStudio

#' set arguments if R is running interactively, otherwise arguments should be passed to the Rscript command
#'  1: working directory of pcvmr (make sure to setwd manually if running interactively)
#'  2: working directory of pcvm
#'  3: name of output folder to store objects (e.g. date, or specific scenario name)
.args = if(interactive()) c(getwd(), "../pcvm", "simulations") else commandArgs(trailingOnly = TRUE)
.args = setNames(.args, c("wd", "metavax_dir", "output_simdir"))
setwd(.args["wd"])

#' Load metavax model and setup initial model_params skeleton
METAVAX_FOLDER = .args["metavax_dir"]
source(sprintf("%s/model/metavax_streppneumo_diamond.R", METAVAX_FOLDER))  

#' Create folder to store output
OUTPUT_FOLDER = setOutputFolder(.args["wd"], "")

#' Set up project specific functions and parameters
source("./functions.R")
source("./data_load_all.R")
source("./model_setup.R")

vaccination_strategies_plot_settings = list(
  list(model_name = "<0", name = "No vaccination", target = T, dosing = F, colour = lshtm_colours$black),
  list(model_name = "<1", name = "<1y", target = T, dosing = T, colour = lshtm_colours$red),
  list(model_name = "<2", name = "<2y", target = T, dosing = T, colour = lshtm_colours$yellow),
  list(model_name = "<5", name = "<5y", target = T, dosing = T, colour = lshtm_colours$blue),
  list(model_name = "<10", name = "<10y", target = T, dosing = T, colour = lshtm_colours$green),
  list(model_name = "<15", name = "<15y", target = T, dosing = T, colour = lshtm_colours$lightgreen)) %>%
  lapply(as.data.table) %>% rbindlist()
vaccination_strategies_plot_settings[, name := factor(name, vaccination_strategies_plot_settings$name)]

plotdata_prevalence = lapply(c("digaale", "bentiu", "bambari", "maiduguri_acute"), function(setting){
  plotdata_prevalence = readRDS(sprintf("./output/%s/simulations/results/plotdata_prevalence.RDS", setting))
  plotdata_prevalence[, setting := setting]
  return(plotdata_prevalence)
}) %>% rbindlist()

plotdata_incidence = lapply(c("digaale", "bentiu", "bambari", "maiduguri_acute"), function(setting){
  plotdata_incidence = readRDS(sprintf("./output/%s/simulations/results/plotdata_incidence.RDS", setting))
  plotdata_incidence[, setting := setting]
  return(plotdata_incidence)
}) %>% rbindlist()

plotdata_cum_incidence = lapply(c("digaale", "bentiu", "bambari", "maiduguri_acute"), function(setting){
  plotdata_cum_incidence = readRDS(sprintf("./output/%s/simulations/results/plotdata_cum_incidence.RDS", setting))
  plotdata_cum_incidence[, setting := setting]
  return(plotdata_cum_incidence)
}) %>% rbindlist()

setting_names = c(digaale = "Digaale (base scenario)",
                  bentiu = "Bentiu",
                  bambari = "Bambari",
                  maiduguri_acute = "Maiduguri")

(plot_prevalence = plotdata_prevalence %>%
  ggplot(aes(x=time, y=med,
             colour=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             fill=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
             group = paste0(vaccine_strategy, age_group)))+
  facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~factor(setting, names(setting_names), setting_names), scales="free")+
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

ggsave(plot_prevalence, filename = sprintf("%s/figure_prevalence_allsettings.png", OUTPUT_FOLDER),
       width = PLOT_WIDTH_ONE, height = PLOT_HEIGHT_ONE, units = "in", bg = "#FFFFFF")

(plot_incidence = plotdata_incidence %>%
    ggplot(aes(x=time, y=med,
               colour=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
               fill=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
               group = paste0(vaccine_strategy, age_group)))+
    facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~factor(setting, names(setting_names), setting_names), scales="free")+
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

ggsave(plot_incidence, filename = sprintf("%s/figure_incidence_allsettings.png", OUTPUT_FOLDER),
       width = PLOT_WIDTH_ONE, height = PLOT_HEIGHT_ONE, units = "in", bg = "#FFFFFF")

(plot_cum_incidence_impact = plotdata_cum_incidence %>%
    ggplot(aes(x=time, y=med,
               colour=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
               fill=factor(vaccine_strategy, vaccination_strategies_plot_settings$model_name, vaccination_strategies_plot_settings$name),
               group = paste0(vaccine_strategy, age_group)))+
    facet_grid(factor(age_group, c("[0y, 120y)", "[0y, 1y)"), c("all ages", "infants"))~factor(setting, names(setting_names), setting_names), scales="free")+
    geom_ribbon(alpha=0.2, aes(ymin=low95, ymax=high95, colour = NULL))+
    geom_line(size=2)+
    scale_colour_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
    scale_fill_manual(values = vaccination_strategies_plot_settings[, c("name", "colour")] %>% as.matrix("name") %>% .[,1])+
    theme_minimal()+labs(x="Years since PCV campaign", y="Cumulative impact on severe\npneumococcal disease cases", colour="Vaccine strategy", fill="Vaccine strategy")+
    theme(legend.position="bottom", panel.grid.minor.y = element_blank(), strip.text = element_text(face="bold", size=14),
          axis.title = element_text(size=14), axis.text = element_text(size = 12), legend.title = element_text(size = 14),
          legend.text = element_text(size = 12), panel.spacing.y=unit(5, "mm"))+
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5) * 365, labels = paste0(c("0", 1, 2, 3, 4, 5), "y"))+
    scale_y_continuous(labels = scales::percent, lim=c(0, NA)))

ggsave(plot_cum_incidence_impact, filename = sprintf("%s/figure_cum_incidence_allsettings.png", OUTPUT_FOLDER),
       width = PLOT_WIDTH_ONE, height = PLOT_HEIGHT_ONE, units = "in", bg = "#FFFFFF")
