#' Get data from Zenodo
digaale_contact_data =
  socialmixr::get_survey("https://zenodo.org/record/5226280")

#' The estimated population size in Digaale (for provided age groups)
#'  can manually be downloaded
digaale_survey_population =
  data.table::fread("https://zenodo.org/record/7071876/files/espicc_somaliland_digaale_survey_population.csv")

#' Note that weekends fall on Fridays and Saturdays in Somaliland.
#' - The dayofweek variable provided in the dataset has been kept
#'   consistent with R defaults (0: Sunday to 6: Saturday)
digaale_contact_data$participants[, c("dayofweek", "dayofweek_name", "weekend")] %>%
  unique %>% setorder(dayofweek) %>% .[]
#' socialmixr currently assumes the weekend to fall on dayofweek
#'  6 (Saturday) and 0 (Sunday)
#' - dayofweek can be manually edited so that Fridays and Saturdays
#'   are taken as the weekend, if you wish to weight contacts by
#'   weekday
digaale_contact_data$participants[, dayofweek := ifelse(dayofweek == 6, 0, dayofweek + 1)]

#' The contact matrix can then be constructed as follows
#' - The provided survey_population can be used to construct a
#'   population representative matrix for Digaale IDP camp
#' - As the sample is not self-weighing (oversampling of young
#'   age groups), it is recommended to apply the survey_weight
#'   as weights
digaale_contact_matrix = digaale_contact_data %>%
  socialmixr::contact_matrix(survey.pop = digaale_survey_population,
                             age.limits = c(seq(0, 20, 5), seq(30, 60, 10)),
                             symmetric = TRUE, weights = "survey_weight", weigh.dayofweek = TRUE)

cmatrix = digaale_contact_matrix$matrix %>% as.data.table()
cnames = copy(colnames(cmatrix))
cmatrix[, contactor_age_group := cnames]
cmatrix = cmatrix %>% melt(measure.vars = cnames, variable.name = c("contactee_age_group"))

age_groups = data.table(name = cnames,
                        from = c(seq(0, 20, 5), seq(30, 60, 10)),
                        to = c(seq(5, 20, 5), seq(30, 60, 10), 100))

plot_cmatrix = cmatrix %>%
  merge(age_groups %>%
          setNames(c("name", "contactor_from", "contactor_to")),
        by.x = "contactor_age_group",
        by.y = "name") %>%
  merge(age_groups %>%
          setNames(c("name", "contactee_from", "contactee_to")),
        by.x = "contactee_age_group",
        by.y = "name") %>%
  #.[, value_binned := cut(value, breaks = seq(0, 4.5, by = 0.5), include.lowest = TRUE, labels = paste0(seq(0, 4, by = 0.5), "–", seq(0.5, 4.5, by = 0.5)))] %>%
  .[, value_binned := cut(value, breaks = seq(0, 4.5, by = 0.5), include.lowest = TRUE)] %>%
  ggplot(aes(xmin = contactor_from, xmax = contactor_to,
             ymin = contactee_from, ymax = contactee_to,
             fill = value))+
  geom_rect()+
  geom_vline(xintercept = seq(5, 55, 5), colour = rgb(1,1,1,0.05))+
  geom_hline(yintercept = seq(5, 55, 5), colour = rgb(1,1,1,0.05))+
  coord_cartesian(xlim = c(0, 60), ylim = c(0, 60), expand = FALSE)+
  scale_x_continuous(breaks = seq(0, 60, 10), labels = c(seq(0, 50, 10), "60+"))+
  scale_y_continuous(breaks = seq(10, 60, 10), labels = c(seq(10, 50, 10), "60+"))+
  labs(x = "Age", y = "Age of contact", fill = "Contacts per day", title = "I. Social contacts")+
  theme_minimal()+
  theme(legend.position = "right",
        legend.direction = "vertical",
        legend.margin = margin(),
        axis.text.y = element_text(hjust = 0.5),
        plot.title = element_text(size = 10, face = "bold"),
        plot.title.position = "plot",
        plot.margin = margin(t = "1"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6.5),
        legend.title = element_text(size = 8, face = "italic"),
        legend.text = element_text(size = 6.5),
        legend.box.spacing = unit(2, "mm"))+
  #scale_fill_manual(values = colorRampPalette(c("#102e21", "#00BF6F"))(9),
  #                  guide = guide_colorsteps(title.position = "left", direction = "vertical",
  #                                          title.theme = element_text(angle = 90)))
  scale_fill_gradientn(colours = c("#102e21", "#013d24", "#00BF6F"),values = c(0, 0.05, 1),
                       guide = guide_colourbar(title.position = "left", direction = "vertical",
                                               title.theme = element_text(angle = 90), barwidth = 0.75, barheight = 6, ticks.colour = NA))

carriage_data = fread("https://raw.githubusercontent.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019-carriage/refs/heads/main/data/participant_data.csv") %>%
  merge(fread("https://raw.githubusercontent.com/kevinvzandvoort/espicc-somaliland-digaale-survey-2019-carriage/refs/heads/main/data/lab_data.csv"), by = "pid")

age_groups_carriage = data.table(low = c(0, 2, 5, 10, 15, 30, 50),
                        high = c(2, 5, 10, 15, 30, 50, 120),
                        label = c("<2", "2-4", "5-9", "10-14", "15-30", "30-49", "50+"))

for(i in 1:nrow(age_groups)){
  carriage_data[participant_age_y >= age_groups_carriage[i, low] & participant_age_y < age_groups_carriage[i, high], age_group := age_groups_carriage[i, label]]  
}
carriage_data[, age_group := factor(age_group, age_groups_carriage$label)]
plot_carriage = carriage_data[!is.na(age_group) & !is.na(pneu_carr_final), .(total = sum(pneu_carr_final == 1), pneusil = sum(pneusil_carr == 1, na.rm = TRUE), N = .N), by = "age_group"] %>%
  .[, c("All serotypes", "Vaccine types") := .(total/N, pneusil/N)] %>%
  melt(measure.vars = c("All serotypes", "Vaccine types"), id.vars = "age_group") %>%
  ggplot(aes(x=age_group, y = value, fill=variable))+
  geom_col(position = "dodge")+
  theme_minimal()+
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.8))+
  geom_hline(yintercept = 0)+
  theme(legend.position = "inside",
        legend.position.inside =  c(1, 1),
        legend.margin = margin(l = 1, t = 1, b = 1, r = 1),
        legend.title = element_blank(),
        legend.box.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"),
        legend.key.size = unit(1, "lines"),
        legend.background = element_rect(fill = "#FFFFFF",
                                         colour = "#000000",
                                         linewidth = 0.25),
        legend.justification = c("right", "top"),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(colour = "#DEDEDE"),
        plot.title = element_text(size = 10, face = "bold"),
        plot.title.position = "plot",
        plot.margin = margin(t = "1"),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6.5),
        legend.text = element_text(size = 6.5))+
  labs(x = "Age", y = "Prevalence", fill = "Serotypes", title = "II. Pneumococcal carriage prevalence")+
  scale_fill_manual(values = c("All serotypes" = "#0D5257",
                               "Vaccine types" = "#00BF6F"))+
  coord_cartesian(expand = FALSE)

library(patchwork)
plot_cmatrix/free(plot_carriage)

plot_cmatrix+
  plot_spacer()+
  free(plot_carriage)+ 
  plot_layout(widths = c(1 - 1/(16/22 + 1), 0.025, 1/(16/22 + 1)))
ggsave("modelling_crises_data.png", width = 6, height = 2.5, units = "in", dpi = 300)
ggsave("modelling_crises_data.eps", width = 6, height = 2.5, units = "in")
ggsave("modelling_crises_data.pdf", width = 6, height = 2.5, units = "in")

plot_title = ggplot()+
  labs(title = "Primary data collection")+
  geom_text(label = "We conducted a cross-sectional survey in\nDigaale IDP camp (Somaliland) to collect key\ndata to parameterize a pneumococcal\ntransmission model, including:",
            aes(x=0, y=0),
            hjust = 0, vjust = 1,
            size = 3.5)+
  coord_cartesian(xlim=c(0, 0.5),
                  ylim = c(-0.0125, 0),
                  expand = FALSE, clip = "off")+
  theme_void()+
  theme(plot.title = element_text(size = 12, face = "bold"), plot.background = element_rect(fill="red"), panel.background = element_rect(fill = "orange"),
        axis.text = element_blank(), axis.title = element_blank(), plot.margin = unit(c(0,0,0,0), "mm"))
free(plot_title)/plot_cmatrix/free(plot_carriage)
