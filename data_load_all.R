getData = function(name){
  if(file.exists(sprintf("%s_processed.RDS", name))){
    return(readRDS(sprintf("%s_processed.RDS", name)))
  } else {
    value = source(sprintf("%s.R", name))$value
    saveRDS(value, sprintf("%s_processed.RDS", name))
    return(value)
  }
}

data = list(
  epidemiology = list(
    prevalence = list(
      digaale_pneumosil = list(
        source = "van Zandvoort K, Hassan AI, Bobe MO, Pell CL, Ahmed MS, Ortika BD, Ibrahim S, Abdi MI, Karim MA, Eggo RM, Ali SY. Pre-vaccination carriage prevalence of Streptococcus pneumoniae serotypes among internally displaced people in Somaliland: a cross-sectional study. Pneumonia. 2024 Dec 5;16(1):25. https://doi.org/10.1186/s41479-024-00148-6",
        data = getData("./data/epidemiology/prevalence/digaale_pneumosil")),
      digaale_pneumosil_alternative = list(
        source = "van Zandvoort K, Hassan AI, Bobe MO, Pell CL, Ahmed MS, Ortika BD, Ibrahim S, Abdi MI, Karim MA, Eggo RM, Ali SY. Pre-vaccination carriage prevalence of Streptococcus pneumoniae serotypes among internally displaced people in Somaliland: a cross-sectional study. Pneumonia. 2024 Dec 5;16(1):25. https://doi.org/10.1186/s41479-024-00148-6",
        data = getData("./data/epidemiology/prevalence/digaale_pneumosil_alternative")),
      digaale_synflorix = list(
        source = "van Zandvoort K, Hassan AI, Bobe MO, Pell CL, Ahmed MS, Ortika BD, Ibrahim S, Abdi MI, Karim MA, Eggo RM, Ali SY. Pre-vaccination carriage prevalence of Streptococcus pneumoniae serotypes among internally displaced people in Somaliland: a cross-sectional study. Pneumonia. 2024 Dec 5;16(1):25. https://doi.org/10.1186/s41479-024-00148-6",
        data = getData("./data/epidemiology/prevalence/digaale_synflorix")),
      digaale_prevenar = list(
        source = "van Zandvoort K, Hassan AI, Bobe MO, Pell CL, Ahmed MS, Ortika BD, Ibrahim S, Abdi MI, Karim MA, Eggo RM, Ali SY. Pre-vaccination carriage prevalence of Streptococcus pneumoniae serotypes among internally displaced people in Somaliland: a cross-sectional study. Pneumonia. 2024 Dec 5;16(1):25. https://doi.org/10.1186/s41479-024-00148-6",
        data = getData("./data/epidemiology/prevalence/digaale_prevenar"))),
    clearance_rate = list(
      source = "Lipsitch M, Abdullahi O, D'Amour A, Xie W, Weinberger DM, Tchetgen ET, Scott JA. Estimating rates of carriage acquisition and clearance and competitive ability for pneumococcal serotypes in Kenya with a Markov transition model. Epidemiology. 2012 Jul 1;23(4):510-9. https://doi.org/10.1097/ede.0b013e31824f2f32",
      data = getData("./data/epidemiology/clearance_rate/kenya_kilifi")),
    case_carrier_ratio = list(
      source = "Flasche S, Ojal J, Le Polain de Waroux O, Otiende M, O’Brien KL, Kiti M, Nokes DJ, Edmunds WJ, Scott JA. Assessing the efficiency of catch-up campaigns for the introduction of pneumococcal conjugate vaccine: a modelling study based on data from PCV10 introduction in Kilifi, Kenya. BMC medicine. 2017 Jun 7;15(1):113. https://doi.org/10.1186/s12916-017-0882-9",
      data = getData("./data/epidemiology/case_carrier_ratio/kenya_kilifi")),
    malnourished_RR_transmission = list(
      source = "Verhagen LM, Gómez-Castellano K, Snelders E, Rivera-Olivero I, Pocaterra L, Melchers WJ, de Waard JH, Hermans PW. Respiratory infections in Enepa Amerindians are related to malnutrition and Streptococcus pneumoniae carriage. Journal of infection. 2013 Oct 1;67(4):273-81.
      Verhagen LM, Hermsen M, Rivera‐Olivero IA, Sisco MC, de Jonge MI, Hermans PW, de Waard JH. Nasopharyngeal carriage of respiratory pathogens in Warao Amerindians: significant relationship with stunting. Tropical Medicine & International Health. 2017 Apr;22(4):407-14
      Gebre T, Tadesse M, Aragaw D, Feye D, Beyene HB, Seyoum D, Mekonnen M. Nasopharyngeal carriage and antimicrobial susceptibility patterns of Streptococcus pneumoniae among children under five in Southwest Ethiopia. Children. 2017 Apr 19;4(4):27.",
      data = setAgeBreaks(0) %>% .[, value := 1.2]),
    malnourished_RR_disease = list(
      source = "https://www.jhsph.edu/ivac/wp-content/uploads/2018/05/9_ISPPD-511.pdf#:~:text=In%20the%20PERCH%20study%2C%20children%20with%20different%20markers,to%20be%20admitted%20to%20hospital%20with%20pneumococcal%20pneumonia.
      Von Mollendorf C, Von Gottberg A, Tempia S, Meiring S, De Gouveia L, Quan V, Lengana S, Avenant T, Du Plessis N, Eley B, Finlayson H. Increased risk for and mortality from invasive pneumococcal disease in HIV-exposed but uninfected infants aged< 1 year in South Africa, 2009–2013. Clinical Infectious Diseases. 2015 May 1;60(9):1346-56. https://academic.oup.com/cid/article/60/9/1346/405394.
      von Mollendorf C, Cohen C, de Gouveia L, Naidoo N, Meiring S, Quan V, Lindani S, Moore DP, Reubenson G, Moshe M, Eley B. Risk factors for invasive pneumococcal disease among children less than 5 years of age in a high HIV prevalence setting, South Africa, 2010 to 2012. The Pediatric infectious disease journal. 2015 Jan 1;34(1):27-34. https://doi.org/10.1097/INF.0000000000000484.",
      data = setAgeBreaks(0) %>% .[, value := 2])),
  vaccine = list(
    efficacy_transmission = list(
      singledose_u1 = list(
        source = "Temple B, Tran HP, Dai VT, Smith-Vaughan H, Balloch A, Beissbarth J, Bright K, Higgins RA, Hinds J, Hoan PT, Nation ML. Efficacy against pneumococcal carriage and the immunogenicity of reduced-dose (0+ 1 and 1+ 1) PCV10 and PCV13 schedules in Ho Chi Minh City, Viet Nam: a parallel, single-blind, randomised controlled trial. The Lancet Infectious diseases. 2023 Aug 1;23(8):933-44. https://doi.org/10.1016/S1473-3099(23)00061-0.",
        data = setAgeBreaks(c(0, 1)) %>% .[, value := c(0.25, 0.5)]),
      twodose_u1 = list(
        source = "Temple B, Tran HP, Dai VT, Smith-Vaughan H, Balloch A, Beissbarth J, Bright K, Higgins RA, Hinds J, Hoan PT, Nation ML. Efficacy against pneumococcal carriage and the immunogenicity of reduced-dose (0+ 1 and 1+ 1) PCV10 and PCV13 schedules in Ho Chi Minh City, Viet Nam: a parallel, single-blind, randomised controlled trial. The Lancet Infectious diseases. 2023 Aug 1;23(8):933-44. https://doi.org/10.1016/S1473-3099(23)00061-0.",
        data = setAgeBreaks(0) %>% .[, value := 0.5])
    ),
    efficacy_total = list(
      singledose_u1 = list(
        source = "Berman-Rosa M, O’Donnell S, Barker M, Quach C. Efficacy and effectiveness of the PCV-10 and PCV-13 vaccines against invasive pneumococcal disease. Pediatrics. 2020 Apr 1;145(4):e20190377. https://doi.org/10.1542/peds.2019-0377.
        Similar for those malnourished: Cohen C, Von Mollendorf C, De Gouveia L, Lengana S, Meiring S, Quan V, Nguweneza A, Moore DP, Reubenson G, Moshe M, Madhi SA. Effectiveness of the 13-valent pneumococcal conjugate vaccine against invasive pneumococcal disease in South African children: a case-control study. The Lancet Global Health. 2017 Mar 1;5(3):e359-69. https://doi.org/10.1016/S2214-109X(17)30043-8.",
        data = setAgeBreaks(c(0, 1)) %>% .[, value := c(0.56, 0.80)]),
      twodose_u1 = list(
        source = "Berman-Rosa M, O’Donnell S, Barker M, Quach C. Efficacy and effectiveness of the PCV-10 and PCV-13 vaccines against invasive pneumococcal disease. Pediatrics. 2020 Apr 1;145(4):e20190377. https://doi.org/10.1542/peds.2019-0377.
        Similar for those malnourished: Cohen C, Von Mollendorf C, De Gouveia L, Lengana S, Meiring S, Quan V, Nguweneza A, Moore DP, Reubenson G, Moshe M, Madhi SA. Effectiveness of the 13-valent pneumococcal conjugate vaccine against invasive pneumococcal disease in South African children: a case-control study. The Lancet Global Health. 2017 Mar 1;5(3):e359-69. https://doi.org/10.1016/S2214-109X(17)30043-8.",
        data = setAgeBreaks(0) %>% .[, value := 0.80])
    ),
    efficacy_duration = list(
      singledose_u1 = list(
        source = "De Waroux OL, Flasche S, Prieto-Merino D, Goldblatt D, Edmunds WJ. The efficacy and duration of protection of pneumococcal conjugate vaccines against nasopharyngeal carriage: a meta-regression model. The Pediatric infectious disease journal. 2015 Aug 1;34(8):858-64. https://doi.org/10.1097/inf.0000000000000717.",
        data = setAgeBreaks(c(0, 1)) %>% .[, value := c(0.5*6, 6)]),
      twodose_u1 = list(
        source = "De Waroux OL, Flasche S, Prieto-Merino D, Goldblatt D, Edmunds WJ. The efficacy and duration of protection of pneumococcal conjugate vaccines against nasopharyngeal carriage: a meta-regression model. The Pediatric infectious disease journal. 2015 Aug 1;34(8):858-64. https://doi.org/10.1097/inf.0000000000000717.",
        data = setAgeBreaks(0) %>% .[, value := c(6)])
    )),
  demography = list(
    digaale = list(
      displaced = list(
        population_data = list(
          source = "van Zandvoort K, Bobe MO, Hassan AI, Abdi MI, Ahmed MS, Soleman SM, Warsame MY, Wais MA, Diggle E, McGowan CR, Satzke C. Social contacts and other risk factors for respiratory infections among internally displaced people in Somaliland. Epidemics. 2022 Dec 1;41:100625. https://doi.org/10.1016/j.epidem.2022.100625",
          data = getData("./data/demography/digaale/population_size_displaced")),
        contact_data = list(
          source = "van Zandvoort K, Bobe MO, Hassan AI, Abdi MI, Ahmed MS, Soleman SM, Warsame MY, Wais MA, Diggle E, McGowan CR, Satzke C. Social contacts and other risk factors for respiratory infections among internally displaced people in Somaliland. Epidemics. 2022 Dec 1;41:100625. https://doi.org/10.1016/j.epidem.2022.100625",
          data = getData("./data/demography/digaale/contact_matrix_displaced")),
        contacts_extra_host = list(
          source = "van Zandvoort K, Bobe MO, Hassan AI, Abdi MI, Ahmed MS, Soleman SM, Warsame MY, Wais MA, Diggle E, McGowan CR, Satzke C. Social contacts and other risk factors for respiratory infections among internally displaced people in Somaliland. Epidemics. 2022 Dec 1;41:100625. https://doi.org/10.1016/j.epidem.2022.100625",
          data = setAgeBreaks(0) %>% .[, value := 0.02658532]),
        household_size = list(
          source = "van Zandvoort K, Bobe MO, Hassan AI, Abdi MI, Ahmed MS, Soleman SM, Warsame MY, Wais MA, Diggle E, McGowan CR, Satzke C. Social contacts and other risk factors for respiratory infections among internally displaced people in Somaliland. Epidemics. 2022 Dec 1;41:100625. https://doi.org/10.1016/j.epidem.2022.100625",
          data = setAgeBreaks(0) %>% .[, value := 4.5]),
        malnourished = list(
          source = "van Zandvoort K, Bobe MO, Hassan AI, Abdi MI, Ahmed MS, Soleman SM, Warsame MY, Wais MA, Diggle E, McGowan CR, Satzke C. Social contacts and other risk factors for respiratory infections among internally displaced people in Somaliland. Epidemics. 2022 Dec 1;41:100625. https://doi.org/10.1016/j.epidem.2022.100625",
          data = setAgeBreaks(0) %>% .[, value := (31+16)/(124+31+16)]),
        migration_rate = list(
          source = "van Zandvoort K, Bobe MO, Hassan AI, Abdi MI, Ahmed MS, Soleman SM, Warsame MY, Wais MA, Diggle E, McGowan CR, Satzke C. Social contacts and other risk factors for respiratory infections among internally displaced people in Somaliland. Epidemics. 2022 Dec 1;41:100625. https://doi.org/10.1016/j.epidem.2022.100625",
          data = setAgeBreaks(0) %>% .[, value := 129.3/1000/365])),
      host = list(
        population_data = list(
          source = "Central Statistics Department, Ministry of Planning and National Development, Somaliland Government, Ministry of Planning and National Development, Somaliland Government. The Somaliland Health and Demographic Survey 2020. 2020.",
          data = getData("./data/demography/digaale/population_size_host")),
        contact_data = list(
          source = "Prem K, Zandvoort KV, Klepac P, Eggo RM, Davies NG, Centre for the Mathematical Modelling of Infectious Diseases COVID-19 Working Group, Cook AR, Jit M. Projecting contact matrices in 177 geographical regions: an update and comparison with empirical data for the COVID-19 era. PLoS computational biology. 2021 Jul 26;17(7):e1009098. https://doi.org/10.1371/journal.pcbi.1009098.",
          data = getData("./data/demography/digaale/contact_matrix_host")),
        malnourished = list(
          source = "Central Statistics Department, Ministry of Planning and National Development, Somaliland Government, Ministry of Planning and National Development, Somaliland Government. The Somaliland Health and Demographic Survey 2020. 2020.",
          data = setAgeBreaks(0) %>% .[, value := 0.2]))),
    bentiu = list(
      displaced = list(
        population_data = list(
          source = "South Sudan: Cumulative Registration Data, January 2014 - August 2015 - South Sudan | ReliefWeb. 2015; published online Aug 31. https://reliefweb.int/report/southsudan/south-sudan-cumulative-registration-data-january-2014-august-2015 (accessed April 16, 2024).",
          data = getData("./data/demography/bentiu/population_size_displaced")),
        contacts_extra_host = list(
          source = "Assume the same as Digaale.",
          data = setAgeBreaks(0) %>% .[, value := 0.02658532]),
        household_size = list(
          source = "South Sudan National Bureau of Statistics (NBS). National Baseline Household Survey 2009 Report. 2012.",
          data = setAgeBreaks(0) %>% .[, value := 7.8]),
        malnourished = list(
          source = "Concern. Bentiu PoC - Nutrition Anthropometry & Retrospective Mortality Survey Report. 2015. https://info.undp.org/docs/pdc/Documents/SSD/2016012447 180405_CWW%20-%20Bentiu%20PoC%20-%20Smart%20Survey%20%20August%202015.pdf (accessed April 16, 2024).",
          data = setAgeBreaks(0) %>% .[, value := 0.179]),
        migration_rate = list(
          source = "Concern. Bentiu PoC - Nutrition Anthropometry & Retrospective Mortality Survey Report. 2015. https://info.undp.org/docs/pdc/Documents/SSD/2016012447 180405_CWW%20-%20Bentiu%20PoC%20-%20Smart%20Survey%20%20August%202015.pdf (accessed April 16, 2024).",
          data = setAgeBreaks(0) %>% .[, value := 12.57/10000])),
      host = list(
        population_data = list(
          source = "Checchi F, Testa A, Warsame A, Bs LQ, Burns R. Estimates of crisis-attributable mortality in South Sudan, December 2013-April 2018. 2018.",
          data = getData("./data/demography/bentiu/population_size_host")),
        contact_data = list(
          source = "Prem K, Zandvoort KV, Klepac P, Eggo RM, Davies NG, Centre for the Mathematical Modelling of Infectious Diseases COVID-19 Working Group, Cook AR, Jit M. Projecting contact matrices in 177 geographical regions: an update and comparison with empirical data for the COVID-19 era. PLoS computational biology. 2021 Jul 26;17(7):e1009098. https://doi.org/10.1371/journal.pcbi.1009098.",
          data = getData("./data/demography/bentiu/contact_matrix_host")),
        malnourished = list(
          source = "Concern. Bentiu PoC - Nutrition Anthropometry & Retrospective Mortality Survey Report. 2015. https://info.undp.org/docs/pdc/Documents/SSD/2016012447 180405_CWW%20-%20Bentiu%20PoC%20-%20Smart%20Survey%20%20August%202015.pdf (accessed April 16, 2024).",
          data = setAgeBreaks(0) %>% .[, value := 0.176]))),
    maiduguri_acute = list(
      displaced = list(
        population_data = list(
          source = "UN OCHA. Northeast Nigeria: Humanitarian emergency - Situation Report No. 1. 2016; published online Nov 28. https://www.unocha.org/publications/report/nigeria/northeast-nigeria-humanitarianemergency-situation-report-no-1-28-november-2016 (accessed April 16, 2024).",
          data = getData("./data/demography/maiduguri/population_size_displaced_acute")),
        contacts_extra_host = list(
          source = "Assume random mixing between host and displaced outside of the household",
          data = setAgeBreaks(0) %>% .[, value := 1000000/(1000000+800000)]),
        household_size = list(
          source = "IOM. Displacement Tracking Matrix. https://dtm.iom.int/ (accessed Dec 18, 2023).",
          data = setAgeBreaks(0) %>% .[, value := 5.03]),
        malnourished = list(
          source = "Nigeria National Bureau of Statistics. Nutrition and food security surveillance: north east Nigeria - Emergency survey. 2019. https://fscluster.org/sites/default/files/documents/nfss_round_8_final_report_november_20 19.pdf (accessed April 16, 2024).",
          data = setAgeBreaks(0) %>% .[, value := 0.192]),
        migration_rate = list(
          source = "Action Against Hunger. Report of Small-Scale SMART Survey in MMC, Jere LGAs, Borno. 2016.",
          data = setAgeBreaks(0) %>% .[, value := 13.69863/10000])),
      host = list(
        population_data = list(
          source = "UN OCHA. Northeast Nigeria: Humanitarian emergency - Situation Report No. 1. 2016; published online Nov 28. https://www.unocha.org/publications/report/nigeria/northeast-nigeria-humanitarianemergency-situation-report-no-1-28-november-2016 (accessed April 16, 2024).",
          data = getData("./data/demography/maiduguri/population_size_host_acute")),
        contact_data = list(
          source = "Prem K, Zandvoort KV, Klepac P, Eggo RM, Davies NG, Centre for the Mathematical Modelling of Infectious Diseases COVID-19 Working Group, Cook AR, Jit M. Projecting contact matrices in 177 geographical regions: an update and comparison with empirical data for the COVID-19 era. PLoS computational biology. 2021 Jul 26;17(7):e1009098. https://doi.org/10.1371/journal.pcbi.1009098.",
          data = getData("./data/demography/maiduguri/contact_matrix_host")),
        malnourished = list(
          source = "Nigeria National Bureau of Statistics. Nutrition and food security surveillance: north east Nigeria - Emergency survey. 2019. https://fscluster.org/sites/default/files/documents/nfss_round_8_final_report_november_20 19.pdf (accessed April 16, 2024).",
          data = setAgeBreaks(0) %>% .[, value := 0.192]))),
    bambari = list(
      displaced = list(
        population_data = list(
          source = "IOM. Displacement Tracking Matrix. https://dtm.iom.int/ (accessed Dec 18, 2023).",
          data = getData("./data/demography/bambari/population_size_displaced")),
        contact_data = list(
          source = "Prem K, Zandvoort KV, Klepac P, Eggo RM, Davies NG, Centre for the Mathematical Modelling of Infectious Diseases COVID-19 Working Group, Cook AR, Jit M. Projecting contact matrices in 177 geographical regions: an update and comparison with empirical data for the COVID-19 era. PLoS computational biology. 2021 Jul 26;17(7):e1009098. https://doi.org/10.1371/journal.pcbi.1009098.",
          data = getData("./data/demography/bambari/contact_matrix_displaced")),
        contacts_extra_host = list(
          source = "N/A",
          data = setAgeBreaks(0) %>% .[, value := NA]),
        household_size = list(
          source = "Not used",
          data = setAgeBreaks(0) %>% .[, value := NA]),
        malnourished = list(
          source = "Unicef. Multiple Indicator Cluster Surveys. https://mics.unicef.org/surveys (accessed April 16, 2024).",
          data = setAgeBreaks(0) %>% .[, value := 0.042]),
        migration_rate = list(
          source = "N/A",
          data = setAgeBreaks(0) %>% .[, value := 0]))),
    kilifi = list(
      host = list(
        malnourished = list(
          source = "County Department of Health. Kilifi County SMART Survey Report. 2016. https://www.nutritionhealth.or.ke/wpcontent/uploads/SMART%20Survey%20Reports/Kilifi%20County%20SMART%20Survey% 20Report%20November2016.pdf.",
          data = setAgeBreaks(0) %>% .[, value := 0.0182])))))

extrapolateContactMatrix = function(popsize_y, household_size_y,
                                    popsize_x = getData("./data/demography/digaale/contact_matrix_displaced")$population$population,
                                    household_size_x = data$demography$digaale$displaced$household_size$data$value,
                                    cm_x_intra = getData("./data/demography/digaale/contact_matrix_displaced")$contact_matrix_intra,
                                    cm_x_extra = getData("./data/demography/digaale/contact_matrix_displaced")$contact_matrix_extra,
                                    cm_agegroups = setAgeBreaks(getData("./data/demography/digaale/contact_matrix_displaced")$population$lower.age.limit),
                                    only_matrix = TRUE){
  
  popsize_y = combineAgeBreaks(cm_agegroups, popsize_y, method="sum")[, value]
  
  expected_matrix_x = createExpectedMatrix(popsize_x)
  expected_matrix_y = createExpectedMatrix(popsize_y)
  
  cm_y_intra = cm_x_intra * (expected_matrix_y/expected_matrix_x) * (household_size_y/household_size_x)
  cm_y_extra = cm_x_extra * (expected_matrix_y/expected_matrix_x)
  
  cm_y = cm_y_intra + cm_y_extra
  
  cm_y = cm_y %>% as.data.table()
  colnames(cm_y) = cm_agegroups[, name] %>% levels
  cm_y[, contactor_age_group := cm_agegroups[, name]]
  cm_y = melt(cm_y, id.vars = "contactor_age_group", variable.name = "contactee_age_group")
  cm_y[, c("contactor_age_group", "contactee_age_group") :=
         .(factor(contactor_age_group, cm_agegroups[, name]),
           factor(contactee_age_group, cm_agegroups[, name]))] %>% .[]
  if(only_matrix){
    return(cm_y)
  } else {
    return(list(cm_intra = cm_y_intra, cm_extra = cm_y_extra, cm_full = cm_y, prop_extra = mean(rowSums(cm_y_extra)/rowSums(cm_y_intra + cm_y_extra))))
  }
}

#' overwrite digaale contact data with actual contact data, and extrapolate to other settings
for(s in c("digaale", "bentiu", "maiduguri_acute")){
  cm = extrapolateContactMatrix(data$demography[[s]]$displaced$population_data$data,
                                data$demography[[s]]$displaced$household_size$data$value, only_matrix = FALSE)  
  
  data$demography[[s]]$displaced$contact_data$data = cm$cm_full
  # proportion of contacts that are extra-household contacts
  data$demography[[s]]$displaced$contacts_prop_extra = list(
    source = "Calculated from extrapolated contact matrix",
    data = setAgeBreaks(0) %>% .[, value := cm$prop_extra])
}
