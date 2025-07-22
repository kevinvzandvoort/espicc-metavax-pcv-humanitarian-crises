#' Get synthetic contact data from Prem et al
if(!file.exists("./data/demography/prem_synthetic_contact_matrices.RDS")){
  #' do this in its own environment
  (function(){
    x = tempfile(fileext = ".rdata")
    download.file("https://github.com/kieshaprem/synthetic-contact-matrices/raw/master/output/syntheticmatrices/contact_all.rdata",
                  x, mode = ifelse(Sys.info()["sysname"] == "Windows", "wb", "w"))
    load("C:/Users/Kevin/Downloads/contact_all.rdata")
    saveRDS(contact_all[c("ETH", "SSD", "NGA", "CAF")], "./data/demography/prem_synthetic_contact_matrices.RDS")
  })() 
}

readRDS("./data/demography/prem_synthetic_contact_matrices.RDS") %>%
  .[["SSD"]] %>%
  as.data.table %>%
  setNames(as.character(setAgeBreaks(seq(0, 75, 5))[, name])) %>%
  .[, contactor_age_group := setAgeBreaks(seq(0, 75, 5))[, name]] %>%
  melt(id.vars = "contactor_age_group", variable.name = "contactee_age_group") %>%
  .[, c("contactor_age_group", "contactee_age_group") :=
      .(factor(contactor_age_group, setAgeBreaks(seq(0, 75, 5))[, name]),
        factor(contactee_age_group, setAgeBreaks(seq(0, 75, 5))[, name]))] %>%
  .[]