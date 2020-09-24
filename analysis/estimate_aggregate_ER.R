
if(!require("here")) {
  install.packages("here")
  library("here")
}
if(!require("readr")) {
  install.packages("readr")
  library("readr")
}

install.packages("tidyverse")
library("tidyverse")
install.packages("magrittr")
library("magrittr")

datadir <- here("data")

AEQ_harvest_SF <-readRDS(paste(datadir,"d_r62_tamm_lmtstkcmplt_no_m_nat_1992_2016.rds",sep = "/")) %>% 
  select(stk_type:fisheryID, contains("Skagit SF")) %>% 
  mutate(tot_aeq = select(., contains("_aeq")) %>% rowSums(na.rm = T)) %>%
  filter(fisheryID != 74) %>% #catch only, leave out for total abundance
  group_by(year) %>% summarise(tot_aeq = sum(tot_aeq, na.rm = T))


write_csv(AEQ_harvest_SF, file.path(datadir, "skagit_sf_chinook_aeq.csv"))


