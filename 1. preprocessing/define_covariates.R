# This file constructs data frame with number of rows as in Betula connectivity data  ('FC_Tanya.mat')
# with the variable "Inclusion" that equals to 1 for rows to keep and 0 for rows to ignore in the analyses.
# The file "Age_Sex_FD_stand_cov.xlsx" with all the covariates to be used in the analyses is also constructed here.

library(xlsx)
library(tidyverse)

# reading cognitive data --------------------------------------------------
file.cognition <- "Betula behavioral data to Tanya update8.xlsx"
file.block.design <- "Block Design ImAGen T1-T6 2015-12-14.xlsx"
file.inclusions <- "Connectome info till Tanya.xlsx"
file.fd <- "FD across the full sample(T5T6T7).xlsx"

data.cognition <- read.xlsx2(
  file = file.cognition,
  1, startRow = 2, endRow = 378, stringsAsFactors = F
)

data.block <- xlsx::read.xlsx(file.block.design, sheetIndex = 1) %>%
  select(Unique.number = Unique.no, T5.Block)
data.cognition <- merge(data.cognition, data.block, by = "Unique.number")
data.fd <- readxl::read_excel(file.fd)

# reading indicators of inclusions ----------------------------------------
inclusion <- readxl::read_excel(file.inclusions)
inclusion <- merge(inclusion, data.fd)

inclusion.long <- inclusion %>%
  gather(id, score, c(FD, Age, includeExclude)) %>%
  unite(temp1, id, Time, sep = ".") %>%
  spread(temp1, score) %>% # transformed to wide
  mutate(
    include.10 = 0,
    include.0 = as.numeric(includeExclude.0 & !Subject_Code_Unique == "B2502.65" & FD.0 < 0.5), # subject B2502.65 has missing cognition at T5
    include.5 = as.numeric(include.0 & includeExclude.5 & !is.na(includeExclude.0) & FD.5 < 0.5)
  ) %>%
  select(Subject_Code_Unique, starts_with("include.")) %>%
  gather(Variable, Inclusion, starts_with("include.")) %>% # transform to long
  separate(Variable, c("idd", "Time")) %>%
  select(-idd)
inclusion <- merge(inclusion, inclusion.long, by = c("Subject_Code_Unique", "Time"), all.x = T, all.y = F)

# subsetting cognitive data -----------------------------------------------
d <- data.cognition %>%
  select(
    MR.subject.ID, T5.SPT, T5.SPT.crc, T5.VT, T5.VT.crc, T5.Wrk00,
    T5.Fluency.A, T5.Fluency.M5, T5.Fluency.Prof.B,
    T5.Lett.comp, T5.Fig.comp, T5.LetterDigit,
    T5.Block,
    Age.cohort.T5, Sex
  ) %>%
  filter(MR.subject.ID %in% inclusion$Subject_Code_Unique) %>%
  mutate_at(vars(T5.SPT:Sex), as.numeric) %>%
  mutate(Sex = Sex - 1) # select cognitive data for those in the subsetted Alireza's fil
d <- merge(inclusion %>%
  select(MR.subject.ID = Subject_Code_Unique, Time, includeExclude, Inclusion, FD), d)
d[d$Inclusion==0, c("T5.SPT", "T5.SPT.crc", "T5.VT", "T5.VT.crc", "T5.Wrk00","T5.Fluency.A",
                    "T5.Fluency.M5", "T5.Fluency.Prof.B", "T5.Lett.comp",
                    "T5.Fig.comp", "T5.LetterDigit", "Sex", "Age.cohort.T5", "T5.Block", "FD")] <- NA

# define T5 and T6 scores -------------------------------------------------
# episodic memory
mean.em <- colMeans(d %>%filter(Inclusion ==1 & Time==0)%>%select(T5.SPT, T5.SPT.crc, T5.VT, T5.VT.crc, T5.Wrk00))
sd.em <- apply(d %>%filter(Inclusion ==1 & Time==0)%>%select(T5.SPT, T5.SPT.crc, T5.VT, T5.VT.crc, T5.Wrk00), 2, sd)
d$T5.EM1.comp <- rowSums(scale(d[, c("T5.SPT", "T5.SPT.crc", "T5.VT", "T5.VT.crc", "T5.Wrk00")], 
                               center = mean.em, scale=sd.em))
# word fluency
mean.fl <- colMeans(d %>%filter(Inclusion ==1 & Time==0)%>%select(T5.Fluency.A, T5.Fluency.M5, T5.Fluency.Prof.B))
sd.fl <- apply(d %>%filter(Inclusion ==1 & Time==0)%>%select(T5.Fluency.A, T5.Fluency.M5, T5.Fluency.Prof.B), 2, sd)
d$T5.Fl1.comp <- rowSums(scale(d[, c("T5.Fluency.A", "T5.Fluency.M5", "T5.Fluency.Prof.B")], 
                               center= mean.fl, scale = sd.fl))

# processing speed
mean.ps <- colMeans(d %>%filter(Inclusion ==1 & Time==0)%>%select(T5.Lett.comp, T5.Fig.comp, T5.LetterDigit))
sd.ps <- apply(d %>%filter(Inclusion ==1 & Time==0)%>%select(T5.Lett.comp, T5.Fig.comp, T5.LetterDigit), 2, sd)
d$T5.PS1.comp <- rowSums(scale(d[, c("T5.Lett.comp", "T5.Fig.comp", "T5.LetterDigit")],
                               center = mean.ps, scale = sd.ps))

# saving ------------------------------------------------------------------
# write.xlsx2(d %>% select(MR.subject.ID, Time, Inclusion, Age.cohort.T5, Sex, T5.EM1.comp, T5.Fl1.comp, T5.PS1.comp, T5.Block),
#   file = "ID_Inclusion_Age_Sex_FD_stand_cov.xlsx",
#   row.names = F
# )
# 
# write.xlsx2(d %>%filter(Inclusion ==1 & Time==0)%>% select(Age.cohort.T5, Sex, FD, T5.EM1.comp, T5.Fl1.comp, T5.PS1.comp, T5.Block),
#   file = "Age_Sex_FD_stand_cov.xlsx",
#   row.names = F, col.names = F
# )
