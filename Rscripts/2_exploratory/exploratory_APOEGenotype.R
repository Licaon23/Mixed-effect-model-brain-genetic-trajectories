source("Rscripts/libraries.R")
load(here(cleandataDir,"root.RData"))

load(here(cleandataDir,"dataForAnalysis_3T.RData"))
paramT <- "Dataset_3T"

# load(here(cleandataDir,"dataForAnalysis_1T.RData"))
# paramT <- "Dataset_1T"

#============================================================================
# B) DESCRIPTIVA CARACTERISTIQUES BASALS DELS PACIENTS
#============================================================================
sampleSubjects <- basal$PTID

genotype <- demo[c("RID","genotype","DX.bl")]
genotype$PTID <- sprintf("%04d", genotype$RID)

genotype <- subset(genotype,PTID %in% sampleSubjects)

tab <- with(genotype,table(genotype,DX.bl))
count <- addmargins(tab)
freqRel <- apply(prop.table(tab,1),2,round,3)


library(kableExtra)
count %>% kbl(caption = "APOE genotype by DX.bl. ADNI-1 1.5T sample") %>% 
  kable_classic(full_width=F,html_font = "Cambria") %>%
  footnote(general = "Number of cases for a given diagnostic group and
           APOE genotype")

freqRel %>% kbl(caption = "APOE genotype by DX.bl. ADNI-1 1.5T sample") %>% 
  kable_classic(full_width=F,html_font = "Cambria") %>%
  footnote(general = "Proportion of individuals with different diagnostic group for each genotype")

