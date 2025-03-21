###############################################################################
### IMPORT AND CLEAN SCRIPT                                                   #
###   ADNI participants Genetic Scores                                        #
###   Version 1                                                               #
###   December 2022                                                           #  
###   Albert Rodrigo Parés                                                    #
###---------------------------------------------------------------------------#
###   Tasks performed: 
###    * Import genetic PRS data from independent .txt files.
###    * Data cleaning and merge all of them in a single dataset
###                                                                           #
###############################################################################
rm(list=ls())
source("Rscripts/libraries.R")
source(here(RscriptsDir,"1_importAndClean","importAndClean_auxFunctions.R"))


#=======================================================================
# A) IMPORT PRS DATA FOR AD AND LONGEVITY
#=======================================================================

### Files names and paths----------------------------------
path <- here(rawdataDir,"geneticData/scores_genetics")
aux <- grep(list.files(path), pattern="^(.*)\\.txt") # select .txt files
files <- list.files(path,full.names = T)[aux]

### FUNCTION: import data-------------------------------
read_geneticScore <- function(file,threshold){
  #--------------------------------------
  # READ DATA
  #--------------------------------------
  score <- read.table(file, header=T)
  
  #-----------------------------------------------------------
  # Replace PTID current value by only its last 4-digit code 
  #-----------------------------------------------------------
  regex_PTID <- "^(\\d{3})_(\\w)_(\\d{4})$" # expected format
  
  #any(!grepl(regex_PTID,score$IID)) #is there any PTID with a different format?
  # we are good;all PTID are consistent format-wise.
  score$PTID <- as.factor(gsub(regex_PTID,"\\3",score$IID))
  
  #------------------------------------------------------
  # Keep Pt_5e.06 column and PTID
  #------------------------------------------------------
  aux <- names(score)[grep(paste0("(*.)",threshold,"$"),names(score))]
  score <- score[,c("PTID",aux)]
  
  #------------------------------------------------------
  # Get target name
  #------------------------------------------------------
  aux <- tail(strsplit(file,split="_")[[1]],1)
  targetName <- sub("(*.)\\.txt","\\1",aux)
  colnames(score) <- c("PTID",paste0(targetName,"_Pt_",threshold))
  
  # RETURN
  return(score)
}

### Apply import function to all files, and merge them by PTID -------------
ptid_scores <- lapply(files[-1],read_geneticScore,threshold="5e.06")
ptid_scores <- Reduce(function(x,y){merge(x,y,by="PTID")},ptid_scores)

### Apply import for AD files (differnet format) ---------------------
ptid_ADScores <- read_geneticScore(files[1],threshold="5e.06")[1:2]
names(ptid_ADScores) <- c("PTID","ADnoAPOE_Pt_5e.06")
ptid_scores <- merge(ptid_scores,ptid_ADScores,by="PTID")


### Rename variables and reorder ----------------------------------------
prs_names <- sub("(*.)(_Pt_5e.06)$",replacement = "\\1",
                 names(ptid_scores)[-1])
names(ptid_scores)[-1] <- prs_names
ptid_scores <- ptid_scores[,c("PTID","AD","ADnoAPOE","LONGEVITY","LONGEVITYnoAPOE",
                              "FTD","PKSON","EEAA","IEAA","FRAILTY","LEXP","TELOM")]

names(ptid_scores)[-1] <- paste("PRS",names(ptid_scores)[-1],sep="_")

#=======================================================================
# B) IMPORT FUNCTIONAL PRS for APOE - Plasma pQTL's
#=======================================================================
### Files names and paths----------------------------------
path <- here(rawdataDir,"geneticData/PRS-APOE-pQTLs-plasma")
aux <- grep(list.files(path), pattern="^(.*)\\.txt") # select .txt files
files <- list.files(path,full.names = T)[aux]

aux <- lapply(files,read_geneticScore,threshold="5e.06")
aux <- lapply(aux,function(x){names(x)<-c("PTID","PRS_APOE.pQTLPlasma");x})
plasma_pQTL_APOE <- do.call(rbind,aux)
plasma_pQTL_APOE <- plasma_pQTL_APOE[order(plasma_pQTL_APOE$PTID),]
rownames(plasma_pQTL_APOE) <- NULL

#=======================================================================
# C) IMPORT FUNCTIONAL PRS for APOE - Brain eQTL's
#=======================================================================
### Files names and paths----------------------------------
path <- here(rawdataDir,"geneticData/PRS-APOE-eQTLs-brain")
aux <- grep(list.files(path), pattern="^(.*)\\.txt") # select .txt files
files <- list.files(path,full.names = T)[aux]

aux <- lapply(files,read_geneticScore,threshold="1")
aux <- lapply(aux,function(x){names(x)<-c("PTID","PRS_APOE.eQTLBrain");x})
brain_eQTL_APOE <- do.call(rbind,aux)
brain_eQTL_APOE <- brain_eQTL_APOE[order(brain_eQTL_APOE$PTID),]
rownames(brain_eQTL_APOE) <- NULL

### MERGE ALL SCORES
geneticData <- Reduce(function(x,y) merge(x,y,by="PTID"),
                      list(ptid_scores,plasma_pQTL_APOE,brain_eQTL_APOE))



#=======================================================================
# D) IMPORT INDIVIDUAL SNPs WITH EFFECT ON BRAIN TRANSCRIPTOME
#=======================================================================
### Files names and paths----------------------------------
path <- here(rawdataDir,"geneticData/Individual-eQTLs-brain")
aux <- grep(list.files(path), pattern="^(.*)\\.txt") # select .txt files
file <- list.files(path,full.names = T)[aux]

### Read data
eQTL.APOE.brain.snps  <- read.table(file, header=T)

### GET 4 digit ID code
regex_PTID <- "^(\\d{3})_(\\w)_(\\d{4})$" # expected format
eQTL.APOE.brain.snps$PTID <- as.factor(gsub(regex_PTID,"\\3",eQTL.APOE.brain.snps$IID))

### Keep useful info and order by PTID
snps <- names(eQTL.APOE.brain.snps)[grep("^rs",names(eQTL.APOE.brain.snps))]
eQTL.APOE.brain.snps <- eQTL.APOE.brain.snps[c("PTID",snps)]
eQTL.APOE.brain.snps <- eQTL.APOE.brain.snps[order(eQTL.APOE.brain.snps$PTID),]
rownames(eQTL.APOE.brain.snps) <- NULL

### Add prefix 'eQTLBrain.'
pos <- names(eQTL.APOE.brain.snps) %in% snps
names(eQTL.APOE.brain.snps)[pos] <- paste("eQTLBrain",names(eQTL.APOE.brain.snps)[pos],sep=".")

### SNPs codified as integer (additive model)
eQTL.APOE.brain.snps[-1] <- lapply(eQTL.APOE.brain.snps[-1],as.integer)


#=======================================================================
# E) IMPORT INDIVIDUAL SNPs WITH EFFECT ON CSF 
#=======================================================================
### Files names and paths----------------------------------
path <- here(rawdataDir,"geneticData/Individual-pQTLs-CSF")
aux <- grep(list.files(path), pattern="^(.*)\\.txt") # select .txt files
file <- list.files(path,full.names = T)[aux]

### Read data
pQTL.APOE.CSF.snps  <- read.table(file, header=T)

### GET 4 digit ID code
regex_PTID <- "^(\\d{3})_(\\w)_(\\d{4})$" # expected format
pQTL.APOE.CSF.snps$PTID <- as.factor(gsub(regex_PTID,"\\3",pQTL.APOE.CSF.snps$IID))

### Keep useful info and order by PTID
snps <- names(pQTL.APOE.CSF.snps)[grep("^rs",names(pQTL.APOE.CSF.snps))]
pQTL.APOE.CSF.snps <- pQTL.APOE.CSF.snps[c("PTID",snps)]
pQTL.APOE.CSF.snps <- pQTL.APOE.CSF.snps[order(pQTL.APOE.CSF.snps$PTID),]
rownames(pQTL.APOE.CSF.snps) <- NULL

### Add prefix 'pQTLCSF.'
pos <- names(pQTL.APOE.CSF.snps) %in% snps
names(pQTL.APOE.CSF.snps)[pos] <- paste("pQTLCSF",names(pQTL.APOE.CSF.snps)[pos],sep=".")

### SNPs codified as integer (additive model)
pQTL.APOE.CSF.snps[-1] <- lapply(pQTL.APOE.CSF.snps[-1],as.integer)


### MERGE INDIVIDUAL SNPs with the rest of genetic data
geneticData <- Reduce(function(x,y) merge(x,y,by="PTID"),
                      list(geneticData,eQTL.APOE.brain.snps,pQTL.APOE.CSF.snps))

head(geneticData)
#=======================================================================
# D) Normalize score value: Z-score center+scale
#=======================================================================
geneticData[2:14] <- lapply(geneticData[2:14],scale)
geneticData[2:14] <- lapply(geneticData[2:14],as.numeric)

#=======================================================================
# E) Histograms
#=======================================================================

# PRS
file <- here("results","Others","PRS_histograms.png")
png(file,width = 800,height = 480)
par(mfrow=c(2,4))
aux <- geneticData[c(2:5,12:14)]
lr <- lapply(1:ncol(aux),function(i){
  h <- hist(aux[,i], xlab="Score value",ylab="",main=colnames(aux)[i],freq=F)
})
dev.off()

# SNPs
file <- here("results","Others","SNPs_barplots.png")
png(file,width = 800,height = 480)
par(mfrow=c(2,5))
aux <- geneticData[15:24]
lr <- lapply(1:ncol(aux),function(i){
  tab <- table(aux[,i])
  b <- barplot(table(aux[,i]), xlab="Number of alleles",
               main = colnames(aux[i]),ylim = c(0,max(table(aux[,i]))+10))
})
dev.off()


#=======================================================================
# C) SAVE clean dataset
#=======================================================================
saveRDS(geneticData,
        file = here(cleandataDir,"geneticData.RDS"))

rm(list=ls())
