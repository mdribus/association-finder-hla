# DESCRIPTION:

# This script is the second step in the disease association pipeline. It accomplishes two tasks:
# 1: performing FDR correction for the P values
# 2: identifying the suggestive HLA variants
# It produces output in the /FDR and /Summary directories.

if(!require("fdrtool"))
  install.packages("fdrtool")
if(!require("dplyr"))
  install.packages("dplyr")

library(fdrtool)
library(dplyr)

# Section 1: perform FDR correction

columns <- read.csv('./disease_columns.csv') # read the disease columns file
columns$analysis <- paste0(columns$locus,'_',columns$type) # combine "locus" and "type" in the disease columns file
offset <- length(which(columns$type=="demo")) + 1 # the offset is where the demographic data ends
ANALYSIS_NAMES <- as.character(unique(columns$analysis[offset:nrow(columns)])) # uncomment this line to run on Cypress

# define command line arguments for disease and population
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0) {
  thisDisease <- args[1]
  pop <- args[2]
}

# set up print statements for the arguments:
print(paste0("This disease: ",thisDisease))
POPULATIONS = as.list(strsplit(pop, ",")[[1]]) # formats the population argument
print(paste0("Population(s): ",POPULATIONS))

for (pop in 1:length(POPULATIONS)){
  for (analysis in 1:length(ANALYSIS_NAMES)) {

    #### load the association files that were generated previously
    thisAnalysis <- ANALYSIS_NAMES[analysis]
    popName<-POPULATIONS[pop]

    orName=paste0(getwd(),"/output/",thisDisease,"/Multi/",popName,"_",thisAnalysis,"_Multi_ORstats.csv")

    run=TRUE
    
    if (file.exists(orName)){
      orData<-read.csv(file=orName,header = T)
      #print(orData)
      if (nrow(orData)<1){
        run=FALSE
      }
    } else {
      run=FALSE
      print(paste0("Odds ratio (OR) file missing: ",orName))
    }
    
    if (run) {

      orData=orData[which(orData$numTested==5),]
      
      if(nrow(orData)>0){
        
        chosenAlleles=as.character(unique(orData$X))
        pvalues=orData$P.value
        
        sortedTable=orData[order(pvalues),]
        
        I=c(1:nrow(sortedTable))
        Q=(I*0.05)/nrow(sortedTable)
        
        sortedTable=cbind(sortedTable,Q)
        
        significant=character()
        
        check= FALSE
        
        adjLCL_OR=numeric()
        adjUCL_OR=numeric()
        
        for (i in 1:nrow(sortedTable)){
          
          p=sortedTable$P.value[i]
          q=sortedTable$Q[i]
          
          adjLCL_OR[i]=exp(log(sortedTable$OR[i]) - qnorm(1-q/2)*sortedTable$SE[i])
          adjUCL_OR[i]=exp(log(sortedTable$OR[i]) + qnorm(1-q/2)*sortedTable$SE[i])
          
          if (is.na(p)){
            significant=c(significant,"-")
          } else {
            if (p<q){
              significant=c(significant,"***")
              check=TRUE
            } else {
              significant=c(significant,"-")
            }
          }# end na check
          
        }# end i loop
        
        adjustedP=t(t(p.adjust(sortedTable$P.value, method="fdr", n=length(sortedTable$P.value))))
        
        if (check){
          maxSig=max(which(significant=="***"))
          significant[1:maxSig]="***"
          
          sortedTable=cbind(sortedTable,significant,adjustedP,adjLCL_OR,adjUCL_OR)
          
          thisSignificant=sortedTable[1:maxSig,]
          
          if (!exists("allSignificant")){
            allSignificant=thisSignificant
          } else {
            allSignificant=rbind(allSignificant, thisSignificant)
          }
          
          write.csv(sortedTable,file=paste0(getwd(),"/output/",thisDisease,"/FDR/",popName,"_",thisAnalysis,"_FDR.csv"))
      
        } else {
          
          sortedTable=cbind(sortedTable,significant,adjustedP,adjLCL_OR,adjUCL_OR)
          write.csv(sortedTable,file=paste0(getwd(),"/output/",thisDisease,"/FDR/",popName,"_",thisAnalysis,"_FDR.csv"))
          }
        }## check
      }##row check
    }# end run check
  
  print(paste0(popName, " allSignificant:"))
  print(allSignificant)
  
  if (exists("allSignificant")){
    write.csv(allSignificant,file=paste0(getwd(),"/output/",thisDisease,"/FDR/",popName,"_FDR_allSignificant.csv")) 
    } else {
      print("allSignificant doesn't exist") 
    }
}
#######################################################################################

# find the suggestive variants

for (pop in 1:length(POPULATIONS)){
  multi_sig_list = list() # this list will be filled with significant (P.value <= 0.05) variants from Multi_ORstats.csv files
  fdr_all_variants_list = list() # this list will contain all the variants and their accompanying statistics after FDR correction
  
  for (analysis in 1:length(ANALYSIS_NAMES)) { # evaluate each analysis separately before combining them in the list

    #### load the files
    thisAnalysis <- ANALYSIS_NAMES[analysis]
    popName<-POPULATIONS[pop]
    
    orName=paste0(getwd(),"/output/",thisDisease,"/Multi/",popName,"_",thisAnalysis,"_Multi_ORstats.csv")
    
    run=TRUE
    
    if (file.exists(orName)){
      orData<-read.csv(file=orName,header = T) # read the Multi_ORstats.csv files
      orData<-filter(orData, P.value <=0.05) # only keep the significant variants
      multi_sig_list[[analysis]] <- orData # add the significant variants to a list
      
      if (nrow(orData)<1){
        run=FALSE
      }
    } else {
      run=FALSE
      print(paste0("Odds ratio (OR) file missing: ",orName))
    }

    fdrName=paste0(getwd(),"/output/",thisDisease,"/FDR/",popName,"_",thisAnalysis,"_FDR.csv")
    
    run=TRUE
    
    if (file.exists(fdrName)){
      fdrData<-read.csv(file=fdrName,header = T) # read the Multi_ORstats.csv files
      fdrData$X.1 <- NULL
      fdr_all_variants_list[[analysis]] <- fdrData
      
      if (nrow(fdrData)<1){
        run=FALSE
      }
    } else {
      run=FALSE
      print(paste0("FDR file missing: ",fdrName))
    }
  }
  
  multi_sig_df = do.call(rbind, multi_sig_list) # convert the list into a data frame for the variants that are significant before FDR
  fdr_all_variants_df = do.call(rbind, fdr_all_variants_list) # converts the list into a data frame for all the FDR metrics
  
  fdrName_all_sig=paste0(getwd(),"/output/",thisDisease,"/FDR/",popName,"_FDR_allSignificant.csv") # identify the FDR_allsignificant.csv files

  if (file.exists(fdrName_all_sig)){
    fdrData_all_sig<-read.csv(file=fdrName_all_sig,header=T) # read the FDR_allsignificant.csv files
    # filter the multi_sig_df; variants that remain significant after FDR correction are removed; identifies suggestive variants
    multi_sig_df_filtered <- filter(multi_sig_df, !(multi_sig_df$X %in% fdrData_all_sig$X))
    multi_sig_df_filtered_final <- filter(fdr_all_variants_df, (fdr_all_variants_df$X %in% multi_sig_df_filtered$X))
    
    names(multi_sig_df_filtered_final)[names(multi_sig_df_filtered_final) == "X"] <- "Variant" # changes the column name
    
    # write the csv file to the summary folder
    write.csv(multi_sig_df_filtered_final,file=paste0(getwd(),"/output/",thisDisease,"/Summary/",popName,"_suggestive_variants.csv"), row.names=FALSE)
  } else {
    run=FALSE
    print(paste0("FDR file missing: ",fdrName_all_sig))
  }
}