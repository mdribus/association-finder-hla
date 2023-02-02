# DESCRIPTION:

# This script is the first step in the disease association pipeline. 
# It accomplishes two tasks:
# 1: performing logistic regression to identify associations
# 2: calculating summary statistics for each analysis category 
# It produces output in the /Multi and /Summary subdirectories.

# SECTION 1: Prepare the environment, load the analysis categories, and define the arguments       

if(!require("data.table"))
  install.packages("data.table")
if(!require("dplyr"))
  install.packages("dplyr")
if(!require("broom"))
  install.packages("broom")
if(!require("purrr"))
  install.packages("purrr")
if(!require("tibble"))
  install.packages("tibble")

library("tibble")
library("purrr")

# these two lines prevent the error files from displaying this warning: 
# `summarise()` ungrouping output (override with `.groups` argument)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)

# clear out the environment and run garbage collection
rm(list=ls(all=TRUE)) 
invisible(gc()) # invisible hides gc output

# read disease_columns.csv to obtain analysis categories
columns <- read.csv('./disease_columns.csv')

# create new analysis column by concatenating locus and type
columns$analysis <- paste0(columns$locus,'_',columns$type)

# specifies which columns are HLA categories to be tested: the offset is where the demographic data ends
offset <- length(which(columns$type=="demo")) + 1
ANALYSIS_NAMES <- as.character(unique(columns$analysis[offset:nrow(columns)])) # uncomment this line to run on Cypress

# define command line arguments for disease, population, imputation number, and top limit
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0) {
  thisDisease <- args[1] # the format will depend on the name of your input file
  pop <- args[2] # the format will depend on the name of your input file
  imputation_replicates <- args[3] # tells the script how many imputation replicates you're using
  top_limit <- args[4] # tells the script how many HLA variants to include in each analysis category
}

# set up print statements for the arguments:
print(paste0("This disease: ",thisDisease))
POPULATIONS = as.list(strsplit(pop, ",")[[1]]) # formats the population argument
print(paste0("Population(s): ",POPULATIONS))
print(paste0("Number of imputation replicates to be analyzed: ", imputation_replicates))
print(paste0("Number of HLA variants to be analyzed in each category: ", top_limit))

# SECTION 2: Run analysis loops                                  

for (pop in 1:length(POPULATIONS)){ # loops through every population
  popName<-POPULATIONS[pop]
  data <- NULL
  for(f in 1:imputation_replicates) { # loops through every imputation file
    fileName=paste0(getwd(),"/HLA_input_",thisDisease,"_",popName,"_",as.character(f),".txt")
    repdata <-fread(fileName,header = T, stringsAsFactors=FALSE,data.table = TRUE)
    repdata$replicate <- f # add replicate column
    data <- bind_rows(data, repdata)
    print(paste("Data Import of", fileName,"completed."))
  }

  print("Finding the associations...")

  # begin analysis loop
  for (analysis in 1:length(ANALYSIS_NAMES)) {
    thisAnalysis <- ANALYSIS_NAMES[analysis] # defines all the analysis categories
    glm_output_df = NULL # define the empty dataframe to fill with glm function output

    # multiple imputation replicates
    for(f in 1:imputation_replicates) { # loops through the different imputation replicates
      thisrepdata <-data[which(data$replicate==f),] # subset replicate
      modelColumns <- thisrepdata[,1:6,with=FALSE]

      ### select analysis columns
      columnInds <- which(columns$analysis==thisAnalysis)
      thisType <- as.character(columns$type[columnInds][1])
      numColumns <- length(columnInds) # number of analysis columns (1 if genotype, 2 if allele or haplotype)
      tempColumns<-thisrepdata[,.SD,.SDcols=columnInds] # analysis columns for this Analysis  - data.table
      analysisColumns <- cbind(modelColumns,tempColumns) # model columns combined with 1 or 2 analysis columns
      if(numColumns==2){inds=7:8} else {inds=7} # inds give column numbers for analysis
      names(analysisColumns)[inds]<- paste0('Allele',1:numColumns) # give column heading for analysis columns - Allele1, (Allele2)
      
      # unique variant names - /Multi/ files key off of this - not impacted by topLimit
      uniAlleles <- as.character(unique(as.vector(t(tempColumns))))

      ### set limit of alleles if needed
      
      # number of HLA variants per analysis
      topLimit<-as.numeric(top_limit)

      fullInds<-c(1:nrow(analysisColumns))

      ### frequency counts
      if (numColumns==2){ # alleles and haplotypes
        fullAllele<-c(as.character(analysisColumns$Allele1),as.character(analysisColumns$Allele2)) # array of all allele names
        fullIndicator<-c(analysisColumns$Disease_Ind,analysisColumns$Disease_Ind) # array of all disease indicators
      } else if (numColumns==1) { # genotypes
        fullAllele<-as.character(analysisColumns$Allele1)
        fullIndicator<-analysisColumns$Disease_Ind
      }
      
      frequency<-table(fullAllele) # overall frequency
      frequency<-frequency[order(frequency,decreasing=T)] # sort the frequency
      frequency2<-table(fullAllele,fullIndicator) # frequency in cases and controls

      fullcount<-nrow(frequency2) # number of HLA variants
      rank<-c(1:fullcount) 
      margin=margin.table(frequency2, 1) # A frequencies (summed over B) - # sum of case and control freqs
      
      prop=as.numeric()
  
      totalCounts=cbind(sum(frequency2[,1]),sum(frequency2[,2])) # case and control total counts
      
      # calculate case and control frequencies
      for (propLoop in 1:nrow(frequency2)){
        thisProp=frequency2[propLoop,]/totalCounts
        prop=rbind(prop,thisProp)
      }
    
      completeFreq=cbind(frequency2,prop,margin) # table with case and control counts, frequencies, total count

      sorted.completeFreq<-completeFreq[order(margin,decreasing=T),] # sort table by total count
      sorted.completeFreq=cbind(sorted.completeFreq,rank) # add rank column
      colnames(sorted.completeFreq)=c('ControlCount','CaseCount','ControlFreq','CaseFreq','Total','TotalRank') # rename columns

      completeFreq2<-as.data.frame(completeFreq) # coverts completeFreq into a data frame
      completeFreq2 <- tibble::rownames_to_column(completeFreq2, "Allele") # moves alleles from index to first column

      # topsorted are most frequent alleles below topLimit threshold
      if (nrow(sorted.completeFreq)>topLimit){
        topsorted.completeFreq=sorted.completeFreq[1:topLimit,]
      } else {
        topsorted.completeFreq=sorted.completeFreq
      }
      
      topsorted.completeFreq<-cbind(rownames(topsorted.completeFreq),topsorted.completeFreq)
      
      if (f==1) {
        allFrequency<-topsorted.completeFreq
      } else {
        allFrequency<-rbind(allFrequency, topsorted.completeFreq) # adds rows for multiple imputation replicate
      }
      
      ### set analysis alleles based on first file alleles
      allAllele=rownames(topsorted.completeFreq)
      
      if (f==1){ # first file only
        chosenAllele=sort(allAllele)
      }
      
      fileCoeffs={}
      fileVcovs={}
      
      if (length(chosenAllele)!=0) { # chosenAllele is set to this
        
        for (i in 1:length(chosenAllele)){ # chosenAllele length depends on the top limit

          currAllele<-as.character(chosenAllele[i]) # which HLA variant - e.g. A*01:01g
          rowforcurrentallele<-filter(completeFreq2, Allele == currAllele)
          controlfrequencyforcurrentallele<-rowforcurrentallele$V3
          
          ALLELE<-numeric(nrow(analysisColumns)) # indicator for which rows have allele - logical OR if two columns
          
          if (numColumns==2){
            b<-analysisColumns[,.I[Allele1==chosenAllele[i]|Allele2==chosenAllele[i]]]
            ALLELE[b]<-1
          } else if (numColumns==1){
            b<-analysisColumns[,.I[Allele1==chosenAllele[i]]]
            ALLELE[b]<-1
          }
          fileNum<-f
          
          ### cut off checks for minimum number of cases
          
          xtabsAllele<-analysisColumns[, .N, by = .(ALLELE,Disease_Ind)]
          
          checkAllele<-xtabsAllele[,min(N)]
          checkAllele2<-xtabsAllele[ ,sum(N>0)]
          
          ### if there are fewer than 5 in cell for allele
          if (checkAllele<5){
            check=0
          } else {
            check=1
          }

          ### if there are fewer than 4 rows in disease check
          if (checkAllele2<4){
            check=0
          } else {
            check=1
          }
          
          ### run glm function to perform logistic regression for variant if checks pass
          if (check>0){
            fullModel<-glm(analysisColumns$Disease_Ind~ALLELE,family=binomial) # this is the logistic regression formula
            OR <- as.data.frame(exp(cbind(OR = coef(fullModel), suppressMessages(confint(fullModel))))) # calculates the odds ratios for terms in fullModel
            OR <- format(OR, scientific = FALSE) # disables scientific notation for the OR df
            colnames(OR) <- c("OR","LCI","UCI") # defines the column names in the OR df
            OR <- lapply(OR, as.numeric) # makes the OR df numeric
            OR$LCI <- round(OR$LCI, digits=6) # rounds the LCI column to 6 digits
            OR$UCI <- round(OR$UCI, digits=6) # rounds the UCI column to 6 digits
            tidy_glm <- tidy(fullModel) # tidies the results
            tidy_glm$OR <- OR$OR # adds the OR column to tidy_glm
            tidy_glm$LCI <- OR$LCI # adds the LCI column to tidy_glm
            tidy_glm$UCI <- OR$UCI # adds the UCI column to tidy_glm
            tidy_glm[2,1] = currAllele # adds the current allele to tidy_glm
            tidy_glm <- tidy_glm[-1,] # gets rid of the (Intercept) row
            glm_output_df <- rbind(glm_output_df, data.frame(tidy_glm)) # fills dataframe with glm results from all 5 imputations for each category
          }
        }
      }
    }

    if (is.null(glm_output_df)){
      next
      
    } else {
      
      ###### create output files for association results
      summarized_results <- glm_output_df %>% # calculates means for all the terms in tidy_glm_list
        group_by(term) %>% # groups by the allele/haplotype/category
        summarize(estimates=mean(estimate), # calculates mean for each column
                  std.error=mean(std.error), 
                  z.value=mean(statistic),
                  p.value=mean(p.value),
                  OR=mean(OR), 
                  LCI=mean(LCI), 
                  UCI=mean(UCI)) %>% as.data.frame
      summarized_results_filtered = subset(summarized_results, select = -c(estimates, z.value) ) # drops unwanted columns
      colnames(summarized_results_filtered) <- c("Variant","SE","P value","OR","LCI","UCI") # defines column names
      summarized_results_filtered$numTested = imputation_replicates # adds numTested column, which is "5" because there are 5 imputations
      summarized_results_filtered <- summarized_results_filtered[, c("Variant","OR","LCI","UCI","P value","numTested","SE")] # defines column order
      colnames(summarized_results_filtered) <- c("","OR","LCI","UCI","P value","numTested","SE") # removes "Variant" column name
      write.csv(summarized_results_filtered, file = paste0(getwd(),"/output/",thisDisease,"/Multi/",popName,"_",thisAnalysis,"_Multi_ORstats.csv"), row.names=FALSE)
      print(summarized_results_filtered)
    }

    if (exists("fileCoeffs")){
      aveCount1<-numeric()
      aveCount0<-numeric()
      aveFreq1<-numeric()
      aveFreq0<-numeric()
      allAlleles<-unique(allFrequency[,1]) # allAlleles is limited to topLimit alleles, unlike /Multi/ which is allCoeffs, uniAlleles
      
      ### get averages and frequencies of each allele across replicates
      for (w in 1:length(allAlleles)){
        subFreq<-allFrequency[which(allFrequency[,1]==allAlleles[w]),]
        rowLen<-rownames(allFrequency)[which(allFrequency[,1]==allAlleles[w])]
        if (length(rowLen)>1){
          aveCount1<-rbind(aveCount1, mean(as.numeric(subFreq[,3])))
          aveCount0<-rbind(aveCount0, mean(as.numeric(subFreq[,2])))
          aveFreq1<-rbind(aveFreq1, mean(as.numeric(subFreq[,5])))
          aveFreq0<-rbind(aveFreq0, mean(as.numeric(subFreq[,4])))
        }else if(length(rowLen)==1){
          aveCount1<-rbind(aveCount1, mean(as.numeric(subFreq[3])))
          aveCount0<-rbind(aveCount0, mean(as.numeric(subFreq[2])))
          aveFreq1<-rbind(aveFreq1, mean(as.numeric(subFreq[5])))
          aveFreq0<-rbind(aveFreq0, mean(as.numeric(subFreq[4])))
        }
      }
      rownames(aveCount0)<-allAlleles
      rownames(aveCount1)<-allAlleles
      rownames(aveFreq1)<-allAlleles
      rownames(aveFreq0)<-allAlleles
      
      allCounts<-as.data.frame(cbind(aveCount1,aveFreq1,aveCount0,aveFreq0))
      names(allCounts)<-c("aveCount1","aveFreq1","aveCount0","aveFreq0")
      
      #### sort and rank by the 0 and 1 counts to compare and write to file
      
      sorted.Count0<-allCounts[order(allCounts$aveCount0,decreasing=T),]
      rank0=c(1:length(allAlleles))
      sorted.Count0.Rank0=cbind(sorted.Count0,rank0)
      
      sorted.Count0.Rank0<-sorted.Count0.Rank0[order(sorted.Count0.Rank0$aveCount1,decreasing=T),]
      rank1=c(1:length(allAlleles))
      sorted.Count0.Rank01=cbind(sorted.Count0.Rank0[,c(1:2)],rank1,sorted.Count0.Rank0[,c(3:5)])
      write.csv(sorted.Count0.Rank01,file=paste0(getwd(),"/output/",thisDisease,"/Summary/",popName,"_",thisAnalysis,"_summaryFile.csv"))
    }
  }
}
