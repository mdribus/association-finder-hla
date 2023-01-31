# DESCRIPTION:

# This script is the fourth step in the disease association pipeline:
# It organizes the factor analysis output into master files in /Masters directory

columns <- read.csv('./disease_columns.csv') # read the disease columns file
columns$analysis <- paste0(columns$locus,'_',columns$type) # combine "locus" and "type" in the disease columns file
offset <- length(which(columns$type=="demo")) + 1 # the offset is where the demographic data ends
ANALYSIS_NAMES <- as.character(unique(columns$analysis[offset:nrow(columns)])) # uncomment this line to run on Cypress

# define command line arguments for disease, population, imputation number, and top limit
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
  popName<-POPULATIONS[pop]

  faName1=paste0(getwd(),"/output/",thisDisease,"/Multi/",popName,"_Grouping.csv")
  faName2=paste0(getwd(),"/output/",thisDisease,"/Multi/",popName,"_Geno_Grouping.csv")
  if (file.exists(faName1)){ # if the directory above exists, read the _Grouping.csv file
    faData1<-read.csv(file=faName1,header = T)
    groupCheck1=TRUE
  } else {
    groupCheck1=FALSE
    print("no fa file")
    faData1={}
  }
  if (file.exists(faName2)){ # if the directory above exists, read the _Geno_Grouping.csv file
    faData2<-read.csv(file=faName2,header = T)
    groupCheck2=TRUE
  } else {
    groupCheck2=FALSE
    faData2={}
  }
  
  if (groupCheck1 | groupCheck2){
    groupCheck=TRUE
  } else {
    groupCheck=FALSE
  }
  
  faData=rbind(faData1,faData2)
  
  grouped={}
  groups={}
  loads={}
  fullMaster={} 
  fullAlleles={}

  if(groupCheck) {
    filegrouped=faData$Allele
    
    for (g in 1:length(filegrouped)){
      
      thisGroup=as.character(filegrouped[g])
      thisNum=as.character(faData$Group.Assignment[g])
      thisLoad=faData$Loading[g]
      
      if (length(grep("%%",thisGroup))>0){
        splits=strsplit(as.character(thisGroup)," %% ")
        splits=splits[[1]]
        numsplit=length(splits)
        theseGroup=as.character(rep(thisNum,numsplit))
        theseLoad=as.character(rep(thisLoad,numsplit))
        grouped=c(grouped,splits)
        groups=c(groups,theseGroup)
        loads=c(loads,theseLoad) 
      } else {
        grouped=c(grouped,thisGroup)
        groups=c(groups,thisNum)
        loads=c(loads,thisLoad)
      }
      
    }## end g

  for (analysis in 1:length(ANALYSIS_NAMES)) {
    thisAnalysis <- ANALYSIS_NAMES[analysis]

    #### load the files
    freqName=paste0(getwd(),"/output/",thisDisease,"/Summary/",popName,"_",thisAnalysis,"_summaryFile.csv")
    fdrName=paste0(getwd(),"/output/",thisDisease,"/FDR/",popName,"_",thisAnalysis,"_FDR.csv")
    
    run=TRUE
    
    if (file.exists(freqName)){
      freqData<-read.csv(file=freqName,header = T,colClasses=c("character","numeric","numeric","numeric","numeric","numeric","integer"))
    } else {
      freqData={}
      run=FALSE
    }

    if (file.exists(fdrName)){
      fdrData<-read.csv(file=fdrName,header = T)
    } else {
      fdrData={}
      run=FALSE
    }
    
    if (run) { 
      
      sigData=fdrData[which(fdrData$significant=="***"),]
      chosenAlleles=unique(as.character(sigData$X))

      masterFile=as.data.frame(matrix(0,nrow=length(chosenAlleles), ncol=13)) # changed 16 to 13
      names(masterFile)=c("OR","FDR LCI", "FDR UCI", "FDR p","Uncorrected p","Group","Loading","CaseRank","ControlRank","CaseCount","ControlCount","CaseFreq","ControlFreq")
      
      # masterFile[,0]=chosenAlleles
      masterFile[,1]=sigData$OR # fills the first column
      masterFile[,2]=sigData$adjLCL_OR # fills the second column
      masterFile[,3]=sigData$adjUCL_OR # fills the third column
      masterFile[,4]=sigData$adjustedP # fills the fourth column
      masterFile[,5]=sigData$P.value # fills the fifth column
      
      if (length(chosenAlleles)>0) {
        for (i in 1:length(chosenAlleles)) {
          
          loads=as.numeric(loads)
          
          grpind=match(chosenAlleles[i],grouped)
          
          if (!is.na(grpind)){
            masterFile[i,7]=loads[grpind]
            if (abs(loads[grpind])<0.02){
              masterFile[i,6]=paste(groups[grpind],"*")
            } else {
              masterFile[i,6]=groups[grpind]
            }
            
          } else {
            masterFile[i,6]="NA"
            masterFile[i,7]="NA"
          }
        
          grpind="NA"
          
          ind=match(chosenAlleles[i],freqData$X)
          
          if (!is.na(ind)){
            
            masterFile[i,8]=freqData[ind,4]
            masterFile[i,9]=freqData[ind,7]
            masterFile[i,10]=freqData[ind,2]
            masterFile[i,11]=freqData[ind,5]
            masterFile[i,12]=freqData[ind,3]
            masterFile[i,13]=freqData[ind,6]
            
          } else {
            
            masterFile[i,c(8:13)]=c("NA","NA","NA","NA","NA","NA") # changed 11-16 to 8-13
            
          } #end freq if
          
          ind="NA"
          
          thisLine=masterFile[i,]
          names(masterFile)=c("OR","FDR LCI", "FDR UCI", "FDR p","Uncorrected p","Group","Loading","CaseRank","ControlRank","CaseCount","ControlCount","CaseFreq","ControlFreq")
          
          if (!exists("fullMaster")){
            fullMaster=thisLine
            names(masterFile)=c("OR","FDR LCI", "FDR UCI", "FDR p","Uncorrected p","Group","Loading","CaseRank","ControlRank","CaseCount","ControlCount","CaseFreq","ControlFreq")
          } else {
            fullMaster=rbind(fullMaster,thisLine)
          }
          
        }# end allele loop
        
        names(masterFile)=c("OR","FDR LCI", "FDR UCI", "FDR p","Uncorrected p","Group","Loading","CaseRank","ControlRank","CaseCount","ControlCount","CaseFreq","ControlFreq")
        rownames(masterFile)=chosenAlleles
        
      } #end allele check
      
      if (!exists("fullAlleles")){
        fullAlleles=chosenAlleles
      } else {
        fullAlleles=c(fullAlleles, chosenAlleles)
      }
      print(paste(popName,thisAnalysis,"completed."))
    } else {
      print(paste(popName,thisAnalysis,"SKIPPED."))
    }
    
  }###end analysis
  
  if (exists("fullMaster")){
    
    names(masterFile)=c("OR","FDR LCI", "FDR UCI", "FDR p","Uncorrected p","Group","Loading","CaseRank","ControlRank","CaseCount","ControlCount","CaseFreq","ControlFreq")
    rownames(fullMaster)=fullAlleles

    write.csv(fullMaster,file=paste0(getwd(),"/output/",thisDisease,"/Masters/",popName,"_Full_ALLELE_FA_MasterFile.csv")) 

  }
}###end populations
}



      
    
    
  
  
