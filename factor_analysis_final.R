# DESCRIPTION:

# This script is the third step in the disease association pipeline:
# It performs factor analysis and produces output in the /Multi directory

if(!require("MASS"))
  install.packages("MASS")
if(!require("psych"))
  install.packages("psych")
if(!require("GPArotation"))
  install.packages("GPArotation")
if(!require("corrplot"))
  install.packages("corrplot")
if(!require("parallel"))
  install.packages("parallel")

# this separates the alleles & haplotypes from the genotypes
currResolution=c("2dig","2digGeno") # HR and HRGeno added later

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

# begin analysis
for (pop in 1:length(POPULATIONS)){
  popName<-POPULATIONS[pop]

  for (resInd in 1:length(currResolution)) {
    thisRes=currResolution[resInd]

    interesting={}
    f=1 # FA based on first imputation replicate
    
    resInds={} # columns for this resolution
    
    for (analysis in 1:length(ANALYSIS_NAMES)) {
      thisAnalysis <- ANALYSIS_NAMES[analysis]
      
      # only columns that match this resolution at a time
      columnInds <- which(columns$analysis==thisAnalysis) # since names are unique should return only one columnInd
      HR_Res = ""
      if (thisRes == "2dig") { # Alleles + Haplotypes
        HR_Res = "HR"
      }
      else if (thisRes == "2digGeno") { # Genotypes
        HR_Res = "HRGeno"        
      }
      if ((columns$type[columnInds[1]]!=thisRes) & (columns$type[columnInds[1]]!=HR_Res)) { next }
      
      resInds=c(resInds,columnInds)
      
      # load significant associations from FDR results files
      fileName=paste0(getwd(),"/output/",thisDisease,"/FDR/",popName,"_",thisAnalysis,"_FDR.csv")
      if (file.exists(fileName)){
        fullORstats<-read.csv(file=fileName,header = T)
        fullORstats=as.data.frame(fullORstats)
        rownames(fullORstats)=fullORstats[,2]
        fullORstats=fullORstats[,3:ncol(fullORstats)]
        these=rownames(fullORstats)[which(fullORstats$significant=="***")]
        interesting=c(interesting,these)
      }
      else {
        print (paste0("Missing file: ",fileName))
      }
    } # end for analysis
    
    if (length(interesting)>1){
      
      fileName=paste0("./HLA_input_",thisDisease,"_",popName,"_",as.character(f),".txt")
      data<-read.table(file=fileName, sep=",",  header = TRUE) 
      
      thisData<-cbind(data[,1:ncol(data)])
      
      patthisAlleleData=thisData[which(thisData$Disease_Ind==1),1:ncol(thisData)]

      conthisAlleleData=thisData[which(thisData$Disease_Ind==0),1:ncol(thisData)]
      
      if ((thisRes == "HR") || (thisRes == "2dig")) {
        thisAlleleData=thisData[,3:ncol(thisData)]
        all=c(1:ncol(thisAlleleData))
        
        hap1Inds <- resInds[seq(1, length(resInds), 2)] # gets every other element from resInds
        hap2Inds <- resInds[seq(2, length(resInds), 2)]
        hap1=hap1Inds
        hap2=hap2Inds
        
        constack1=conthisAlleleData[,hap1]
        constack2=conthisAlleleData[,hap2]
        names(constack2)=names(constack1)
        conthisAlleleDataStack=rbind(constack1,constack2)
        
        patstack1=patthisAlleleData[,hap1]
        patstack2=patthisAlleleData[,hap2]
        names(patstack2)=names(patstack1)

        patthisAlleleDataStack=rbind(patstack1,patstack2)

      } else {
        patthisAlleleDataStack=patthisAlleleData[,resInds]
        conthisAlleleDataStack=conthisAlleleData[,resInds]
      }
      
      patfactormatrix=as.data.frame(matrix(0,nrow=nrow(patthisAlleleDataStack),ncol=1))

      confactormatrix=as.data.frame(matrix(0,nrow=nrow(conthisAlleleDataStack),ncol=1))
      
      for (cols in 1:ncol(patthisAlleleDataStack)){
        
        conthisCol=conthisAlleleDataStack[,cols]
        patthisCol=patthisAlleleDataStack[,cols]
        theseUni=unique(conthisCol)
        
        for (saps in 1:length(theseUni)){
          
          thisAllele=theseUni[saps]
          
          if (names(patthisAlleleDataStack)[cols]=="KIR_B1_80"  & thisAllele=="KIR_Ligand_Bw6"){
            thisAllele="KIR_Ligand_Bw6-NA"
          }
          
          if (thisAllele %in% interesting){
            
            contempAlleles=as.data.frame(matrix(0,nrow=nrow(conthisAlleleDataStack),ncol=1))
            pattempAlleles=as.data.frame(matrix(0,nrow=nrow(patthisAlleleDataStack),ncol=1))
            
            names(contempAlleles)=thisAllele
            names(pattempAlleles)=thisAllele
            
            coninds=which(conthisCol==thisAllele)
            patinds=which(patthisCol==thisAllele)
            
            contempAlleles[coninds,]=1
            pattempAlleles[patinds,]=1
            
            confactormatrix=cbind(confactormatrix,contempAlleles)

            patfactormatrix=cbind(patfactormatrix,pattempAlleles)
          }##end if
        } #end saps
      } #end col

      patfullfactormatrix=patfactormatrix[,2:ncol(patfactormatrix)]
      confullfactormatrix=confactormatrix[,2:ncol(confactormatrix)]
  
      patsiginds=match(interesting,names(patfullfactormatrix))
      consiginds=match(interesting,names(confullfactormatrix))
      
      patsigfactormatrix=patfullfactormatrix[,patsiginds]
      consigfactormatrix=confullfactormatrix[,consiginds]
      
      patcorrMat=cor(patsigfactormatrix)
      concorrMat=cor(consigfactormatrix)
      
      if (thisRes == "HR"){
        consigcount1=consigfactormatrix[1:nrow(constack1),]
        consigcount2=consigfactormatrix[(nrow(constack1)+1):nrow(consigfactormatrix),]
        consigcountmatrix=consigcount1+consigcount2
        
        patsigcount1=patsigfactormatrix[1:nrow(patstack1),]
        patsigcount2=patsigfactormatrix[(nrow(patstack1)+1):nrow(patsigfactormatrix),]
        patsigcountmatrix=patsigcount1+patsigcount2
        
      } else { 
        consigcountmatrix=consigfactormatrix
        patsigcountmatrix=patsigfactormatrix
      }

      checkData=as.data.frame(consigcountmatrix)
      patCheckData=as.data.frame(patsigcountmatrix)
      
      checkMat=cor(consigcountmatrix)
      checkCorr=as.matrix(checkMat)
      patCheckCorr=as.matrix(cor(patCheckData))

      xNames=names(checkData)
      contin=TRUE
      
      while (contin){
        diag(patCheckCorr)=0
        absCheckCorr=abs(patCheckCorr)
        allInds=which(absCheckCorr>0.99,arr.ind=TRUE)
        
        if (nrow(allInds)>1){
          contin=TRUE
        }else {
          contin=FALSE
        }# end if 1
        
        if (nrow(allInds)>0){
          rowInd=allInds[1,1]
          colInd=allInds[1,2]
          a=checkData[1:500,rowInd]
          b=checkData[1:500,colInd]
          invertCheck=patCheckCorr[rowInd,colInd]
          xNames=names(patCheckData)
          if (invertCheck<0){
            newName=paste0(xNames[rowInd]," %% *C* ",xNames[colInd])
          } else {
            newName=paste0(xNames[rowInd]," %% ",xNames[colInd])
          }
          newCol=as.data.frame(checkData[,rowInd])
          newPatCol=as.data.frame(patCheckData[,rowInd])
          names(newCol)=newName
          names(newPatCol)=newName
          
          low=min(c(rowInd, colInd))
          hi=max(c(rowInd, colInd))
          checkData[,hi]=NULL
          checkData[,low]=NULL
          checkData=cbind(checkData,newCol)
          checkCorr=as.matrix(cor(checkData))
          
          patCheckData[,hi]=NULL
          patCheckData[,low]=NULL
          patCheckData=cbind(patCheckData,newPatCol)

          patCheckCorr=as.matrix(cor(patCheckData))
          
        }#end if2
        
      }#end while
      
      
      ##################
      ## factor cases ##
      ##################
      
      patCheckData = patCheckData[, colSums(patCheckData !=0) >0] # this is new: deletes columns with only 0s
      
      patCheckCorr=as.matrix(cor(patCheckData))
      checkCorr=as.matrix(cor(checkData))

      if (length(rownames(patCheckCorr))>2) {
        
        png(filename=paste0(getwd(),"/output/",thisDisease,"/Multi/",popName,"_",thisRes,"_PATscreeplot.png"),height=700, width=800,type="cairo")
        patscree=fa.parallel(patCheckData, fm="minres",fa="fa")
        dev.off()
        
        patnumfactors=patscree$nfact
        
      } else {
        
        patnumfactors=0
      }
      
      if (patnumfactors>1){
        
        patsigfactors=fa(patCheckData,nfactors=patnumfactors,fm="minres",covar=FALSE)
        
        patloadings=as.matrix(patsigfactors$loadings[1:ncol(patCheckData),1:patnumfactors])
        #patloadings=eig$vectors[,1:patnumfactors]
        
        groupNum=matrix(0,nrow=nrow(patloadings), ncol=1)
        loadings=matrix(0,nrow=nrow(patloadings), ncol=1)
        ##*** loadings ==1 model of independce could not be rejected
        
        for (loads in 1:nrow(patloadings)){
          num=which.max(abs(patloadings[loads,]))
          groupNum[loads,1]=num
          loadings[loads,1]=patloadings[loads,num]
          
        } # end loadings loop
        
        if ((thisRes == "HR") || (thisRes == "2dig")) {
          groupNum=as.matrix(paste("Group",groupNum))
        } else {
          groupNum=as.matrix(paste("Geno Group",groupNum))
        }
        
        groupTable=as.data.frame(cbind(rownames(patloadings),groupNum,loadings))
        names(groupTable)=c("Allele","Group Assignment","Loading")
        groupTable
        
        if ((thisRes == "HR") || (thisRes == "2dig")) {
          write.csv(groupTable,file=paste0(getwd(),"/output/",thisDisease,"/Multi/",popName,"_Grouping.csv"))
        } else {
          write.csv(groupTable,file=paste0(getwd(),"/output/",thisDisease,"/Multi/",popName,"_Geno_Grouping.csv"))
        }
        
      }## end if
    } ## end interesting check
  } # end for resolution
} # end for populations