# this script organizes all of the master files into an Excel spreadsheet

# install the required packages and open the required libraries
if(!require("openxlsx"))
  install.packages("openxlsx")
if(!require("gtools"))
  install.packages("gtools")

library(openxlsx)
library(gtools)

# define command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 0) {
  thisDisease <- args[1]
  pop <- args[2]
}

# set up print statements for the arguments:
print(paste0("This disease: ",thisDisease))
POPULATIONS = as.list(strsplit(pop, ",")[[1]]) # formats the population argument
print(paste0("Population(s): ",POPULATIONS))

# identify the name of the final tables file
final_tables <- paste0(getwd(),"/output/",thisDisease,"/Masters/","Final_Tables.xlsx")

# create an empty dataframe that has the column names of the final tables file
empty_df <- setNames(data.frame(matrix(ncol=13, nrow=0)),c("X","OR","FDR LCI","FDR_p","Uncorrected p","Group","Loading","CaseRank","ControlRank","CaseCount","ControlCount","CaseFreq","ControlFreq"))
write.xlsx(empty_df, final_tables) # write the empty dataframe to an Excel spreadsheet
empty_wb <- loadWorkbook(final_tables) # load the workbook

# loop through all of the populations 
for (pop in POPULATIONS){
    master_file <- paste0(getwd(),"/output/",thisDisease,"/Masters/",pop,"_Full_ALLELE_FA_MasterFile.csv") # identify the name of the master file
    master_file <- read.csv(master_file) # read the master file
    names(master_file)[names(master_file)== "X"] <- "HLA.variant" # rename the column that contains the HLA variants
    master_file$value <- gsub("[^0-9.-]", "", master_file$Group) # make a "value" column with the "Group" number
    master_file$remainder <- as.numeric(as.character(master_file$value))%%2 # assign labels of 1 or 0 for odd or even numbers
    master_file$color <- ifelse((master_file["remainder"] != 0), "White","Gray") # assign color labels based on odd or even numbers
    
    master_file <- master_file[order(master_file$FDR.p),] # sorts the dataframe first based on p value
    master_file <- master_file[mixedorder(master_file$Group),] # sorts the dataframe second based on group

    addWorksheet(empty_wb, pop) # add a worksheet for each population
    header_style <- createStyle(textDecoration = "Bold") # format the header on each worksheet to bold
    writeData(empty_wb,pop,master_file, headerStyle=header_style) # write each worksheet and apply the formatting
    
    columns = 1:which(colnames(master_file)=="ControlFreq") # gives you a column range that excludes the extra columns
    rows = 1:(nrow(master_file)+1) # gives you the range of rows in the master_file

    color <- createStyle(fgFill="#CCCCCC") # gray color
    addStyle(wb = empty_wb, sheet=pop, style=color, rows=which(grepl(master_file$color, pattern="White"))+1, col=columns, gridExpand=TRUE) # add the style
    deleteData(empty_wb, sheet=pop, rows=rows, cols=15:17, gridExpand=TRUE) # delete the columns you don't need: value, remainder, and color

    saveWorkbook(empty_wb,final_tables,overwrite=TRUE) # save the workbook

}

# remove the original empty sheet from the workbook
filled_wb <- loadWorkbook(final_tables)
removeWorksheet(filled_wb, sheet="Sheet 1")
saveWorkbook(filled_wb, final_tables, overwrite=TRUE)
