#' Check that sumstats has correct columns and that they are in the correct order for MAGMA and LDSC, pure R version (can be very slow)
#'
#' @return col_headers The new column headers for the sumstats file
#'
#' @examples
#' format_sumstats_for_magma_crossplatform(path)
#'
#' @import data.table
#' @import stringr
#' @export



format_sumstats_for_magma_crossplatform <- function(path) {
  # Checking if the file exists should happen first
  if (!file.exists(path)) {stop("Path to GWAS sumstats is not valid")}
  
  # Make another file in order to not overwrite the original.
  processed_file = paste0(path, "_processed.sumstats")
  if (file.exists(processed_file)) {
    print(paste0(processed_file, " already exists."))
    if (toupper(readline("Do you want to overwrite? [ y/n ]: ")) %in% c("YES", "Y")) {
      file.remove(processed_file)
      print("Copying into intermediary file (in order to not modify the original).")
      file.copy(path, processed_file, overwrite = FALSE)
      path <- processed_file
    } else {
      stop("User does not want overwriting of currently existing processed sumstats file.")
    }
  }
  
  # Read the sumstats file
  print("Loading sumstats file.")
  sumstats_file <- readLines(path)
  
  
  print("Ensuring that FS are tabs, not spaces.")
  row_of_data <- strsplit(sumstats_file[2], "\t")
  if (length(row_of_data) == 1) {
    if (grep(" ", row_of_data) == 1) {
      print("This GWAS sumstat file has 'space' FS instead of tabs (unusual, not proper input for MAGMA). Correcting in output file.")
      sumstats_file <- gsub(pattern = " ", replace = "\t", x = sumstats_file)
    }
  }
  
  # This standardises the sumstats column headers
  sumstats_file[1] = standardise.sumstats.column.headers.crossplatform(sumstats_file[1])
  
  # Check if there are CHR and BP columns
  if(!sum(c("SNP","BP") %in% sumstats_file[1])==2){
    # If not, see if there is a column storing the data in a format like: CHR:BP:A2:A1 or just CHR:BP
    # - UKBB data from Ben Neale has a Variant column with CHR:POS:REF:ALT where ALT allele is the effect allele in the model [NB: the ALT allele is NOT always the minor allele]
    # -- For input to LDSC, A1 is effect allele, A2 is non-effect allele
    # - DIAGRAM diabetes data has a Chr:Position column with CHR:BP
    # - BMI adjusted for smoking has markername with CHR:BP (with the chromosome name having 'chr' preceeding)
    # - Agression [EAGLE] just doesn't have any CHR or BP data
    print("Summary statistics file does not have obvious CHR or BP columns. Checking to see if they are joined in another column")
    
    # Obtain a row of the actual data
    row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
    
    # Check if there is a column of data with CHR:BP:A2:A1 format
    fourStepCol = grep(".*:.*:\\w:\\w",row_of_data) # This grabs the column that contains {numeric (.) + anything (*) + ":" + numeric (.) + anything (*) + ":" + character + ":" + character}
    
    if (length(fourStepCol)) { 
      
      header <- strsplit(sumstats_file[1], "\t")[[1]]
      
      print(sprintf("Column %s is being replaced with CHR BP A2 A1", header[fourStepCol]))
      
      header[fourStepCol] <- "CHR\tBP\tA2\tA1"
      
      sumstats_file[1] <- paste(header, collapse = "\t")
      
      # Now that the header is fixed, replace the first three ":" with tabs
      sumstats_file <- sub(pattern = ":", replace = "\t", x = sub(pattern = ":", replace = "\t", x = sub(pattern = ":", replace = "\t", x = sumstats_file)))
      
      #print(sumstats_file[1])
      row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
    }
    
    # Check if there is a column of data with CHR:BP format
    twoStepCol = grep(".*:.*", row_of_data)
    if (length(twoStepCol)) {
      
      header <- strsplit(sumstats_file[1], "\t")[[1]]
      
      print(sprintf("Column %s is being replaced with CHR BP", header[twoStepCol]))
      
      header[twoStepCol] <- "CHR\tBP"
      
      sumstats_file[1] <- paste(header, collapse = "\t")
      
      # Now that the header is fixed, replace the first ":" with tabs
      sumstats_file <- sub(pattern = ":", replace = "\t", x = sumstats_file)
      
      #print(sumstats_file[1])
      row_of_data <- strsplit(sumstats_file[2], "\t")[[1]]
    }
    
    # Restandardise in case the joined column headers were unusual
    sumstats_file[1] = standardise.sumstats.column.headers.crossplatform(sumstats_file[1])
  }
  
  
  # OPTIONAL CODE, WHICH WILL ELIMINATE THE ROWS WHICH HAVE MORE COLUMNS THAN THERE ARE COLUMN HEADERS.
  
  print(paste0("Eliminating any rows (SNPs) that have more columns than the header. Current number of rows: ", length(sumstats_file)))
  expected_num_of_columns <- length(strsplit(sumstats_file[1], "\t")[[1]])
  rows_to_eliminate <- c()
  how_many_rows_to_eliminate <- 0
  tab_freq <- function(x){
   count <- table(strsplit(x, '')[[1]])
   return(count[names(count) == "\t"])
  }
  for (row in 2:length(sumstats_file)) {
   if (row %% 100000 == 0) {print(paste0(row, " checked rows."))}
   if (tab_freq(sumstats_file[row]) != (expected_num_of_columns-1)) {
     rows_to_eliminate[how_many_rows_to_eliminate+1] <- row
     how_many_rows_to_eliminate <- how_many_rows_to_eliminate + 1
   }
  }
  if (how_many_rows_to_eliminate > 0) {sumstats_file <- sumstats_file[-rows_to_eliminate]}
  print(paste0("Rows left after processing: ", length(sumstats_file)))
  
  
  # If SNP is present... BUT not CHR or BP then need to find the relevant locations
  col_headers = strsplit(sumstats_file[1], "\t")[[1]]; writeLines(sumstats_file, con = path)
  if(sum(c("CHR","BP") %in% col_headers)==0 & sum("SNP" %in% col_headers)==1){
    print("SNP is present, but not CHR or BP. Finding the relevant genome locations.")
    library(data.table)
    sumstats = fread(path)
    data("SNP_LOC_DATA")
    SNP_LOC_DATA_2 = SNP_LOC_DATA[,1:3]
    sumstats2 = merge(sumstats,SNP_LOC_DATA_2,by="SNP")
    sumstats3 = data.frame(sumstats2)[,c("SNP","CHR","BP",setdiff(colnames(sumstats2),c("SNP","CHR","BP")))]
    fwrite(sumstats3,file=path,sep="\t"); sumstats_file <- readLines(path)
  }
  
  # If CHR and BP are present... BUT not SNP then need to find the relevant SNP ids
  col_headers = strsplit(sumstats_file[1], "\t")[[1]]
  if(sum(c("CHR","BP") %in% col_headers)==2 & sum("SNP" %in% col_headers)==0){
    print("There is no SNP column found within the data. It must be inferred from CHR and BP information.")
    print("Note: this process drops any SNPs which are not from Hapmap.")
    genomebuild <- as.numeric(readline("Which genome build is the data from? 1 for GRCh37, 2 for GRCh38"))
    if(!genomebuild %in% c(1,2)){stop("Genome build must be entered as either 1 (for GRCh37) or 2 (for GRCh38)")}
    data("SNP_LOC_DATA")
    if(genomebuild==1){genomebuild="GRCh37"}else{genomebuild="GRCh38"}
    snpLocDat = SNP_LOC_DATA[SNP_LOC_DATA$Build==genomebuild,][,-4]
    library(data.table)
    sumstats = fread(path)
    sumstats$CHR = as.factor(sumstats$CHR)
    if(length(grep("chr",sumstats$CHR[1]))!=0){sumstats$CHR = gsub("chr","",sumstats$CHR)}
    sumstats2 = merge(sumstats,snpLocDat,by=c("CHR","BP"))
    fwrite(sumstats2,file=path,sep="\t"); sumstats_file <- readLines(path)
  }
  
  
  print("Checking that all the vital columns are present (SNP, CHR, BP, P, A1, A2).")
  rows_of_data <- c(sumstats_file[1], sumstats_file[2]); col_headers = strsplit(rows_of_data[1], "\t")[[1]]
  for(key_column in c("SNP","CHR","BP","P","A1","A2")){
    code_example = "sed -i '' '1s/p_value/P/' IQ.Sniekers.2017.txt"
    if(!key_column %in% col_headers){
      print("Header of file:")
      print(rows_of_data)
      stop(sprintf("Cannot find a %s column in GWAS sumstats file. \nUse code such as '%s' to fix",key_column,code_example))
    }
  }
  
  
  # Check there is at least one signed sumstats column
  print("Checking that there is at least one signed sumstats column (eg.: Z, OR, BETA, LOG_ODDS, SIGNED_SUMSTAT)")
  signed_stat_column_names = c("Z","OR","BETA","LOG_ODDS","SIGNED_SUMSTAT")
  if(sum(signed_stat_column_names %in% col_headers)<1 %in% col_headers){
    print("Header of file:")
    print(rows_of_data)
    stop("ERROR: cannot find a column name representing signed statistic in GWAS sumstats file. I.e. Z, OR, BETA")
  }
  print("All vital columns present.")
  
  
  # Check that first three column headers are SNP, CHR, BP (in that order)
  print("Checking that the first three column headers are SNP, CHR and BP in this order.")
  if(!sum(col_headers[1:3]==c("SNP","CHR","BP"))==3){
    print("Reordering the columns (to be SNP, CHR and BP et al.). This may take 30 seconds - 15 mins depending on hardware and sumstats file size.")
    whichSNP = which(col_headers=="SNP")[1]
    whichCHR = which(col_headers=="CHR")[1]
    whichBP = which(col_headers=="BP")[1]
    otherCols = setdiff(1:length(col_headers),c(whichSNP,whichCHR,whichBP))
    
    write.table(x=read.table(path)[,c(whichSNP,whichCHR,whichBP,otherCols)], file=path, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  
  # MAGMA cannot handle P-values as low as 3e-400... so convert them to zeros
  print("The following operation may take 30 seconds - 15 minutes depending on hardware and sumstats file size.")
  print("MAGMA cannot handle P-values as low as 3e-400.")
  if (toupper(readline("Do you want MAGMA.celltyping to convert any (if) existing ones to zeroes? [ y/n ]: ")) %in% c("YES", "Y")) {
    rows_of_data <- c(sumstats_file[1], sumstats_file[2]); col_headers = strsplit(rows_of_data[1], "\t")[[1]]
    sumstats <- read.table(path, stringsAsFactors = FALSE)
    for (i in seq_along(sumstats[,which(col_headers=="P")])) {
      if (sumstats[i,which(col_headers=="P")]=="P") {next} # To skip the header.
      sumstats[i,which(col_headers=="P")] <- as.numeric(as.character(sumstats[i,which(col_headers=="P")])) # This converts anything under 3e-400 to zeros.
    }
    write.table(x=sumstats, file=path, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE); sumstats_file <- readLines(path)
  }
  
  
  # Sometimes the N column is not all integers... so round it up
  if("N" %in% col_headers) {
    #whichN = which(col_headers %in% "N")
    print("The following operation may take 30 seconds - 15 minutes depending on hardware and sumstats file size.")
    print("Sometimes the N column is not all integers.")
    if (toupper(readline("Do you want MAGMA.celltyping to round them up if existing? [ y/n ]: ")) %in% c("YES", "Y")) {
      rows_of_data <- c(sumstats_file[1], sumstats_file[2]); col_headers = strsplit(rows_of_data[1], "\t")[[1]]
      sumstats <- read.table(path, stringsAsFactors = FALSE)
      for (i in seq_along(sumstats[,which(col_headers=="N")])) {
        if (sumstats[i,which(col_headers=="N")]=="N") {next} # To skip the header.
        sumstats[i,which(col_headers=="N")] <- round(as.numeric(as.character(sumstats[i,which(col_headers=="N")]))) # This converts anything under 3e-400 to zeros.
      }
      write.table(x=sumstats, file=path, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE); sumstats_file <- readLines(path)
    }
  }
  
  # All rows should start with either SNP or rs... if they don't drop them
  print("Dropping all rows that don't start with 'rs'")
  sumstats_file <- readLines(path)
  sumstats_file <- c(sumstats_file[1],sumstats_file[grepl("^rs",sumstats_file)])
  writeLines(text=sumstats_file, con = path)
  
  print("Removing duplicated RSIDs.")
  sumstats <- read.table(path)
  if(sum(duplicated(sumstats[,1]))>0){
    notDup = which(!duplicated(sumstats[,1]))
    write.table(x=read.table(path)[notDup,], file=path, sep="\t", quote=FALSE, row.names = FALSE, col.names = FALSE); sumstats_file <- readLines(path)
  }
  
  # Show how the data now looks
  print("Succesfully finished preparing sumstats file:")
  print("Header of file:")
  con <- file(path,"r") ; rows_of_data <- readLines(con,n=2) ; close(con)
  print(rows_of_data)
  col_headers = strsplit(rows_of_data[1],"\t")[[1]]
  print(col_headers)
  return(path) # Returns address of modified file
}
