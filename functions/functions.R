#functions.R

#reads csv file into dataframe
morph1.morph2 <- read.csv("data/Transcription_Morph-Morph.csv")

# if the publication is from Mack, replace publication with, "4"
morph1.morph2$Publication[morph1.morph2$Publication == "Mack_et_al_2020"] <- rep("4", sum(morph1.morph2$Publication == "Mack_et_al_2020"))
# if the publication is from McGaugh, replace with, "5"
morph1.morph2$Publication[morph1.morph2$Publicaton == "McGaugh_et_al_2020"] <- rep("5", sum(morph1.morph2$Publication == "McGaugh_et_al_2020"))

#reads csv file into GeneToGO dataframe, blank fields are added where rows are unequal length
GeneToGO <- read.csv("data/AMexGOTerms.csv", fill = T)

TranscTable <- function(morph1, morph2, condition, direction, tr.stat, tr.thresh, GOTable){
  wrnings <- c("Errors: ")
  # If condition is NOT "Between morph"...
  if(condition != "Between morph"){
    # Use transcription data of morph-control comparisons
    # empty dataframe for control condition
    condition_control <- data.frame()
    in_table <- condition_control
    # Set morph2 to control
    morph2 <- "Control"
    # Set grep pattern to morph1
    search_cond <- morph1
    # If searching between morphs, set search condition to morphs-of-comparison
  }else if(condition == "Between morph"){
    # Use transcription data for between-morph comparisons
    in_table <- morph1.morph2
    # Set grep pattern to comparison
    search_cond <- paste(c(morph1, morph2), collapse = "-")
  }
  # Store comparison for output. Note: "morph2" is control if "Between morph"
  # was not entered
  comp <- paste(c(morph1, morph2), collapse = "-")
  
  # If genes above threshold were requested...
  if(direction == "Above"){
    # Find all rows-of-interest (ROIs) for morph(s)-of-interest where logFC is
    # above specified value, and condition matches the input specification
    if(tr.stat == "logFC"){
      ROIs <- in_table[(grepl(search_cond, in_table$Comparison)
                        & (in_table$logFC > tr.thresh) &
                          (in_table$Condition == condition)), ]
      # Sort candidate rows with highest logFC values on top
      ROIs <- ROIs[order(ROIs$logFC, decreasing = T),]
      # Find all rows-of-interest (ROIs) for morph(s)-of-interest where logFC is
      # below specified value, and condition matches the input specification
    }else{
      ROIs <- in_table[(grepl(search_cond, in_table$Comparison)
                        & (in_table$PValue > tr.thresh) &
                          (in_table$Condition == condition)), ]
      # Sort candidate rows with highest p values on top
      ROIs <- ROIs[order(ROIs$PValue, decreasing = T),]
    }
    
    # If genes whose value is below the provided stat were requested...
  }else if(direction == "Below"){
    # Find all rows-of-interest (ROIs) for morph(s)-of-interest where logFC is
    # BELOW specified value, and condition matches the input specification
    if(tr.stat == "logFC"){
      ROIs <- in_table[(grepl(search_cond, in_table$Comparison)
                        & (in_table$logFC < tr.thresh) &
                          (in_table$Condition == condition)), ]
      # Sort candidate rows with SMALLEST logFC values on top
      ROIs <- ROIs[order(ROIs$logFC, decreasing = F),]
    }else{
      ROIs <- in_table[(grepl(search_cond, in_table$Comparison)
                        & (in_table$PValue < tr.thresh) &
                          (in_table$Condition == condition)), ]
      ROIs <- ROIs[order(ROIs$PValue, decreasing = F),]
    }
  }
  # Check if any genes were found for specified conditions
  if(nrow(ROIs) == 0){
    wrnings <- append(wrnings, "No genes found matching given parameters.")
    output.df <- as.data.frame(matrix(rep(NA,10), ncol = 10))
    names(output.df) <- c("Gene Name","Gene Stable ID","GO Term(s)","Comparison",
                          "logFC","p-value","Age at Sampling","Tissue","Ensembl Family Description",
                          "Publication"
    )
    return(list(output.df, wrnings))
  }
  # Obtain GO terms for ROIs
  GOTerms <- character(length = nrow(ROIs))
  for(i in 1:length(GOTerms)){
    if(length(grep(ROIs$Gene_stable_ID[i], GOTable$Ensembl)) != 0){
      GOTerms[i] = paste(GOTable$Gene.ontology.IDs[GOTable$Ensembl == ROIs$Gene_stable_ID[i]], collapse = " ")
    }else{
      GOTerms[i] = NA
    }
  }
  # Output gene names, gene stable IDs, GO terms, morph-of-comparison, logFC,
  # p-value, Ensembl family information (or study specific details, if Between
  # morph), and publication name to a dataframe
  if(condition == "Between morph"){
    special_info <- ROIs$study_specific_gene_details
    special_name <- "Study-Specific Details"
  }else{
    special_info <- ROIs$Ensembl_Family_Description
    special_name <- "Ensembl Family Description"
  }
  output.df <- data.frame(
    tolower(ROIs$Gene_name),
    ROIs$Gene_stable_ID,
    GOTerms,
    rep(comp, nrow(ROIs)),
    ROIs$logFC,
    ROIs$PValue,
    ROIs$Age_at_Sampling,
    ROIs$Tissue,
    special_info,
    ROIs$Publication
  )
  names(output.df) <- c(
    "Gene Name",
    "Gene Stable ID",
    "GO Term(s)",
    "Comparison",
    "logFC",
    "p-value",
    "Age at Sampling",
    "Tissue",
    special_name,
    "Publication"
  )
  return(list(output.df,wrnings))
}
