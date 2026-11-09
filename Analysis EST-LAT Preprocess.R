#Set working directory to source file location

#Libraries
library(qs)         
library(data.table) 
library(readxl)     
library(janitor)    
library(tidyverse)  
library(SRS)        
library(Biostrings) 

# install.packages(c("qs2", "data.table", "readxl", "janitor", "tidyverse", "SRS"))
# BiocManager::install("Biostrings")

#Data
qs::qload("../data/EstLat.qs")
ft <- read_excel("../data/13225_2020_466_MOESM4_ESM.xlsx") #Fungal traits

#Metadata
setDT(metadata)                
metadata <- metadata[!(control %chin% "MOLD")]  #Remove MOLD controls

#Filter low-depth samples (sum of Abundance >= 1223)
TAB_flt <- TAB_smr[SampleID %in% metadata$SampleID, .SD[sum(Abundance) >= 1223], by = SampleID]


#Filter metadata for vegetation type
metaford <- metadata[
  `dataset#4:_ectomycorrhizal_fungi` == "yes" &
    EcM_plants_sampled %between% c(1, 9) &
    !is.na(`vegetation_type_(for_habitat_type_analysis)`) &
    SampleID %in% TAB_flt$SampleID,
  .(SampleID, BIOME, EcM_plants_sampled, EcM_OTUs, `vegetation_type_(for_habitat_type_analysis)`)
]

#Keep vegetation types with >=14 samples
metaford <- metaford[, if(.N >= 14) .SD, by = `vegetation_type_(for_habitat_type_analysis)`]
metaford <- metaford[order(BIOME, EcM_plants_sampled, -EcM_OTUs)]


#Subsampling per vegetation type
subveg <- metaford[, {
  n_categories <- uniqueN(EcM_plants_sampled)           #Number of categories
  n_per_category <- floor(25 / n_categories)            #Base samples per category
  remainder <- 25 %% n_categories                       #Extra samples to distribute
  
  #Sort by EcM_plants_sampled (ascending), EcM_OTUs (descending), SampleID
  sorted_SD <- .SD[order(EcM_plants_sampled, -EcM_OTUs, SampleID)]
  
  #Take base number per category
  selected <- sorted_SD[, head(.SD, min(.N, n_per_category)), by = EcM_plants_sampled]
  
  #Distribute remainder to highest OTU samples
  if (remainder > 0) {
    extra_samples <- sorted_SD[!SampleID %in% selected$SampleID][order(-EcM_OTUs)][1:remainder]
    selected <- rbind(selected, extra_samples)
  }
  
  #Limit to max 25 samples
  if (nrow(selected) > 25) {
    selected <- selected[order(-EcM_OTUs, SampleID)][1:25]
  }
  
  #Complete if fewer than 25
  missing <- max(0, 25 - nrow(selected))
  if (missing > 0) {
    extra_samples <- sorted_SD[!SampleID %in% selected$SampleID][order(-EcM_OTUs, SampleID)][1:missing]
    selected <- rbind(selected, extra_samples)
  }
  
  selected[!is.na(SampleID)]
}, by = `vegetation_type_(for_habitat_type_analysis)`]

#Rename column
setnames(subveg, "vegetation_type_(for_habitat_type_analysis)", "vegetation_type")

#Filtered and full metadata for vegetation
metaveg_r <- subveg
metaveg_c <- metadata[SampleID %in% subveg$SampleID]
setnames(metaveg_c, "vegetation_type_(for_habitat_type_analysis)", "vegetation_type")


#OTU table for vegetation
otu <- dcast(TAB_flt, ID_Clustered ~ SampleID, value.var = "Abundance", fun.aggregate = sum, fill = 0)
otu <- as.data.frame(otu)
otu <- column_to_rownames(otu, "ID_Clustered")
otu_ids <- rownames(otu)

#Subset OTUs for selected vegetation samples
otuveg <- otu[, intersect(colnames(otu), subveg$SampleID), drop = FALSE]


#Normalize OTU table with SRS
#Calculate minimum sequencing depth per sample
vegsam <- merge(
  TAB_smr[, .(abundance = sum(Abundance)), by = SampleID], 
  subveg[, .SD, .SDcols = c("SampleID", "vegetation_type")], 
  by = "SampleID"
)
Cmin <- min(vegsam$abundance)

#Apply SRS normalization
otuvnorm <- SRS(otuveg, Cmin, set_seed = TRUE, seed = 1)

#Restore OTU row names and remove OTUs absent in all samples
rownames(otuvnorm) <- otu_ids
otuvnorm <- rownames_to_column(otuvnorm, var = "ID_Clustered")
otuvnorm <- otuvnorm[rowSums(otuvnorm[, -1] == 0, na.rm = TRUE) != (ncol(otuvnorm) - 1), ]
setDT(otuvnorm)


#Prepare taxonomy
TAXbiom <- as.data.frame(TAXf)
TAXbiom2 <- TAXbiom %>% 
  pivot_longer(cols = -c(QueryName, Seq), names_to = c("Number", ".value"), 
               names_pattern = "(\\d+)_(.*)")
TAXbiom2$Genus <- gsub("\\(Fungi\\)|\\(|\\)|\\.s\\.str\\.?$", "", as.character(TAXbiom2$Genus))
TAXbiom3 <- subset(TAXbiom2, Number == 1)


#Assign guilds
ft[ft == ""] <- NA

assign_guilds <- function(tax) {
  genus_index <- which(names(ft) == "GENUS")
  columns_after_genus <- ft[, (genus_index + 2):(ncol(ft) - 1)]
  
  guild_a <- merge(tax, 
                   ft[, c("GENUS", names(columns_after_genus))], 
                   by.x = "Genus", 
                   by.y = "GENUS", 
                   all.x = TRUE)
  return(guild_a)
}

guild <- assign_guilds(TAXbiom3)
guildsub <- guild[primary_lifestyle == "ectomycorrhizal"]


#Filter taxonomy for vegetation
taxveg <- guildsub[QueryName %chin% otuvnorm$ID_Clustered]
setcolorder(taxveg, c(setdiff(names(taxveg), "Genus"), "Genus"))
setcolorder(taxveg, c(
  names(taxveg)[1:which(names(taxveg) == "Family")],
  "Genus",
  setdiff(names(taxveg), c(names(taxveg)[1:which(names(taxveg) == "Family")], "Genus"))
))

#Final filtered OTU table
otuvst <- otuvnorm[ID_Clustered %in% taxveg$QueryName]

#Subset sequences for vegetation
smrveg <- TAB_smr[ID_Clustered %in% taxveg$QueryName & SampleID %in% subveg$SampleID]
seqveg <- SEQS[names(SEQS) %in% smrveg$ID_Denoised]

#Write sequences to fasta
# writeXStringSet(seqveg, filepath = "../data/seqveg.fasta")

#Save all objects for vegetation subsample

qs::qsave(list(
  otuvst = otuvst,
  metaveg_c = metaveg_c,
  metaveg_r = metaveg_r,
  taxveg = taxveg,
  seqveg = seqveg
), "vegetation.qs")
