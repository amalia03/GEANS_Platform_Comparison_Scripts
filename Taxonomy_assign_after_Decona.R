library(taxonomizr)
library(dplyr)
library(worrms)
library(stringr)
library(data.table)
library(tidyr)

# Load the final output file from Decona. This file contains all consensus sequences and the blast results.
dataset <- read.csv("OTUcalled_finalblastoutput.txt", header=FALSE, sep="\t")

#Get the taxonomy of the accession numbers from the accessionTaxa.sql file
#res$accessionNumber <- sapply(strsplit(as.character(res[,16]),'\\|'),'[',4)
dataset$accessionNumber <- gsub(".*\\|([^|]+)\\|", "\\1", dataset[,16])
dataset$accessionNumber <- gsub("_\\(2\\)", "", dataset$accessionNumber)
dataset$accessionNumber <- ifelse(grepl("\\.[0-9]+$", dataset[, 18]), dataset[, 18], paste0(dataset[, 18], ".1"))
dataset$accession <- accessionToTaxa(dataset$accessionNumber, "/home/reindert/Blast_database/Accession_db_13122023/accessionTaxa.sql" , version = c("version", "base"))
dataset$Taxonomy <-getTaxonomy(dataset$accession , "/home/reindert/Blast_database/Accession_db_13122023/accessionTaxa.sql")

### Filter data and extract read count 
#Set column headers
names(dataset) <- c("Barcode", "Total_readcount", "SeqID", "PercIdent", "Alignmentlength", "Mismatch", "Gapopen", "QueryStart", "QueryEnd", "SubjectStart", "SubjectEnd", "Evalue", "Bitscore", "Querylength", "Subjectlength", "SseqID", "OTU", "accessionNumber", "accession", "Taxonomy") 
# Make sure sequence ID's are unique
dataset <- unite(dataset, uniqueID, c(Barcode, SeqID), remove=FALSE)
# Sort data based on Evalue
dataset <- dataset[order(dataset$Evalue),] 
# keep the top 5 hits per sequence
dataset <- Reduce(rbind, by(dataset, dataset["uniqueID"], head, n=5))
# Extract read count sequences from sequence name
#dataset$Clust_nr <- str_extract(dataset$SeqID, "\\d+-\\d+$")
#dataset$Read_Count <- as.numeric(str_extract(dataset$Clust_nr, "\\d+"))\
dataset <- cbind(dataset, reshape2::colsplit(dataset$SeqID, pattern="-", c("a", "b", "c", "d")))
dataset$Read_Count<- ifelse(grepl("*subsampled*", dataset$b), dataset$c, dataset$b)

#sort by SequenceID
dataset <- dataset[order(dataset$uniqueID),]
#put taxononomy in new columns
dataset$phylum <- dataset$Taxonomy[,"phylum"]
dataset$class <- dataset$Taxonomy[,"class"]
dataset$order <- dataset$Taxonomy[,"order"]
dataset$family <- dataset$Taxonomy[,"family"]
dataset$genus <- dataset$Taxonomy[,"genus"]
dataset$species <- dataset$Taxonomy[,"species"]


### All that have no hit with the database is names as unclassified
#first only percentage ID
#Change dataset$PercIdent that contain NA in the column to 0
dataset <- dataset %>% mutate(PercIdent = ifelse(is.na(PercIdent), 0, PercIdent))

#but than also taxonomy needs to be renamed to unclassified:
dataset$phylum[dataset$PercIdent==0] <- "unclassified"
dataset$class[dataset$PercIdent==0] <- "unclassified" 
dataset$order[dataset$PercIdent==0] <- "unclassified"
dataset$family[dataset$PercIdent==0] <- "unclassified"
dataset$genus[dataset$PercIdent==0] <- "unclassified"
dataset$species[dataset$PercIdent==0] <- "unclassified"

################################################################################################################################################################################################
#If multiple options for a single sequence are still present, and one or more of those options have 100% Percentage identity, 
#then drop other options that have less than 100% identity

#first remove the duplicates so only hits that are truely different are kept
dataset <- dataset %>% 
  distinct(uniqueID, Evalue, PercIdent, Alignmentlength, species, .keep_all=TRUE)

# Then based on all selected columns rename species that occur more then once for "unclassified"
#and put the original species name in a new column
#When only one species was present, then it shows "no problem" in the warning column

dataset <- dataset %>% 
group_by(uniqueID, Evalue, PercIdent, Alignmentlength) %>%
  mutate(warning = paste( species, collapse = " + ")) %>%
  mutate(phylum = ifelse(n_distinct(phylum) > 1, "unclassified", phylum)) %>%
  mutate(class = ifelse(n_distinct(class) > 1, "unclassified", class)) %>%
  mutate(order = ifelse(n_distinct(order) > 1, "unclassified", order)) %>%
  mutate(family = ifelse(n_distinct(family) > 1, "unclassified", family)) %>%
  mutate(genus = ifelse(n_distinct(genus) > 1, "unclassified", genus)) %>%
  mutate(species = ifelse(n_distinct(species) > 1, "unclassified", species)) %>%
  mutate(warning = ifelse(grepl("\\+", warning), warning, "No problem")) %>%
  ungroup()

# warning, but this just puts the whole "origninal_values list" so not perfect but ok
if (any(!is.na(dataset$warning))) {
  warning(paste("some taxa where not identified to species level, check column 'warning' in the R object "))
}

# Keep the highest scoring hit
dataset <- Reduce(rbind, by(dataset, dataset["uniqueID"], head, n=1))

## Add consensus sequences to the dataset

Sequences <- lapply(paste0(list.files(pattern=".*_concatenated.fasta")), function(filename){
  barcode = str_extract(filename, "(barcode)[0-9]+")
  c(paste0(barcode), read.delim(filename, sep= ">", header=FALSE))
})

Sequences <- rbindlist(lapply(Sequences, as.data.table), use.names=TRUE, fill=TRUE, idcol=TRUE)[,-1]
#Change column names to something more appropriate
colnames(Sequences) <- c("Barcode", "Sequence", "SeqID")
#Tidy up data
Sequences <- Sequences %>% mutate(SeqID = lag(SeqID)) ### Wat doet deze? En kan hij weg?
Sequences <- Sequences[!(Sequences$Sequence == ""),] ### Wat doet deze? En kan hij weg?
Sequences$SeqID <- gsub(".fasta", "", Sequences$SeqID)

#Make sure seqID's are unique
Sequences <- unite(Sequences, uniqueID, c(Barcode, SeqID), remove=FALSE)

#Merge two datasets (dataset and Sequences)
Sequencedataset <- merge(dataset, Sequences, by="uniqueID")
#Remove redundant columns and rename some
Sequencedataset <- Sequencedataset %>% select(-Barcode.y, -a, -b, -c, -d)
Sequencedataset <- setnames(Sequencedataset, old = c("Barcode.x"), new=c("Barcode"))


################################################################################################################################################################################################
#Save and extract the final dataset
################################################################################################################################################################################################

save(Sequencedataset, file="Decona_R_Object.RData")

getwd()
