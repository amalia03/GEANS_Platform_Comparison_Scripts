{library(tidyverse)
library(phyloseq)}


my_data_decona  <- Sequencedataset
#remove"barcode and keep only the number
my_data_decon2$Barcode <- gsub('barcode', '', my_data_decona$Barcode )
head(my_data_decon2)

#add a list with general metadat!!! FOR PHYLOSEQ MAKE SURE THAT THE NAME OF THIS COLUMN IS THE SAME AS IN YOUR METADATA SHEET SO IT CAN COMMUNICATE
#PHYLOSEQ IS THE WHOLE REASON WHY YOU PUT THE LIST HERE
list_rep <- as.data.frame(R_list_2023_ERIK_blast)
names(list_rep)[1] <- c("Barcode")
correction_data <- merge(my_data_decona,list_rep,by="Barcode")



#           Taxonomy assignment correction          

correction_data$phylum[correction_data$Alignmentlength <= 250] <- ""
correction_data$species[correction_data$Mismatch >=4] <- "unclassified"
correction_data$species[correction_data$PercIdent <= 98] <- "unclassified"
correction_data$genus[correction_data$PercIdent <= 97] <- "unclassified"
correction_data$family[correction_data$PercIdent <= 95] <- "unclassified"
correction_data$order[correction_data$PercIdent <= 93] <- "unclassified"
correction_data$phylum[correction_data$PercIdent <= 80] <- ""
#correction_data <- correction_data[!(is.na(correction_data$Phylum) | correction_data$Phylum==""), ]

my_data_decona <- as.data.frame(correction_data)


#           Rarefaction     

```
head(df_merge)

df<- dcast(my_data_decona, Sample_data ~ uniqueID, value.var = "Read_Count", fun.aggregate=sum)
df<- arrange(df,Sample_data)
rownames(df) <- df$Sample_data
df$Sample_data <-NULL 
matrix <- as.matrix(df)

#Plotting rarecurve 

raremax <- min(rowSums(matrix))
raremax
col <- c("black", "darkred", "forestgreen", "orange", "blue", "yellow", "hotpink")
lty <- c("solid", "dashed", "longdash", "dotdash")
pars <- expand.grid(col = col, lty = lty, stringsAsFactors = FALSE)
head(pars)
 
with(pars[1:12, ],
     rarecurve(matrix, step = 100, col = col,
               lty = lty, label= TRUE))




# Matrix  To make  OTU table. 
Matrix_data_frame <- data.frame(t(matrix))



#   adding the  matrix to the rest of the seq info
uniqueID <- rownames(Matrix_data_frame)
rownames(Matrix_data_frame) <- NULL
data <- cbind(uniqueID,Matrix_data_frame)

#here make sure you select and thus name the right column; only keep the ones that you need for your analysis
Tidy_data <- select(my_data_decona, uniqueID, Sample_data, Total_readcount, Alignmentlength , Read_Count,  phylum, class, order, family, genus, species)
merged_data <-left_join(Tidy_data, data, by= 'uniqueID')
merged_data = merged_data[!duplicated(merged_data$uniqueID),]



#   change your frame for OTU table                          

#deselect columns that are non- numeric except for uniqueID

OTU.frame <- merged_data  %>% ungroup() %>% select( -Sample_data, -Total_readcount, -Alignmentlength , -Read_Count, -phylum, -class, -order, -family, -genus, -species) 
head(OTU.frame)
OTU.frame <- as.data.frame(OTU.frame)
OTU.frame = OTU.frame[!duplicated(OTU.frame$uniqueID),]
OTU_frame <- OTU.frame[,-1]
rownames(OTU_frame) <- OTU.frame[,1]


# change your frame for TAXA table                         


#select the columns that are important 

TAX.table <- select(merged_data, uniqueID, phylum, order, family, genus, species)
head(TAX.table)
colnames(TAX.table) <- c("uniqueID", "phylum", "order", "family", "genus", "species")
TAX.matrix  <- as.matrix(TAX.table)
TAX_matrix <- TAX.matrix[,-1]
rownames(TAX_matrix) <- TAX.matrix[,1]

# Create your sample data                     


Samples <- as.data.frame(Metadata_ERIK_23_2kb_BLAST)
rownames(Samples) <- paste0(Samples$Sample_data)
Sample = sample_data(Samples)

# Produce your phyloseq object                        
#produce the correct table inputs

OTU = otu_table(OTU_frame, taxa_are_rows = TRUE)
TAX = tax_table(TAX_matrix)
SAMPLE = sample_data(Sample)

#phyloseq 

physeq_data = phyloseq(OTU, TAX, SAMPLE)
save(physeq_data, file="C:/Path/to/directory/physeq_data.RData")

 
