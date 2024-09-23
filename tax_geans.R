#Locate the BLAST files (here it assumes it is in the same directory and they have the names written as such "X_reads.bl")
bl.files <- list.files("./","reads.bl$")

#Import the files in a list
bl.l <- lapply(1:length(bl.files), function(x){
    read.delim(paste0("../contig_vs_GEANS_dtst/",bl.files[x]), header=F,sep="\t", stringsAsFactors=FALSE, quote= "" )
})

#Use the file name as an ID
sample.names <- sub("([0-9A-Z]+).*", "\\1",bl.files)
sample.names 

names(bl.l) <- sample.names

#Name BLAST columns
for(x in 1:length(bl.l)){    
     colnames(bl.l[[x]]) <- c('query', 'subject', 'title', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'ssend','evalue', 'score', 'length', 'pident', 'nident', 'gapopen', 'gaps', 'qcovs')
}

#Enforce alignment quality filters (alignment length and/or percent similarity)
bl.l <- lapply(bl.l, function(x){x[x$length>=100,]})

#For each query, remove matches that aligned to more than one phylum (high taxonomic ambiguity)
m.ph <- sapply(unique(bl.l$query), function(x){length(unique(bl.l[bl.l$query==x,"phylum"]))})
print(paste0("Proportion of single phyla per query sequence ", round((sum(m.ph[m.ph==1])/length(m.ph))*100, 2), "%"))
s.ph <- m.ph[m.ph==1]
tax.info.sph <- bl.l[bl.l$query %in% names(s.ph), ]

#Keep only the top scoring (bitscore) match per contig
bl.s.top <- lapply(1:length(bl.l),function(x){
    bl.s.i <- bl.l[[x]]
        bl.s.sc <- tapply(1:nrow(bl.s.i), bl.s.i$query, function(i){
        i[ order(bl.s.i[i, 'score']) ]
        })
    bl.s.sc <- sapply(bl.s.sc, function(x){ x[1] })
    bl.s.ts <- bl.s.i[bl.s.sc,]
    bl.s.ts
})

#Import the mapped reads per contig from the bbmap analysis
quant.l <- read.delim("file_fpkm.tsv", header=F,sep= "\t", stringsAsFactors=FALSE, quote= "" )
colnames(quant.l) <- c("id", "read_length",  "bases", "coverage", "reads", "RPKM", "frags", "FPKM")}

#Extract the contig from the fpkm file identifier. 
quant.id <- strsplit(quant.l$id, split = " ")
quant.id <- lapply(all.id, function(y){(y[1])})
quant.l.id <- cbind(quant.l, quant.id)

#Then merge the taxonomically assigned dataframe with the read mapping dataframe. 
bl.top.quant <- merge(bl.top.sp[[x]], quant.l, by="query"))

#Make a genus and species column based on the match identifiers (which usually contain the complete taxonomy of each match)
bl.tax <- lapply(bl.s.top, function(x){
    x.tax <- strsplit(x$subject, ";")
    x$genus <- unlist(lapply(x.tax, function(y){
        y[5]}))
    x$species <- unlist(lapply(x.tax, function(y){
        y[6]}))
    x$species <- paste(x$species, x$title) 
    x
})

#Rename certain species with defunct/incorrect appelations 
bl.tax <- lapply(bl.tax, function(x){
    x[x$species=="Limecola balthica", "species"] <- "Macoma balthica"
    x[x$species=="Owenia fisiformis", "species"] <- "Owenia fusiformis"
    x[x$species=="Sagartia elegans", "species"] <- "Cylista elegans"
    x[x$species=="Tritia reticulatus", "species"] <- "Tritia reticulata" 
    x[x$species=="Ensis directus", "species"] <- "Ensis leei"
    x[x$species=="Venerupis senegalensis", "species"] <- "Venerupis corrugata"
    x[x$species=="Haliplanella lineata", "species"] <- "Diadumene lineata"
    x
})

##Set a limit to the lowest number of reads matching to each species (here set to 5 matches per species)
bl.tax.t <- lapply(bl.tax, function(x){table(x$species)})
bl.tax.t <- lapply(bl.tax.t, function(x){x[x>=5]})
bl.tax.t
unlist(lapply(bl.tax.t, sum))

##Create an abundance table for each species per location
all.tax <- unique(unlist(lapply(bl.tax.t, names)))
gr.matrix <- matrix(0, nrow=length(all.tax), ncol=12)
rownames(gr.matrix) <- sort(all.tax)
colnames(gr.matrix) <- sample.names

for(i in 1:length(sample.names)){
    for(j in names(bl.tax.t[[i]])){
        if(which(rownames(gr.matrix)==j)>0){
            gr.matrix[which(rownames(gr.matrix)==j),i] <- bl.tax.t[[i]][names(bl.tax.t[[i]])==j]
        }}
}

#Export the matrix in a table.
write.table(gr.matrix, "geans_reads.tsv", sep="\t", quote=F) 
