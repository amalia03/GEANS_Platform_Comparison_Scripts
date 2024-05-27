library(vegan)
library(MASS)
source("tax_viz_functions.R")

bl.files <- list.files("../contig_vs_GEANS_dtst/","reads.bl$")

bl.l <- lapply(1:length(bl.files), function(x){
    read.delim(paste0("../contig_vs_GEANS_dtst/",bl.files[x]), header=F,sep="\t", stringsAsFactors=FALSE, quote= "" )
})

sample.names <- sub("([0-9A-Z]+).*", "\\1",bl.files)
sample.names 

names(bl.l) <- sample.names
##Remove renegade column and rename the files

for(x in 1:length(bl.l)){    
     colnames(bl.l[[x]]) <- c('query', 'subject', 'title', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'ssend','evalue', 'score', 'length', 'pident', 'nident', 'gapopen', 'gaps', 'qcovs')
}
unlist(lapply(bl.l, function(x){length(unique(x$query))}))


#bl.l <- lapply(bl.l, function(x){x[x$qlen>150& x$length>=100,]})
bl.l <- lapply(bl.l, function(x){x[x$length>=100,]})
unlist(lapply(bl.l, function(x){length(unique(x$query))}))

m.ph <- sapply(unique(bl.l$query), function(x){length(unique(bl.l[bl.l$query==x,"phylum"]))})
print(paste0("Proportion of single phyla per query sequence ", round((sum(m.ph[m.ph==1])/length(m.ph))*100, 2), "%"))
s.ph <- m.ph[m.ph==1]
tax.info.sph <- bl.l[bl.l$query %in% names(s.ph), ]


bl.s.top <- lapply(1:length(bl.l),function(x){
    bl.s.i <- bl.l[[x]]
        bl.s.sc <- tapply(1:nrow(bl.s.i), bl.s.i$query, function(i){
        i[ order(bl.s.i[i, 'score']) ]
        })
    bl.s.sc <- sapply(bl.s.sc, function(x){ x[1] })
    bl.s.ts <- bl.s.i[bl.s.sc,]
    bl.s.ts
})

bl.tax <- lapply(bl.s.top, function(x){
    x.tax <- strsplit(x$subject, ";")
    x$genus <- unlist(lapply(x.tax, function(y){
        y[5]}))
    x$species <- unlist(lapply(x.tax, function(y){
        y[6]}))
    x$species <- paste(x$species, x$title) 
    x
})

bl.tax <- lapply(bl.tax, function(x){
    x[x$species=="Limecola balthica", "species"] <- "Macoma balthica"
    x[x$species=="Owenia fisiformis", "species"] <- "Owenia fusiformis"
    x[x$species=="Sagartia elegans", "species"] <- "Cylista elegans"
    x
})

bl.tax.t <- lapply(bl.tax, function(x){table(x$species)})
bl.tax.t <- lapply(bl.tax.t, function(x){x[x>=5]})
bl.tax.t
unlist(lapply(bl.tax.t, sum))

unique(names(unlist(lapply(1:3, function(x){bl.tax.t[[x]]}))))

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

write.table(gr.matrix, "geans_reads.tsv", sep="\t", quote=F) 
