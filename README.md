## Shotgun Sequencing data processing 
This directory serves as a script depository for *Paper X*. It includes scripts that process the NovaSeq sequences and their subsequent taxonomic assignment. 

### _In silico_ steps before taxonomic assignment processing. 
Starting with the direct output from NovaSeq:
1. **Quality filtering**: trimmed using **cutadapt** (--minimum-length=100, q=30)
2. **Sequence pairing**: Reads were merged with **PEAR** v 1.7.2 (Zhang et al., 2014) using the default parameters.
3. **Contig assembly**: Reads were assembled into contigs using **idba_ud** (Peng, 2012) using the default parameters.
4. **Read quantification**: Reads were mapped back to the assembled contigs using **BBMap** using minid=0.90.
5. **Taxonomic assignment**: Contigs were aligned to a customized version of NCBI-nt database using BLAST (-pid 97, -evalue 1E-10, -length 100)

The script **tax_geans.R** was then used to merge and process the taxonomically assigned contigs and their respective mapped reads for the data analysis. 

## Taxonomic assignment 
Below is the script content for the file tax_script.R:

1. Import and setup taxonomic assignment dataframe (lines 1-18).
2. Apply additional alignment filtering steps (line 21).
3. Remove taxonomic assignmnets to more than one phyla per contig (lines 23-27).
4. Keep best scoring match per query (lines 29-38).
5. Import and setup read mapping dataframe (lines 40-42).
6. Merge the read mapping and taxonomic assignment dataframes together (lines 44-50)
7. Set up a genus and a species column based on the database match (lines 52-61)
8. Rename species with defunct/incorrect appelations (lines 63-73)
9. Set a limit to the lowest number of reads matching per species (lines 75-79)
10. Create and export a matrix with the read counts per species and location (lines 81-95)
