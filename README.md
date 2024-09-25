In this depository, you will find the scripts and steps that were used for the bioinformatic analysis for our paper, *Paper X*. 
For more information on the metabarcoding steps, please refer to https://github.com/karlijn-doorenspleet/decona-postprocessing/tree/main. 

## Metabarcoding data processing
For Illumina Miseq data we used  
1. Seqtab – to excel
2. Excel to phyloseq

## Nanopore Sequencing Analysis

### Pre- processing and taxonomic linieage retreaval
**Taxonomy_assign_after_Decona.R** to download OTUcalled_finalblastoutput.txt. Also the accessionTaxa.sql file is needed (to be downloaded from NCBI).
Script content:
1.	load libraries (line 1-6)
2.	Read and load the necassary files (line 7- 15)
3.	Get the taxonomy and add to the dataframe (line 16-17)
4.	Filters data and p hits with highest e-values (line 18-27)
5.	Extract read count information (31-32)
6.	Add columns taxonomy and rename to unclassified if % identity =0 (line 33-56)
7.	Rename name in taxonomy to "unclassified" when there are multiple hits with the same e-value (line 57- 80)
8.	Store taxonomy names (species) in new column and keep highest hit (line 81-89)
9.	Add consensus sequence (line 90-96)
10.	Filter and rename data (line 97-114)
11.	Save final dataset (line 115-121)

### Processing and phyloseq object obtaining

**From_decona_to_phyloseq.R**
To try and check out file structures download and load R_list_github.xlxs, metadata_github.xlxs, and Decona_R_Object_github.Rdata to your environment.
Script content:
1.	Load packages and load objects (line 1-7)
2.	Rename some columns and rows (line 8-10)
3.	Add R_list with Sample_data to dataframe (line 11-16)
4.	Taxonomy assignment correction: allow identification on certain taxonomic level; based on subsettable criteria (line 20-30).
5.	Rarefaction curve of the data (line 31-60)
6.	Prepare OTU frame, Taxa table and Metadata for phyloseq object (line 61 -107)
7.	Produce phyloseq object (line 108-119)

Additionally perform Barcode leakage correction:

1.	Take 0.1% reads of each species (line 1-13)
2.	Substract from each sample (line 14-20)

## Shotgun Sequencing data processing 
Below you will find the the scripts that processed the NovaSeq sequences and their subsequent taxonomic assignment. 

### _In silico_ steps before taxonomic assignment processing. 
Starting with the direct output from NovaSeq:
1. **Quality filtering**: trimmed using **cutadapt** (--minimum-length=100, q=30)
2. **Sequence pairing**: Reads were merged with **PEAR** v 1.7.2 (Zhang et al., 2014) using the default parameters.
3. **Contig assembly**: Reads were assembled into contigs using **idba_ud** (Peng, 2012) using the default parameters.
4. **Read quantification**: Reads were mapped back to the assembled contigs using **BBMap** using minid=0.90.
5. **Taxonomic assignment**: Contigs were aligned to a customized version of NCBI-nt database using BLAST (-pid 97, -evalue 1E-10, -length 100)

The script **tax_geans.R** was then used to merge and process the taxonomically assigned contigs and their respective mapped reads for the data analysis. 

### Taxonomic assignment 
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
