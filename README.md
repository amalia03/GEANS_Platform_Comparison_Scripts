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
