# sysndd_gene_table.xlsx, sysndd_raw_to_rdata.R, sysndd_genes.RData
Sheet was downloaded from https://sysndd.dbmr.unibe.ch/Genes on Wednesday
12 March 2025. Website version v0.1.0.
The R script sysndd_raw_to_rdata.R was used to convert the raw data from
xlsx to RData format. The RData file contains the data in a df format.


# DDG2P_2025-01-28.csv.gz, ddg2p_raw_to_rdata.R, ddg2p_genes.RData
https://ftp.ebi.ac.uk/pub/databases/gene2phenotype/G2P_data_downloads_archive/2025_01_28/DDG2P_2025-01-28.csv.gz
Accessed on Thursday 13MAR25. Accessing the previous one 2025-01-28 (not the latest,
2025-02-28) because they haven't assigned organ to the latest. Merging together
definitive and strong to keep for the seed genes. Then keeping as moderate and
limited genes as it is.
The R script ddg2p_raw_to_rdata.R was used to convert the raw data from
.csv.gz to RData format. The RData file contains the data in a df format.


# gel_panelapp_raw.R, gel_panelapp_genes.csv
R script to download the raw data from GEL panelapp. 13MAR25.
The only modification to the file is a final column with the panel number.
These raw files are not used by the functions, but are provided for transparency.


# aus_panelapp_raw.R, aus_panelapp_genes.csv
R script to download the raw data from Australia panelapp. 13MAR25.
The only modification to the file is a final column with the panel number.
These raw files are not used by the functions, but are provided for transparency.


# raw_protein_coding_genes.R, protein_coding_genes.csv.gz
R script to download the raw protein coding genes file from HGNC. 13MAR25.
These raw files are not used by the functions, but are provided for transparency.



# set_internal_data.R, sysdata.rda
When all RData files are ready, the script set_internal_data.R should be run to
set the internal data. This is found in R/
