# GeneAnnotationsR

Package for getting gene-level annotations in R.

##### Firstly, all the protein coding genes are needed

-   A function called 'get_protein_coding_genes()' gets a df of all protein coding genes HGNC IDs.

#### List of gene annotations:

##### Constraint metrics

-   Alphamissense

-   Gnomad constraint metrics

-   Shet metrics

-   DOMINO metric

-   SCoNes metric

##### Cell line screens gene essentiality

-   Cancer cell line (DepMap) screen gene essentiality

-   hPSC screen gene essentiality

##### Mouse model gene essentiality

-   Mouse (MGI + IMPC) gene viability

##### Phenotype

-   Mouse and zebrafish Phenodigm scores (only available for NDD and cardiac rare disease genes)

##### Interactions

-   Chemical-gene interactions

##### Gene sequence

-   Gene expression (human and mouse, evodevo group)

-   Gene sequence annotations (gene length, coding sequence length, 3' and 5' UTR length)

-   *De novo* mutation rates (synonymous, missense, loss-of-function, splice site, nonsense, and all) (<https://pmc.ncbi.nlm.nih.gov/articles/PMC4222185/>)

##### Protein

-   Protein expression

-   Sub-cellular location

-   GO terms

#### Other functions available:

-   for summary stats, plotting

### Requirements:

-   BiocManager

-   biomaRt

-   dplyr

-   readr

-   tibble

-   data.table

-   magrittr

-   tidyr

-   org.Hs.eg.db

-   AnnotationDbi

-   GO.db

### Example:

Check rmarkdown.
