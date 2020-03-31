# Single Cell Human Kidney Organoid Analysis

<img src="https://github.com/jpromeror/SC_KidneyOrganoid_ACE2/blob/master/SuppFigure2.png?raw=true" width="760" height="544">

This repository includes the script for the analysis of the single cell kidney organoid for the paper 
*"Inhibition of SARS-CoV-2 infections in engineered human tissues using clinical-grade soluble human ACE2"*.

# Citation
Add paper citation

# Data
The repository includes an .rds file with a data.frame that includes the metadata of the Seurat object and the UMAP coordinates.

```{r, eval=FALSE}
KidneyOrganoid_MD<-readRDS("KidneyOrganoid_MetaData.rds")
```

The raw data can be retrieved from GEO under the accession number [GSEXXXX](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108349)
