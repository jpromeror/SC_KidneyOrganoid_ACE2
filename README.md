# Single Cell Kidney Organoid Analysis

![**Figure 1.** SingleCell Human Kidney Organoid ](https://github.com/jpromeror/EventPointer/blob/master/vignettes/Figure1.png)

This repository includes the script for the analysis of the single cell kidney organoid for the paper "Inhibition of SARS-CoV-2 infections in engineered human tissues using clinical-grade soluble human ACE2".

# Citation
Add paper citation

# Data
The repository includes an .rda file with a data.frame that includes the metadata of the Seurat object and the UMAP coordinates.

```{r, eval=FALSE}
KidneyOrganoid<-readRDS("KidneyOrganoid_MetaData.rds")
```

The raw data can be retrieved from GEO under the accession number GSEXXXX
