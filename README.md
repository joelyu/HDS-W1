# Healthcare Data Science Module 1 Assignment

## Effects of age and gender on using gene expression levels as prognostic markers for diffuse large B cell lymphoma (DLBCL)

This is the submission GitHub repository for the Module 1 Assignment for the MSt program in Healthcare Data Science. The datasets are derived from, and are therefore fully attribute to Reddy et. al for their work on Genetic and Functional Drivers of Diffuse Large B Cell Lymphoma ([ScienceDirect link here](https://www.sciencedirect.com/science/article/pii/S0092867417311212)), namely the [S1 Table](https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311212-mmc1.xlsx) and [S2 Table](https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311212-mmc2.xlsx) from their Supplementary Materials. See [Data Attribution](#data-attribution) for details.

---

### Data Attribution

#### Main citation

Anupama Reddy, Jenny Zhang, Nicholas S. Davis, Andrea B. Moffitt, Cassandra L. Love, Alexander Waldrop, Sirpa Leppa, Annika Pasanen, Leo Meriranta, Marja-Liisa Karjalainen-Lindsberg, Peter Nørgaard, Mette Pedersen, Anne O. Gang, Estrid Høgdall, Tayla B. Heavican, Waseem Lone, Javeed Iqbal, Qiu Qin, Guojie Li, So Young Kim, Jane Healy, Kristy L. Richards, Yuri Fedoriw, Leon Bernal-Mizrachi, Jean L. Koff, Ashley D. Staton, Christopher R. Flowers, Ora Paltiel, Neta Goldschmidt, Maria Calaminici, Andrew Clear, John Gribben, Evelyn Nguyen, Magdalena B. Czader, Sarah L. Ondrejka, Angela Collie, Eric D. Hsi, Eric Tse, Rex K.H. Au-Yeung, Yok-Lam Kwong, Gopesh Srivastava, William W.L. Choi, Andrew M. Evens, Monika Pilichowska, Manju Sengar, Nishitha Reddy, Shaoying Li, Amy Chadburn, Leo I. Gordon, Elaine S. Jaffe, Shawn Levy, Rachel Rempel, Tiffany Tzeng, Lanie E. Happ, Tushar Dave, Deepthi Rajagopalan, Jyotishka Datta, David B. Dunson, Sandeep S. Dave,

Genetic and Functional Drivers of Diffuse Large B Cell Lymphoma, Cell, Volume 171, Issue 2, 2017, Pages 481-494.e15, ISSN 0092-8674, [https://doi.org/10.1016/j.cell.2017.09.027](https://www.sciencedirect.com/science/article/pii/S0092867417311212)

Abstract: Summary

Diffuse large B cell lymphoma (DLBCL) is the most common form of blood cancer and is characterized by a striking degree of genetic and clinical heterogeneity. This heterogeneity poses a major barrier to understanding the genetic basis of the disease and its response to therapy. Here, we performed an integrative analysis of whole-exome sequencing and transcriptome sequencing in a cohort of 1,001 DLBCL patients to comprehensively define the landscape of 150 genetic drivers of the disease. We characterized the functional impact of these genes using an unbiased CRISPR screen of DLBCL cell lines to define oncogenes that promote cell growth. A prognostic model comprising these genetic alterations outperformed current established methods: cell of origin, the International Prognostic Index comprising clinical variables, and dual MYC and BCL2 expression. These results comprehensively define the genetic drivers and their functional roles in DLBCL to identify new therapeutic opportunities in the disease.

Keywords: exome sequencing; genetic mutations; diffuse large B cell lymphoma; DLBCL; TCGA; The Cancer Genome Atlas

#### Data Source

- [Supplementary Table S1](https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311212-mmc1.xlsx): Clinical Information and Genetic Alteration Data for 1,001 DLBCL Samples
- [Supplementary Table S2](https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311212-mmc2.xlsx): ABC/GCB Classification Using Gene Expression Data

---

### Repository file structure

```
├── data
│   ├── mmc1-ClinicalInformation.csv
│   └── mmc2-GeneExpression.csv
├── supplementary-materials/
│   ├── 01-data-preprocessing.R
│   └── 02-data-processing.R
│   ├── 03-data-exploration.R
│   ├── 04-data-analysis.R
│   ├── 05-response-analysis.R
│   ├── 06-effect-analysis.R
│   ├── 07-effect-visualisation.R
├── LICENSE.txt
├── Module1-assignment.html
├── Module1-assignment.Rmd
├── README.md
├── references.bib
└── Source.txt
```

### Supplementary Analysis Scripts

The `supplementary-materials/` directory contains reproducible scripts for data processing and quality assessment:

- **[01-data-preprocessing.R](supplementary-materials/01-data-preprocessing.R)** - Downloads and preprocesses raw Excel data from ScienceDirect to CSV format
- **[02-data-processing.R](supplementary-materials/02-data-processing.R)** - Combines datasets and performs comprehensive data quality analysis
- **[03-data-exploration.R](supplementary-materials/03-data-exploration.R)** - Initial visual inspection and assumption checking
- **[04-data-analysis.R](supplementary-materials/04-data-analysis.R)** - ANOVA and t-tests

### Data cleaning

Both tables ([Supplementary Table S1](https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311212-mmc1.xlsx) & [Supplementary Table S2](https://ars.els-cdn.com/content/image/1-s2.0-S0092867417311212-mmc2.xlsx)) were downloaded, and the first sheets named `Clinical Information` and `Gene Expression` respectively were exported as `mmc1-clinicalInformation.csv` and `mmc2-GeneExpression.csv`, both with the first 3 rows removed. The [`01-data-preprocessing.R`](supplementary-materials/01-data-preprocessing.R) script will pull the datasets from source, extract the right sheets from the excel files, clean headers (3 rows) and export to cleaned `.csv` files. If all else fails, this simple cleaning could be done via Microsoft Excel or any other major spreadsheet software.

### License and Usage

This dataset is used under [Elsevier User License](https://www.elsevier.com/about/policies-and-standards/open-access-licenses/elsevier-user) for
non-commercial academic purposes only. See [LICENSE.txt](LICENSE.txt) for details.
