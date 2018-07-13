## DWLS- Gene expression deconvolution using dampened weighted least squares

Dampened weighted least squares (DWLS) is an estimation method for gene expression deconvolution, in which the cell-type composition of a bulk RNA-seq data set is computationally inferred.  This method corrects common biases towards cell types that are characterized by highly expressed genes and/or are highly prevalent, to provide accurate detection across diverse cell types.  To begin, the user must input a bulk RNA-seq data set, along with a labeled representative single-cell data set that will serve to generate cell-type-specific gene expression profiles.  Ideally, the single-cell data set will contain cells from all cell types that may be found in the bulk data.  DWLS will return the cell-type composition of the bulk data.

DWLS is written in the R programming language.  To apply DWLS to any RNA-seq data set, please refer to the manual, **Manual.pdf**.  This includes code to run DWLS on a simple intestinal stem cell example.  This project includes a script containing the core functions, under the **Deconvolution_functions.R**, as well as three applications of DWLS to various data sets.  These are each contained in a separate folder, and include:
 
 - **Simulation_Schelker**: This application uses data from a previous deconvolution analysis by Schelker et al. [1], where bulk data is simulated by adding together gene expression profiles from single-cell data that is derived from human donor peripheral blood mononuclear cells (PBMCs), tumor-derived melanoma patient samples, and ovarian cancer ascites samples.
 - **ISC**: This analysis deconvolves bulk intestinal stem cell RNA-seq data under three conditions, using a gene expression signature created from intestinal stem cell single-cell data sets under multiple conditions.  Both single-cell and bulk data sets are taken from Yan et al. [2].
 - **MCA**: We generate bulk RNA-seq data from four healthy mouse tissues, and deconvolve each using a comprehensive set of single-cell data sets from the Mouse Cell Atlas from Han et al. [3].
 
 ###### Data sources:<br />
 <sub> [1] Schelker M, Feau S, Du J, Ranu N, Klipp E, MacBeath G, Schoeberl B, Raue A: Estimation of immune cell content in tumour tissue using single-cell RNA-seq data. Nat Commun 2017, 8:2032. <br />
 [2] Yan KS, Janda CY, Chang J, Zheng GXY, Larkin KA, Luca VC, Chia LA, Mah AT, Han A, Terry JM, et al: Non-equivalence of Wnt and R-spondin ligands during Lgr5. Nature 2017, 545:238-242. <br />
 [3] Han X, Wang R, Zhou Y, Fei L, Sun H, Lai S, Saadatpour A, Zhou Z, Chen H, Ye F, et al: Mapping the Mouse Cell Atlas by Microwell-Seq. Cell 2018, 172:1091-1107.e1017.
