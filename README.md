# FLASH-MM: Fast and Scalable Single-Cell Differential Expression Analysis Using Linear Mixed-Effects Models

**FLASH-MM** is a fast and scalable algorithm for differential expression (DE) analysis in large-scale single-cell RNA-seq (scRNA-seq) datasets. It addresses challenges such as intra-subject correlation, inter-subject variability, and the computational demands of analyzing millions of cells.

## Key Features
- **Efficient and Scalable**: Precomputes summary statistics to handle large datasets while maintaining single-cell resolution.
- **Accurate DE Analysis**: Controls type-I error rates and maintains high statistical power.
- **Broad Applications**: Supports case-control comparisons, cell-type-specific analyses, and multi-subject studies.
- **Simulation Tool**: Includes `simuRNAseq` for generating realistic scRNA-seq datasets.

## Applications
FLASH-MM has been applied to:
- Case-control comparisons in tuberculosis immune atlases.
- Cell-type-specific sex comparisons in kidney datasets.

With its speed, accuracy, and flexibility, FLASH-MM enables robust DE analysis for large-scale single-cell studies across diverse biological contexts.
