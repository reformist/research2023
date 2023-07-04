CMap LINCS 2020 Data Release (Beta)
December 2020
============================

This file describes the CMap LINCS resource released
in December 2020. This release is an expansion upon
the previous 2017 data release and contains ~3M gene
expression profiles and ~1M replicate-collapsed signatures.

Data are provided in GCTX format. There is one file
for each data level / perturbagen type combination, in
order to facilitate easier downloading. Metadata for the
experiments (columns) and genes (rows) of the matrices
are provided as tab-delimited text files. 

Note that while these static files are availabe for download,
the underlying data are also made availabe through Google
BigQuery, which enables extracting arbitrary subsets of
data. This data can be accessed at: 
https://console.cloud.google.com/bigquery?p=cmap-big-table&d=cmap_lincs_public_views&page=dataset 

Please also note that these data are released as a
beta version and will likely be subject to change as
issues are corrected and/or updates are made. File names
will include version numbers that will be updated as
the associated data changes.



Files
============================

siginfo.txt	Metadata for level 5 signatures
instinfo.txt	Metadata for levels 3 and 4 profiles
geneinfo.txt	Metadata for genes (the features of the data matrices)
cellinfo.txt	Metadata for cell lines
compoundinfo.txt	Metadata for cell lines (note there is one entry per compound/moa/target/structure combination, so some compounds appear more than once)

LINCS2020 Release Metadata Field Definitions.xlsx	Excel spreadsheet with definitions for the fields in each metadata file

These treatment abbreviations indicate which perturbagen
type(s) are included in the given GCTX file and the
conventions apply to all data levels.

trt_cp	compounds
trt_sh	shRNA
trt_oe	over-expression
trt_xpr	CRISPR
trt_misc	other treatments (ex: ligands, siRNA, etc)
ctl	negative controls (ex: DMSO, untreated)

These files may be updated periodically as new data
are generated and/or annotations are updated.

For reference
============================

GCTX format: https://clue.io/connectopedia/gctx_format
Data levels: https://clue.io/connectopedia/data_levels
