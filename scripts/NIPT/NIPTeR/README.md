INTRODUCTION
------------

Non-Invasive Prenatal Testing (NIPT) is a method based on analysis of cell-free fetal (cff) DNA in maternal blood using next-generation sequencing. The basic approach is relatively straightforward: cell-free DNA is isolated from the mothers blood plasma and the obtained DNA fragments are sequenced. Subsequently the number of DNA fragments originating from the different chromosomes is counted, since in case of a fetal trisomy a relative increase of the fraction of reads of the trisomic chromosome is expected.

Various analysis methods for trisomy prediction, variation reduction and quality control have been described in the literature. These methods can be combined to improve the sensitivity of the analysis. The NIPTeR package makes these methods available and allows for flexibility in workflow design.

NIPTeR uses BAM files as input and calculates chromosomal fractions, which are used for trisomy prediction. These fractions are based on read counts of the different chromosomes. The chromosomes are divided in bins to allow for read count correction.

INSTALLATION INSTRUCTIONS
-------------------------

NIPTeR is dependent on three packages: sets, RSamtools and S4Vectors. For installation use:

install.packages("sets") source("<http://bioconductor.org/biocLite.R>") biocLite("Rsamtools") install.packages("NIPTeR")

EXAMPLES AND OVERVIEW
---------------------

Workflow examples and components overview are included in the NIPTeR vignette.
