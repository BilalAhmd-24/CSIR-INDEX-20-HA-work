
# the essentials ----
# this chunk contains the minimal essential code from this script. Simply uncomment the lines below and run the code.
library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package

#Change your path where you have .tsv files
setwd("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/alignment")
#Load the studydesign file
targets <- read_tsv("studydesign_1.tsv")# read in your study design

#Load the path where all files are located 
path <- file.path(targets$sample, "abundance.tsv") # set file paths to your mapped data
path

#Annotating, retrieving the transcripts from Ensemble Human sapiens
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")

# now check to make sure this path is correct by seeing if the files exist
(file.exists(path)) 

Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = Tx,
                     #If txOut is TRUE, it will give transcript level data
                     txOut = FALSE, #determines whether your data represented at transcript or gene level
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)
#take a look at the type of object you just created
class(Txi_gene)
names(Txi_gene)

print("Step 1 complete!")

# if you want to write your counts or TPM to a file on your harddrive
capture.output(summary(Txi_gene), file = "My New File.txt")
capture.output(Txi_gene, file = "Txi_list_Index.csv")
