"Let's try studyng Dexa LA (08) versus Dexa HA (07)"

library(tidyverse) # provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tximport) # package for getting Kallisto results into R
library(ensembldb) #helps deal with ensembl
library(EnsDb.Hsapiens.v86) #replace with your organism-specific database package

#Change your path where you have .tsv files
setwd("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/alignment/")

#Load the studydesign file
targets <- read_tsv("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/studydesign_group3.tsv")# read in your study design

#Load the path where all files are located 
path <- file.path(targets$sample, "abundance.tsv") # set file paths to your mapped data
path

#Annotating, retrieving the transcripts from Ensemble Human sapiens
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))
Tx <- as_tibble(Tx)
Tx <- dplyr::rename(Tx, target_id = tx_id)
Tx <- dplyr::select(Tx, "target_id", "gene_name")


#Checking the annotations, transcript_mapped_to_Gene
Tx
# now check to make sure this path is correct by seeing if the files exist
(file.exists(path)) 


"Kallisto estimates transcript abundances in terms of TPM (Transcripts Per Million), 
which is a normalized measure that accounts for both the length of the transcript 
and the total number of reads in the sample. However, TPM is not a raw count and 
is not suitable for use as input to DESeq2 or other tools that expect count data.
To convert TPMs to raw counts, tximport uses the lengthScaledTPM method, 
which scales the TPM values by the length of the transcript and the total 
number of mapped reads, and then rounds the result to the nearest integer count value. 
This method is a simple way to convert TPMs to counts that can be used as input to 
differential expression analysis tools"

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

# if you want to write your summary counts or TPM to a file on your harddrive
#capture.output(summary(Txi_gene), file = "C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_vs_HA_samples(16_14)/TPM_counts_30_samples.txt")
# if you want to write your counts or TPM to a file on your harddrive
#capture.output(Txi_gene, file = "C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_vs_HA_samples(16_14)/TPM_counts_30_samples.xlsx")
"This is not a perfect way to store files"


View(Txi_gene)


"______________________STEP-2_:Lets's visualise data and filter it, normalize it__________________________________________________________________"
# the essentials ----
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)

# Examine your data up to this point ----
#Putting Txi_gene list (abundance) to myTPM object
#Converting the list object to array for operations
myTPM <- Txi_gene$abundance
View(myTPM)

#This is a matrix array
class(myTPM)

#Lets capture the output and store it to  hard disk as myTPM.csv
"Probably next time?"
capture.output(myTPM, file = "C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/myTPM.csv")
" What a poor way to save files?"

#Putting the counts data into myCounts
myCounts <- Txi_gene$counts

#Lets check myCOunts
View(myCounts)

#See the sum of columns of all those counts // whether its equal?
colSums(myTPM)

#See the sum of columns of all those counts // whether its equal?
colSums(myCounts)

"A count is simply that, a count of reads on some feature"
"An abundance is a more biologically meaningful (though not necessarily statistically useful)"
"quantification of expression of a gene or transcript that is normalized in some way. "

#Loading the  sample labels 
sampleLabels <- targets$sample
sampleLabels

"we will be working with TPM onwards"

# Generate summary stats for your data ----
# 1st, calculate summary stats for each transcript or gene, and add these to your data matrix
# then use the base R function 'transform' to modify the data matrix (equivalent of Excel's '=')
# then we use the 'rowSds', 'rowMeans' and 'rowMedians' functions from the matrixStats package
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM),
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))



# look at what you created
"AVerage and Median should be somewhat similar in ideal case?"
head(myTPM.stats)
"You added three columns of SD, AVG, MED in the last to TPM Object"

# Create your first plot using ggplot2 ----
# produce a scatter plot of the transformed data
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape =16, size=2)+
  geom_hex(show.legend = FALSE) +
  theme_half_open()

#Try with a better one?
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED, ) +
  geom_hex(shape=23, size=3)+
  theme_classic() +
  theme_dark() + 
  theme_bw()

"Remove this comment next time"
ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure1_Scatterplot_TPM_stats",  width = 4, height = 4, device='tiff', dpi=500)

# Experiment with point shape and size in the plot above
# Experiment with other plot types (e.g. 'geom_hex' instead of 'geom_point')
# Add a theme to your ggplot code above.  Try 'theme_bw()'
# How would these graphs change if you log2 converted the data?

# Let's expand on the plot above a bit more and take a look at each 'layer' of the ggplot code
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=16, size=2) +
  geom_smooth(method=lm) +
  geom_hex(show.legend = FALSE) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized LA versus HA data",
       caption="INDEX - CMREF 2020") +
  theme_classic() +
  theme_dark() + 
  theme_bw()
ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure2_Scatterplot_TPM_stats_Line",  width = 4, height = 4, device='tiff', dpi=500)


"Crucial step: Now, lets deal with DGE data"

"#1Figure out, Why we are using counts here again?: Resolved comment above"

#Creating a DGE list with function DGE list from EdgeR package
"Function used in DGEList byEdgeR"
#see, counts again here
myDGEList <- DGEList(Txi_gene$counts)
View(myDGEList)
class(myDGEList)

#DEGList objects are a good R data file to consider saving to you working directory
#save(myDGEList, file = "C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_vs_HA_samples(16_14)/myDGEList")

#Compute counts per million (CPM) or reads per kilobase per million (RPKM)
log2.cpm <- cpm(myDGEList, log=TRUE)
View(log2.cpm)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
colnames(log2.cpm.df) <- c("geneID", sampleLabels)

log2.cpm.df.pivot <- pivot_longer(log2.cpm.df, # dataframe to be pivoted
                                  cols = sampleLabels, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

p1 <- ggplot(log2.cpm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM) LA versus HA",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  coord_flip()

ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure3_Violin_Plot_Unfiltered_Non_Normalized",  width = 10, height = 4, device='tiff', dpi=500)

#Alternate graphs without samples information
#__________________________________________________________
"For Plot-1"
Alternate_df <- log2.cpm.df
new_sampleLabels <-  c( "Sample-1", "Sample-2", "Sample-3", "Sample-4", "Sample-5",
                          "Sample-6", "Sample-7", "Sample-8", "Sample-9", "Sample-10", "Sample-11", "Sample-12", "Sample-13", "Sample-14", "Sample-15")
colnames(Alternate_df) <- c("geneID", new_sampleLabels)
Alternate_df.pivot <- pivot_longer(Alternate_df, # dataframe to be pivoted
                                  cols = new_sampleLabels, # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)
pq1 <- ggplot(Alternate_df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM) LA versus HA",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  coord_flip()
ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure3_BlindedViolin_Plot_Unfiltered_Non_Normalized",  width = 10, height = 4, device='tiff', dpi=500)

"For Plot-2"
Alternate1_df<-log2.cpm.filtered
Alternate1.cpm.filtered.df <- as_tibble(Alternate1_df, rownames = "geneID")
colnames(Alternate1.cpm.filtered.df) <- c("geneID", new_sampleLabels)
Alternate1.cpm.filtered.df.pivot <- pivot_longer(Alternate1.cpm.filtered.df, # dataframe to be pivoted
                                           cols = new_sampleLabels, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

pq2 <- ggplot(Alternate1.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  coord_flip()
ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure4_BlindedViolin_Plot_filtered_Non_Normalized",  width = 10, height = 4, device='tiff', dpi=500)

"For third one"
Alternate2_df <-log2.cpm.filtered.norm.df
colnames(Alternate2_df) <- c("geneID", new_sampleLabels)


Alternate2_df.pivot <- pivot_longer(Alternate2_df, # dataframe to be pivoted
                                                cols = new_sampleLabels, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

pq3 <- ggplot(Alternate2_df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  coord_flip()
ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure5_BlindedViolin_Plot_filtered_TMM_Normalized",  width = 10, height = 4, device='tiff', dpi=500)


plot_grid(pq1, pq2, pq3, labels = c('A', 'B', 'C'), label_size = 12)

"For Cluster Dendogram"
plot(clusters, labels=new_sampleLabels)
#_________________________________________________________________________________________
#first, take a look at how many genes or transcripts have no read counts at all
"Focus on this, it was 10 first"
table(rowSums(myDGEList$counts==0)==1)
#1022 genes were absent, and have no counts
#344226 genes have some data


#Compute counts per million (CPM) or reads per kilobase per million (RPKM).
cpm <- cpm(myDGEList)

# now set some cut-off to get rid of genes/transcripts with low counts
# again using rowSums to tally up the 'TRUE' results of a simple evaluation
# how many genes had more than 1 CPM (TRUE) in at least 3 samples
#How many rows/genes have counts per million greater than one in more than 10 samples
#Focus on this, it was 10 first
keepers <- rowSums(cpm>1)>=7 #user defined

#now use base R's simple subsetting method to filter your DGEList based on the logical produced above
myDGEList.filtered <- myDGEList[keepers,]

#Making filtered data_________________
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
colnames(log2.cpm.filtered.df) <- c("geneID", sampleLabels)
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = sampleLabels, # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

p2 <- ggplot(log2.cpm.filtered.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  coord_flip()
ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure4_Violin_Plot_filtered_Non_Normalized",  width = 10, height = 4, device='tiff', dpi=500)


#Making the filtered data using Trimmed trimmed mean of M values #EdgeR
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")
colnames(log2.cpm.filtered.norm.df) <- c("geneID", sampleLabels)


log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = sampleLabels, # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

p3 <- ggplot(log2.cpm.filtered.norm.df.pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()+
  coord_flip()
ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure5_Violin_Plot_filtered_TMM_Normalized",  width = 10, height = 4, device='tiff', dpi=500)


plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)
print("Step 2 complete!")

"________________________STEP:3_______________________________________________________"

# the essentials ----
library(tidyverse)
library(DT)
library(gt)
library(plotly)

# Identify variables of interest in study design file ----
group <- targets$Group
group <- factor(group)


# Hierarchical clustering ---------------
#hierarchical clustering can only work on a data matrix, not a data frame
#try using filtered and unfiltered data...how does this change the result?
#try other distance methods (e.g. switch from 'maximum' to 'euclidean')...how does this change the result?
distance <- dist(t(log2.cpm.filtered.norm), method = "euclidean") #other distance methods are "euclidean", maximum", "manhattan", "canberra", "binary" or "minkowski"
clusters <- hclust(distance, method = "average") #other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
plot(clusters, labels=sampleLabels)
plot(clusters, labels=new_sampleLabels)
#ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_vs_HA_samples(16_14)/Figure6_BlindedEuclidean_Cluster_Dendogram",  width = 10, height = 4, device='tiff', dpi=500)




# Principal component analysis (PCA) -------------
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)

#look at the PCA result (pca.res) that you just created
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
screeplot(pca.res)# A screeplot is a standard way to view eigenvalues for each PCA

pc.var <- pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per <- round(pc.var/sum(pc.var)*100, 1) 
pca.res.df <- as_tibble(pca.res$x)

#Visualize your PCA result ------------------
#lets first plot any two PCs against each other
#We know how much each sample contributes to each PC (loadings), so let's plot
pca.res.df <- as_tibble(pca.res$x)
pca.plot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, label=new_sampleLabels, color = group) +
  geom_point(size=4) +
  stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure7_PCA_plot_HA_vs_LA_42",  width = 10, height = 4, device='tiff', dpi=500)
ggplotly(pca.plot)
 
#Create a PCA 'small multiples' chart ----
  # this is another way to view PCA laodings to understand impact of each sample on each pricipal component
pca.res.df <- pca.res$x[,1:4] %>% # note that this is the first time you've seen the 'pipe' operator from the magrittr package
as_tibble() %>%
add_column(sample = new_sampleLabels,
           group = group)

pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=group) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()
ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure9Blinded_PCA_Small_multiple_plots",  width = 10, height = 4, device='tiff', dpi=500)


#Calculating gene expression data 
view(log2.cpm.filtered.norm.df)
#Assigning te data into another variable
mydata.df <- log2.cpm.filtered.norm.df

#Changing the names of columns into character vector
colnames(mydata.df) <- c('GeneID', 'La1','La2','La3','La4','La5','La6','La7','La8',
                         'Ha1','Ha2','Ha3','Ha4','Ha5','Ha6','Ha7')
view(mydata.df)

#Calculating the average gene expression data
mydata.df <- mydata.df %>% 
  mutate(Low_altitude.AVG = (La1+ La2 + La3 + La4 + La5 + La6 + La7 + La8)/8,
         High_altitude.AVG =(Ha1+ Ha2 + Ha3 + Ha4 + Ha5 + Ha6 + Ha7)/7,
         #now make columns comparing each of the averages above that you're interested in
         #LogFC would be, Genes expressed higher in Low altitude would be positive
         #Whilw genes overexpressed at High-altitude would be negative
         #now make columns comparing each of the averages above that you're interested in
         LogFC = (High_altitude.AVG - Low_altitude.AVG)) %>% #note that this is the first time you've seen the 'pipe' operator
  mutate_if(is.numeric, round, 2)


#now look at this modified data table
mydata.df

# Use dplyr 'arrange' and 'select' to sort your dataframe based on any variable
# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
mydata.sort <- mydata.df %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(GeneID, LogFC)

view(mydata.sort)
# Use dplyr "filter" and "select" functions to pick out genes of interest 
# ways to tweak the 'select' function:
# use ':' between two column names to select all columns between
# use 'contains', 'starts_with' or 'ends_with' to modify how you select
# can refer to columns using exact name or numerical indicator
# use boolean operators such as '&' (and), '|' (or), '==' (equal to), '!' (not)
library(dplyr)
mydata.filter <- mydata.df %>%
  dplyr::filter(GeneID=="NOS3" | geneID=="EGLN1" | geneID=="EPAS1" | geneID=="ET1" | %>%
                  dplyr::select(GeneID, Low_altitude.AVG, High_altitude.AVG) %>%
                  dplyr::arrange(desc(LogFC))
                
view(mydata.filter)   

# you can also filter based on any regular expression
mydata.grep <- mydata.df %>%
  dplyr::filter(grepl('EPAS1|EGLN1|NOS3', GeneID)) %>%
  dplyr::select(GeneID, Low_altitude.AVG, High_altitude.AVG, LogFC) %>%
  dplyr::arrange(desc(GeneID))

# Produce publication-quality tables using the gt package ----
gt(mydata.grep)

# Make an interactive scatter plot with plotly -----
# begin by storing your ggplot object
myplot <- ggplot(mydata.df) + 
  aes(x=Low_altitude.AVG, y=High_altitude.AVG) +
  geom_point(shape=16, size=1) +
  ggtitle("Low_altitude.AVG vs. High_altitude.AVG") +
  theme_bw()

ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure10_HighaltiudevsLowaltitude_plot",  width = 10, height = 4, device='tiff', dpi=500)

p1 <- ggplotly(myplot)


#let's customize this graphic by adding a more informative mouseover tooltip
myplot <- ggplot(mydata.df) +
  aes(x=Low_altitude.AVG, y=High_altitude.AVG, 
      text = paste("Symbol:", GeneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("Low_altitude.AVG vs. High_altitude.AVG") +
  theme_bw()

ggplotly(myplot)


"________________________STEP:4_______________________________________________________"
# the goal of this script is to identify differentially expressed genes (DEGs) and differential transcript usage (DTU)

# Load packages -----
library(tidyverse) # you know it well by now!
library(limma) # venerable package for differential gene expression using linear modeling
library(edgeR)
library(gt)
library(DT)
library(plotly)

# Set up your design matrix AGAIN----
group <- factor(targets$Group)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Model mean-variance trend and fit linear model to data ----
# Use VOOM function from Limma package to model the mean-variance relationship
#Transform RNA-Seq Data Ready for Linear Modelling
v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = FALSE)
fit <- lmFit(v.DEGList.filtered.norm, design)

# Contrast matrix ----
contrast.matrix <- makeContrasts(Exposure = High_altitude - Low_altitude,
                                 levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)

# TopTable to view DEGs -----
"First_Main_Table_to_view"
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC")
view(myTopHits)


# convert to a tibble
myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")


gt(myTopHits.df)
#TopTable (from Limma) outputs a few different stats:
# logFC, AveExpr, and P.Value should be self-explanatory
# adj.P.Val is your adjusted P value, also known as an FDR (if BH method was used for multiple testing correction)
# B statistic is the log-odds that that gene is differentially expressed. If B = 1.5, then log odds is e^1.5, where e is euler's constant (approx. 2.718).  So, the odds of differential expression os about 4.8 to 1
# t statistic is ratio of the logFC to the standard error (where the error has been moderated across all genes...because of Bayesian approach)
  
# Volcano Plots ----
# in topTable function above, set 'number=40000' to capture all genes

# now plot
vplot <- ggplot(myTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#BE684D") +
  annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
  labs(title="Volcano plot",
       subtitle = "INDEX-20 Exposure to Hypobaric-Hypoixa",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
ggsave("C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/Figure11_Volcano_plot_LogversusLogFC_21_April",  width = 10, height = 4, device='tiff', dpi=500)

# Now make the volcano plot above interactive with plotly
ggplotly(vplot)

#decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=2)


# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include="up")

# retrieve expression data for your DEGs ----
head(v.DEGList.filtered.norm$E)
colnames(v.DEGList.filtered.norm$E) <- sampleLabels

#selecting only differentally expressed genes
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes)

#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")
# create interactive tables to display your DEGs ----
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in COmparison 1: High-altitude versus Low-altitude',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

#write your DEGs to a file
"Saving those 282_DEGs, with 0.05 p-value"
write_tsv(diffGenes.df,"C:/Users/mohin/Desktop/Medgenome_data/Analysis_Bilal/Comparisons/LA_Dexa_vs_HA_Dexa(8_7)/DiffGenes_2.txt") #NOTE: this .txt file can be directly used for input into other clustering or network analysis tools (e.g., String, Clust (https://github.com/BaselAbujamous/clust, etc.)

"Step: 5 complete"

#___________________________________________________STEP:6_____________________________________
"this script creates heatmaps from your differentially expressed genes or transcripts
#and selects modules of co-expressed genes based on pearson correlations"

# the essentials ----
library(tidyverse)
library(gplots)
library(RColorBrewer)

# Choose your color pallette ----
myheatcolors <- rev(brewer.pal(name="RdBu", n=11))

# Data----
"diffGenes created in previous step"
# you can make a heatmap out of any datamatrix
# we'll use our 'diffgenes' datamatrix that was produced at the end of the last class in the Step 5 script
# as a reminder, this was produced as follows

# Cluster DEGs ----
"#cluster rows by pearson correlation"
#begin by clustering the genes (rows) in each set of differentially expressed genes
# we use the 'cor' function and the pearson method for finding all pairwise correlations of genes
# '1-cor' converts this to a 0-2 scale for each of these correlations, which can then be used to calculate a distance matrix using 'as.dist'
clustRows <- hclust(as.dist(1-cor(t(diffGenes), method="pearson")), method="complete")


# hierarchical clustering is a type of unsupervised clustering. Related methods include K-means, SOM, etc 
# unsupervised methods are blind to sample/group identity
# in contrast, supervised methods 'train' on a set of labeled data.  
# supervised clustering methods include random forest, and artificial neural networks
"cluster columns by spearman correlation"
#now cluster your samples (columns)
#we may not acutally use this clustering result, but it's good to have just in case
clustColumns <- hclust(as.dist(1-cor(diffGenes, method="spearman")), method="complete")

#Cut the resulting tree and create color vector for clusters.  
#Vary the cut height (h =) to give more or fewer clusters, or use force k= number of clusters
#we'll look at these clusters in more detail later
module.assign <- cutree(clustRows, k=2)


#now assign a color to each module (makes it easy to identify and manipulate)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 



# Produce a static heatmap of DEGs ----
#plot the hclust results as a heatmap
heatmap.2(diffGenes, 
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          RowSideColors=module.color,
          col=myheatcolors, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(10,20))

#choose a cluster(s) of interest by selecting the corresponding number based on the previous graph
modulePick <- 2 #use 'c()' to grab more than one cluster from the heatmap.  e.g., c(1,2)

#now we pull out the genes from this module using a fancy subsetting operation on a named vector
myModule_up <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_up <- hclust(as.dist(1-cor(t(myModule_up), method="pearson")), method="complete") 

# Create heatmap for chosen sub-cluster.
heatmap.2(myModule_up, 
          Rowv=as.dendrogram(hrsub_up), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

modulePick <- 1 
myModule_down <- diffGenes[names(module.assign[module.assign %in% modulePick]),] 
hrsub_down <- hclust(as.dist(1-cor(t(myModule_down), method="pearson")), method="complete") 

heatmap.2(myModule_down, 
          Rowv=as.dendrogram(hrsub_down), 
          Colv=NA, 
          labRow = NA,
          col=myheatcolors, scale="row", 
          density.info="none", trace="none", 
          RowSideColors=module.color[module.assign%in%modulePick], margins=c(8,20))

"Step6: Finished"
"However, uncomplete yet!!"
#__________________________STEP7: Enrirchment Analysis___________________________________
# the essentials ----
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots


# Carry out GO enrichment using gProfiler2 ----
"Tool is gProfiler2"
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=300, sort.by="logFC")

View(myTopHits)
#use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
gost.res_up <- gost(rownames(myModule_up), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_up, interactive = T, capped = T) #set interactive=FALSE to get plot for publications

#For low expressed genes
gost.res_down <- gost(rownames(myModule_down), organism = "hsapiens", correction_method = "fdr")
gostplot(gost.res_down, interactive = T, capped = T) #set interactive=FALSE to get plot for publications

hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 

# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(mydata.df, GeneID, LogFC)
mydata.gsea <- mydata.df.sub$LogFC
names(mydata.gsea) <- as.character(mydata.df.sub$GeneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)
myGSEA.df <- as_tibble(myGSEA.res@result)

#view results as an interactive table
datatable(myGSEA.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in exposure to Hypobaric-Hypoxia',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(3:10), digits=2)

# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res, 
          geneSetID = 47, #can choose multiple signatures to overlay in this plot
          pvalue_table = TRUE, #can set this to FALSE for a cleaner plot
          title = myGSEA.res$Description[47]) #can also turn off this title

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "HYPOXIA",
    NES < 0 ~ "NORMOXIA"))

# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df[1:20,], aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()
