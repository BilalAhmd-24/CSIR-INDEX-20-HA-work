# the essentials ----
library(tidyverse)
library(edgeR)
library(matrixStats)
library(cowplot)

# Examine your data up to this point ----
myTPM <- Txi_gene$abundance
View(myTPM)
capture.output(myTPM, file = "myTPM.csv")
myCounts <- Txi_gene$counts
colSums(myTPM)
colSums(myCounts)

#Loading the  sample labels 
sampleLabels <- targets$sample
sampleLabels


# Generate summary stats for your data ----
# 1st, calculate summary stats for each transcript or gene, and add these to your data matrix
# then use the base R function 'transform' to modify the data matrix (equivalent of Excel's '=')
# then we use the 'rowSds', 'rowMeans' and 'rowMedians' functions from the matrixStats package
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))

# look at what you created
head(myTPM.stats)

# Create your first plot using ggplot2 ----
# produce a scatter plot of the transformed data
ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=23, size=3)+
  theme_half_open()
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
       subtitle="unfiltered, non-normalized data",
       caption="INDEX - CMREF 2020") +
  theme_classic() +
  theme_dark() + 
  theme_bw()

#Creating a DGE list with function DGE list
myDGEList <- DGEList(Txi_gene$counts)
View(myDGEList)

#DEGList objects are a good R data file to consider saving to you working directory
save(myDGEList, file = "myDGEList")

log2.cpm <- cpm(myDGEList, log=TRUE)

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
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
  #coord_flip()

#first, take a look at how many genes or transcripts have no read counts at all
table(rowSums(myDGEList$counts==0)==10)
#181 genes were absent, and have no counts

cpm <- cpm(myDGEList)
# now set some cut-off to get rid of genes/transcripts with low counts
# again using rowSums to tally up the 'TRUE' results of a simple evaluation
# how many genes had more than 1 CPM (TRUE) in at least 3 samples
#How many rows/genes have counts per million greater than one in more than 10 samples
keepers <- rowSums(cpm>1)>=10 #user defined

#now use base R's simple subsetting method to filter your DGEList based on the logical produced above
myDGEList.filtered <- myDGEList[keepers,]

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
  theme_bw()#+
  #coord_flip()

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
  theme_bw()#+
  #coord_flip()

plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 12)
print("Step 2 complete!")
