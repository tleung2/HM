  ## SET WORKING DIRECTORY

  ## Install Phyloseq if needed (for R ver 4.0+)
if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

  ## Alternative for installing phyloseq (from Phyloseq website)
  ## I always have problem going this route
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

  ## ------- LOAD LIBRARIES -----------
  #install.packages(c("vegan","tidyverse","scales","gridGraphics", "reshape2"))
library(tidyverse) ## Package includes ggplot and dplyr
library(scales)
library(gridGraphics) 
library(reshape2)
library(vegan)
library(phyloseq)


############################################################################################
  ## ------- ASSIGNING VARIABLES FOR PHYLOSEQ-----------
  ## Make sure that the tax and shared files are in the working directory
  ## Copy shared, tax, and map file names with extension to corresponding values
  ## Assign variables for data that will be imported
sharedfile_gg = "HM_gg.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
taxfile_gg = "HM_gg.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"
  ## Repeat for rdp ref database
sharedfile_rdp = "HM_rdp.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.shared"
taxfile_rdp = "HM_rdp.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.0.03.cons.taxonomy"


  ## ------- IMPORT MOTHUR DATA -----------
  ## Combines the shared and taxonomy files
gg_data<-import_mothur(mothur_shared_file = sharedfile_gg,
                           mothur_constaxonomy_file = taxfile_gg)

  ## Repeat for rdp ref database
rdp_data<-import_mothur(mothur_shared_file = sharedfile_rdp,
                       mothur_constaxonomy_file = taxfile_rdp)
  ## ------- EXAMINING MOTHUR DATA-------
  ## check the tax names, which are organized as ranks in columns
head(tax_table(rdp_data)) 
head(tax_table(gg_data)) 

############################################################################################
  ## ------- IMPORT METADATA ---------
  ##  Can also import as excel, does not have to be csv
map<-read.csv("HM_mapfile.csv", stringsAsFactors = FALSE)
head(map) ## check headings

  ## Convert the metadata into phyloseq format
  ## Make sre that rownames must match the 
  ## sample names in your shared and taxonomy file
map2 <- sample_data(map)

  ## Assign ID created by Mothur as rownames (NOT SAMPLE ID!)
  ## This will allow phyloseq to merge Mothur output with mapfile
rownames(map2) <- map2$Mothur_ID 
map2$Depth <- factor(map2$Depth, levels=c('0','1', '2', '3', '4', '5', '5.1', '5.8','6', '6.8','7', '8', '9', '9.5','10'))
  ## Merge mothurdata with sample metadata
  ## Phyloseq will merge both datasets using the Mothur_ID
gg_merge <- merge_phyloseq(gg_data, map2)
rdp_merge <- merge_phyloseq(rdp_data, map2)


colnames(tax_table(gg_merge))  ## Check to see how many tax ranks there are
colnames(tax_table(gg_merge)) <- c("Kingdom", "Phylum", "Class", 
                                     "Order", "Family", "Genus", "Species") ## assigns names

  ## Repeat for rdp ref database
colnames(tax_table(rdp_merge))  ## Check to see how many tax ranks there are
colnames(tax_table(rdp_merge)) <- c("Kingdom", "Phylum", "Class", 
                                   "Order", "Family", "Genus") ## assigns names

############################################################################################
  ## ------- PRUNING THE DATASET -------------
  ## filter out samples we don't want to include in our analysis 
  ## such as OTU with 0 counts
  ## Note that some OTU have 0.0001 abundance
gg.prune.data <- gg_merge %>%
  prune_taxa(taxa_sums(.) > 1e-5, .) ## Value can be changed accordingly
gg.prune.data@tax_table  ## lists the taxonomy


  ## Repeat for rdp ref database
rdp.prune.data <- rdp_merge %>%
  prune_taxa(taxa_sums(.) > 1e-5, .) ## Value can be changed accordingly
rdp.prune.data@tax_table  ## lists the taxonomy


############################################################################################
  ## ------- CHECKING COMMUNITY IN DATA --------------------
  ## Use %>% to pass from left to right operator (chain functions together)
  ## epa.data %>% tax_glom() = tax_glom(epa.data)
  ## Change the rank to something you want
  ## --- What organisms are present? ---
  ## First, we need to reshape the data in a table format to see what 
  ## organisms were present
gg.species <- gg.prune.data %>%
  tax_glom(taxrank = "Species") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

   ## Write a csv file of the gg table
write.table(gg.species, "HM_gg_species.csv", row.names = FALSE, sep = ";")


  ## Repeat for rdp ref database
rdp.species <- rdp.prune.data %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

  ## Write a csv file of the rdp table
write.table(rdp.species, "HM_rdp_species.csv", row.names = FALSE, sep = ";")


############################################################################################
  ## ------- SUBSETTING IN DATA IN PHYLOSEQ -----------
  ## Use %in% and subset_taxa
  ## Species table shows both non-photosynthetic and photosynthetic eukaryotes
  ## We wanted to include only photosynthetic aquatic eukaryotes and exclude all others
  ## *Note that some reference database use Dinophyta whereas others use Dinoflagellata
gg.cyano.data <- subset_taxa(gg.prune.data, Phylum == "p__Cyanobacteria")

  ## Repeat for rdp ref database
rdp.cyano.data <- subset_taxa(rdp.prune.data, Phylum == "Cyanobacteria")

############################################################################################
  ## ------- LOOK AT DISTRIBUTION OF READ COUNTS  -------------
  ## Make a data frame with a column for the read counts of each sample
sample_sum_df <- data.frame(sum = sample_sums(gg.cyano.data))

sample_sum_df2 <- data.frame(sum = sample_sums(rdp.cyano.data))

  ## Histogram of sample read counts
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())


############################################################################################
  ## -------------------------- CHECKING DIVERSITY ----------------------------------- ##
  #######################################################################################
  
  ## This function estimates a number of alpha-diversity metrics and returns 
  ## a ggplot plotting object
  ## You must use untrimmed, non-normalized count data for meaningful results,
  ## as many of these estimates are highly dependent on the number of 
  ## singletons (single reads that were not stitched). 
  ## In the 'measures()' you can add different types of metrics. 
  ## See the R help for what kinds of Alpha diversity metrics can be used
  ## --- How does diversity vary with depth over time? ---
  ## Plot Depth as x-axis, and identify Month by different shapes
P1 = plot_richness(gg.prune.data, x="Depth", shape="Month", measures=c("Shannon")) +
  geom_point(size=5) + labs(x = "Depth (m)") + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        text = element_text(size = 16),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16), 
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.key=element_rect(fill='white'))
P1

#########################################################################################
  ## ----------------------- RAREFRACTION CURVE: CYANOBACTERIA ---------------------- ##
  ######################################################################################  

  ## The rarecurve function is part of "vegan" package, not phyloseq
  ## Make rarefaction curve
rarecurve(t(otu_table(gg_merge)), step=50, cex=0.5)
  ## rarefy without replacement
ps.rarefied = rarefy_even_depth(gg.cyano.data, rngseed=1, sample.size=0.9*min(sample_sums(gg.cyano.data)), replace=F)
ps.rarefied
rarecurve(t(otu_table(ps.rarefied)), step=50, cex=0.5)
  ## rarefied dataset
plot_bar(ps.rarefied, fill="Order")
#rarefied dataset by month
plot_bar(ps.rarefied, fill="Order") + facet_wrap(Month~., scales="free_x", nrow=1)
ps.phylum = tax_glom(ps.rarefied, taxrank="Rank2", NArm=FALSE)            
ps.phylum
plot_bar(ps.phylum, fill="Rank2") + facet_wrap(~Month, scales= "free_x", nrow=1)


#########################################################################################
  ## --------------------BARPLOTS USING PHYLOSEQ ----------------------------------- ##
  ####################################################################################
p <- plot_bar(gg.prune.data, fill="Phylum") + facet_wrap(~Month, scales= "free_x", ncol=1)
p + geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack") +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        panel.background = element_blank(),
        plot.title = element_text(size = 22),
        panel.grid.major = element_blank(),   ## Change this to element_line if you want lines
        axis.text = element_text(size = 14, color = "black"), 
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 22))



############################################################################################
  ## ------- RESHAPE DATA FOR BARPLOT (RAW COUNTS)---------
  ## Data needs to be reshaped for ggplot
  ## This will produce relative abundance 
  ## We break it down to "Species" tax but this can be changed
gg.cyano.counts <- gg.cyano.data %>%
  tax_glom(taxrank = "Species") %>%                     # agglomerate at Order level
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

  ## Check Order (or any other taxa)
levels(gg.cyano.species$Order)
  ## If the command: levels(epa_order$Order or any tax column) shows NULL 
  ## Need to convert taxonomy columns to factor, R ver. 4.0.2 identifies these as character
  ## SKip this if your R ver identifies taxonomy as factor
gg.cyano.counts[,c(8:15)] <- lapply(gg.cyano.counts[,c(8:15)], as.factor) #convert phyla to factor
levels(gg.cyano.counts$Phylum) # should see levels (taxa), check that haptophyte exists
  ## ------- EXPORT SUBSETTED DATA ---- skip if not needed---------
  ##  Exports the transformed subset dataframe if needed
write.table(photo.counts, "name_file.csv", row.names = FALSE, sep = ";")


############################################################################################
  ## ------- COMMUNITY COMPOSITION : BARPLOTS (RAW COUNTS)------------
  ##  Using stacked bar plots to look at community composition
  ##  Set colors for plotting
  ##  For this project, interested in 5 phyla (subgroups)
  ##  Assign more colors if you want more than 5 subgroups
phylum_colors <- c("#669900", "#CCCFFF", "#CC9933","#663300", "#FFCC66")

  ## Plot stacked bargraph
p2<-ggplot(data = gg.cyano.counts, aes(x = Abundance, y = reorder(Depth, desc(Depth)), fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  labs(x = "Abundance", y = "Depth(m)") +
  #ggtitle("Phylum Composition of Photosynthetic Eukaryotic Community") +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        panel.background = element_blank(),
        plot.title = element_text(size = 22),
        panel.grid.major = element_blank(),   ## Change this to element_line if you want lines
        axis.text = element_text(size = 14, color = "black"), 
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 22)) +
        facet_wrap(~Month, ncol = 1, scales = "free")   ## ncol = # columns
p2
############################################################################################
  ## ------- RESHAPE and TRANSFORM DATA FOR BARPLOT (RELATIVE ABUNDANCE)---------
  ## Organizes, formats, and transforms subsetted dataframe
  ## This will produce relative abundance, using ggplot
  ## We break it down to "Species" tax but this can be changed
gg.cyano.trans <- gg.cyano.data %>%
  tax_glom(taxrank = "Species") %>%                     # agglomerate at Order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  #filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa
  arrange(Phylum)                                      # Sort data frame alphabetically by phylum

  ## Check Order (or any other taxa)
levels(gg.cyano.trans$Order)
  ## If the command: levels(epa_order$Order or any tax column) shows NULL 
  ## Need to convert taxonomy columns to factor, R ver. 4.0.2 identifies these as character
  ## SKip this if your R ver identifies taxonomy as factor
gg.cyano.trans[,c(8:15)] <- lapply(gg.cyano.trans[,c(8:15)], as.factor) #convert phyla to factor
levels(gg.cyano.counts$Phylum) # should see levels (taxa), check that haptophyte exists
  ## ------- EXPORT SUBSETTED DATA ---- skip if not needed---------
  ##  Exports the transformed subset dataframe if needed
write.table(photo.trans, "name_file.csv", row.names = FALSE, sep = ";")

############################################################################################
  ## ------- COMMUNITY COMPOSITION : BARPLOTS (RELATIVE ABUNDANCE)------------
  ##  Using stacked bar plots to look at community composition
  ##  Set colors for plotting
  ##  For this project, interested in 5 phyla (subgroups)
  ##  Assign more colors if you want more than 5 subgroups
phylum_colors <- c("#669900", "#CCCFFF", "#CC9933","#663300", "#FFCC66")

  ## Plot stacked bargraph
p3<-ggplot(data = gg.cyano.trans, aes(x = Depth, y = Abundance, fill = Phylum)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  guides(fill = guide_legend(reverse = TRUE, keywidth = 1, keyheight = 1)) +
  labs(x = "Depth(m)", y = "Relative Abundance") +
  #ggtitle("Phylum Composition of Photosynthetic Eukaryotic Community") +
  theme(axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14), 
        panel.background = element_blank(),
        plot.title = element_text(size = 22),
        panel.grid.major = element_line(color = "black"),
        axis.text = element_text(size = 14, color = "black"), 
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color="black", fill=NA),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 22)) +
        facet_wrap(~Month, ncol = 3, scales = "free")
p3

#############################################################################################
  ## ------- PREPARING DATA FOR ORDINATION-------
  ## ** Note that ordination in phyloseq uses only phyloseq objects
  ## ** Therefore, all dataframe should have "phyloseq"
  ## R sees "Month" in the phyloseq data as characters not factor or levels
  ## Convernt month to factor and assign levels using sample_data()
sample_data(photo.data)$Month <- factor(
  sample_data(photo.data)$Month, 
  levels = c("February", "May", "September")
)

  ## Converts Year to factor type and assign levels
sample_data(photo.data)$Year <- factor(
  sample_data(photo.data)$Year,
  levels = c("2017", "2018"))

  ## Converts Depth to factor type and assign levels
sample_data(photo.data)$Depth <- factor(
  sample_data(photo.data)$Depth,
  levels = c("1", "2", "3", "5", "7", "11", "12",
             "13","13Dp","14","15"))

  ## Normalize number of reads in each sample using median sequencing depth.
  ## Nick's tutorial first prunes the data and then rarefies it 
  ## but literature advises using rarefy with caution because data is lost
total = median(sample_sums(photo.data))
standf = function(x, t=total) round(t * (x / sum(x)))
  ## assigned as photo.data2 to keep orignial intact
photo.data2 = transform_sample_counts(photo.data, standf)  


##############################################################################################
  ## ------- PLOT ALL ORDINATIONS ----------------
  ## Setting pipline for plotting all ordinations
  ## Need to load plyr and ape package for this
library(plyr)
library(ape)
  ## First need to make a raondon phylum tree
random_tree = rtree(ntaxa(photo.data), rooted=TRUE, tip.label=taxa_names(photo.data))
plot(random_tree)
  ## Then add tree to photo.data
photo.data2 = merge_phyloseq(photo.data, random_tree)
photo.data2
  ## Set R pipeline for distance type and orindation methods
dist = "bray"
ord_meths = c("DCA", "CCA", "RDA", "DPCoA", "NMDS", "MDS", "PCoA")

  ## Loops through different method parameter options to the plot_ordination function, 
  ## and stores the plot results in a list, and then plot these results in a combined 
  ## graphic using ggplot2. 
  ## You will not see any plots, it is a function and takes a hot moment to run
  ## ** Note: Can change "taxa" to "samples" to look at sample, see photo.data2 dataframe 
  ## for options, however color must variable must correspond with corresponding selection
plist = llply(as.list(ord_meths), function(i, photo.data2, dist){
  ordi = ordinate(photo.data2, method=i, distance=dist)
  plot_ordination(photo.data2, ordi, "taxa", color="Phylum")
}, photo.data2, dist)
  ## Assigns ordination method types to plist
names(plist) <- ord_meths
  ## extract the data from each of those individual plots from plist,
  ## and put it back together in one big data.frame 
pdataframe = ldply(plist, function(x){
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "method"

  ## Plot the faceted scatterplot
  ## Change the color and shape according to "taxa" or "samples" in the plist
p4 = ggplot(pdataframe, aes(Axis_1, Axis_2, color=Phylum, shape=Phylum)) + 
  geom_point(size=4) + 
  facet_wrap(~method, scales="free") +
  scale_fill_brewer(type="qual", palette="Set1") +
  scale_colour_brewer(type="qual", palette="Set1") +
  theme_classic() + theme(axis.text.y.left = element_text(size=12, color = "black"), 
                          axis.text.x.bottom = element_text(size=12, color = "black"),
                          legend.text = element_text(size = 12),
                          legend.title = element_text(size=12),
                          text = element_text(size = 12),
                          axis.title.x = element_text(size=15, face="bold"),
                          axis.title.y = element_text(size=15, face="bold"))
p4

###########################################################################################
  ## ------- UNCONSTRAINED ORDINATIONS--------
  ## Reminder: the dataset needs to be a phyloseq object
  ## the Environment panel will contain datasets that have the word "phyloseq"
  ## Run the PCoA --Note:can change method to NMDS,CCA,RDA,DCA,CAP, see R documentaion
photo_pcoa <- ordinate(
  physeq = photo.data, 
  method = "NMDS", 
  distance = "bray"
)

  ## ------- Plot the PCoA ------------
  ## Note: "split" type result in a combined plot with both taxa and samples
  ## supported options are "samples", "sites", "taxa", "biplot", "split", "scree"
p5<-plot_ordination(
  physeq = photo.data,
  ordination = photo_pcoa,
  type = "Split", # change this to taxa to see just taxa graph or sample graph
  shape = "Depth", #v ariables are with repect to type, eg. "taxa" type has phylum while "Sample" type has others (see object in environment)
  color = "Phylum", # same comment as above
  title = "NMDS of Photosynthetic Eukaryote Communities") + 
  geom_point(size =3) +
  geom_text(mapping = aes(label = Month), size = 3, vjust = 1.5) + ## Labels sample
  scale_shape_manual(values=c(13, 16, 17, 18, 19, 0, 2, 4, 3, 8, 11, 7)) + ## Assign symbols for shape
  #scale_color_manual(values = c("#A65628", "red", "#FFAE19", "#4DAF4A", "#1919FF",
                                #"darkorchid3", "magenta", "#FF9900", "#00CC33", "#FF3366",
                                #"#CC00CC", "#6633CC")) + ## Assign colors
  theme_classic() + theme(axis.text.y.left = element_text(size=12, color = "black"), 
                          axis.text.x.bottom = element_text(size=12, color = "black"),
                          legend.text = element_text(size = 12),
                          legend.title = element_text(size=12),
                          text = element_text(size = 12),
                          axis.title.x = element_text(size=15, face="bold"),
                          axis.title.y = element_text(size=15, face="bold"))
p5

  ## Facet if needed but keep in mind that "sample" and "Taxa" are separate column in the 
  ## phyloseq object so it may look confusing
  ## Reminder: Facet by will depend on "type" chosen
p6<-p4+facet_wrap(~Phylum, scales = "free") ## can removed scales to keep axes the same
p6


#############################################################################################
  ## ------- SUBSET BY MONTH AND ORDINATE PCoA --------
  ## Regardless of the ordinate method, it seems that there is a temporal effect
  ## September samples cluster together and more so for May samples
  ## 1) Subset photo.data into months : May and September (February ony had 2 samples)
photo.may <- subset_samples(photo.data, Month == "May")

  ## Run the PCoA --Note:can change method to NMDS,CCA,RDA,DCA,CAP, see R documentaion
  ## We have unsufficient data but running it anyway to see if there is anything worthwhile
may_nmds <- ordinate(
  physeq = photo.may, 
  method = "NMDS", 
  distance = "bray"
)

  ## ------- Plot the PCoA --------------
  ## Note: "split" type result in a combined plot with both taxa and samples
  ## supported options are "samples", "sites", "taxa", "biplot", "split", "scree"
p7<-plot_ordination(
  physeq = photo.may,
  ordination = may_nmds,
  type = "Split", # change this to taxa to see just taxa graph or sample graph
  shape = "Depth", #v ariables are with repect to type, eg. "taxa" type has phylum while "Sample" type has others (see object in environment)
  color = "Phylum", # same comment as above
  title = "NMDS of Photosynthetic Eukaryote Communities in May") + 
  geom_point(size =3) +
  geom_text(mapping = aes(label = Month), size = 3, vjust = 1.5) + ## Labels sample
  scale_shape_manual(values=c(13, 16, 17, 18, 19, 0, 3)) + ## Assign symbols for shape
  #scale_color_manual(values = c("#A65628", "red", "#FFAE19", "#4DAF4A", "#1919FF",
  #"darkorchid3", "magenta", "#FF9900", "#00CC33", "#FF3366",
  #"#CC00CC", "#6633CC")) + ## Assign colors
  theme_classic() + theme(axis.text.y.left = element_text(size=12, color = "black"), 
                          axis.text.x.bottom = element_text(size=12, color = "black"),
                          legend.text = element_text(size = 12),
                          legend.title = element_text(size=12),
                          text = element_text(size = 12),
                          axis.title.x = element_text(size=15, face="bold"),
                          axis.title.y = element_text(size=15, face="bold"))
p7

  ## Facet if needed but keep in mind that "sample" and "Taxa" are separate column in the 
  ## phyloseq object so it may look confusing
  ## Reminder: Facet by will depend on "type" chosen
  ## Can set scales free by scales = "free", "free_x" or "free_y"
p8<-p7+facet_wrap(~Phylum) 
p8


#############################################################################################
  ## ------- RESHAPE PHYLOSEQ DATA FOR VEGAN  ------------
  ## CCA compares variables within datasets and between 2 matrices/datasets
  ## Convert the phyloseq data into a dataframe if you haven't done so
  ## For this pipline, we have converted it in the bar plot sections
  ## You can use raw or relative counts, change the data accordingly
  ## Turn on the vegan package to do this
library(vegan)
  ## ------- Create a new dataframe and reshaping it ------------
  ## It is less messy to create a new dataframe then reshape it because there are too many 
  ## variables (too many taxa) and the table will be very wide
  ## For this study, we are interested in eukaryotic phyla
photo.trans2 <- photo.data %>%
  tax_glom(taxrank = "Phylum") %>%                     # agglomerate at Phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>%
  psmelt()                                          # Melt to long format
  
  ## To simplify things, include only columns of interest
  ## For this study, we are removing columns Doman and Kingdom
photo.trans3 <- subset(photo.trans2, select = -c(8,9))

  ## Next, pivot the table so that levels in Phyla becomes it's own column
  ## Check the output in the environment, you should seea dataframe
  ##**Note that the number of observations will increase if agglomerating at Genus or Species level
photo.trans4 <- photo.trans3 %>%
  pivot_wider(names_from = Phylum, values_from = Abundance, id_cols = Mothur_ID)

  ## Add the metadata (in csv format) to the new dataframe
  ## Now we have a dataframe that consists of both environment and response variables
photo.trans4 <- merge(photo.trans4, map, by = "Mothur_ID")


##############################################################################################
  ## ------- TESTING DATASET FOR CCA -----------
  ## CCA has 2 assumptions:
  ## 1) Some linearity between independent and response variables
  ## 2) Response variables are unimodal (even one!)
  ## Will take a while of the dataframe is large
  ## Test your data by using GGally package
library(GGally)
ggpairs(photo.trans4[,c(2:6)])  ## Select columns 2 to 6 for this study
  ## ------- RUN THE CCA USING VEGAN ------------
  ## Define the environmental and response variables
  ## code explanation: cca(species data,environmental data)
ccamodel1<-cca(photo.trans4[,c(2:6)],photo.trans4[,c(9:10)])
ccamodel1 ## Axis 1 and 2 explain 53% of the variation in OTU
ccamodel1$CCA  ## This gives all the data that will be used for plotting


##############################################################################################
  ## ------- PREPARING FOR PLOT WITH GGPLOT ----------------
  ## ------- Installing ggvegan package-----------
  ## Need ggvegan package to do this
  ## This worked for R version 4.0.2
install.packages("remotes")
remotes::install_github("gavinsimpson/ggvegan", force = TRUE)
library(ggvegan)
  ## ------- Convert the CCA results into a readable format-----------
  ## Convert the CCA results into a readable format
  ## The fortify command will recognize the CCA
cca.res<-fortify(ccamodel1, axes = 1:2)  ## CCA results in a dataframe format, see environment


  ## ------- Preparing for plotting -------------
  ## subset sites/samples
site.data<-subset(cca.res, Score == "sites") 
species.data<-subset(cca.res, Score == "species")

  ## Add Depth, sample ID, Month to species/site subset data:
  ## binds_col() is from dplyr package (within tidyverse)
site.data<-bind_cols(site.data, map[,c(1:5)])  ## column 1 to 5

  ## Scale environmental arrows to the plot 
  ## subset environmental variables -- these will plot as arrows later
arrows<-subset(cca.res, Score == "biplot")

  ## multiply the environmental variables (arrows) to scale it to the plot
scores<-c('CCA1','CCA2')
mul<-ggvegan:::arrowMul(arrows[,scores],
                        subset(cca.res, select = scores, Score == "sites"))

  ## Scales the biplot arrows
arrows[,scores] <-arrows[,scores] * mul


  ## ------- Plot CCA using ggplot---------

p9<-ggplot() +
  geom_point(site.data, mapping = aes(x = CCA1, y = CCA2, color = Month,
                                      shape = Depth), size = 3) + #leave out color not wanted
  geom_segment(arrows, mapping = aes(x = 0, xend = CCA1, y = 0, yend = CCA2),
               arrow = arrow(length = unit(0.03, "npc")), # unit = arrow end
               color = "blue", size = 1.5) +
  geom_text(arrows, mapping = aes(label = Label, x = CCA1*1.1, y = CCA2*1.1),
            size = 5) +
  geom_text(species.data, mapping = aes(label = Label, x = CCA1*1.1, y = CCA2*1.1),
            size = 5) +
  scale_shape_manual(values=c(13, 16, 17, 18, 19, 0, 2, 4, 3, 8, 11, 7)) +
  coord_fixed() +
  theme(legend.background = element_rect(fill="white", size=0.3, 
                                         linetype="solid", colour="black"),
        panel.background = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=16), axis.title = element_text(size = 16),
        legend.position ="right", legend.key = element_rect(fill = "white"),
        legend.title = element_text(size = 16), legend.text = element_text(size = 16))
p9  

