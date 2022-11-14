
### METAGENOMICS
##
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", version = "3.10")

library(dada2); packageVersion("dada2") ## Dada2 v.1.18.0
library(ShortRead); packageVersion("ShortRead") ## ShortRead 1.48.0
library(phyloseq); packageVersion("phyloseq") ## phyloseq 1.34.0
library(gridExtra)
require(digest)

### METAGENOMIC ANALYSIS 
### ALBA BURILLO 
### SAMPLES: VILLAVICIOSA 16S (4 SAMPLES)


#### Set Working directory to the root directory where you have your data code
code_path<-("C:/Users/Albaburillo/Desktop/fastq_16S")
code_path
setwd(code_path)

#### Download RawData from

#### Download Training Data from
RawDataPath <- paste(code_path,"/fastq_16S/",sep="")
TrainingPath<-paste0(code_path,"/Training/")


# Sort ensures forward/reverse reads are in same order
fnFs <- sort(list.files(pattern="_R1_001.fastq.gz"))
fnRs <- sort(list.files(pattern="_R2_001.fastq.gz"))

# Infer sample names from file names
sample.names <- sapply(strsplit(fnFs, "_S"), `[`, 1)

# #metadata <-sapply(json_data, `[[`, "metadata")
# metadata<-read.table(paste(data_path,"Metadata/metadata.csv",sep=""))

##DBModule$insertLogs(sprintf("Processing %d samples",length(sample.names)),job_id,1)
fnFs <- file.path(fnFs)
fnRs <- file.path(fnRs)
plotQualityProfile(fnFs[1:4])
plotQualityProfile(fnRs[1:4])
filt_path <- file.path(code_path, "DADA2/filtered") # Place filtered files in filtered/ subdirectory

filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

##filterAndTrim
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,250),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE,verbose=T)

### How many sequences are passing the quality Filters
head(out)

plotQualityProfile(filtFs[1:2])
plotQualityProfile(filtRs[1:2])


### From filtered reads, we iteratively learn an error profile
### Should take less than 5 minutes on macboookpro
errF <- learnErrors(filtFs, multithread=TRUE,nbases=50000000,randomize=T,verbose=T)
errR <- learnErrors(filtRs, multithread=TRUE,nbases=50000000,randomize=T,verbose=T)

###
plotErrors(errF, nominalQ=TRUE)


### Using the error model that we have learnt from the data
### We can infer which reads are a result of true existence or a results from error from a true existing sequence.
### DERREPLICATION
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

##DBModule$insertLogs("Removing sequencing errors",job_id,1)
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]

##DBModule$insertLogs("Merging denoised forward and reverse reads.",job_id,1)
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])

seqtab <- makeSequenceTable(mergers)
dim(seqtab)


##DBModule$insertLogs("Removing chimeras",job_id,1)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim.pooled <- removeBimeraDenovo(seqtab, method="pooled", multithread=TRUE, verbose=TRUE)
### The following takes too long.
#seqtab.nochim.persample <- removeBimeraDenovo(seqtab, method="per-sample", multithread=TRUE, verbose=TRUE)


dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN),rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged",  "nonchim")
rownames(track) <- sample.names
head(track)


library(openxlsx)
write.xlsx(track,"table16S.xlsx")




list.files(paste(code_path,"/Training/",sep=""))

taxa.rdp<- assignTaxonomy(seqtab.nochim, paste(code_path,"/Training/rdp_train_set_16.fa.gz",sep=""), multithread=TRUE)
taxa.silva<-readRDS("taxa.silva.rds")
taxa.silva<- assignTaxonomy(seqtab.nochim, paste(code_path,"/Training/silva_nr_v132_train_set.fa.gz",sep=""), multithread=TRUE)


### Assign species if possible
### Since exact sequences are assigned we can try to assign these sequences to species level in some cases (specific genus, not generalizable) ### Assignment to species highly depends on the database of choice.
### RDP
taxa.rdp <- addSpecies(taxa.silva, paste0(code_path,"/Training/rdp_species_assignment_16.fa.gz",""),verbose=T)
taxa.silva <- addSpecies(taxa.silva, paste0(code_path,"/Training/silva_species_assignment_v132.fa.gz",""),verbose=T)

seqs <- getSequences(seqtab.nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA,verbose=FALSE)


phangAlign <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phangAlign)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
detach("package:phangorn", unload=TRUE)



### DFForSeq<-colnames(seqtab.nochim)
### DADA_seq.fas <- DFForSeq
### for(seq in DADA_seq.fas){
   ### write(paste(">",digest(seq),"\n",seq,sep=""),file="DADA2/DADA_seq.fas",append=T)}

#### Taxa rownames are now the sequence of the ASV itself. This is not manageable downstream
#### We transform sequences to shorter strings using md5-algorithm so we can keep track
#rownames(taxa.silva)<-sapply(rownames(taxa.silva),digest,algo="md5")
#colnames(seqtab.nochim)<-sapply(colnames(seqtab.nochim),digest,algo="md5")

### We have the sequences for each ASV in DADA_seq.fas
### We can use it in any othter downstream analysis, but right now we want
### to create a phylogenetic tree to attach it to our data. This tree will allow
### calculation of different classes of phylogneetic-based distances (Unifrac, Faith's index,...)
### We build the alignment with mafft, mafft should be in your path (conda install -c bioconda mafft)
### Need to know where mafft executable is (typically in ~/miniconda3/bin or similar)

#system("~/Users/Albaburillo/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/albaburillo/anaconda3/envs/qiime2-2020.11/bin/mafft --auto --thread -1 DADA2/DADA_seq.fas > DADA2/DADA_seq.afa")

### Based on the alignmnent (afa fill  calculate the tree
### We use FastTree for that, using a generalized time-reversible model
### need to know where fastree executable is (conda install -c bioconda fasttree

#system("~/Users/Albaburillo/AppData/Local/Packages/CanonicalGroupLimited.Ubuntu20.04onWindows_79rhkp1fndgsc/LocalState/rootfs/home/albaburillo/anaconda3/envs/qiime2-2020.11/bin/fasttree -nt -gtr DADA2/DADA_seq.afa > DADA2/DADA_seq.nwk")

### Now we read the phylogenetic tree into a phylogenetic object (ape package)
# library(ape)
#treeSilva <- read_tree("DADA2/DADA_seq.nwk")
#treeSilva <-read_tree(system.file("DADA2/DADA_seq.fas", package="phyloseq"))

### How did we do for our mock positive control samples?
#unqs.mock <- seqtab.nochim["Sample-Cpos1",]
#unqs.mock <- seqtab.nochim["Sample-Cpos2",]
#unqs.mock <- seqtab.nochim["Sample-Cpos3-1rep",
#unqs.mock <- sort(unqs.mock[unqs.mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
#cat("DADA2 inferred", length(unqs.mock), "sample sequences present in the Mock community.\n")
#rownames(unqs.mock)


#### Read metadata
metadata<-read.csv(paste0(code_path,"/metadata.csv"))
rownames(metadata)<-metadata$SampleID

#### Create PhyloSeq Object
#### File named metadata.csv must exist and have SampleID variable for sample name
ps_silva<-phyloseq(otu_table(seqtab,taxa_are_rows=FALSE),tax_table(taxa.silva),phy_tree(fitGTR$tree))


### Clean phyloseq object from mock controls
sample_data(ps_silva)<-metadata

### Clean workspace
###rm(DFForSeq,filtFs,filtRs,fnRs,fnFs,seq,track,mergers,seqtab.nochim.pooled)
# #asiggn metadata
# ##DBModule$insertLogs("Assign metadata if exists",job_id,1)
# if(!any(sapply(metadata,function(list) length(list)==0))) {
#   sampledata<-metModule$obtainSampleMet(metadata,sample.names)
#   otu_samples<-rownames(otu_table(ps_silva))
#   meta_samples<-rownames(sampledata)
#   if (all(length(otu_samples)==length(meta_samples)) && all(otu_samples==meta_samples)) {
#     sample_data(ps_silva)<-sampledata
#     sample_names(ps_silva)
#   }
# }
# ps_silva

### We keep the session
#save.image(file=paste0(code_path,"/DADA2/DADA2_Rsession.RData"))

rank_names(ps_silva)


table(tax_table(ps_silva)[, "Phylum"], exclude = NULL)


# Compute prevalence of each feature, store as data.frame
prevdf <- apply(X = otu_table(ps_silva),
                MARGIN = ifelse(taxa_are_rows(ps_silva), yes = 1, no = 2),
                FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf <- data.frame(Prevalence = prevdf,
                     TotalAbundance = taxa_sums(ps_silva),
                     tax_table(ps_silva))

#Then compute the total and average prevalences of each feature:
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})


# Define phyla to filter
filterPhyla <- c("Fusobacteria","WPS-2")

# Filter entries with unidentified Phylum.
ps <- subset_taxa(ps_silva, !Phylum %in% filterPhyla)
ps


samples.out <- rownames(seqtab.nochim) 
T1 <- samples.out[1]
T2 <- samples.out[2]
T3 <-samples.out[3]
T4 <-samples.out[4]


prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(ps, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +  geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")


#Plot the species richness

plot_ordination(ps, ord.nmds.bray, title="ORDINATION PLOT")
plot_richness(ps, measures=c("Shannon", "Simpson")) + theme_bw()
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")

# Subset to the remaining phyla
prevdf1 <- subset(prevdf, Phylum %in% get_taxa_unique(ps_silva, "Phylum"))

# Plot prevalenced vs abundances. What do we see?
require(ggplot2)
library(ggplot2)
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(ps_silva),color=Phylum)) +
  # Include a guess for parameter
  geom_point(size = 2, alpha = 0.7) +
  geom_hline(yintercept = 0.05, alpha = 0.5, linetype = 2) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")


## Let's keep frequent ASVs
# Define prevalence threshold as 5% of total samples
prevalenceThreshold <- 0.05 * nsamples(ps_silva)

# Execute prevalence filter, using `prune_taxa()` function
keepTaxa <- rownames(prevdf1)[(prevdf1$Prevalence >= prevalenceThreshold)]
ps1 <- prune_taxa(keepTaxa, ps_silva)


### Diversity impact on phylogenetic tree
#Despite having many ASVs we have very little number of bacterial genus
# How many genera would be present after filtering?
length(get_taxa_unique(ps1, taxonomic.rank = "Genus"))


#######   We will use the following packages: phyloseq, vegan, ggplot2 and DeSeq2
suppressPackageStartupMessages(library(phyloseq))### Bioconductors
suppressPackageStartupMessages(library(vegan))### Bioconductors
suppressPackageStartupMessages(library(ggplot2))### CRAN
suppressPackageStartupMessages(library(DESeq2)) ### Bioconductors
suppressPackageStartupMessages(require(doBy))### CRAN
suppressPackageStartupMessages(require(RColorBrewer))### CRAN
#######   gdata packages allows for microsofft office filetype read-in
suppressPackageStartupMessages(library(gdata))### CRAN
packageVersion("phyloseq")
packageVersion("vegan")
packageVersion("ggplot2")
packageVersion("DESeq2")
packageVersion("RColorBrewer")
packageVersion("doBy")


### Let's plot the number of counts for the remaining samples
barplot(colSums(otu_table(ps_silva)),las=2,cex.names=0.6,main="#Counts/Sample")


#### Start Alpha - Diversity Analysis
### We will characterize alpha-diversity indices using a rarefied subset of 5000 counts
### in order to balance sampling between samples
### First we try to apply a prevalence filter
x.2.0 <- ps_silva
wh0=genefilter_sample(x.2.0,filterfun_sample(function(x) x>1), A=0.01*nsamples(x.2.0))
x.2.0<-prune_taxa(wh0,x.2.0)

# Compute prevalence of each feature, store as data.frame
prevdf = apply(X = otu_table(x.2.0),
               MARGIN = ifelse(taxa_are_rows(x.2.0), yes = 1, no = 2),
               FUN = function(x){sum(x > 0)})
# Add taxonomy and total read counts to this data.frame
prevdf = data.frame(Prevalence = prevdf,
                    TotalAbundance = taxa_sums(x.2.0),
                    tax_table(x.2.0))
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
prevdf1 = subset(prevdf, Phylum %in% get_taxa_unique(x.2.0, "Phylum"))
ggplot(prevdf1, aes(TotalAbundance, Prevalence / nsamples(x.2.0),color=Phylum)) +
  # Include a guess for parameter
  geom_hline(yintercept = 0.01, alpha = 0.5, linetype = 2) + geom_point(size = 2, alpha = 0.7) +
  scale_x_log10() +  xlab("Total Abundance") + ylab("Prevalence [Frac. Samples]") +
  facet_wrap(~Phylum) + theme(legend.position="none")

### Let's go back
x.2.0 <- ps_silva

##########Generate Diversity Plots
p<-plot_richness(ps_silva,measures=c("Shannon","Simpson","InvSimpson"))
p+geom_boxplot()

#####################################################################################################################
######################################### Taxonomical analysis ######################################################
### Let's make an initial taxonomy description of the data
### Cumulative Stacked Barplots at Phylum, Genus and Species
### 3.0
x.3.0<-ps_silva
### We apply a more stringent filter, basically we are not interested in rare OTUs but on
### general main trends of taxonomical composition
wh0=genefilter_sample(x.3.0,filterfun_sample(function(x) x>5), A=0.1*nsamples(x.3.0))
x.3.0<-prune_taxa(wh0,x.3.0)
### tax_glom function agglomerates/collapses all OTU belonging to the same taxonomical level
x.3.0.phylum<-tax_glom(x.3.0,taxrank="Phylum")
x.3.0.genus<-tax_glom(x.3.0,taxrank="Genus")
x.3.0.genus<-subset_taxa(x.3.0.genus,Genus != "unclassified")
x.3.0.genus<-subset_taxa(x.3.0.genus,Genus != "Incertae_Sedis")
### Convert to data.frame for easier manipulation
ps.melt.x.3.0.phylum<-psmelt(x.3.0.phylum)
ps.melt.x.3.0.genus<-psmelt(x.3.0.genus)

### Let's look at the phylum data.frame
summary(ps.melt.x.3.0.phylum)

### We will calculate relative abundances as proportions for further analysis
ps.melt.x.3.0.phylum$AbundanceProportion <- ave(ps.melt.x.3.0.phylum$Abundance,list(ps.melt.x.3.0.phylum[,"Sample"]), FUN=function(L) L/sum(L))
ps.melt.x.3.0.genus$AbundanceProportion <- ave(ps.melt.x.3.0.genus$Abundance,list(ps.melt.x.3.0.genus[,"Sample"]), FUN=function(L) L/sum(L))

reverse.levels <- function(x) {
  if(is.factor(x)) {
    x <- factor(as.character(x), levels=rev(levels(x)), ordered=TRUE)
  } else if(is.data.frame(x)) {
    for(i in seq_along(x)) {
      if(is.factor(x[,i])) {
        x[,i] <- factor(as.character(x[,i]), levels=rev(levels(x[,i])), ordered=TRUE)
      } else {
        warning(paste0('Column ', i, ' is not a factor.'))
      }
    }
  } else {
    stop(paste0('Unsupported format: ', class(x)))
  }
  return(x)
}

### Barplots at phylum. We use decreasing phylum abundances as X order
### and the phylum abundance of first sample as y order
### This careful use of X/Y order will already reveal patterns
### Define Level Order for X axis (SampleID)
my.levels=orderBy(~-AbundanceProportion,data=ps.melt.x.3.0.phylum[ps.melt.x.3.0.phylum$Phylum=="Fusobacteria",])$Sample
ps.melt.x.3.0.phylum$Sample<-factor(ps.melt.x.3.0.phylum$Sample,
                                      levels=my.levels,ordered=T)

### Define Level Order for Y axis (phylum,genus,species)
my.levels<-orderBy(~-AbundanceProportion,data=ps.melt.x.3.0.phylum[ps.melt.x.3.0.phylum$Sample== as.character(my.levels[1]),])$Phylum
ps.melt.x.3.0.phylum$Phylum<-factor(ps.melt.x.3.0.phylum$Phylum,levels=my.levels,ordered=T)

### We also take care of colors
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorCount=length(unique(ps.melt.x.3.0.phylum$Phylum))
getPalette=colorRampPalette(brewer.pal(12,"Set3"))
p<-ggplot(ps.melt.x.3.0.phylum,aes(x=Sample,y=AbundanceProportion,order=ps.melt.x.3.0.phylum$Phylum))
p+geom_bar(stat="identity",aes(fill=Phylum))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(title="Phylum Level")+
  scale_fill_manual(values=getPalette(colorCount))
### Let's save the plot
#ggsave("x.3.0.phylum.barplot.ordered.pdf")

#Barplots at genus
# Define Level Order for X axis (SampleID)
my.levels=orderBy(~-AbundanceProportion,data=ps.melt.x.3.0.genus[ps.melt.x.3.0.genus$Genus=="Alistipes",])$Sample
ps.melt.x.3.0.genus$Sample<-factor(ps.melt.x.3.0.genus$Sample,
                                     levels=my.levels,ordered=T)
#Define Level Order for Y axis (genus,genus,species)
my.levels<-orderBy(~-AbundanceProportion,data=ps.melt.x.3.0.genus[ps.melt.x.3.0.genus$Sample == as.character(my.levels[1]),])$Genus
ps.melt.x.3.0.genus$Genus<-factor(ps.melt.x.3.0.genus$Genus,levels=my.levels,ordered=T)

ps.melt.x.3.0.genus<-subset(ps.melt.x.3.0.genus,ps.melt.x.3.0.genus$AbundanceProportion>0.02)
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colorCount=length(unique(ps.melt.x.3.0.genus$Genus))
getPalette=colorRampPalette(brewer.pal(12,"Set3"))
p<-ggplot(ps.melt.x.3.0.genus,aes(x=Sample,y=AbundanceProportion,order=Genus))
p+geom_bar(stat="identity",aes(fill=Genus))+ theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(title="genus Level")+
  scale_fill_manual(values=getPalette(colorCount))+guides(fill=guide_legend(ncol=2))
#ggsave("x.3.0.genus.barplot.ordered.pdf")


### We can assess Differential Abundance, related to any variable in our metadata
### 3.0
x.3.0<-ps_silva
wh0=genefilter_sample(x.3.0,filterfun_sample(function(x) x>1), A=0.1*nsamples(x.3.0))
x.3.0<-prune_taxa(wh0,x.3.0)
### tax_glom function agglomerates/collapses all OTU belonging to the same taxonomical level
x.3.0.phylum<-tax_glom(x.3.0,taxrank="Phylum")
x.3.0.genus<-tax_glom(x.3.0,taxrank="Genus")

ps.melt.x.3.0.phylum<-psmelt(x.3.0.phylum)
ps.melt.x.3.0.genus<-psmelt(x.3.0.genus)

### Let's transform our data in different ways
ps.melt.x.3.0.phylum$Abundance<-floor(ps.melt.x.3.0.phylum$Abundance)
### Add a pseudocount to allow for log transformation
ps.melt.x.3.0.phylum$Abundance<-ps.melt.x.3.0.phylum$Abundance+1
ps.melt.x.3.0.phylum$AbundanceProportion <- ave(ps.melt.x.3.0.phylum$Abundance,list(ps.melt.x.3.0.phylum[,"Sample"]), FUN=function(L) L/sum(L))
ps.melt.x.3.0.phylum$AbundanceProportionLog <- log10(ps.melt.x.3.0.phylum$AbundanceProportion)
ps.melt.x.3.0.phylum$AbundanceNorm<-floor(10000*ps.melt.x.3.0.phylum$AbundanceProportion)

ggplot(ps.melt.x.3.0.phylum,aes(x=Phylum,y=AbundanceProportionLog,fill=Phylum))+geom_boxplot()+facet_wrap(~Sample)+ theme(axis.text.x = element_text(angle = 90))

ps.melt.x.3.0.genus$Abundance<-floor(ps.melt.x.3.0.genus$Abundance)
ps.melt.x.3.0.genus$Abundance<-ps.melt.x.3.0.genus$Abundance+1
ps.melt.x.3.0.genus$AbundanceProportion <- ave(ps.melt.x.3.0.genus$Abundance,list(ps.melt.x.3.0.genus[,"Sample"]), FUN=function(L) L/sum(L))
ps.melt.x.3.0.genus$AbundanceProportionLog <- log10(ps.melt.x.3.0.genus$AbundanceProportion)
ps.melt.x.3.0.genus$AbundanceNorm<-floor(10000*ps.melt.x.3.0.genus$AbundanceProportion)

### To plot at Genus level we first need to select the most abundant genus
### plotting all of them is a terrible mess
levels(ps.melt.x.3.0.genus$Genus)
meanByGenus<-aggregate(ps.melt.x.3.0.genus$AbundanceProportion,list(ps.melt.x.3.0.genus$Genus),mean)
colnames(meanByGenus)<-c("Genus","meanAbundance")
meanByGenus<-meanByGenus[with(meanByGenus, order(-meanAbundance, Genus)), ]


#### Plot relative logProportions of most 10 abundant genus
ggplot(ps.melt.x.3.0.genus[ps.melt.x.3.0.genus$Genus %in% meanByGenus[1:10,1],],aes(x=Genus,y=AbundanceProportionLog,fill=Genus))+geom_boxplot()+facet_wrap(~Sample)+ theme(axis.text.x = element_text(angle = 90))

########################################################################
############ Start Genus comparison barplots #########################
########################################################################
### We have seen that there are differences but we need to see which these differences are
###Barplots for Significant Genus difference
x.4.0<-ps_silva
wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>1), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))
x.4.0<-tax_glom(x.4.0,taxrank="Genus")
psmelt.x.4.0.genus<-psmelt(x.4.0)

###
statisticValue<-function(value=NULL){
  if(value>0.1){
    return (paste("=",round(value,digits=2)))
  }
  else if(value>=0.05){
    return (paste("=",round(value,digits=2)))
  }
  else if(value>=0.01){
    return (paste("<0.05"))
  }
  else if(value>=0.001){
    return (paste("<0.01"))
  }
  else if(value>=0.0001){
    return(paste("<0.001"))
  }
  else{return (paste("<0.0001"))}
}
boxplotNumericByGroup<-function(mydata,category,variable,nbvariable,test,Rank=NULL){
  # if(is.null(Rank)){
  #    title<-paste(as.character(variable), " by ",category)
  #    fileForOutput<-paste(as.character(variable),"by",category,sep="_")
  #  }else{
  Phylum<-unique(as.vector(mydata[,Rank]))
  title<-paste(as.character(variable), " of\n ",Phylum,as.character(Rank)," by ",category)
  fileForOutput<-paste(as.character(variable),"of",Phylum,as.character(Rank),"by",category,sep="_")
  #  }
  numberOfLevels<-length(unique(mydata[,category]))
  colorsPlotRight<-c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
  require(gridExtra)
  if(numberOfLevels==2){
    test<-t.test(mydata[,variable]~mydata[,category])
    testString<-paste("Students t-test. p-value",statisticValue(test$p.value))
  }
  if(numberOfLevels>=3){
    require(MASS)
    glm.nb.model<-glm.nb(mydata[,nbvariable]~mydata[,category],method="glm.fit")
    glm.nb.aov<-aov(glm.nb.model)
    test<-aov(mydata[,variable]~mydata[,category])
    tukey.test<-TukeyHSD(test)
    print(tukey.test$mydata)
    mymatrix<-tukey.test$mydata
    mymatrix<-mymatrix[,c("diff","p adj")]
    mymatrix<-as.matrix(mymatrix)
    #mymatrix<-round(mymatrix,3)
    for(i in 1:nrow(mymatrix)){
      print(mymatrix[i,"p adj"])
      mymatrix[i,"p adj"]<-statisticValue(as.numeric(mymatrix[i,"p adj"]))
    }
    text2.df<-as.table(mymatrix)
    testString<-paste("ANOVA PR(>F)", statisticValue(summary(test)[[1]][["Pr(>F)"]][[1]]))
    testString<-paste(testString,"\n","NegBin ANOVA PR(>F)",statisticValue(summary(glm.nb.aov)[[1]][["Pr(>F)"]][[1]]))
  }
  #  if(test=="lm"){
  #    model<-lm(mydata[,variable]~mydata[,category])
  #  }
  mydata$xaxis<-mydata[,category]
  mydata$yvalue<-mydata[,variable]
  p<-ggplot(mydata,aes(x=xaxis,y=yvalue,fill=as.factor(xaxis)))+
    geom_boxplot()+geom_jitter(color="DarkRed")+
    #ggtitle(title)+
    xlab(category)+ylab(variable)+
    ylim(min(mydata$yvalue-1),1)+
    scale_fill_manual(values=colorsPlotRight[1:numberOfLevels])+
    theme(legend.position=c(1,1),legend.justification=c(1,1))+
    #annotate("text",x=numberOfLevels/2.5,y=0.5,label=testString,size=3)+
    annotate("text",x=1,y=0,label=testString,size=3)+
    annotate("text",y=0.4,x=1.5,label=title,size=4)
  plotRight<-ggplot(mydata,aes(yvalue,fill=xaxis))+geom_density(alpha=.5)+
    coord_flip()+scale_fill_manual(values=colorsPlotRight[1:numberOfLevels])+
    theme(legend.position="none")+
    xlim(min(mydata$yvalue-1),max(mydata$yvalue+1))
  if(numberOfLevels>=3){
    p<-p+annotation_custom(tableGrob(text2.df), ymin=min(mydata$yvalue)-1, ymax=min(mydata$yvalue), xmax=numberOfLevels/1.2, xmin=numberOfLevels/2)
  }  #p2<-tableGrob(text2.df)
  #grid.arrange(p2,p,main="prova",ncol=2)
  fileForPlot <- paste(fileForOutput,".pdf")
  pdf(fileForPlot,paper="a4r")
  grid.arrange(p,plotRight,nrow=1,ncol=2,widths=c(4,1),heights=c(4))
  #p2<-ggplot(p2)
  #ggsave(filename = fileForPlot,dpi=600,width=11, height=8.5)
  ##dev.off()
  grid.arrange(p,plotRight,nrow=1,ncol=2,widths=c(4,1),heights=c(4))
  #print(p2)
  #return(p2)
}
compareTwoGroups<-function(mydata=NULL,variable=NULL,category1=NULL,category2=NULL,fileForPlot=NULL,minCounts=500,maxAlpha=0.01, design=NULL){
  fileForPlot=paste("NegBin_DiffTest_",as.character(deparse(substitute(mydata))),"_",variable,"_",category1,"vs",category2,".pdf",sep="")
  print(paste("Output file for Plots: ",fileForPlot))
  fileForTable=paste("NegBin_DiffTest_",as.character(deparse(substitute(mydata))),"_",eval(variable),"_",eval(category1),"vs",eval(category2),".txt",sep="")
  print(paste("Output file for Table: ",fileForTable))
  stringForTitle=paste(variable,"/",category1,"vs",category2)
  kostic <- mydata
  LastRank<-rank_names(mydata)[length(rank_names(mydata))][[1]]
  require(phyloseq)
  require(DESeq2)
  kostic <- prune_samples(sample_sums(kostic) > minCounts, kostic)
  diagdds = phyloseq_to_deseq2(kostic, as.formula(design))
  #print(diagdds)
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
  geoMeans = apply(counts(diagdds), 1, gm_mean)
  diagdds = estimateSizeFactors(diagdds, geoMeans = geoMeans)
  #colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c("untreated","treated"))
  #colData(diagdds)$condition<-factor(colData(diagdds)$condition,levels=c(category1,category2))
  diagdds = DESeq(diagdds, fitType="parametric",test="Wald")
  res=results(diagdds,contrast=c(variable,category1,category2))
  #res=results(diagdds)
  print(res)
  res = res[order(res$padj, na.last=NA), ]
  sigtab = res[(res$padj < maxAlpha), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kostic)[rownames(sigtab), ], "matrix"))
  sigtab$OtuID<-rownames(sigtab)
  head(sigtab)
  write.table(sigtab,file=fileForTable,sep="\t")
  #Cleanup for Positive enrichment in csigtabarcinoma
  posigtab=sigtab
  #posigtab = sigtab[sigtab[, "log2FoldChange"] > 0, ]
  #posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", "Phylum", "Class", "Family", "Genus")]
  posigtab = posigtab[, c("baseMean", "log2FoldChange", "lfcSE", "padj", rank_names(mydata))]
  library("ggplot2")
  theme_set(theme_bw())
  scale_fill_discrete <- function(palname = "Set1", ...) {
    scale_fill_brewer(palette = palname, ...)
  }
  sigtabgen=sigtab
  
  #sigtabgen = subset(sigtab, !is.na(Genus))
  #sigtabgen = subset(sigtabgen, sigtabgen$Genus != "unclassified")
  # Phylum order
  x = tapply(sigtabgen$log2FoldChange, sigtabgen$Phylum, function(x) max(x))
  x = sort(x, TRUE)
  sigtabgen$Phylum = factor(as.character(sigtabgen$Phylum), levels=names(x))
  # Genus order
  if(as.character(LastRank)== "Genus"){
    sigtabgen$LastRank<-sigtabgen[,"Genus"]
  }else{
    sigtabgen$LastRank<-paste(sigtabgen[,"Genus"]," ",sigtabgen[,as.character(LastRank)])
  }
  x = tapply(sigtabgen$log2FoldChange, sigtabgen$LastRank, function(x) max(x))
  x = sort(x, TRUE)
  
  sigtabgen$LastRank = factor(as.character(sigtabgen$LastRank), levels=names(x))
  sigtabgen$log2Counts<-log2(sigtabgen$baseMean)
  sigtabgen$alpha<- 1 - sigtabgen$padj
  #pdf(fileForPlot)
  p<-ggplot(sigtabgen,aes(x=LastRank,y=log2FoldChange))
  p<-p+geom_point(aes(colour=Phylum,size=log2Counts,alpha=alpha))
  p<-p+scale_size_continuous(range=c(1,20))
  #+geom_point(aes(size=sigtabgen$log2Counts))+scale_size_area()
  p<-p+theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=0.5,size=10))
  p<-p+theme(legend.key.size=unit(1,"cm"))
  p<-p+ ggtitle(paste(stringForTitle," Data:",as.character(deparse(substitute(mydata))))) +
    theme(plot.title = element_text(lineheight=.7, face="bold"))
  print(p)
  #ggplot(sigtabgen, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=6)+scale_size(range=c(1,5))+
  # theme(axis.text.x = element_text(angle = -90, hjust = 0,size=3, vjust=0.5), legend.key.size=unit(0.5,"cm"),legend.text=element_text(size=3))
  ggsave(filename = fileForPlot,dpi=600,width=11, height=8.5)
  return(sigtab)
  #dev.off()
}


my.subset<-subset(ps.melt.x.3.0.phylum,ps.melt.x.3.0.phylum$Phylum=="Firmicutes")
p1<-boxplotNumericByGroup(my.subset,category="RiskGroup",variable="AbundanceProportionLog",nbvariable="Abundance",Rank="Phylum")
my.subset<-subset(ps.melt.x.3.0.genus,ps.melt.x.3.0.genus$Genus=="Prevotella")
p1<-boxplotNumericByGroup(my.subset,category="RiskGroup",variable="AbundanceProportionLog",nbvariable="Abundance",Rank="Genus")
p1<-boxplotNumericByGroup(my.subset,category="HIVStatus",variable="AbundanceProportionLog",nbvariable="Abundance",Rank="Genus")


### Let's move one step further and screen all genus for significant difference by using
### a more complex (and adequate) statistical framework using DESeq2 package and negative binomial
### distribution fits to detect over(under)-represented genus in a dichotomic condition

compareTwoGroups(mydata=x.3.0.genus,variable="HIVStatus",category1="HIVpos",category2="HIVneg",design=~HIVStatus,maxAlpha=0.01)
compareTwoGroups(mydata=x.3.0.genus,variable="SexualPractice",category1="MSM",category2="HET",design=~SexualPractice,maxAlpha=0.01)

#### End of Bacterial taxonomical analysis

### Start Ordination Analysis
### For this we need to use ecological distance. We will be using both Bray-Curtis distance
### and WUnifrac phylogenetic distance for this tutorial

#Let put our data aside
x.4.0<-ps_silva

### Remember that we have more than 16k OTUs in our data. Most of them are unique to a single sample
### thus indicating that they are probably artifact. Usually, these OTU don't have a great impact on distance
### calculations for this particular reason

#####Intended for de novo OTU-picking only, we keep OTUs that appear at least in two different samples
wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>0), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
#x.4.0<-tax_glom(x.4.0,taxrank="Genus")

### We transform the data to proportion abundances
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))

### We will be using NMDS which used an initial random seed for iteration
### Since we want to be able to repeat the same analysis we fix a seed for analysis
set.seed(12345)


### NMDS Ordination with Bray-Curtis distance
x.4.0.ord<-ordinate(x.4.0,"NMDS",distance="bray",trymax=200)
capture.output(file="NMDS_Bray_proportions_ordinfo.txt",x.4.0.ord)
#pdf("NMDS_Bray_proportions_stressplot.pdf")
stressplot(x.4.0.ord)
dev.off()
### Let's take a look at stressplot. Is it good?
stressplot(x.4.0.ord)
### Simple Plot
p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples

### Let's map some metadata with coloured and ellipses
p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples  + geom_point(aes(color=SexualPractice,shape=HIVStatus,size=Age)) +ggtitle("Unconstrained NMDS(Bray-Curtis)")+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14))+
  theme(legend.text=element_text(size=12,face="italic"),legend.title=element_text(size=12))+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=SexualPractice),level=0.95)+
  theme(plot.title=element_text(lineheight=1,face="bold",size=15))
#ggsave("NMDS_Bray_Proportion_Colour_withClusters_SizeByAge_BCN.pdf",) ################### <----  Save Pdf

### Now using HIV STatus
p.4.0.samples  + geom_point(size = 5.5,aes(color=HIVStatus,fill=HIVStatus,shape=SexualPractice)) +ggtitle("Unconstrained NMDS(Bray-Curtis)")+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=HIVStatus),level=0.95)+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))
ggsave("NMDS_Bray_Proportion_Colour_withClusters_byHIVStatus.pdf")

### Let's do a similar analysis but using phylogenetic OTU-OTU relationship
### using unifrac distances
x.4.0.ord<-ordinate(x.4.0,"NMDS",distance="wunifrac")
capture.output(file="NMDS_WUnifrac_proportions_ordinfo.txt",x.4.0.ord)
pdf("NMDS_WUnifrac_proportions_stressplot.pdf")
stressplot(x.4.0.ord)
dev.off()

### Simple Plot
p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples

### Now with metadata mapping
p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples  + geom_point(size = 5.5,aes(color=SexualPractice,fill=RiskGroup,shape=HIVStatus,size=Age)) +ggtitle("Unconstrained NMDS(Weighted Unifrac)")+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=SexualPractice),level=0.95)+
  ggtitle("Genus level")+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))
ggsave("NMDS_Wunifrac_Proportion_Colour_withClusters.pdf")

p.4.0.samples  + geom_point(size = 5.5,aes(color=HIVStatus,fill=HIVStatus,shape=SexualPractice,size=Age)) +ggtitle("Unconstrained NMDS(Weighted UniFrac)")+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=HIVStatus),level=0.95)+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))
ggsave("NMDS_WUnifrac_Proportion_Colour_withClusters_byHIVStatus.pdf")


############## Start Ordination Analysis Multiple Distance measures ##################################
require(plyr)
x.4.0<-xB
#####Intended for de novo OTU-picking only
wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>1), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))
set.seed(12345)


mymethod="NMDS"
dist=c("bray","canberra","manhattan", "euclidean", "kulczynski", "jaccard", "gower", "altGower", "horn", "mountford" , "binomial")
plist= llply(as.list(dist),function(i,physeq,mymethod){
  set.seed(1)
  ordi=ordinate(physeq,method=mymethod,distance=i,trymax=100)
  plot_ordination(physeq,ordi,"samples",color="RiskGroup")
},x.4.0,mymethod)
names(plist)<-dist

pdataframe = ldply(plist, function(x) {
  df = x$data[, 1:2]
  colnames(df) = c("Axis_1", "Axis_2")
  return(cbind(df, x$data))
})
names(pdataframe)[1] = "distance"


p = ggplot(pdataframe, aes(Axis_1, Axis_2)) +
  geom_point(size = 2, aes(shape=HIVStatus,color=SexualPractice,fill=SexualPractice)) +
  facet_wrap(~distance, scales = "free")+
  scale_colour_manual(values=c("darkred","darkolivegreen","dodgerblue2")) +
  theme_bw() +
  theme(panel.border=element_blank(),strip.text = element_text(size=12,face="bold"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),
        axis.text.x=element_text(size=4),
        axis.text.y=element_text(size=4))+
  theme(axis.title.x=element_text(size=4),axis.title.y=element_text(size=4))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  stat_ellipse(geom="polygon",aes(fill=SexualPractice),alpha=0.25,level=0.95)+
  ggtitle("Genus level")+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))
p
ggsave("NMDS_MultipleDistances_WithGenusProportion_WithClusters_byRiskGroup.pdf")
############## End Ordination Analysis Multiple Distance measures ####################################
################ End ordination ########################################


########################################################################
############## Start NMDS cluster Analysis #############################
### Until now we have just mapped our metadata into ordination space and have seen
### that there is an apparently clear association between some of our metadata
### and bacterial composition. Can we further confirm that?
### Let's do a formal clustering analysis

require(vegan)
require(cluster)
x.4.0<-xps_silva
#####Intended for de novo OTU-picking only
wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>0), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
#x.4.0<-tax_glom(x.4.0,taxrank="Genus")
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))
set.seed(12345)
x.4.0.ord<-ordinate(x.4.0,"NMDS",distance="bray",trymax=200)
x.4.0.ord

# Distance matrix (Bray-Curtis distance)
Dist<-vegdist(t(otu_table(x.4.0)),method="bray")
# Compute cluster analysis from 2 to 10 groups
################+++++++
Clusterization<-list()
for (i in 1:9){
  Clusterization[[i]]<-pam(Dist,k=i+1)
}

# Silhouette coefficients
Sil<-vector()
for (i in 1:9){Sil[i]<-Clusterization[[i]]$silinfo$avg.width}
# Plot Silhouette coefficient
plot(2:10,Sil,xlab="Number of clusters",ylab="Silhouette coefficient",
     pch=16,lwd=2,type="o",main="Silhouette coefficient")

# So, (by previous analysis) 2 clusters is the best number of groups according to
# silhouette coefficient
Cl<-pam(Dist,k=2)


# Cluster assignation, we define a cluster-label according to PAM analysis
Vector.K2<-as.factor(Cl$clustering)
my.sample_data<-data.frame(sample_data(x.4.0))
my.sample_data$Cluster<-Vector.K2
my.sample_data$SexualPractice<-as.factor(my.sample_data$SexualPractice)
sample_data(xB)<-my.sample_data
data<-data.frame(sample_data(xB))

### Statistics for association between PAM clustering and SexualPractice, in this case
xtabs(~data$Cluster+data$SexualPractice)
fisher.test(xtabs(~data$Cluster+data$SexualPractice))

### Statistics for association between PAM clustering and HIVStatus, in this case
xtabs(~data$Cluster+data$HIV_Status)
fisher.test(xtabs(~data$Cluster+data$HIV_Status))

### Let's use cluster labels for plotting
x.4.0<-ps_silva

wh0=genefilter_sample(x.4.0,filterfun_sample(function(x) x>0), A=0.01*nsamples(x.4.0))
x.4.0<-prune_taxa(wh0,x.4.0)
#x.4.0<-tax_glom(x.4.0,taxrank="Genus")
x.4.0 = transform_sample_counts(x.4.0, function(x) ((x/sum(x))))
set.seed(12345)
x.4.0.ord<-ordinate(x.4.0,"NMDS",distance="bray",trymax=200)

p.4.0.samples=plot_ordination(x.4.0,x.4.0.ord)
p.4.0.samples  + geom_point(size = 5,aes(color=Cluster,shape=SexualPractice)) + ggtitle("Unconstrained PCoA(Bray-Curtis)")+
  scale_colour_manual(values=c("darkgreen","darkorange"))+
  theme_bw()+
  theme(panel.border=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
        axis.line=element_line(colour="black"),axis.text.x=element_text(size=16),axis.text.y=element_text(size=16))+
  theme(axis.title.x=element_text(size=18),axis.title.y=element_text(size=18))+
  theme(legend.text=element_text(size=16,face="italic"),legend.title=element_text(size=16))+
  stat_ellipse(geom="polygon",alpha=0.25,aes(fill=Cluster),level=0.95)+
  ggtitle("Genus level")+
  theme(plot.title=element_text(lineheight=1,face="bold",size=19))
ggsave("NMDS_Bray_byNMDSCluster_color.pdf")
capture.output(file="NMDS_Bray_byNMDSCluster_FisherTest.txt",fisher.test(ftable(xtabs(~sample_data(x.4.0)$SexualPractice+sample_data(x.4.0)$Cluster))))

#######Adonis TESTs##########################################
### Until now we have addressed association

x.6.0<-ps_silva
wh0=genefilter_sample(x.6.0,filterfun_sample(function(x) x>1), A=0.01*nsamples(x.6.0))
x.6.0<-prune_taxa(wh0,x.6.0)
x.6.0.genus<-tax_glom(x.6.0,taxrank="Genus")
x.6.0<-x.6.0.genus
x.6.0 = transform_sample_counts(x.6.0, function(x) (((1* (x))/sum(x))))
#### Non-parametric analysis of variance using NPMANOVA through adonis function {vegan}package

x.6.0.response<-t(otu_table(x.6.0))
metadata<-data.frame(sample_data(x.6.0))

x.6.0.explanatory<-metadata[,c("Season","Sample.site","Sample.type")]
for(var in colnames(x.6.0.explanatory)){
  explanatory=data.frame(x.6.0.explanatory[,eval(var)])
  myadonis<-adonis(x.6.0.response~.,data=explanatory)
  capture.output(paste("Adonis on: ",var),file="adonis_tests.txt",append=T)
  capture.output(myadonis,file="adonis_tests.txt",append=T)
}

adonis(x.6.0.response~.,data=x.6.0.explanatory)
x.6.0.adonis<-adonis(x.6.0.response~.,data=x.6.0.explanatory)
simper(x.6.0.response,x.6.0.explanatory[,"Season"])


