## Making a dotplot of LefSe-identified OTUs
# 1.4.16
# Anna M. Seekatz

# adapted from http://polisci.msu.edu/jacoby/research/dotplots/tpm/Creating%20figures/Creating%20Figure%204.R
library(plyr)
library(reshape2)
library(Hmisc)

library(lattice)
library(dplyr)

# files used:
	# suberin_all.alpha.phylo.txt: contains summary and phylotype information collated previously
	# suberin.relOTUs.txt: relative abundance of OTUs, filtered
	# erinfmt.new.taxonomy.names.txt: edited taxonomy files with OTU info
	# erinsubset_all.lefse.results.txt: lefse results, edited from mothur

# meta and OTU files
all<-read.table(file="suberin_all.alpha.phylo.txt", header=TRUE)
all<-all[, c(1, 104:142)]				
	# this had phylotypes, and we want OTUs; thus we subset
sub.otus<-read.table(file="suberin.relOTUs.txt", header=TRUE)
tax<-read.table(file="erinfmt.new.taxonomy.names.txt", header=TRUE)
keep<-as.character(colnames(sub.otus[1:506]))
filtered.tax<-tax[tax$OTU %in% keep, ]
colnames(sub.otus)<-filtered.tax$taxname
sub.otus$sampleID<-rownames(sub.otus)

# merge meta and OTUs:
sub.all<-merge(all, sub.otus, by.x="sampleID", by.y="sampleID")
rownames(sub.all)<-sub.all$sampleID

# lefse results (filtered file to include only the significant ones)
lefse<-read.table(file="erinsubset_all.lefse.results.txt", header=TRUE)
	# get taxnames for lefse OTU file:
keep<-as.character(lefse$OTU)
filtered.tax<-tax[tax$OTU %in% keep, ]
lefse$otuname<-filtered.tax$taxname
lefse <- lefse[order(lefse$clin_Class, lefse$group_Class, lefse$otuname),]
lefse.otus<-as.character(lefse$otuname)

#### Fig. 5B: plotting lefse OTUs by group
# these are the OTUs significant between positive and negative samples
# let's limit our files to only the lefse results, and whatever category you are using: 
otus.df<-sub.all[, which(colnames(sub.all) %in% lefse.otus)]
otus.df$sampleID<-rownames(otus.df)
groups<-sub.all[, c("sampleID", "group_reinfection")]
otus<-merge(groups, otus.df, by="sampleID")
otus<-otus[,2:14]
	# now you have your list of significant otus in a dataframe!
	
# create your stats summary dataframe:
# column names: Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria
# column names, with "": "Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"
x<-otus	
means<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) mean =mean(x) )
means2<-melt(means, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.means<-means2$value

sds<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) sd =sd(x) )
sds2<-melt(sds, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.sds<-sds2$value
otu.highsd<-otu.means+sds2$value
otu.lowsd<-otu.means-sds2$value

medians<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) median =median(x) )
medians2<-melt(medians, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.medians<-medians2$value

ses<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) se=sd(x)/sqrt(length(x)) )
ses2<-melt(ses, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.ses<-ses2$value
otu.highse<-otu.means+ses2$value
otu.lowse<-otu.means-ses2$value

highq<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) highq=quantile(x, 0.75) )
highq2<-melt(ses, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.highq<-otu.medians+highq2$value

lowq<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) lowq=quantile(x, 0.75) )
lowq2<-melt(ses, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.lowq<-otu.medians-lowq2$value

# add all the data together:
labels<-paste(as.character(means2$group_reinfection), as.character(means2$variable))
summary<-data.frame(label=labels, 
			mean = otu.means, 
			sd = otu.sds, 
			se = otu.ses, 
			median = otu.medians, 
			highq = otu.highq,
			lowq = otu.lowq,
			highse = otu.highse,
			lowse = otu.lowse,
			highsd = otu.highsd,
			lowsd = otu.lowsd
			)
summary$sequence <- seq(1, length(summary$label))
summary$label <- reorder(summary$label, summary$sequence)
summary$otu <- means2$variable
summary$group <-means2$group_reinfection

# plot it:
labels <- summary$otu
labels <- gsub("_", ": ", labels)
summary$otu<-gsub("_", ": ", summary$otu)
averages <- summary$mean
ranges <- summary$se
groupnames<-as.character(unique(summary$group))

# option 1 (used as manuscript Supplemental Figure S3:
dotchart(averages, labels=labels, xlab='relative abundance (mean + se)',  pch=20, col=c("chartreuse3", "orange", "darkgoldenrod"),
		groups=summary$group, cex=0.7,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(29, 15, 1, 30, 16, 2, 31, 17, 3, 32, 18, 4, 33, 19, 5, 34, 20, 6, 35, 21, 7, 36, 22, 8, 37, 23, 9, 38, 24, 10, 39, 25, 11, 40, 26, 12), 
	averages+ranges, c(29, 15, 1, 30, 16, 2, 31, 17, 3, 32, 18, 4, 33, 19, 5, 34, 20, 6, 35, 21, 7, 36, 22, 8, 37, 23, 9, 38, 24, 10, 39, 25, 11, 40, 26, 12), 
	col=c("chartreuse3", "orange", "darkgoldenrod"), lwd=2)
# note: adding the bars on your graph is a bit confusing
# you may have to play around with the ordering, since this is a grouped dotplot
# future ways of making this step simpler are appreciated!
	
# option 2:  
dotchart(averages, labels="", xlab='relative abundance (mean + se)',  pch=20, col=c("chartreuse3", "orange", "darkgoldenrod"),
		groups=as.factor(summary$otu), cex=0.8,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(56:58, 51:53, 46:48, 41:43, 36:38, 31:33, 26:28, 21:23, 16:18, 11:13, 6:8, 1:3), 
	averages+ranges, c(56:58, 51:53, 46:48, 41:43, 36:38, 31:33, 26:28, 21:23, 16:18, 11:13, 6:8, 1:3), 
	col=c("chartreuse3", "orange", "darkgoldenrod"), lwd=2)
legend("bottomright", groupnames, col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")
	#Looks BEAUTIFUL
	
#### Fig. 5A: plotting lefse OTUs by clinical status (negative vs. positive)
# these are the OTUs significant between positive and negative samples
otus.df<-sub.all[, which(colnames(sub.all) %in% lefse.otus)]
otus.df$sampleID<-rownames(otus.df)
groups<-sub.all[, c("sampleID", "clinical_result")]
otus.clinical<-merge(groups, otus.df, by="sampleID")
otus.clinical<-otus.clinical[,2:14]
	
# create your stats summary dataframe:
# column names: Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria
# column names, with "": "Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"
x<-otus.clinical
means<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) mean =mean(x) )
means2<-melt(means, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.means<-means2$value

sds<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) sd =sd(x) )
sds2<-melt(sds, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.sds<-sds2$value
otu.highsd<-otu.means+sds2$value
otu.lowsd<-otu.means-sds2$value

medians<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) median =median(x) )
medians2<-melt(medians, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.medians<-medians2$value

ses<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) se=sd(x)/sqrt(length(x)) )
ses2<-melt(ses, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.ses<-ses2$value
otu.highse<-otu.means+ses2$value
otu.lowse<-otu.means-ses2$value

highq<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) highq=quantile(x, 0.75) )
highq2<-melt(ses, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.highq<-otu.medians+highq2$value

lowq<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) lowq=quantile(x, 0.75) )
lowq2<-melt(ses, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.lowq<-otu.medians-lowq2$value

# add all the data together:
labels<-paste(as.character(means2$clinical_result), as.character(means2$variable))
summary<-data.frame(label=labels, 
			mean = otu.means, 
			sd = otu.sds, 
			se = otu.ses, 
			median = otu.medians, 
			highq = otu.highq,
			lowq = otu.lowq,
			highse = otu.highse,
			lowse = otu.lowse,
			highsd = otu.highsd,
			lowsd = otu.lowsd
			)
summary$sequence <- seq(1, length(summary$label))
summary$label <- reorder(summary$label, summary$sequence)
summary$otu <- means2$variable
summary$group <-means2$clinical_result

# plot it:
labels <- summary$otu
labels <- gsub("_", ": ", labels)
labels <- gsub("000", "", labels)
summary$otu<-gsub("_", ": ", summary$otu)
summary$otu<-gsub("000", "", summary$otu)
averages <- summary$mean
ranges <- summary$se
groupnames<-as.character(unique(summary$group))

# option 1 (used as manuscript Figure 5):
dotchart(averages, labels="", xlab='relative abundance (mean + se)',  pch=20, col=c("grey47", "magenta"),
		groups=as.factor(summary$otu), cex=0.8, cex.lab=0.8,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(45:46, 41:42, 37:38, 33:34, 29:30, 25:26, 21:22, 17:18, 13:14, 9:10, 5:6, 1:2), 
	averages+ranges, c(45:46, 41:42, 37:38, 33:34, 29:30, 25:26, 21:22, 17:18, 13:14, 9:10, 5:6, 1:2), 
	col=c("grey47", "magenta"), lwd=2)


# option 2:
dotchart(averages, labels=labels, xlab='relative abundance (mean + se)',  pch=20, col=c("grey47", "magenta"),
		groups=summary$group, cex=0.7,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(15, 1, 16, 2, 17, 3, 18, 4, 19, 5, 20, 6, 21, 7, 22, 8, 23, 9, 24, 10, 25, 11, 26, 12), 
	averages+ranges, c(15, 1, 16, 2, 17, 3, 18, 4, 19, 5, 20, 6, 21, 7, 22, 8, 23, 9, 24, 10, 25, 11, 26, 12), 
	col=c("grey47", "magenta"), lwd=2)
# note: adding the bars on your graph is a bit confusing
# you may have to play around with the ordering, since this is a grouped dotplot
# future ways of making this step simpler are appreciated!


#####################
# for double graphs of both, copy/paste below:

par(mai=c(1.02,0.82,0.82,0.22))
otus.df<-sub.all[, which(colnames(sub.all) %in% lefse.otus)]
otus.df$sampleID<-rownames(otus.df)
groups<-sub.all[, c("sampleID", "clinical_result")]
otus.clinical<-merge(groups, otus.df, by="sampleID")
otus.clinical<-otus.clinical[,2:14]
# create your stats summary dataframe:
# column names: Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria
# column names, with "": "Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"
x<-otus.clinical
means<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) mean =mean(x) )
means2<-melt(means, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.means<-means2$value

sds<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) sd =sd(x) )
sds2<-melt(sds, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.sds<-sds2$value
otu.highsd<-otu.means+sds2$value
otu.lowsd<-otu.means-sds2$value

medians<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) median =median(x) )
medians2<-melt(medians, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.medians<-medians2$value

ses<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) se=sd(x)/sqrt(length(x)) )
ses2<-melt(ses, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.ses<-ses2$value
otu.highse<-otu.means+ses2$value
otu.lowse<-otu.means-ses2$value

highq<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) highq=quantile(x, 0.75) )
highq2<-melt(ses, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.highq<-otu.medians+highq2$value

lowq<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ clinical_result,data = x,FUN=function(x) lowq=quantile(x, 0.75) )
lowq2<-melt(ses, id.vars = "clinical_result", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.lowq<-otu.medians-lowq2$value

# add all the data together:
labels<-paste(as.character(means2$clinical_result), as.character(means2$variable))
summary<-data.frame(label=labels, 
			mean = otu.means, 
			sd = otu.sds, 
			se = otu.ses, 
			median = otu.medians, 
			highq = otu.highq,
			lowq = otu.lowq,
			highse = otu.highse,
			lowse = otu.lowse,
			highsd = otu.highsd,
			lowsd = otu.lowsd
			)
summary$sequence <- seq(1, length(summary$label))
summary$label <- reorder(summary$label, summary$sequence)
summary$otu <- means2$variable
summary$group <-means2$clinical_result

# plot it:
labels <- summary$otu
labels <- gsub("_", ": ", labels)
summary$otu<-gsub("_", ": ", summary$otu)
averages <- summary$mean
ranges <- summary$se
groupnames<-as.character(unique(summary$group))

# option 2:
dotchart(averages, labels=labels, xlab='relative abundance (mean + se)',  pch=20, col=c("grey47", "magenta"),
		groups=summary$group, cex=0.7,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(15, 1, 16, 2, 17, 3, 18, 4, 19, 5, 20, 6, 21, 7, 22, 8, 23, 9, 24, 10, 25, 11, 26, 12), 
	averages+ranges, c(15, 1, 16, 2, 17, 3, 18, 4, 19, 5, 20, 6, 21, 7, 22, 8, 23, 9, 24, 10, 25, 11, 26, 12), 
	col=c("grey47", "magenta"), lwd=2)
# note: adding the bars on your graph is a bit confusing
# you may have to play around with the ordering, since this is a grouped dotplot
# future ways of making this step simpler are appreciated!

## 2nd graph:
otus.df<-sub.all[, which(colnames(sub.all) %in% lefse.otus)]
otus.df$sampleID<-rownames(otus.df)
groups<-sub.all[, c("sampleID", "group_reinfection")]
otus<-merge(groups, otus.df, by="sampleID")
otus<-otus[,2:14]
	# now you have your list of significant otus in a dataframe!
# create your stats summary dataframe:
# column names: Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria
# column names, with "": "Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"
x<-otus	
means<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) mean =mean(x) )
means2<-melt(means, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.means<-means2$value

sds<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) sd =sd(x) )
sds2<-melt(sds, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.sds<-sds2$value
otu.highsd<-otu.means+sds2$value
otu.lowsd<-otu.means-sds2$value

medians<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) median =median(x) )
medians2<-melt(medians, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.medians<-medians2$value

ses<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) se=sd(x)/sqrt(length(x)) )
ses2<-melt(ses, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.ses<-ses2$value
otu.highse<-otu.means+ses2$value
otu.lowse<-otu.means-ses2$value

highq<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) highq=quantile(x, 0.75) )
highq2<-melt(ses, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.highq<-otu.medians+highq2$value

lowq<-aggregate(cbind(Otu00001_Enterobacteriaceae, Otu00002_Enterococcus, Otu00003_Bacteroides, Otu00006_Bacteroides, Otu00010_Bacteroides, Otu00016_Clostridium_XlVa, Otu00023_Roseburia, Otu00030_Blautia, Otu00034_Clostridium_XVIII, Otu00044_Flavonifractor, Otu00048_Subdoligranulum, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) lowq=quantile(x, 0.75) )
lowq2<-melt(ses, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae", "Otu00002_Enterococcus", "Otu00003_Bacteroides", "Otu00006_Bacteroides", "Otu00010_Bacteroides", "Otu00016_Clostridium_XlVa", "Otu00023_Roseburia", "Otu00030_Blautia", "Otu00034_Clostridium_XVIII", "Otu00044_Flavonifractor", "Otu00048_Subdoligranulum", "Otu00050_Bacteria"))
otu.lowq<-otu.medians-lowq2$value

# add all the data together:
labels<-paste(as.character(means2$group_reinfection), as.character(means2$variable))
summary<-data.frame(label=labels, 
			mean = otu.means, 
			sd = otu.sds, 
			se = otu.ses, 
			median = otu.medians, 
			highq = otu.highq,
			lowq = otu.lowq,
			highse = otu.highse,
			lowse = otu.lowse,
			highsd = otu.highsd,
			lowsd = otu.lowsd
			)
summary$sequence <- seq(1, length(summary$label))
summary$label <- reorder(summary$label, summary$sequence)
summary$otu <- means2$variable
summary$group <-means2$group_reinfection

# plot it:
labels <- summary$otu
labels <- gsub("_", ": ", labels)
summary$otu<-gsub("_", ": ", summary$otu)
averages <- summary$mean
ranges <- summary$se
groupnames<-as.character(unique(summary$group))

dotchart(averages, labels="", xlab='relative abundance (mean + se)',  pch=20, col=c("chartreuse3", "orange", "darkgoldenrod"),
		groups=as.factor(summary$otu), cex=0.7,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(56:58, 51:53, 46:48, 41:43, 36:38, 31:33, 26:28, 21:23, 16:18, 11:13, 6:8, 1:3), 
	averages+ranges, c(56:58, 51:53, 46:48, 41:43, 36:38, 31:33, 26:28, 21:23, 16:18, 11:13, 6:8, 1:3), 
	col=c("chartreuse3", "orange", "darkgoldenrod"), lwd=2)
legend("bottomright", groupnames, col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")