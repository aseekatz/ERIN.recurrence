## Making a dotplot of LefSe-identified OTUs--UPDATED
# 3.18.16
# Anna M. Seekatz

# adapted from http://polisci.msu.edu/jacoby/research/dotplots/tpm/Creating%20figures/Creating%20Figure%204.R
library(plyr)
library(reshape2)
library(Hmisc)

library(lattice)
library(dplyr)

# files used:
	# updated_meta_fromkrishna_03.08.16.txt: updated meta
	# suberin.relOTUs.txt: relative abundance of OTUs, filtered (same info applicable as before)
	# erinfmt.new.taxonomy.names.txt: edited taxonomy files with OTU info (same info applicable as before)
	# updated_erinsubset_all.lefse.results.txt: new lefse results, edited from mothur
	


# meta and OTU files
meta<-read.table(file="updated_meta_fromkrishna_03.08.16.txt", header=TRUE)
sub.otus<-read.table(file="../suberin.relOTUs.txt", header=TRUE)
meta2<-meta[-which(meta$seqID %in% c("DA3240", "DA3260_P", "DA3298", "DA3299", "DA3376")), ]
meta<-meta2
tax<-read.table(file="../erinfmt.new.taxonomy.names.txt", header=TRUE)
keep<-as.character(colnames(sub.otus[1:506]))
filtered.tax<-tax[tax$OTU %in% keep, ]
colnames(sub.otus)<-filtered.tax$taxname
sub.otus$sampleID<-rownames(sub.otus)

# merge meta and OTUs:
sub.otus$sampleID<-as.factor(sub.otus$sampleID)
sub.all<-merge(meta, sub.otus, by.x="seqID", by.y="sampleID")
rownames(sub.all)<-sub.all$seqID

# lefse results (filtered file to include only the significant ones)
lefse<-read.table(file="updated_erinsubset_all.lefse.results.txt", header=TRUE)
	# get taxnames for lefse OTU file:
keep<-as.character(lefse$OTU)
filtered.tax<-tax[tax$OTU %in% keep, ]
lefse$otuname<-filtered.tax$taxname
lefse <- lefse[order(lefse$clinical_Class, lefse$otuname),]
lefse.otus<-as.character(lefse$otuname[1:5])
# for severe ones:
#lefse <- lefse[order(lefse$severe_Class, lefse$otuname),]
#lefse.otus<-as.character(lefse$otuname[1:7])
	# these were all the sign. OTUs in the index samples of recurrent vs. nonrecurrent patients

#### Fig. 4B: plotting lefse OTUs by group
# these are the OTUs significant between positive and negative samples
# let's limit our files to only the lefse results, and whatever category you are using: 
otus.df<-sub.all[, which(colnames(sub.all) %in% lefse.otus)]
otus.df$sampleID<-rownames(otus.df)
groups<-sub.all[, c("seqID", "group_reinfection")]
otus<-merge(groups, otus.df, by.x="seqID", by.y="sampleID")
otus<-otus[,2:7]
	# now you have your list of significant otus in a dataframe!
	
# create your stats summary dataframe:
# column names: Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria
# column names, with "": "Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"
x<-otus	
means<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) mean =mean(x) )
means2<-melt(means, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.means<-means2$value

sds<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) sd =sd(x) )
sds2<-melt(sds, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.sds<-sds2$value
otu.highsd<-otu.means+sds2$value
otu.lowsd<-otu.means-sds2$value

medians<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) median =median(x) )
medians2<-melt(medians, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.medians<-medians2$value

ses<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) se=sd(x)/sqrt(length(x)) )
ses2<-melt(ses, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.ses<-ses2$value
otu.highse<-otu.means+ses2$value
otu.lowse<-otu.means-ses2$value

highq<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) highq=quantile(x, 0.75) )
highq2<-melt(ses, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.highq<-otu.medians+highq2$value

lowq<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ group_reinfection,data = x,FUN=function(x) lowq=quantile(x, 0.75) )
lowq2<-melt(ses, id.vars = "group_reinfection", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
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

# option 1 (used as manuscript Supplemental Figure S1:
dotchart(averages, labels=labels, xlab='relative abundance (mean + se)',  pch=20, col=c("chartreuse3", "orange", "darkgoldenrod"),
		groups=summary$group, cex=0.7,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(15, 8, 1, 16, 9, 2, 17, 10, 3, 18, 11, 4, 19, 12, 5), 
	averages+ranges, c(15, 8, 1, 16, 9, 2, 17, 10, 3, 18, 11, 4, 19, 12, 5), 
	col=c("chartreuse3", "orange", "darkgoldenrod"), lwd=2)
# note: adding the bars on your graph is a bit confusing
# you may have to play around with the ordering, since this is a grouped dotplot
# future ways of making this step simpler are appreciated!
	
# option 2:
## this one was used for Fig 4:  
dotchart(averages, labels="", xlab='relative abundance (mean + se)',  pch=20, col=c("chartreuse3", "orange", "darkgoldenrod"),
		groups=as.factor(summary$otu), cex=0.8,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(21:23, 16:18, 11:13, 6:8, 1:3), 
	averages+ranges, c(21:23, 16:18, 11:13, 6:8, 1:3), 
	col=c("chartreuse3", "orange", "darkgoldenrod"), lwd=2)
legend("bottomright", groupnames, col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")
	#Looks BEAUTIFUL
	

#### Fig. 4A: plotting lefse OTUs by clinical status (negative vs. positive)
# these are the OTUs significant between positive and negative samples
otus.df<-sub.all[, which(colnames(sub.all) %in% lefse.otus)]
otus.df$sampleID<-rownames(otus.df)
groups<-sub.all[, c("seqID", "POS_NEG")]
otus.clinical<-merge(groups, otus.df, by="sampleID")
otus.clinical<-otus.clinical[,2:7]

otus.df<-sub.all[, which(colnames(sub.all) %in% lefse.otus)]
otus.df$sampleID<-rownames(otus.df)
groups<-sub.all[, c("seqID", "POS_NEG")]
otus.clinical<-merge(groups, otus.df, by.x="seqID", by.y="sampleID")
otus.clinical<-otus.clinical[,2:7]
	
# create your stats summary dataframe:
# column names: Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria
# column names, with "": "Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"
x<-otus.clinical
means<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ POS_NEG,data = x,FUN=function(x) mean =mean(x) )
means2<-melt(means, id.vars = "POS_NEG", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.means<-means2$value

sds<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ POS_NEG,data = x,FUN=function(x) sd =sd(x) )
sds2<-melt(sds, id.vars = "POS_NEG", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.sds<-sds2$value
otu.highsd<-otu.means+sds2$value
otu.lowsd<-otu.means-sds2$value

medians<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ POS_NEG,data = x,FUN=function(x) median =median(x) )
medians2<-melt(medians, id.vars = "POS_NEG", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.medians<-medians2$value

ses<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ POS_NEG,data = x,FUN=function(x) se=sd(x)/sqrt(length(x)) )
ses2<-melt(ses, id.vars = "POS_NEG", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.ses<-ses2$value
otu.highse<-otu.means+ses2$value
otu.lowse<-otu.means-ses2$value

highq<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ POS_NEG,data = x,FUN=function(x) highq=quantile(x, 0.75) )
highq2<-melt(ses, id.vars = "POS_NEG", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.highq<-otu.medians+highq2$value

lowq<-aggregate(cbind(Otu00001_Enterobacteriaceae,Otu00003_Bacteroides,Otu00012_Clostridium_XI,Otu00014_Streptococcus, Otu00050_Bacteria) ~ POS_NEG,data = x,FUN=function(x) lowq=quantile(x, 0.75) )
lowq2<-melt(ses, id.vars = "POS_NEG", measure.vars = c("Otu00001_Enterobacteriaceae","Otu00003_Bacteroides","Otu00012_Clostridium_XI","Otu00014_Streptococcus", "Otu00050_Bacteria"))
otu.lowq<-otu.medians-lowq2$value

# add all the data together:
labels<-paste(as.character(means2$POS_NEG), as.character(means2$variable))
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
summary$group <-means2$POS_NEG

# plot it:
labels <- summary$otu
labels <- gsub("_", ": ", labels)
labels <- gsub("000", "", labels)
summary$otu<-gsub("_", ": ", summary$otu)
summary$otu<-gsub("000", "", summary$otu)
averages <- summary$mean
ranges <- summary$se
groupnames<-as.character(unique(summary$group))

# option 1 (used as manuscript Figure 4):
dotchart(averages, labels="", xlab='relative abundance (mean + se)',  pch=20, col=c("grey47", "magenta"),
		groups=as.factor(summary$otu), cex=0.8, cex.lab=0.8,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(17:18, 13:14, 9:10, 5:6, 1:2), 
	averages+ranges, c(17:18, 13:14, 9:10, 5:6, 1:2), 
	col=c("grey47", "magenta"), lwd=2)
legend("bottomright", groupnames, col=c("grey47", "magenta"), pch=19, cex=0.6, bg="white")


# option 2:
dotchart(averages, labels=labels, xlab='relative abundance (mean + se)',  pch=20, col=c("grey47", "magenta"),
		groups=summary$group, cex=0.7,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(8, 1, 9, 2, 10, 3, 11, 4, 12, 5), 
	averages+ranges, c(8, 1, 9, 2, 10, 3, 11, 4, 12, 5), 
	col=c("grey47", "magenta"), lwd=2)
# note: adding the bars on your graph is a bit confusing
# you may have to play around with the ordering, since this is a grouped dotplot
# future ways of making this step simpler are appreciated!


###
#---
###

# for severity:
# used the same steps above to get 'lefse.otus' that were different with severity

otus.df<-sub.all[, which(colnames(sub.all) %in% lefse.otus)]
otus.df$sampleID<-rownames(otus.df)
groups<-sub.all[, c("seqID", "severeidsa")]
otus.severe<-merge(groups, otus.df, by.x="seqID", by.y="sampleID")
otus.severe<-otus.severe[,2:9]
	
# create your stats summary dataframe:
# column names: Otu00004_Escherichia, Otu00015_Lachnospiraceae, Otu00017_Bifidobacterium, Otu00030_Blautia, Otu00044_Flavonifractor, Otu00050_Bacteria, Otu00097_Lactococcus
# column names, with "": "Otu00004_Escherichia", "Otu00015_Lachnospiraceae", "Otu00017_Bifidobacterium", "Otu00030_Blautia", "Otu00044_Flavonifractor", "Otu00050_Bacteria", "Otu00097_Lactococcus"
x<-otus.severe
means<-aggregate(cbind(Otu00004_Escherichia, Otu00015_Lachnospiraceae, Otu00017_Bifidobacterium, Otu00030_Blautia, Otu00044_Flavonifractor, Otu00050_Bacteria, Otu00097_Lactococcus) ~ severeidsa,data = x,FUN=function(x) mean =mean(x) )
means2<-melt(means, id.vars = "severeidsa", measure.vars = c("Otu00004_Escherichia", "Otu00015_Lachnospiraceae", "Otu00017_Bifidobacterium", "Otu00030_Blautia", "Otu00044_Flavonifractor", "Otu00050_Bacteria", "Otu00097_Lactococcus"))
otu.means<-means2$value

sds<-aggregate(cbind(Otu00004_Escherichia, Otu00015_Lachnospiraceae, Otu00017_Bifidobacterium, Otu00030_Blautia, Otu00044_Flavonifractor, Otu00050_Bacteria, Otu00097_Lactococcus) ~ severeidsa,data = x,FUN=function(x) sd =sd(x) )
sds2<-melt(sds, id.vars = "severeidsa", measure.vars = c("Otu00004_Escherichia", "Otu00015_Lachnospiraceae", "Otu00017_Bifidobacterium", "Otu00030_Blautia", "Otu00044_Flavonifractor", "Otu00050_Bacteria", "Otu00097_Lactococcus"))
otu.sds<-sds2$value
otu.highsd<-otu.means+sds2$value
otu.lowsd<-otu.means-sds2$value

medians<-aggregate(cbind(Otu00004_Escherichia, Otu00015_Lachnospiraceae, Otu00017_Bifidobacterium, Otu00030_Blautia, Otu00044_Flavonifractor, Otu00050_Bacteria, Otu00097_Lactococcus) ~ severeidsa,data = x,FUN=function(x) median =median(x) )
medians2<-melt(medians, id.vars = "severeidsa", measure.vars = c("Otu00004_Escherichia", "Otu00015_Lachnospiraceae", "Otu00017_Bifidobacterium", "Otu00030_Blautia", "Otu00044_Flavonifractor", "Otu00050_Bacteria", "Otu00097_Lactococcus"))
otu.medians<-medians2$value

ses<-aggregate(cbind(Otu00004_Escherichia, Otu00015_Lachnospiraceae, Otu00017_Bifidobacterium, Otu00030_Blautia, Otu00044_Flavonifractor, Otu00050_Bacteria, Otu00097_Lactococcus) ~ severeidsa,data = x,FUN=function(x) se=sd(x)/sqrt(length(x)) )
ses2<-melt(ses, id.vars = "severeidsa", measure.vars = c("Otu00004_Escherichia", "Otu00015_Lachnospiraceae", "Otu00017_Bifidobacterium", "Otu00030_Blautia", "Otu00044_Flavonifractor", "Otu00050_Bacteria", "Otu00097_Lactococcus"))
otu.ses<-ses2$value
otu.highse<-otu.means+ses2$value
otu.lowse<-otu.means-ses2$value

highq<-aggregate(cbind(Otu00004_Escherichia, Otu00015_Lachnospiraceae, Otu00017_Bifidobacterium, Otu00030_Blautia, Otu00044_Flavonifractor, Otu00050_Bacteria, Otu00097_Lactococcus) ~ severeidsa,data = x,FUN=function(x) highq=quantile(x, 0.75) )
highq2<-melt(ses, id.vars = "severeidsa", measure.vars = c("Otu00004_Escherichia", "Otu00015_Lachnospiraceae", "Otu00017_Bifidobacterium", "Otu00030_Blautia", "Otu00044_Flavonifractor", "Otu00050_Bacteria", "Otu00097_Lactococcus"))
otu.highq<-otu.medians+highq2$value

lowq<-aggregate(cbind(Otu00004_Escherichia, Otu00015_Lachnospiraceae, Otu00017_Bifidobacterium, Otu00030_Blautia, Otu00044_Flavonifractor, Otu00050_Bacteria, Otu00097_Lactococcus) ~ severeidsa,data = x,FUN=function(x) lowq=quantile(x, 0.75) )
lowq2<-melt(ses, id.vars = "severeidsa", measure.vars = c("Otu00004_Escherichia", "Otu00015_Lachnospiraceae", "Otu00017_Bifidobacterium", "Otu00030_Blautia", "Otu00044_Flavonifractor", "Otu00050_Bacteria", "Otu00097_Lactococcus"))
otu.lowq<-otu.medians-lowq2$value

# add all the data together:
labels<-paste(as.character(means2$severeidsa), as.character(means2$variable))
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
summary$group <-means2$severeidsa

# plot it:
labels <- summary$otu
labels <- gsub("_", ": ", labels)
labels <- gsub("000", "", labels)
summary$otu<-gsub("_", ": ", summary$otu)
summary$otu<-gsub("000", "", summary$otu)
averages <- summary$mean
ranges <- summary$se
groupnames<-as.character(unique(summary$group))

# option 1 (used as manuscript Figure 4):
# since the first OTU is significantly more abundant than others, let's split these up into 2 quadrants to emphasize the others...
summary1<-summary[1:2, ]
labels1 <- summary1$otu
labels1 <- gsub("_", ": ", labels)
labels1 <- gsub("000", "", labels)
summary1$otu<-gsub("_", ": ", summary1$otu)
summary1$otu<-gsub("000", "", summary1$otu)
averages1 <- summary1$mean
ranges1 <- summary1$se
groupnames1<-as.character(unique(summary1$group))

summary2<-summary[3:14, ]
labels2 <- summary2$otu
labels2 <- gsub("_", ": ", labels)
labels2 <- gsub("000", "", labels)
summary2$otu<-gsub("_", ": ", summary2$otu)
summary2$otu<-gsub("000", "", summary2$otu)
averages2 <- summary2$mean
ranges2 <- summary2$se
groupnames2<-as.character(unique(summary2$group))

# summary 1:
par(mfrow=c(2,1))
dotchart(averages1, labels="", xlab='relative abundance (mean + se)',  pch=20, col=c("blue", "red"),
		groups=as.factor(summary1$otu), cex=0.8, cex.lab=0.8,
         xlim=c(min(averages1-ranges1)-1, max(averages1+ranges1)+1))
segments(averages1-ranges1, 
	c(1:2), 
	averages1+ranges1, c(1:2), 
	col=c("blue", "red"), lwd=2)
legend("bottomright", c("not severe", "severe"), col=c("blue", "red"), pch=19, cex=0.6, bg="white")
# summary 2:
dotchart(averages2, labels="", xlab='relative abundance (mean + se)',  pch=20, col=c("blue", "red"),
		groups=as.factor(summary2$otu), cex=0.8, cex.lab=0.8,
         xlim=c(min(averages2-ranges2)-1, max(averages2+ranges2)+1))
segments(averages2-ranges2, 
	c(21:22, 17:18, 13:14, 9:10, 5:6, 1:2), 
	averages2+ranges2, c(21:22, 17:18, 13:14, 9:10, 5:6, 1:2), 
	col=c("blue", "red"), lwd=2)
legend("bottomright", c("not severe", "severe"), col=c("blue", "red"), pch=19, cex=0.6, bg="white")
### this still doesn't give even graphs...just go with summary, and chop it down:

# summary (all):
dotchart(averages, labels="", xlab='relative abundance (mean + se)',  pch=20, col=c("blue", "red"),
		groups=as.factor(summary$otu), cex=0.8, cex.lab=0.8,
         xlim=c(0, max(averages+ranges)+1))
segments(averages-ranges, 
	c(25:26, 21:22, 17:18, 13:14, 9:10, 5:6, 1:2), 
	averages+ranges, c(25:26, 21:22, 17:18, 13:14, 9:10, 5:6, 1:2), 
	col=c("blue", "red"), lwd=2)
legend("bottomright", c("not severe", "severe"), col=c("blue", "red"), pch=19, cex=0.6, bg="white")


# option 2:
dotchart(averages, labels=labels, xlab='relative abundance (mean + se)',  pch=20, col=c("blue", "red"),
		groups=summary$group, cex=0.7,
         xlim=c(min(averages-ranges)-1, max(averages+ranges)+1))
segments(averages-ranges, 
	c(8, 1, 9, 2, 10, 3, 11, 4, 12, 5), 
	averages+ranges, c(8, 1, 9, 2, 10, 3, 11, 4, 12, 5), 
	col=c("blue", "red"), lwd=2)
# note: adding the bars on your graph is a bit confusing
# you may have to play around with the ordering, since this is a grouped dotplot
# future ways of making this step simpler are appreciated!

