## Making a dotplot of LefSe-identified OTUs
# 1.7.15
# Anna M. Seekatz

library(shape)

# used files:
	# erinsubset.thetayc0.03.pcoa.axes
	# erinsubset.thetayc0.03.nmds.axes
	# erinfmt.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.summary
	# ERIN.UMFMT_metadata_filtered.txt
	# erinfmt_extra.meta.txt
	# erinsubset_summary.txt #created from the mothur files
	
# read in files and merge together:
#meta<-read.table(file="ERIN.UMFMT_metadata.txt", header=TRUE)
meta<-read.table(file="../ERIN.UMFMT_metadata_filtered.txt", header=TRUE)
	erinsub.meta<-meta[meta$group %in% c("recurrent", "nonrecurrent"), ]
pcoa<-read.table(file="erinsubset_mothurfiles/erinsubset.thetayc0.03.pcoa.axes", header=TRUE)
	colnames(pcoa)[1:4] <- c("sampleID", "pcoa_axis1", "pcoa_axis2", "pcoa_axis3")
nmds<-read.table(file="erinsubset_mothurfiles/erinsubset.thetayc0.03.nmds.axes", header=TRUE)
	colnames(nmds)[1:4] <- c("sampleID", "nmds_axis1", "nmds_axis2", "nmds_axis3")
sum<-read.table(file="../mothur.files_3.18.15/erinfmt.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.summary", header=TRUE)

combined.pcoa<-merge(erinsub.meta, pcoa, by.x=c("seqID"), by.y=c("sampleID"))
combined.nmds<-merge(combined.pcoa, nmds, by.x=c("seqID"), by.y=c("sampleID"))
combined.sum<-merge(combined.nmds, sum, by.x=c("seqID"), by.y=c("group"))
combined.all<-combined.sum[,c(1:26, 254:256, 259:273)]  #this picks the specific columns
#write.table(combined.all, 'erinsubset_summary.txt',quote=FALSE,sep="\t", col.names=NA)

combined.all<-read.table(file="erinsubset_summary.txt", header=TRUE)
extra<-read.table(file="erinfmt_extra.meta.txt", header=TRUE)
filtered.extra<-extra[, c("seqID", "abx_prior", "abx_typeprior", "ppi", "final_plating", "severeidsa")]

community<-merge(combined.all, filtered.extra, by.x="seqID", by.y="seqID", all.x=TRUE)
#write.table(community, file="erinsub_all.extrameta.summary.txt", sep="\t", quote=FALSE, col.names=NA)

# color vectors
group.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "recurrent" ) {
	colorvec[i] = "orange"
	}
	if ( n[i] == "nonrecurrent" ) {
	colorvec[i] = "chartreuse3"
	}
	if ( n[i] == "reinfection" ) {
	colorvec[i] = "darkgoldenrod"
	}
	}
	c(colorvec)
	}

clin.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "positive" ) {
	colorvec[i] = "magenta"
	}
	if ( n[i] == "negative" ) {
	colorvec[i] = "grey87"
	}
	}
	c(colorvec)
	}

par(mfrow=c(1,2))
# PCOA:
df<-community
par(mfrow=c(1,2))
plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.55,0.55), xlim=c(-0.55,0.55), main="PCOA: by patient group")
legend("topright",legend=c("recurrent", "nonrecurrent", "reinfection"), col="black", pt.bg=c("orange", "chartreuse3", "darkgoldenrod"), cex=0.7, pch=21)
text(0.38,-0.55, labels="AMOVA, p=0.107", cex=0.7)

plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=clin.col(df$clinical_result), cex=1, xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.55,0.55), xlim=c(-0.55,0.55), main="PCOA: by sample clinical result")
legend("topright",legend=c("positive", "negative"), col="black", pt.bg=c("magenta", "grey87"), cex=0.7, pch=21)
text(0.38,-0.55, labels="**AMOVA, p=0.002", cex=0.7)

#NMDS:
df<-community
plot(df$nmds_axis2~df$nmds_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, xlab="NMDS 1", ylab="NMDS 2", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.65,0.65), xlim=c(-0.65,0.65), main="nmds: by patient group")
legend("topleft",legend=c("recurrent", "nonrecurrent", "reinfection"), col="black", pt.bg=c("orange", "chartreuse3", "darkgoldenrod"), cex=0.7, pch=21)
text(0.42,-0.65, labels="AMOVA, p=0.107", cex=0.7)

plot(df$nmds_axis2~df$nmds_axis1, pch=21, col="black", bg=clin.col(df$clinical_result), cex=1, xlab="NMDS 1", ylab="NMDS 2", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.65,0.65), xlim=c(-0.65,0.65), main="nmds: by sample clinical result")
legend("topleft",legend=c("positive", "negative"), col="black", pt.bg=c("magenta", "grey87"), cex=0.7, pch=21)
text(0.42,-0.65, labels="**AMOVA, p=0.002", cex=0.7)

####
# However, let's do an AMOVA just to see if it separates statistically:
erinsubset.group.design<-combined.all[, c("seqID", "group")]
erinsubset.group.design<-rename(erinsubset.group.design, replace=c("seqID"="group", "group"="CDI"))
#write.table(erinsubset.group.design, file="erinsubset.group.design", sep="\t", quote=FALSE, col.names=NA)

# AMOVA results: NOPE
#nonrecurrent-recurrent	Among	Within	Total
#SS	0.576944	92.9168	93.4937
#df	1	228	229
#MS	0.576944	0.40753
#Fs:	1.41571
#p-value: 0.143

# how about by sample type:
erinsubset.clinical.design<-combined.all[, c("seqID", "clinical_result")]
erinsubset.clinical.design<-rename(erinsubset.clinical.design, replace=c("seqID"="group"))
write.table(erinsubset.clinical.design, file="erinsubset.clinical.design", sep="\t", quote=FALSE, col.names=NA)
	# this actually DID separate significantly...

####

## corr.axes to pcoa or nmds plots:
# produced corr.axes files based on erinsub samples only:
	# erinsubset.0.03.spearman.pcoa.corr.axes
		# filter: p<0.0005 for both axes, length>0.5 = 12 OTUs
	# erinsubset.0.03.spearman.nmds.corr.axes

library(shape)

corr.pcoa<-read.table(file="erinsubset_mothurfiles/pcoa.nmds/erinsubset.0.03.spearman.pcoa.corr.axes", header=TRUE)
tax<-read.table(file="erinfmt.new.taxonomy.names.txt", header=TRUE)
	corr.pcoa<-merge(corr.pcoa, tax, by="OTU", all.y=TRUE)
	
corr.nmds<-read.table(file="erinsubset_mothurfiles/pcoa.nmds/erinsubset.0.03.spearman.nmds.corr.axes", header=TRUE)
tax<-read.table(file="erinfmt.new.taxonomy.names.txt", header=TRUE)
	corr.nmds<-merge(corr.nmds, tax, by="OTU", all.y=TRUE)

	
# subset to get significant OTUs:
corr.pcoa1<-subset(corr.pcoa, p.value < 0.001)
dim(corr.pcoa1)										#78 OTUs
corr.pcoa2<-subset(corr.pcoa1, p.value.1 < 0.001)			#50 OTUs
corr.pcoa2.abundant<-subset(corr.pcoa2, Size > 100000 & length > 0.4)		# 4 OTUs!

corr.nmds1<-subset(corr.nmds, p.value < 0.001)
dim(corr.nmds1)										#142 OTUs
corr.nmds2<-subset(corr.nmds1, p.value.1 < 0.001)			#2 OTUs
corr.nmds2.abundant<-subset(corr.nmds2, Size > 100000 & length > 0.4)		# 2 OTUs!

# PCOA:
par(mfrow=c(1,2))
plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, xlab="NMDS 1", ylab="NMDS 2", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.7,0.7), xlim=c(-0.7,0.7), main="nmds: by patient group")
Arrows(0, 0, x1=0.637056, y1=0.348986, lty=1, arr.length=0.3, arr.type="triangle")
text(0.637056, 0.348986, label="Otu2: \nEnterococcus", cex=.6, pos=2)

Arrows(0, 0, x1=-0.560707, y1=-0.479283, lty=1, arr.length=0.3, arr.type="triangle")
text(-0.560707, -0.479283, label="Otu3 & Otu6: \nBacteroides", cex=.6, pos=1)

Arrows(0, 0, x1=-0.516695, y1=0.570671, lty=1, arr.length=0.3, arr.type="triangle")
text(-0.516695, 0.570671, label="Otu4: \nEscherichia", cex=.6, pos=4)

Arrows(0, 0, x1=-0.414691, y1=-0.349732, lty=1, arr.length=0.3, arr.type="triangle")
#text(-0.414691, -0.349732, label="Otu6: \nBacteroides", cex=.6, pos=1)

# NMDS:
par(mfrow=c(1,2))
plot(df$nmds_axis2~df$nmds_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, xlab="NMDS 1", ylab="NMDS 2", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.75,0.75), xlim=c(-0.75,0.75), main="nmds: by patient group")
Arrows(0, 0, x1=0.404481, y1=0.277719, lty=1, arr.length=0.3, arr.type="triangle")
text(0.404481, 0.277719, label="Otu1: \nEnterobacteriaceae", cex=.6, pos=4)

Arrows(0, 0, x1=0.6412387, y1=-0.412713, lty=1, arr.length=0.3, arr.type="triangle")
text(0.641238, -0.412713, label="Otu2: \nEnterococcus", cex=.6, pos=1)


#### Used for Figure 4A,B:
# both PCOA w/axes:
par(mfrow=c(1,2))
plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.7,0.7), xlim=c(-0.7,0.7), main="PCOA: by patient group")
Arrows(0, 0, x1=0.637056, y1=0.348986, lty=1, arr.length=0.3, arr.type="triangle")
	text(0.637056, 0.348986, label="Otu2: \nEnterococcus", cex=.6, pos=2)
Arrows(0, 0, x1=-0.560707, y1=-0.479283, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.560707, -0.479283, label="Otu3 & Otu6: \nBacteroides", cex=.6, pos=1)
Arrows(0, 0, x1=-0.516695, y1=0.570671, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.516695, 0.570671, label="Otu4: \nEscherichia", cex=.6, pos=4)
Arrows(0, 0, x1=-0.414691, y1=-0.349732, lty=1, arr.length=0.3, arr.type="triangle")
legend("topright",legend=c("recurrent", "nonrecurrent", "reinfection"), col="black", pt.bg=c("orange", "chartreuse3", "darkgoldenrod"), cex=0.7, pch=21)
	text(0.38,-0.55, labels="AMOVA, p=0.107", cex=0.7)

plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=clin.col(df$clinical_result), cex=1, xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.7,0.7), xlim=c(-0.7,0.7), main="PCOA: by sample clinical status")
Arrows(0, 0, x1=0.637056, y1=0.348986, lty=1, arr.length=0.3, arr.type="triangle")
	text(0.637056, 0.348986, label="Otu2: \nEnterococcus", cex=.6, pos=2)
Arrows(0, 0, x1=-0.560707, y1=-0.479283, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.560707, -0.479283, label="Otu3 & Otu6: \nBacteroides", cex=.6, pos=1)
Arrows(0, 0, x1=-0.516695, y1=0.570671, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.516695, 0.570671, label="Otu4: \nEscherichia", cex=.6, pos=4)
Arrows(0, 0, x1=-0.414691, y1=-0.349732, lty=1, arr.length=0.3, arr.type="triangle")
legend("topright",legend=c("positive", "negative"), col="black", pt.bg=c("magenta", "grey87"), cex=0.7, pch=21)
	text(0.38,-0.55, labels="**AMOVA, p=0.002", cex=0.7)



# both NMDS w/axes:
par(mfrow=c(1,2))
plot(df$nmds_axis2~df$nmds_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, xlab="NMDS 1", ylab="NMDS 2", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.7,0.7), xlim=c(-0.7,0.7), main="NMDS: by patient group")
Arrows(0, 0, x1=0.404481, y1=0.277719, lty=1, arr.length=0.3, arr.type="triangle")
	text(0.404481, 0.277719, label="Otu1: \nEnterobacteriaceae", cex=.6, pos=4)
Arrows(0, 0, x1=0.6412387, y1=-0.412713, lty=1, arr.length=0.3, arr.type="triangle")
	text(0.641238, -0.412713, label="Otu2: \nEnterococcus", cex=.6, pos=1)
legend("topright",legend=c("recurrent", "nonrecurrent", "reinfection"), col="black", pt.bg=c("orange", "chartreuse3", "darkgoldenrod"), cex=0.7, pch=21)
	text(0.38,-0.65, labels="AMOVA, p=0.107", cex=0.7)

plot(df$nmds_axis2~df$nmds_axis1, pch=21, col="black", bg=clin.col(df$clinical_result), cex=1, xlab="NMDS 1", ylab="NMDS 2", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.7,0.7), xlim=c(-0.7,0.7), main="NMDS: by sample clinical status")
Arrows(0, 0, x1=0.404481, y1=0.277719, lty=1, arr.length=0.3, arr.type="triangle")
	text(0.404481, 0.277719, label="Otu1: \nEnterobacteriaceae", cex=.6, pos=4)
Arrows(0, 0, x1=0.6412387, y1=-0.412713, lty=1, arr.length=0.3, arr.type="triangle")
	text(0.641238, -0.412713, label="Otu2: \nEnterococcus", cex=.6, pos=1)
legend("topright",legend=c("positive", "negative"), col="black", pt.bg=c("magenta", "grey87"), cex=0.7, pch=21)
	text(0.38,-0.65, labels="**AMOVA, p=0.002", cex=0.7)