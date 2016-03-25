## Making a dotplot of LefSe-identified OTUs
# 3.16.16
# Anna M. Seekatz

library(shape)

# used files:
	# erinsubset.thetayc0.03.pcoa.axes		# data
	# erinsubset.thetayc0.03.nmds.axes		# data
	# erinfmt.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.summary
	# ERIN.UMFMT_metadata_filtered.txt		# OLD metadata
	# updated_suberin_summary.phylotyped.txt
	# updated_erinsubset_summary.txt 		#created from the mothur files
	
# read in files and merge together:
#meta<-read.table(file="ERIN.UMFMT_metadata.txt", header=TRUE)
meta<-read.table(file="updated_meta_fromkrishna_03.08.16.txt", header=TRUE)
	erinsub.meta<-meta[meta$group_reinfection %in% c("recurrent", "nonrecurrent", "reinfection), ]
pcoa<-read.table(file="erinsubset_mothurfiles/erinsubset.thetayc0.03.pcoa.axes", header=TRUE)
	colnames(pcoa)[1:4] <- c("sampleID", "pcoa_axis1", "pcoa_axis2", "pcoa_axis3")
nmds<-read.table(file="erinsubset_mothurfiles/erinsubset.thetayc0.03.nmds.axes", header=TRUE)
	colnames(nmds)[1:4] <- c("sampleID", "nmds_axis1", "nmds_axis2", "nmds_axis3")
sum<-read.table(file="../mothur.files_3.18.15/erinfmt.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.groups.summary", header=TRUE)

combined.pcoa<-merge(erinsub.meta, pcoa, by.x=c("seqID"), by.y=c("sampleID"))
combined.nmds<-merge(combined.pcoa, nmds, by.x=c("seqID"), by.y=c("sampleID"))
combined.sum<-merge(combined.nmds, sum, by.x=c("seqID"), by.y=c("group"))
combined.all<-combined.sum[,c(1:26, 254:256, 259:273)]  #this picks the specific columns
#write.table(combined.all, 'updated_erinsubset_summary.txt',quote=FALSE,sep="\t", col.names=NA)

#combined.all<-read.table(file="updated_erinsubset_summary.txt", header=TRUE)
#extra<-read.table(file="erinfmt_extra.meta.txt", header=TRUE)
#filtered.extra<-extra[, c("seqID", "abx_prior", "abx_typeprior", "ppi", "final_plating", "severeidsa")]

#community<-merge(combined.all, filtered.extra, by.x="seqID", by.y="seqID", all.x=TRUE)
#write.table(community, file="erinsub_all.extrameta.summary.txt", sep="\t", quote=FALSE, col.names=NA)
# already created from updated data:
### if these steps were already done previously, read in the file:
community<-read.table(file='updated_erinsubset_extra.summary.txt', header=TRUE)
	
	
	# adjust some levels of unknowns:
levels(community$POS_NEG) <- c(levels(community$POS_NEG), "not_tested")
community$POS_NEG[is.na(community$POS_NEG)]<-'not_tested'
# let's also replace some others in the final_plating
levels(community$final_plating) <- c(levels(community$final_plating), "unknown")
community$final_plating[is.na(community$final_plating)]<-'unknown'
community$final_plating[community$final_plating==c("positive_maybe")]<-'unknown'
community$final_plating[community$final_plating==c("No_sample")]<-'unknown'


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
	if ( n[i] == "not_tested" ) {
	colorvec[i] = "grey27"
	}
	}
	c(colorvec)
	}

par(mfrow=c(1,2))
# PCOA:
df<-community
par(mfrow=c(1,2))
plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.55,0.55), xlim=c(-0.55,0.55), main="")
legend("topright",legend=c("recurrent", "nonrecurrent", "reinfection"), col="black", pt.bg=c("orange", "chartreuse3", "darkgoldenrod"), cex=0.7, pch=21)
text(0.38,-0.55, labels="AMOVA, p=0.107", cex=0.7)

plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=clin.col(df$POS_NEG), cex=1, xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.55,0.55), xlim=c(-0.55,0.55), main="")
legend("topright",legend=c("positive", "negative", "unknown"), col="black", pt.bg=c("magenta", "grey87", "grey27"), cex=0.7, pch=21)
text(0.38,-0.55, labels="**AMOVA, p=0.002", cex=0.7)

# PCOA, but ONLY first index samples:
df<-community[community$index_sample_n==1, ]
par(mfrow=c(1,2))
plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.55,0.55), xlim=c(-0.55,0.55), main="")
legend("topright",legend=c("recurrent", "nonrecurrent", "reinfection"), col="black", pt.bg=c("orange", "chartreuse3", "darkgoldenrod"), cex=0.7, pch=21)
text(0.38,-0.55, labels="AMOVA, p=", cex=0.7)

plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=clin.col(df$POS_NEG), cex=1, xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.55,0.55), xlim=c(-0.55,0.55), main="")
legend("topright",legend=c("positive", "negative", "unknown"), col="black", pt.bg=c("magenta", "grey87", "grey27"), cex=0.7, pch=21)
text(0.38,-0.55, labels="**AMOVA, p=", cex=0.7)


#NMDS:
df<-community
plot(df$nmds_axis2~df$nmds_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, xlab="NMDS 1", ylab="NMDS 2", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.65,0.65), xlim=c(-0.65,0.65), main="nmds: by patient group")
legend("topleft",legend=c("recurrent", "nonrecurrent", "reinfection"), col="black", pt.bg=c("orange", "chartreuse3", "darkgoldenrod"), cex=0.7, pch=21)
text(0.42,-0.65, labels="AMOVA, p=0.107", cex=0.7)

plot(df$nmds_axis2~df$nmds_axis1, pch=21, col="black", bg=clin.col(df$clinical_result), cex=1, xlab="NMDS 1", ylab="NMDS 2", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.65,0.65), xlim=c(-0.65,0.65), main="nmds: by sample clinical result")
legend("topleft",legend=c("positive", "negative"), col="black", pt.bg=c("magenta", "grey87"), cex=0.7, pch=21)
text(0.42,-0.65, labels="**AMOVA, p=0.002", cex=0.7)

####
# updated amova: some groups are now significant, woo hoo

# 3.21.16

## comparing groups:
# can use similar files as for lefse:
combined<-read.table(file="updated_erinsubset_summary.txt", header=TRUE)
erinsub<-combined[combined$index_sample_n %in% c(1), ]		# for only index
updated_all.groups.design<-combined[, c("seqID", "group_reinfection")]
updated_all.index.design<-erinsub[, c("seqID", "group_reinfection")]
write.table(updated_all.groups.design, file="UPDATED_amova/updated_all.groups.design", sep="\t", quote=FALSE, col.names=NA)
write.table(updated_all.index.design, file="UPDATED_amova/updated_all.index.design", sep="\t", quote=FALSE, col.names=NA)

$ /Users/annaseekatz/Desktop/umich/mothur_v1.33.3/mothur
	# in mothur:
dist.shared(shared=../UPDATED_lefse/updated_recurrent.shared, calc=thetayc)
system(mv ../UPDATED_lefse/updated_recurrent.thetayc.0.03.lt.dist erinsubset.allgrps.dist)
amova(phylip=erinsubset.allgrps.dist, design=updated_all.groups.design)
				# holy crap! this is now significant... p-value: 0.016*; let's try with only the index samples:
dist.shared(shared=../UPDATED_lefse/recurrent_index.shared, calc=thetayc)
system(mv ../UPDATED_lefse/recurrent_index.thetayc.0.03.lt.dist erinsubset_index.dist)
amova(phylip=erinsubset_index.dist, design=updated_all.index.design)
				# ALMOST significant...0.068, so trending
				
## comparing neg/pos (and plating):
# can use similar files as for lefse:
# let's use the design files we made for lefse:

#(in directory):
$ ln -s ../UPDATED_lefse/updated_plating.design
$ ln -s ../UPDATED_lefse/updated_clinical2.design

$ /Users/annaseekatz/Desktop/umich/mothur_v1.33.3/mothur
	# in mothur:
			# clinical test:
dist.shared(shared=../UPDATED_lefse/updated_clinical2.shared, calc=thetayc)
system(mv ../UPDATED_lefse/updated_clinical2.thetayc.0.03.lt.dist updated_clinical2.thetayc.0.03.lt.dist)
amova(phylip=updated_clinical2.thetayc.0.03.lt.dist, design=updated_clinical2.design)
				# significant!! p-value: 0.015*
			# plating:
dist.shared(shared=../UPDATED_lefse/updated_plating.shared, calc=thetayc)
system(mv ../UPDATED_lefse/updated_plating.thetayc.0.03.lt.dist updated_plating.thetayc.0.03.lt.dist)
amova(phylip=updated_plating.thetayc.0.03.lt.dist, design=updated_plating.design)
				# also significant! p-value: <0.001*

####

## corr.axes to pcoa or nmds plots:
# produced corr.axes files based on erinsub samples only:
	# erinsubset.0.03.spearman.pcoa.corr.axes
		# filter: p<0.0005 for both axes, length>0.5 = 12 OTUs
	# erinsubset.0.03.spearman.nmds.corr.axes

library(shape)

corr.pcoa<-read.table(file="../erinsubset_mothurfiles/pcoa.nmds/erinsubset.0.03.spearman.pcoa.corr.axes", header=TRUE)
tax<-read.table(file="../erinfmt.new.taxonomy.names.txt", header=TRUE)
	corr.pcoa<-merge(corr.pcoa, tax, by="OTU", all.y=TRUE)
	
corr.nmds<-read.table(file="../erinsubset_mothurfiles/pcoa.nmds/erinsubset.0.03.spearman.nmds.corr.axes", header=TRUE)
tax<-read.table(file="../erinfmt.new.taxonomy.names.txt", header=TRUE)
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
plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, xlab="NMDS 1", ylab="NMDS 2", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.7,0.7), xlim=c(-0.7,0.7), main="pcoa: by patient group")
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


#### Used for Figure 3A,B:
# both PCOA w/axes:
par(mfrow=c(2,1))
plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, 
	xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.7,0.7), 
	xlim=c(-0.7,0.7))
Arrows(0, 0, x1=0.637056, y1=0.348986, lty=1, arr.length=0.3, arr.type="triangle")
	text(0.637056, 0.348986, label="Otu2: \nEnterococcus", cex=.6, pos=2)
Arrows(0, 0, x1=-0.560707, y1=-0.479283, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.560707, -0.479283, label="Otu3 & Otu6: \nBacteroides", cex=.6, pos=1)
Arrows(0, 0, x1=-0.516695, y1=0.570671, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.516695, 0.570671, label="Otu4: \nEscherichia", cex=.6, pos=4)
Arrows(0, 0, x1=-0.414691, y1=-0.349732, lty=1, arr.length=0.3, arr.type="triangle")
legend("topright",legend=c("recurrent", "nonrecurrent", "reinfection"), col="black", pt.bg=c("orange", "chartreuse3", "darkgoldenrod"), cex=0.7, pch=21)
	text(0.38,-0.55, labels="*AMOVA, p=0.016", cex=0.7)

plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=clin.col(df$POS_NEG), cex=1, 
	xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.7,0.7), xlim=c(-0.7,0.7))
Arrows(0, 0, x1=0.637056, y1=0.348986, lty=1, arr.length=0.3, arr.type="triangle")
	text(0.637056, 0.348986, label="Otu2: \nEnterococcus", cex=.6, pos=2)
Arrows(0, 0, x1=-0.560707, y1=-0.479283, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.560707, -0.479283, label="Otu3 & Otu6: \nBacteroides", cex=.6, pos=1)
Arrows(0, 0, x1=-0.516695, y1=0.570671, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.516695, 0.570671, label="Otu4: \nEscherichia", cex=.6, pos=4)
Arrows(0, 0, x1=-0.414691, y1=-0.349732, lty=1, arr.length=0.3, arr.type="triangle")
legend("topright",legend=c("positive", "negative", "unknown"), col="black", pt.bg=c("magenta", "grey87", "grey27"), cex=0.7, pch=21)
	text(0.38,-0.55, labels="clinical lab: *AMOVA, p=0.015", cex=0.7)
	#text(0.38,-0.6, labels="cultivation: **AMOVA, p=0.001", cex=0.7)


#### Used for Figure S2A,B:
# both PCOA w/axes:
df<-community[community$index_sample_n==1, ]
par(mfrow=c(1,2))
plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=group.col(df$group_reinfection), cex=1, 
	xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.7,0.7), 
	xlim=c(-0.7,0.7))
Arrows(0, 0, x1=0.637056, y1=0.348986, lty=1, arr.length=0.3, arr.type="triangle")
	text(0.637056, 0.348986, label="Otu2: \nEnterococcus", cex=.6, pos=2)
Arrows(0, 0, x1=-0.560707, y1=-0.479283, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.560707, -0.479283, label="Otu3 & Otu6: \nBacteroides", cex=.6, pos=1)
Arrows(0, 0, x1=-0.516695, y1=0.570671, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.516695, 0.570671, label="Otu4: \nEscherichia", cex=.6, pos=4)
Arrows(0, 0, x1=-0.414691, y1=-0.349732, lty=1, arr.length=0.3, arr.type="triangle")
legend("topright",legend=c("recurrent", "nonrecurrent", "reinfection"), col="black", pt.bg=c("orange", "chartreuse3", "darkgoldenrod"), cex=0.7, pch=21)
	text(0.38,-0.55, labels="AMOVA, p=0.068", cex=0.7)

# also, by plating:
df<-community
plate.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "positive" ) {
	colorvec[i] = "magenta"
	}
	if ( n[i] == "negative" ) {
	colorvec[i] = "grey87"
	}
	if ( n[i] == "unknown" ) {
	colorvec[i] = "grey27"
	}
	}
	c(colorvec)
	}
plot(df$pcoa_axis2~df$pcoa_axis1, pch=21, col="black", bg=plate.col(df$final_plating), cex=1, 
	xlab="PCOA 1 (12.4%)", ylab="PCOA 2 (9.9%)", cex.lab=0.8, cex.axis=0.8, ylim=c(-0.7,0.7), xlim=c(-0.7,0.7))
Arrows(0, 0, x1=0.637056, y1=0.348986, lty=1, arr.length=0.3, arr.type="triangle")
	text(0.637056, 0.348986, label="Otu2: \nEnterococcus", cex=.6, pos=2)
Arrows(0, 0, x1=-0.560707, y1=-0.479283, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.560707, -0.479283, label="Otu3 & Otu6: \nBacteroides", cex=.6, pos=1)
Arrows(0, 0, x1=-0.516695, y1=0.570671, lty=1, arr.length=0.3, arr.type="triangle")
	text(-0.516695, 0.570671, label="Otu4: \nEscherichia", cex=.6, pos=4)
Arrows(0, 0, x1=-0.414691, y1=-0.349732, lty=1, arr.length=0.3, arr.type="triangle")
legend("topright",legend=c("positive", "negative", "unknown"), col="black", pt.bg=c("magenta", "grey87", "grey27"), cex=0.7, pch=21)
	text(0.38,-0.6, labels="cultivation: **AMOVA, p=0.001", cex=0.7)




