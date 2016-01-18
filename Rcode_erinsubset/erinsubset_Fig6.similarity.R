## Making diversity plots (index single-point comparison AND over time)
# 1.7.16
# Anna M. Seekatz

# files used:
	# erinsubset.0.03.summary		#note: you have to fix the columns before reading into R
	# erinsub_all.extrameta.summary.txt
	# erinsub_ind.dist.byindex.txt (created)
	# erinsub_shared.meta.dist.txt (created)
	# erinsub_ind.dist.overtime.txt (created)
	
library(shape)
library(plyr)

# have to recalculate a summary.shared file from the .dist file created earlier:
#(in mothur):
	# summary.shared(shared=/Users/annaseekatz/Desktop/umich/projects/ERIN_clinical/erinfmt_analysis_athome/ERINsubset/erinsubset_mothurfiles/erinsubset.0.03.shared, calc=sharedsobs-braycurtis-spearman-thetayc-jsd-sharednseqs)

## create a shared distance file that incorporates the metadata for both sample columns:

# merge shared files with meta (includes extra meta):
var<-read.table(file="erinsub_all.extrameta.summary.txt", header=TRUE)
mdist<-read.table(file="erinsubset_mothurfiles/erinsubset.0.03.summary", header=TRUE)

#add sample information for 'group1':
m<-merge(var, mdist, by.x=c("seqID"), by.y=c("s1"))				#merge for sample1 metadata 			
m.2<-m[,1:49]													#separate to add a s1 to colnames
m.3<-m[,50:58]
colnames(m.2) <- paste("s1", colnames(m.2), sep = "_")
m<-cbind(m.2,m.3)

#add sample information for 2nd sample: ('group2'):
m2<-merge(var, m, by.x=c("seqID"), by.y=c("s2"))				#merge for sample2 metadata 			
m2.2<-m2[,1:49]													#separate to add a s1 to colnames
m2.3<-m2[,50:106]
colnames(m2.2) <- paste("s2", colnames(m2.2), sep = "_")
m.shared<-cbind(m2.2,m2.3)

#write.table(m.shared, file="erinsub_shared.meta.dist.txt", quote=FALSE, sep="\t")

## next, we need to eliminate some sample comparisons
	# in this case, we are ONLY doing within patient comparisons
	# you do not need comparisons between patient x and patient y, just samples within patient x

# filter distances ONLY within a patient:
# first, we will split the data so that we ONLY keep distances between the same patientID:
ind<-m.shared[(m.shared$s1_patientID==m.shared$s2_patientID), ]

test<-ind[,c("s2_total_sample_n", "s2_patientID", "s1_patientID", "s2_sample_n", "s1_sample_n")]			
	#check if this looks good....
#write.table(test, file="test_dist.within.individual.txt", quote=FALSE, sep="\t", col.names=NA)
length(unique(ind$s1_patientID))
	#[1] 69		#this appears to be right!
#write.table(ind, file="erinsub_dist.within.individual.txt", quote=FALSE, sep="\t", col.names=NA)

# you can separate your samples by total n samples collected if you wanted to concentrate on particular samplings:
	summary(as.factor(ind$s1.total_sample_n))							#how many sampling levels do you have?
	two<-ind[(ind$s1_total_sample_n=="2" & ind$s2_total_sample_n=="2"), ]
	three<-ind[(ind$s1_total_sample_n=="3" & ind$s2_total_sample_n=="3"), ]
	four<-ind[(ind$s1_total_sample_n=="4" & ind$s2_total_sample_n=="4"), ]
	five<-ind[(ind$s1_total_sample_n=="5" & ind$s2_total_sample_n=="5"), ]
	six<-ind[(ind$s1_total_sample_n=="6" & ind$s2_total_sample_n=="6"), ]
	seven<-ind[(ind$s1_total_sample_n=="7" & ind$s2_total_sample_n=="7"), ]
		#this can be used if you need to split each of the samplings
		#you could also work from the large file (ind) itself
		# for instance, you can see how many patients are in each sampling group:
	length(unique(three$s1_patientID))
		# 21 patients with 3 samplings

# now, for each level, let's make sure that we are only doing within-sample comparisons-sample1 compared to sample 2, etc:
#three[,c("s2.seqID", "s1.seqID", "s2.sample_n", "s1.sample_n", "thetayc", "s1.patientID", "s2.patientID")]
	#can subset, but let's just include all of the info just in case we want to do another type of comparison

## now, let's separate the data into ONLY comparisons between consecutive samplings
	# this requires some work (apologies for the long code)
	# we will also label each comparison consecutively ('sampling' column)

# differences between sampling 1 and 2:
ind.1a<-ind[(ind$s1_sample_n=="1" & ind$s2_sample_n=="2"), ]
ind.1b<-ind[(ind$s1_sample_n=="2" & ind$s2_sample_n=="1"), ]		#none this way
ind1<-rbind(ind.1a, ind.1b)
ind1$sampling<- c(1)

# differences between sampling 2 and 3:
ind.2a<-ind[(ind$s1_sample_n=="3" & ind$s2_sample_n=="2"), ]
ind.2b<-ind[(ind$s1_sample_n=="2" & ind$s2_sample_n=="3"), ]		
ind2<-rbind(ind.2a, ind.2b)
ind2$sampling<- c(2)

# differences between sampling 3 and 4:
ind.3a<-ind[(ind$s1_sample_n=="3" & ind$s2_sample_n=="4"), ]
ind.3b<-ind[(ind$s1_sample_n=="4" & ind$s2_sample_n=="3"), ]		
ind3<-rbind(ind.3a, ind.3b)
ind3$sampling<- c(3)

# differences between sampling 4 and 5:
ind.4a<-ind[(ind$s1_sample_n=="5" & ind$s2_sample_n=="4"), ]
ind.4b<-ind[(ind$s1_sample_n=="4" & ind$s2_sample_n=="5"), ]		
ind4<-rbind(ind.4a, ind.4b)
ind4$sampling<- c(4)

# differences between sampling 5 and 6:
ind.5a<-ind[(ind$s1_sample_n=="5" & ind$s2_sample_n=="6"), ]
ind.5b<-ind[(ind$s1_sample_n=="6" & ind$s2_sample_n=="5"), ]		
ind5<-rbind(ind.5a, ind.5b)
ind5$sampling<- c(5)

# differences between sampling 6 and 7:
ind.6a<-ind[(ind$s1_sample_n=="7" & ind$s2_sample_n=="6"), ]
ind.6b<-ind[(ind$s1_sample_n=="6" & ind$s2_sample_n=="7"), ]		
ind6<-rbind(ind.6a, ind.6b)
ind6$sampling<- c(6)

ind.overtime<-rbind(ind1, ind2, ind3, ind4, ind5, ind6)
#write.table(ind.overtime, file="erinsub_ind.dist.overtime.txt", quote=FALSE, sep="\t", col.names=NA)

# plot it:

group.col<-function(n) {
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
# summary stats:
tyc.data <- ddply(ind.overtime, c("sampling", "s1_group_reinfection"), summarise,
               N    = length(thetayc),
               mean = mean(thetayc),
               sd   = sd(thetayc),
               se   = sd / sqrt(N),
               median = median(thetayc),
               lq	 = quantile(thetayc, 0.25, na.rm=TRUE),
			   hq	 = quantile(thetayc, 0.75, na.rm=TRUE) )

	# if you want to eliminate samplings past 4
split.tyc<-split(tyc.data, tyc.data$s1_group_reinfection)
tyc.rec<-split.tyc$'recurrent'
tyc.non<-split.tyc$'nonrecurrent'
tyc.reinf<-split.tyc$'reinfection'

# in looking at the total (n) per sampling per group, we should eliminate all points after 4 since after that there are 2 or less patients per group
tyc.data<-tyc.data[tyc.data$sampling %in% c(1,2,3),]
split.tyc<-split(tyc.data, tyc.data$s1_group_reinfection)
tyc.rec<-split.tyc$'recurrent'
tyc.non<-split.tyc$'nonrecurrent'
tyc.reinf<-split.tyc$'reinfection'

# plot it!
## Fig. 6A:
x<-tyc.rec$sampling
d<-0.15
	# this offsets your points so they don't overlap
plot(x, tyc.rec$mean,
    ylim=c(0,1.2),
    pch=21, xlab="", ylab=expression(paste("", theta, "yc  distance")),
    xaxt='n', col="black", bg="orange", cex=0.8, xlim=c(0.75, 3.25), cex.lab=0.8, cex.axis=0.8)
axis(side=2, labels="(community dissimilarity)", at=0.55, line=0.5, lwd=0, cex.axis=0.8, mgp=c(1,1.8,0))
arrows(x, tyc.rec$mean-tyc.rec$sd, x, tyc.rec$mean+tyc.rec$sd, length=0.05, angle=90, code=3, col="orange")
lines(x, tyc.rec$mean, col="orange", pch=19)
points(x, tyc.rec$mean, col="black", bg="orange", pch=21, cex=0.8)

arrows(tyc.non$sampling-d, tyc.non$mean-tyc.non$sd, tyc.non$sampling-d, 
	tyc.non$mean+tyc.non$sd, length=0.05, angle=90, code=3, col="chartreuse3")
lines(tyc.non$sampling-d, tyc.non$mean, col="chartreuse3", pch=19)
points(tyc.non$sampling-d, tyc.non$mean, col="black", bg="chartreuse3", pch=21, cex=0.8)

arrows(tyc.reinf$sampling-2*d, tyc.reinf$mean-tyc.reinf$sd, tyc.reinf$sampling-2*d, 
	tyc.reinf$mean+tyc.reinf$sd, length=0.05, angle=90, code=3, col="darkgoldenrod")
lines(tyc.reinf$sampling-2*d, tyc.reinf$mean, col="darkgoldenrod", pch=19)
points(tyc.reinf$sampling-2*d, tyc.reinf$mean, col="black", bg="darkgoldenrod", pch=21, cex=0.8)
		   
legend("bottomleft",legend=c("nonrecurrent","recurrent", "reinfected"), 
	col=c("chartreuse3", "orange", "darkgoldenrod"), cex=0.7, pch=19)
axis(side=1, at=c(1-d, 2-d, 3-d), labels=c('s1 vs. s2', 's2 vs. s3', 's3 vs. s4'), 
	line=0.5, lwd=0, mgp=c(1,0,0), cex.axis=0.8)
#text(c(1.5), c(0.6), c('*FMT'), cex=0.8)	# can add text if necessary (such as n?)
#title('Change in community composition between samplings')

### Community dissimilarity during clinical status change:
	# using the same files as before, we will not split the files according to clinical status change

# first, split the samples in the original ind file by group:
ind.all<-split(ind, ind$s1_group_reinfection)
ind.reinf<-ind.all$'reinfection'
ind.non<-ind.all$'nonrecurrent'
ind.rec<-ind.all$'recurrent'

# within group, ID when there is a change in index status:
# note: this is regardless of sampling time, just change between + and - samples WITHIN a patient

# recurrent:
rec0<-ind.rec[(ind.rec$s1_index==ind.rec$s2_index), ]
rec0$delta<-c(0)

rec1<-ind.rec[(ind.rec$s1_index!=ind.rec$s2_index), ]
rec1$delta<-c(1)

rec<-rbind(rec1, rec0)

# you could also get the sim. between positive samples or negative samples only within a patient:
	allrec1<-ind.rec[(ind.rec$s1_index==1 & ind.rec$s2_index==1), ]
	allrec1$STATUS<-c("positive")
	allrec0<-ind.rec[(ind.rec$s1_index==0 & ind.rec$s2_index==0), ]
	allrec0$STATUS<-c("negative")
	allrecs<-ind.rec[(ind.rec$s1_index!=ind.rec$s2_index), ]
	allrecs$STATUS<-c("switch")
	allrec<-rbind(allrec0, allrec1, allrecs)

# nonrecurrent:
non0<-ind.non[(ind.non$s1_index==ind.non$s2_index), ]
non0$delta<-c(0)

non1<-ind.non[(ind.non$s1_index!=ind.non$s2_index), ]
non1$delta<-c(1)

non<-rbind(non1, non0)

# you could also get the sim. between positive samples or negative samples only within a patient:
	allnon1<-ind.non[(ind.non$s1_index==1 & ind.non$s2_index==1), ]
	allnon1$STATUS<-c("positive")
	allnon0<-ind.non[(ind.non$s1_index==0 & ind.non$s2_index==0), ]
	allnon0$STATUS<-c("negative")
	allnons<-ind.non[(ind.non$s1_index!=ind.non$s2_index), ]
	allnons$STATUS<-c("switch")
	allnon<-rbind(allnon0, allnon1, allnons)

# reinfurrent:
reinf0<-ind.reinf[(ind.reinf$s1_index==ind.reinf$s2_index), ]
reinf0$delta<-c(0)

reinf1<-ind.reinf[(ind.reinf$s1_index!=ind.reinf$s2_index), ]
reinf1$delta<-c(1)

reinf<-rbind(reinf1, reinf0)

# you could also get the sim. between positive samples or negative samples only within a patient:
	allreinf1<-ind.reinf[(ind.reinf$s1_index==1 & ind.reinf$s2_index==1), ]
	allreinf1$STATUS<-c("positive")
	allreinf0<-ind.reinf[(ind.reinf$s1_index==0 & ind.reinf$s2_index==0), ]
	allreinf0$STATUS<-c("negative")
	allreinfs<-ind.reinf[(ind.reinf$s1_index!=ind.reinf$s2_index), ]
	allreinfs$STATUS<-c("switch")
	allreinf<-rbind(allreinf0, allreinf1, allreinfs)

ind.byindex<-rbind(rec, non, reinf)
#write.table(ind.byindex, file="erinsub_ind.dist.byindex.txt", quote=FALSE, sep="\t", col.names=NA)

# now plot it:

## Fig. 6B:
x<-ind.byindex$delta
plot(x, ind.byindex$thetayc, col=grouped.col(ind.byindex$s1.group), pch=19)
# ind.byindex$GROUP<-paste("s1_group_reinfection", "delta", sep="_")

data <- ddply(ind.byindex, c("delta", "s1_group_reinfection"), summarise,
               N    = length(thetayc),
               mean.tyc = mean(thetayc),
               mean.jsd = mean(jsd),
               mean.bc = mean(braycurtis),
               mean.sobs = mean(sharedsobs),
               sd   = sd(thetayc),
               se   = sd / sqrt(N) )
# this might actually look good!!

d.split<-split(data, data$s1_group_reinfection)
d.rec<-d.split$'recurrent'
d.non<-d.split$'nonrecurrent'
d.reinf<-d.split$'reinfection'

# let's just plot all of the samples, without separating them into groups:
# this is significant (see below)
boxplot(thetayc~delta, data = ind.byindex, cex=0.6, xaxt='n', cex.lab=0.8, 
	ylab=expression(paste("", theta, "yc distance")), cex.lab=0.8, cex.axis=0.8)
axis(side=1, at=c(1, 2), labels=c('No change \n(CDI status)', 'Index change \n(CDI status)'), 
	line=0.5, lwd=0, cex.axis=0.8, mgp=c(1,0.5,0))
axis(side=3, at=c(1.5), labels=c("Change in community composition \nduring CDI diagnosis \n(all groups together)"), 
	line=0.5, lwd=0, cex.axis=0.8, mgp=c(1,0.5,0))

# making the graph look better, grouping by multiple factors:
# option 1: by major category of change vs. no change:
ind.byindex$s1_group<-droplevels(ind.byindex$s1_group)				#you must drop the ghost levels
boxplot.double = boxplot(thetayc~s1_group_reinfection + delta, data = ind.byindex, at = c(1, 2, 3, 5, 6, 7), 
	ylim = c(0,1), col = c('chartreuse3', 'orange', 'darkgoldenrod'), xlim=c(0,8), 
	ylab=expression(paste("", theta, "yc distance")), las=2, mgp=c(0.5,1,1), xaxt='n', xlab='', cex.axis=0.8)
legend("bottomright", c("nonrecurrent", "recurrent", "reinfected"), col=c('chartreuse3', 'orange', 'darkgoldenrod'), pch=15, cex=0.7)
axis(side=2, labels="(community dissimilarity)", at=0.55, line=0.5, lwd=0, cex.axis=0.8, mgp=c(1,1.8,0))
axis(side=1, at=c(2, 6), labels=c('No change \n(CDI status)', 'Index change \n(CDI status)'), line=0.5, lwd=0, cex.axis=0.8, mgp=c(1,0.5,0))
title('Change in community composition \nduring CDI diagnosis')
	# communities are more similar when there is NO change in clinical status 
	# significant by delta in kruskal test (also wilcox)

# option 2: by major category of groups:
ind.byindex$s1_group<-droplevels(ind.byindex$s1_group)				#you must drop the ghost levels
boxplot.double = boxplot(thetayc~delta+s1_group_reinfection, data = ind.byindex, at = c(1, 2, 4, 5, 7, 8), 
	ylim = c(0,1.2), col = c(rep('chartreuse3', 2), rep('orange', 2), rep('darkgoldenrod', 2)), xlim=c(0,9), 
	ylab=expression(paste("", theta, "yc distance")), las=2, mgp=c(0.5,1,1), xaxt='n', xlab='', cex.axis=0.8, cex=0.6, cex.lab=0.8)
legend("bottomright", c("nonrecurrent", "recurrent", 'reinfection'), col=c('chartreuse3', 'orange', 'darkgoldenrod'), pch=15, cex=0.6)
axis(side=2, labels="(community dissimilarity)", at=0.55, line=0.5, lwd=0, cex.axis=0.8, mgp=c(1,1.8,0))
axis(side=1, at=c(1, 2, 4, 5, 7, 8), labels=rep(c('-', '+'),3), line=0.5, lwd=0, cex.axis=1.5, mgp=c(1,0.5,0))
#title('Change in community composition \nduring CDI diagnosis')
mtext(expression(paste("", Delta, " in")), 1, cex=0.8, at=c(-1,-1), padj=0.6)
mtext("clinical \nstatus:", 1, cex=0.8, at=c(-1,-1), padj=1.5)
#segments(x0=0.8, y0=1.05, x1 = 2.2, y1 = 1.05, lwd=1.5)
#segments(x0=3.8, y0=1.05, x1 = 5.2, y1 = 1.05, lwd=1.5)
#segments(x0=6.8, y0=1.05, x1 = 8.2, y1 = 1.05, lwd=1.5)
	#useful, but we can do better:
Arrows(0.8, 1.05, 2.2, 1.05, code = 3, arr.type = "T", arr.length=0.2)
Arrows(3.8, 1.05, 5.2, 1.05, code = 3, arr.type = "T", arr.length=0.2)
Arrows(6.8, 1.05, 8.2, 1.05, code = 3, arr.type = "T", arr.length=0.2)
text(1.5,1.1, labels="** p=0.0014", cex=0.6)
text(4.5,1.1, labels="ns", cex=0.6)
text(7.5,1.1, labels="ns", cex=0.6)

# stats:
	# couldn't get it to separate by two variables
kruskal.test(thetayc~delta, data=ind.byindex)
	#data:  thetayc by delta
	#Kruskal-Wallis chi-squared = 3.9912, df = 1, p-value = 0.04574
	# but, significant if only looking at change (even if including all groups)
# if using option 2, can also do stats within group:
# had previously done this: ind.byindex<-rbind(rec, non, reinf)
wilcox.test(thetayc~delta, data=non)
	#W = 764, p-value = 0.00232
	# woo hoo!
wilcox.test(thetayc~delta, data=rec)
	#W = 1328, p-value = 0.3492
	# not sign. awesome!
wilcox.test(thetayc~delta, data=reinf)
	#W = 66, p-value = 0.05776
	# not sig--ok, since there are only a couple points anyhow...
#### stats between neg/nonrecurrent and neg/recurrent:
ind.byindex$GROUP<-paste(ind.byindex$s1_group_reinfection, ind.byindex$delta, sep="_")
test<-ind.byindex[ind.byindex$GROUP %in% c("recurrent_1", "nonrecurrent_1"), ]
test<-droplevels(test)
wilcox.test(thetayc~GROUP, data=test)
	#Wilcoxon rank sum test with continuity correction
	#data:  thetayc by GROUP
	#W = 2054, p-value = 0.0194
	#alternative hypothesis: true location shift is not equal to 0
test2<-ind.byindex[ind.byindex$GROUP %in% c("recurrent_0", "nonrecurrent_0"), ]
test2<-droplevels(test2)
wilcox.test(thetayc~GROUP, data=test2)
	#Wilcoxon rank sum test with continuity correction
	#data:  thetayc by GROUP
	#W = 677, p-value = 0.07288
	#alternative hypothesis: true location shift is not equal to 0
test3<-ind.byindex[ind.byindex$GROUP %in% c("recurrent_1", "nonrecurrent_0"), ]
test3<-droplevels(test3)
wilcox.test(thetayc~GROUP, data=test3)
	#Wilcoxon rank sum test with continuity correction
	#data:  thetayc by GROUP
	#W = 859, p-value = 0.2755
	#alternative hypothesis: true location shift is not equal to 0
test4<-ind.byindex[ind.byindex$GROUP %in% c("recurrent_0", "nonrecurrent_1"), ]
test4<-droplevels(test4)
wilcox.test(thetayc~GROUP, data=test4)
	#Wilcoxon rank sum test with continuity correction
	#data:  thetayc by GROUP
	#W = 1599, p-value = 0.359
	#alternative hypothesis: true location shift is not equal to 0



	
### both graphs together:
par(mfrow=c(1,2))
par(mar=c(4,4.1,1.0,1.1))
#layout(matrix(c(1,1,0,2), 2, 2, byrow=TRUE), respect=TRUE)
	# A: over time
x<-tyc.rec$sampling
d<-0.15
plot(x, tyc.rec$mean,
    ylim=c(0,1.2),
    pch=21, xlab="", ylab=expression(paste("", theta, "yc  distance")),
    xaxt='n', col="black", bg="orange", cex=0.8, xlim=c(0.5, 3.25), cex.lab=0.8, cex.axis=0.8)
axis(side=2, labels="(community dissimilarity)", at=0.55, line=0.5, lwd=0, cex.axis=0.8, mgp=c(1,1.8,0))
#axis(side=1, labels="sampling", at=0.55, line=0.5, lwd=0, cex.axis=0.8, mgp=c(0,1.8,2))
arrows(x, tyc.rec$mean-tyc.rec$sd, x, tyc.rec$mean+tyc.rec$sd, length=0.05, angle=90, code=3, col="orange")
lines(x, tyc.rec$mean, col="orange", pch=19)
points(x, tyc.rec$mean, col="black", bg="orange", pch=21, cex=0.8)

arrows(tyc.non$sampling-d, tyc.non$mean-tyc.non$sd, tyc.non$sampling-d, 
	tyc.non$mean+tyc.non$sd, length=0.05, angle=90, code=3, col="chartreuse3")
lines(tyc.non$sampling-d, tyc.non$mean, col="chartreuse3", pch=19)
points(tyc.non$sampling-d, tyc.non$mean, col="black", bg="chartreuse3", pch=21, cex=0.8)

arrows(tyc.reinf$sampling-2*d, tyc.reinf$mean-tyc.reinf$sd, tyc.reinf$sampling-2*d, 
	tyc.reinf$mean+tyc.reinf$sd, length=0.05, angle=90, code=3, col="darkgoldenrod")
lines(tyc.reinf$sampling-2*d, tyc.reinf$mean, col="darkgoldenrod", pch=19)
points(tyc.reinf$sampling-2*d, tyc.reinf$mean, col="black", bg="darkgoldenrod", pch=21, cex=0.8)
legend("bottomleft",legend=c("nonrecurrent","recurrent", "reinfected"), 
	col=c("chartreuse3", "orange", "darkgoldenrod"), cex=0.7, pch=19)
axis(side=1, at=c(1-d, 2-d, 3-d), labels=c('s1 vs. s2', 's2 vs. s3', 's3 vs. s4'), 
	line=0.5, lwd=0, mgp=c(1,0,0), cex.axis=0.8)
	# B: during clinical conversion
ind.byindex$s1_group<-droplevels(ind.byindex$s1_group)				#you must drop the ghost levels
boxplot.double = boxplot(thetayc~delta+s1_group_reinfection, data = ind.byindex, at = c(1, 2, 4, 5, 7, 8), 
	ylim = c(0,1.2), col = c(rep('chartreuse3', 2), rep('orange', 2), rep('darkgoldenrod', 2)), xlim=c(0,9), 
	ylab=expression(paste("", theta, "yc distance")), las=2, mgp=c(0.5,1,1), xaxt='n', xlab='', cex.axis=0.8, cex=0.6, cex.lab=0.8)
#legend("bottomright", c("nonrecurrent", "recurrent", 'reinfection'), col=c('chartreuse3', 'orange', 'darkgoldenrod'), pch=15, cex=0.6)
axis(side=2, labels="(community dissimilarity)", at=0.55, line=0.5, lwd=0, cex.axis=0.8, mgp=c(1,1.8,0))
axis(side=1, at=c(1, 2, 4, 5, 7, 8), labels=rep(c('-', '+'),3), line=0.5, lwd=0, cex.axis=1.5, mgp=c(1,0.5,0))
mtext(expression(paste("", Delta, " in")), 1, cex=0.8, at=c(-1,-1), padj=0.6)
mtext("clinical \nstatus:", 1, cex=0.8, at=c(-1,-1), padj=1.5)
Arrows(0.8, 1.1, 2.2, 1.1, code = 3, arr.type = "T", arr.length=0.2)
Arrows(3.8, 1.1, 5.2, 1.1, code = 3, arr.type = "T", arr.length=0.2)
Arrows(6.8, 1.1, 8.2, 1.1, code = 3, arr.type = "T", arr.length=0.2)
text(1.5,1.15, labels="** p=0.0023", cex=0.6)
text(4.5,1.15, labels="ns", cex=0.6)
text(7.5,1.15, labels="p=0.058", cex=0.6)