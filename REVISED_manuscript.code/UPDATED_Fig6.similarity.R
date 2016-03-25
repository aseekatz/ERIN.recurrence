## Making diversity plots (index single-point comparison AND over time)
# 3.16.16
# Anna M. Seekatz

# files used:
	# ../erinsubset.0.03.summary		#note: you have to fix the columns before reading into R
	# updated_erinsubset_extra.summary.txt
	# updated_erinsub_ind.dist.byindex.txt (created)
	# updated_erinsub_shared.meta.dist.txt (created)
	# updated_erinsub_ind.dist.overtime.txt (created)
	
library(shape)
library(plyr)

# have to recalculate a summary.shared file from the .dist file created earlier:
#(in mothur):
	# summary.shared(shared=/Users/annaseekatz/Desktop/umich/projects/ERIN_clinical/erinfmt_analysis_athome/ERINsubset/erinsubset_mothurfiles/erinsubset.0.03.shared, calc=sharedsobs-braycurtis-spearman-thetayc-jsd-sharednseqs)

## create a shared distance file that incorporates the metadata for both sample columns:

# merge shared files with meta (includes extra meta):
var<-read.table(file="updated_erinsubset_extra.summary.txt", header=TRUE)
mdist<-read.table(file="../erinsubset_mothurfiles/erinsubset.0.03.summary", header=TRUE)

#add sample information for 'group1':
m<-merge(var, mdist, by.x=c("seqID"), by.y=c("s1"))				#merge for sample1 metadata 			
m.2<-m[,1:54]													#separate to add a s1 to colnames
m.3<-m[,57:65]
colnames(m.2) <- paste("s1", colnames(m.2), sep = "_")
m<-cbind(m.2,m.3)

#add sample information for 2nd sample: ('group2'):
m2<-merge(var, m, by.x=c("seqID"), by.y=c("s2"))				#merge for sample2 metadata 			
m2.2<-m2[,1:54]													#separate to add a s1 to colnames
m2.3<-m2[,57:118]
colnames(m2.2) <- paste("s2", colnames(m2.2), sep = "_")
m.shared<-cbind(m2.2,m2.3)

#write.table(m.shared, file="updated_erinsub_shared.meta.dist.txt", quote=FALSE, sep="\t")

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
#write.table(ind, file="updated_erinsub_dist.within.individual.txt", quote=FALSE, sep="\t", col.names=NA)

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
ind.1a<-ind[(ind$s1_index_sample_n=="1" & ind$s2_index_sample_n=="2"), ]
ind.1b<-ind[(ind$s1_index_sample_n=="2" & ind$s2_index_sample_n=="1"), ]		#none this way
ind1<-rbind(ind.1a, ind.1b)
ind1$sampling<- c(1)

# differences between sampling 2 and 3:
ind.2a<-ind[(ind$s1_index_sample_n=="3" & ind$s2_index_sample_n=="2"), ]
ind.2b<-ind[(ind$s1_index_sample_n=="2" & ind$s2_index_sample_n=="3"), ]		
ind2<-rbind(ind.2a, ind.2b)
ind2$sampling<- c(2)

# differences between sampling 3 and 4:
ind.3a<-ind[(ind$s1_index_sample_n=="3" & ind$s2_index_sample_n=="4"), ]
ind.3b<-ind[(ind$s1_index_sample_n=="4" & ind$s2_index_sample_n=="3"), ]		
ind3<-rbind(ind.3a, ind.3b)
ind3$sampling<- c(3)

# differences between sampling 4 and 5:
ind.4a<-ind[(ind$s1_index_sample_n=="5" & ind$s2_index_sample_n=="4"), ]
ind.4b<-ind[(ind$s1_index_sample_n=="4" & ind$s2_index_sample_n=="5"), ]		
ind4<-rbind(ind.4a, ind.4b)
ind4$sampling<- c(4)

# differences between sampling 5 and 6:
ind.5a<-ind[(ind$s1_index_sample_n=="5" & ind$s2_index_sample_n=="6"), ]
ind.5b<-ind[(ind$s1_index_sample_n=="6" & ind$s2_index_sample_n=="5"), ]		
ind5<-rbind(ind.5a, ind.5b)
ind5$sampling<- c(5)

# differences between sampling 6 and 7:
ind.6a<-ind[(ind$s1_index_sample_n=="7" & ind$s2_index_sample_n=="6"), ]
ind.6b<-ind[(ind$s1_index_sample_n=="6" & ind$s2_index_sample_n=="7"), ]		
ind6<-rbind(ind.6a, ind.6b)
ind6$sampling<- c(6)

ind.overtime<-rbind(ind1, ind2, ind3, ind4, ind5, ind6)
#write.table(ind.overtime, file="updated_erinsub_ind.dist.overtime.txt", quote=FALSE, sep="\t", col.names=NA)


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
tyc.data2<-tyc.data[tyc.data$sampling %in% c(1,2),]
split.tyc<-split(tyc.data2, tyc.data2$s1_group_reinfection)
tyc.rec<-split.tyc$'recurrent'
tyc.non<-split.tyc$'nonrecurrent'
tyc.reinf<-split.tyc$'reinfection'

# plot it!
## Fig. S4:
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



### some extra analyses:

# using new category, sample status
ind<-read.table(file="updated_erinsub_dist.within.individual.txt", header=TRUE)

# let's compare the sample statuses of longitudinally sampled individuals...
ind.all<-split(ind, ind$s1_group_reinfection)
ind.reinf<-ind.all$'reinfection'
ind.non<-ind.all$'nonrecurrent'
ind.rec<-ind.all$'recurrent'

# reinfected group:
ind.reinf<-droplevels(ind.reinf)
summary(ind.reinf$s1_sample_status)
	# check out the spread...
	# let's compare index to recovery/treatment, index to reinfection, 
reinf.1a<-ind.reinf[(ind.reinf$s1_sample_status=="index" & ind.reinf$s2_sample_status=="reinfection"), ]
reinf.1b<-ind.reinf[(ind.reinf$s2_sample_status=="index" & ind.reinf$s1_sample_status=="reinfection"), ]
reinf1<-rbind(reinf.1a, reinf.1b)
reinf1$STATUS<-c("index_to_reinfection")

#reinf.2a<-ind.reinf[(ind.reinf$s1_sample_status=="index" & ind.reinf$s2_sample_status %in% c("recovery", "treatment")), ]
#reinf.2b<-ind.reinf[(ind.reinf$s2_sample_status=="index" & ind.reinf$s1_sample_status %in% c("recovery", "treatment")), ]
#reinf2<-rbind(reinf.2a, reinf.2b)
#reinf2$STATUS<-c("index_to_other")

reinf.2a<-ind.reinf[(ind.reinf$s1_sample_status=="index" & ind.reinf$s2_sample_status %in% c("recovery")), ]
reinf.2b<-ind.reinf[(ind.reinf$s2_sample_status=="index" & ind.reinf$s1_sample_status %in% c("recovery")), ]
reinf2<-rbind(reinf.2a, reinf.2b)
reinf2$STATUS<-c("index_to_recovery")

reinf.3a<-ind.reinf[(ind.reinf$s1_sample_status=="index" & ind.reinf$s2_sample_status %in% c("treatment")), ]
reinf.3b<-ind.reinf[(ind.reinf$s2_sample_status=="index" & ind.reinf$s1_sample_status %in% c("treatment")), ]
reinf3<-rbind(reinf.3a, reinf.3b)
reinf3$STATUS<-c("index_to_treatment")

reinf<-rbind(reinf1, reinf2, reinf3)

# nonrecurrent group:
ind.non<-droplevels(ind.non)
summary(ind.non$s1_sample_status)
	# check out the spread...
	# let's compare index to recovery/treatment
non.1a<-ind.non[(ind.non$s1_sample_status=="index" & ind.non$s2_sample_status=="recovery"), ]
non.1b<-ind.non[(ind.non$s2_sample_status=="index" & ind.non$s1_sample_status=="recovery"), ]
non1<-rbind(non.1a, non.1b)
non1$STATUS<-c("index_to_recovery")

non.2a<-ind.non[(ind.non$s1_sample_status=="index" & ind.non$s2_sample_status=="treatment"), ]
non.2b<-ind.non[(ind.non$s2_sample_status=="index" & ind.non$s1_sample_status=="treatment"), ]
non2<-rbind(non.2a, non.2b)
non2$STATUS<-c("index_to_treatment")

#non.2a<-ind.non[(ind.non$s1_sample_status=="index" & ind.non$s2_sample_status %in% c("recovery", "treatment")), ]
#non.2b<-ind.non[(ind.non$s2_sample_status=="index" & ind.non$s1_sample_status %in% c("recovery", "treatment")), ]
#non2<-rbind(non.2a, non.2b)
#non2$STATUS<-c("index_to_other")
#non<-non2

non<-rbind(non1, non2)

# recurrent group:
ind.rec<-droplevels(ind.rec)
summary(ind.rec$s1_sample_status)
	# check out the spread...
	# let's compare index to recovery/treatment, index to reinfection, index to recurrence: 
rec.1a<-ind.rec[(ind.rec$s1_sample_status=="index" & ind.rec$s2_sample_status=="recurrence"), ]
rec.1b<-ind.rec[(ind.rec$s2_sample_status=="index" & ind.rec$s1_sample_status=="recurrence"), ]
rec1<-rbind(rec.1a, rec.1b)
rec1$STATUS<-c("index_to_recurrence")

#rec.2a<-ind.rec[(ind.rec$s1_sample_status=="index" & ind.rec$s2_sample_status %in% c("recovery", "treatment")), ]
#rec.2b<-ind.rec[(ind.rec$s2_sample_status=="index" & ind.rec$s1_sample_status %in% c("recovery", "treatment")), ]
#rec2<-rbind(rec.2a, rec.2b)
#rec2$STATUS<-c("index_to_other")

rec.3a<-ind.rec[(ind.rec$s1_sample_status=="index" & ind.rec$s2_sample_status=="reinfection"), ]
rec.3b<-ind.rec[(ind.rec$s2_sample_status=="index" & ind.rec$s1_sample_status=="reinfection"), ]
rec3<-rbind(rec.3a, rec.3b)
rec3$STATUS<-c("index_to_reinfection")

rec.2a<-ind.rec[(ind.rec$s1_sample_status=="index" & ind.rec$s2_sample_status %in% c("recovery")), ]
rec.2b<-ind.rec[(ind.rec$s2_sample_status=="index" & ind.rec$s1_sample_status %in% c("recovery")), ]
rec2<-rbind(rec.2a, rec.2b)
rec2$STATUS<-c("index_to_recovery")

rec.4a<-ind.rec[(ind.rec$s1_sample_status=="index" & ind.rec$s2_sample_status %in% c("treatment")), ]
rec.4b<-ind.rec[(ind.rec$s2_sample_status=="index" & ind.rec$s1_sample_status %in% c("treatment")), ]
rec4<-rbind(rec.4a, rec.4b)
rec4$STATUS<-c("index_to_treatment")

rec<-rbind(rec1, rec2, rec3, rec4)

# combine all?
all<-rbind(reinf, non, rec)
summary(as.factor(all$STATUS))

# ok...we now have some comparisons.
# let's graph:

## status within group:
boxplot.double = boxplot(thetayc~STATUS+s1_group_reinfection, data = all, 
	ylim = c(0,1), xlim=c(0,14), 
	ylab=expression(paste("", theta, "yc distance")), las=2, mgp=c(0.5,1,1), xlab='', cex.axis=0.8)
# this works!

## trying it with dots:
group2.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "index_to_recovery_recurrent" ) {
	colorvec[i] = "orange"
	}
	if ( n[i] == "index_to_recurrence_recurrent" ) {
	colorvec[i] = "orange"
	}
	if ( n[i] == "index_to_reinfection_recurrent" ) {
	colorvec[i] = "orange"
	}
	if ( n[i] == "index_to_treatment_recurrent" ) {
	colorvec[i] = "orange"
	}
	if ( n[i] == "index_to_recovery_nonrecurrent" ) {
	colorvec[i] = "chartreuse3"
	}
	if ( n[i] == "index_to_treatment_nonrecurrent" ) {
	colorvec[i] = "chartreuse3"
	}
	if ( n[i] == "index_to_recovery_reinfection" ) {
	colorvec[i] = "darkgoldenrod"
	}
	if ( n[i] == "index_to_reinfection_reinfection" ) {
	colorvec[i] = "darkgoldenrod"
	}
	if ( n[i] == "index_to_treatment_reinfection" ) {
	colorvec[i] = "darkgoldenrod"
	}
	}
	c(colorvec)
	}

### Fig. 6A:
## similarity, all samples:	
data<-all		# ALL index
plot<-plot(thetayc ~ as.factor(s1_group_reinfection), data = data, ylab=expression(paste("", theta, "yc distance")), ylim=c(0,1.2),
	xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), cex=0.8, cex.lab=0.8, cex.axis=0.8)
points(thetayc ~ jitter(as.numeric(s1_group_reinfection, factor=0)), data = data, bg=group.col(data$s1_group_reinfection), col="black", pch=21, cex=0.8)
names<-c("nonrecurrent", "recurrent", "reinfection")
text(x =  seq(1,3,by=1), y = par("usr")[3]-0, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=0.8)
kruskal.test(thetayc~s1_group_reinfection, data=data)
text(1.5,19, labels="Kruskal-Wallis test \np=0.025", cex=0.6)


## Figure 6B:
# reorder some things:
data<-all
group<-data$s2_group_reinfection
stat<-as.character(data$STATUS)
data$combo<-paste(stat, group, sep="_")
data$combo<-as.factor(data$combo)
levels(data$combo)
data$combo<- ordered(data$combo, 
	levels = c("index_to_recovery_nonrecurrent", "index_to_recovery_recurrent", "index_to_recovery_reinfection", 
				"index_to_recurrence_recurrent", "index_to_reinfection_recurrent", "index_to_reinfection_reinfection", 
				"index_to_treatment_nonrecurrent", "index_to_treatment_recurrent", "index_to_treatment_reinfection"))
ycat<-data$combo
plot<-plot(thetayc ~ as.factor(ycat), data = data, ylab=expression(paste("", theta, "yc distance")), xlab="", outline=FALSE, 
	ylim=c(0,1.2), cex.lab=1, cex.axis=0.6, las=2)
points(thetayc ~ jitter(as.numeric(ycat, factor=0)), data = data, pch=21,
	bg=group2.col(data$combo), col="black", cex=0.9)
names<-as.character(unique(data$combo))
text(x =  c(1.5, 3.5, 5.5), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
legend("topright", c("nonrecurrent", "recurrent", "reinfected"), col=c("chartreuse3", "orange", "darkgoldenrod"), pch=c(19), cex=0.8)
kruskal.test(thetayc~combo, data=data)

# option 1: by major category of change vs. no change:
boxplot.double = boxplot(thetayc~s1_group_reinfection + STATUS, data = all, #at = c(1, 2, 3, 5, 6, 7), 
	ylim = c(0,1), col = c('chartreuse3', 'orange', 'darkgoldenrod'), xlim=c(0,8), 
	ylab=expression(paste("", theta, "yc distance")), las=2, mgp=c(0.5,1,1), xaxt='n', xlab='', cex.axis=0.8)
legend("bottomright", c("nonrecurrent", "recurrent", "reinfected"), col=c('chartreuse3', 'orange', 'darkgoldenrod'), pch=15, cex=0.7)
axis(side=2, labels="(community dissimilarity)", at=0.55, line=0.5, lwd=0, cex.axis=0.8, mgp=c(1,1.8,0))
axis(side=1, at=c(2, 6), labels=c('No change \n(CDI status)', 'Index change \n(CDI status)'), line=0.5, lwd=0, cex.axis=0.8, mgp=c(1,0.5,0))
title('Change in community composition \nduring CDI diagnosis')






