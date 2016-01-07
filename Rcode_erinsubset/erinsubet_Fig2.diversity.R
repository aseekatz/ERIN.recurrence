## Making diversity plots (index single-point comparison AND over time)
# 1.4.15
# Anna M. Seekatz

# files used:
	# allerin_summary.phylo2.txt
	# erinfmt_summary.txt
	# erinfmt_extra.meta.txt


### For diversity between groups, did the following:

### Fig. 2A, 2C: diversity between patient groups, severity index

extra<-read.table(file="erinfmt_extra.meta.txt", header=TRUE)
	# from this, we want to include severeidsa, abx_prior, abx_allprior, ppi, final_plating
filtered.extra<-extra[, c("seqID", "abx_prior", "abx_typeprior", "ppi", "final_plating", "severeidsa")]

sum<-read.table(file='../erinfmt_summary.txt', header=TRUE)
subsum<-sum[sum$group_reinfection %in% c("recurrent", "nonrecurrent", "reinfection"), ]
subsum<-droplevels(subsum)
subsum2<-merge(subsum, filtered.extra, by.x="seqID", by.y="seqID", all.x=TRUE)
	# finally, take only first sample:
#subsum3<-subsum2[subsum2$sample_n==1, ]

# let's subset only the first sample (as above):
# let's also only take the patients that were sasmples longitudinally:
# we will also add a new column, 'change.simpson' to reflect the change in the invsimpson index from each patient's initial sample

# take first day for each patient, then calculate the change from that...
firstdiv<-subsum2[which(subsum2$sample_n==1), c("patientID","invsimpson","group2", "group_reinfection")]
duplicated(firstdiv$patientID)		# patient EIN191 is duplicated?
	#fix this:
	subsum2[subsum2$patientID==c("EIN00191"),]
	subsum2$sample_n[subsum2$patientID==c("EIN00191") & subsum2$rel_day==7]<-2
# now, create 'change in diversity' column:
	# combined to get a column with the diversity upon first sampling
subsum3<-merge(subsum2, firstdiv, by.x="patientID", by.y="patientID")
	# get difference in diversity:
subsum3$change.simpson<-subsum3$invsimpson.y-subsum3$invsimpson.x
	# if you want to eliminate first point:
#subsum4<-subsum3[-which(subsum3$sample_n==1), ]
subsum2.longi<-subsum2[subsum2$total_sample_n>1, ]		
firstdiv2<-subsum2.longi[which(subsum2.longi$sample_n==1), c("patientID","invsimpson","group2", "group_reinfection", "clinical_result", "final_plating", "severeidsa", "abx_prior", "abx_typeprior", "ppi")]

	# had previously 
boxplot(invsimpson~group_reinfection, data=firstdiv)
	# looks promising: let's only use longitudinal patients, since we know more about them...

# define some color vectors:
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
sev.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "yes" ) {
	colorvec[i] = "hotpink"
	}
	if ( n[i] == "no" ) {
	colorvec[i] = "thistle1"
	}
	}
	c(colorvec)
	}

# plot it:
par(mfrow=c(1,2))
par(mai=c(0.75,0.82,0.42,0.22))
# 2A: diversity between different patient groups (first sample only):
plot<-plot(invsimpson ~ as.factor(group_reinfection), data = firstdiv2, ylab="inverse Simpson", xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,20), cex=0.8, cex.lab=0.8, cex.axis=0.8)
points(invsimpson ~ jitter(as.numeric(group_reinfection, factor=0)), data = firstdiv2, bg=group.col(firstdiv2$group_reinfection), col="black", pch=21, cex=0.8)
names<-as.character(unique(firstdiv2$group_reinfection))
text(x =  seq(1,3,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=0.8)
kruskal.test(invsimpson~group_reinfection, data=firstdiv2)
text(1.5,19, labels="Kruskal-Wallis test \np=0.070", cex=0.6)
	#Kruskal-Wallis rank sum test
	#data:  invsimpson by group_reinfection
	#Kruskal-Wallis chi-squared = 5.3128, df = 2, p-value = 0.0702
	# so close!!!
	summary(as.factor(firstdiv$group_reinfection))		# uneven proportions

# 2B: diversity between severe and non severe samples
# using all samples (not just index sample)
# we need to subset the samples that we have severity info for (about 2/3 samples)
subsum2$severeidsa<-as.factor(subsum2$severeidsa)
subsum3<-subsum2[subsum2$severeidsa %in% c("0", "1"), ]
subsum3<-droplevels(subsum3)
ycat<-subsum3$severeidsa
data<-subsum3
plot<-plot(invsimpson ~ as.factor(ycat), data = data, ylab="inverse Simpson", xlab="", outline=FALSE, xaxt="n", ylim=c(0,20), cex=0.8, cex.lab=0.8, cex.axis=0.8)
points(invsimpson ~ jitter(as.numeric(ycat, factor=0)), data = data, bg=sev.col(severeidsa), cex=0.8, col="black", pch=21, mgp=c(1,2,2))
names<-c("no", "severe")
text(x =  seq(1,3,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=0.8)
wilcox.test(invsimpson~ycat, data=data)
	## significance!!!!!
text(1,19, labels="*Wilcoxon test \np=0.037", cex=0.6)

### Fig. S1A, B, C, D: Other diversity index comparisons

# some color/shape vectors first:
clin.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "positive" ) {
	colorvec[i] = "magenta"
	}
	if ( n[i] == "negative" ) {
	colorvec[i] = "grey47"
	}
	}
	c(colorvec)
	}
abx.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "yes" ) {
	colorvec[i] = "hotpink"
	}
	if ( n[i] == "no" ) {
	colorvec[i] = "thistle1"
	}
	}
	c(colorvec)
	}
sev.pch <- function(n) {
	pchvec <- vector(mode="numeric", length=length(n))
	for (i in 1:length(n)) {
	pchvec[i] = 15
	if ( n[i] == 1 ) {
	pchvec[i] = 24
	}
	if (n[i] == 0 ) {
	pchvec[i] = 21
	}
	}
	(pchvec)
	}
	
# plotting:
par(mfrow=c(1,3))
par(mai=c(0.75,0.45,0.2,0.1))

# S1A: by clinical result (not significant):
# using all samples (subsum2 file)	
ycat<-subsum2$clinical_result
plot<-plot(invsimpson ~ as.factor(ycat), data = subsum2, ylab="inverse Simpson", xlab="", 
	outline=FALSE, xaxt="n", ylim=c(0,20),
	cex.axis=1, cex.lab=1)
points(invsimpson ~ jitter(as.numeric(ycat, factor=0)), data = subsum2, 
	bg=clin.col(ycat), col="black", pch=21, cex=1)
names<-levels(ycat)
text(x =  seq(1,3,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
wilcox.test(invsimpson~ycat, data=subsum2)
text(1,19, labels="Wilcoxon test\np=0.14", cex=0.8)

#s1B: by abx use
# must subset some samples (about 2/3)
subsum3<-subsum2[subsum2$abx_prior %in% c("yes", "no"), ]
subsum3<-droplevels(subsum3)
ycat<-subsum3$abx_prior
data<-subsum3
plot<-plot(invsimpson ~ as.factor(ycat), data = data, ylab="inverse Simpson", xlab="", 
	outline=FALSE, xaxt="n", ylim=c(0,20), cex.axis=1, cex.lab=1)
points(invsimpson ~ jitter(as.numeric(ycat, factor=0)), 
	data = data, bg=abx.col(ycat), col="black", pch=21, cex=1)
names<-levels(unique(ycat))
text(x =  seq(1,3,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
wilcox.test(invsimpson~ycat, data=data)
text(1,19, labels="Wilcoxon test \np=0.67", cex=0.8)

## S1C: slicing and dicing: severity WITHIN patient groups:
prop.table(table(subsum3$group_reinfection, subsum3$severeidsa))
table(subsum3$group_reinfection, subsum3$severeidsa)
	# pretty even distributions of severity within these groups...
subsum2$severeidsa<-as.factor(subsum2$severeidsa)
subsum3<-subsum2[subsum2$severeidsa %in% c("0", "1"), ]
subsum3<-droplevels(subsum3)
ycat<-subsum3$severeidsa
ycat2<-subsum3$group_reinfection
data<-subsum3
# for now, let's do it this way:
group<-as.character(subsum3$group_reinfection)
sev<-as.character(subsum3$severeidsa)
subsum3$combo<-paste(group, sev, sep="_")
subsum3$combo<-as.factor(subsum3$combo)
ycat<-subsum3$combo
data<-subsum3
plot<-plot(invsimpson ~ as.factor(ycat), data = data, ylab="inverse Simpson", xlab="", outline=FALSE, xaxt="n", 
	ylim=c(0,20), cex.lab=1, cex.axis=1)
		# although the at function can be used here to separate the groups, I can't figure out how to add the points when the 'at' is specified...
points(invsimpson ~ jitter(as.numeric(ycat, factor=0)), data = data, 
	bg=group.col(group_reinfection), col="black", pch=sev.pch(severeidsa), cex=0.9)
names<-as.character(unique(subsum3$group_reinfection))
text(x =  c(1.5, 3.5, 5.5), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
legend("topright", c("not severe", "severe"), col=c("black"), pch=c(21,24), cex=0.8)
kruskal.test(invsimpson~ycat, data=data)
text(4.9,15, labels="*Kruskal-Wallis test \np=0.049", cex=0.8)
	## BARELY significant


###
#---
###

### For over time graphing, need to subset groups:

sum<-read.table(file='../erinfmt_summary.txt', header=TRUE)
subsum<-sum[sum$group2 %in% c("recurrent", "nonrecurrent"), ]
levels(as.factor(subsum$total_sample_n))
subsum2<-subsum[subsum$total_sample_n > 1, ]
subsum2<-droplevels(subsum2)

## Fig. 2C: diversity over time

# graphing each group, with a regression line:
recur<-subset(subsum2, (subset=group2 %in% c("recurrent")))				#subset data by group
non<-subset(subsum2, (subset=group2 %in% c("nonrecurrent")))

# let's subset ONLY the first (index) sample collected, which is sample_n=1:
firstdiv2<-subsum2[which(subsum2$sample_n==1), c("patientID","invsimpson","group2", "group_reinfection")]
duplicated(firstdiv$patientID)		# check if there are any mistakes in the meta files

# let's 
# now, create 'change in diversity' column:
	# combined to get a column with the diversity upon first sampling
subsum3<-merge(subsum2, firstdiv, by.x="patientID", by.y="patientID")
	# get difference in diversity:
subsum3$change.simpson<-subsum3$invsimpson.y-subsum3$invsimpson.x

# for some of the graphs, we want to eliminate 
subsum4<-subsum3[-which(subsum3$sample_n==1), ]

# then, subset each group you are graphing:
recur<-subset(subsum4, (subset=group_reinfection.x %in% c("recurrent")))				#subset data by group
non<-subset(subsum4, (subset=group_reinfection.x %in% c("nonrecurrent")))
reinf<-subset(subsum4, (subset=group_reinfection.x %in% c("reinfection")))

### Fig. S2: 
	#A: diversity over time, by sampling period:
	#B: Change in diversity per patient, by sampling period:
recur<-subset(subsum3, (subset=group_reinfection.x %in% c("recurrent")))				#subset data by group
non<-subset(subsum3, (subset=group_reinfection.x %in% c("nonrecurrent")))
reinf<-subset(subsum3, (subset=group_reinfection.x %in% c("reinfection")))

par(mfrow=c(2,1))
par(mai=c(0.75,0.82,0.42,0.22))
# S2A: flat diversity over time
plot(jitter(recur$sample_n, amount=0.05), recur$invsimpson.x, mgp=c(2,1,0),
	pch=21, col="black", bg="orange", cex=0.8, ylim=c(0,50), cex.lab=0.8, cex.axis=0.8,
	main="Diversity over time", ylab="invsimpson", xlab="relative sampling")
points(0.25+jitter(non$sample_n, amount=0.05), non$invsimpson.x, pch=21, xaxt='n', col="black", bg="chartreuse3", cex=0.8)
points(0.5+jitter(reinf$sample_n, amount=0.05), reinf$invsimpson.x, pch=21, xaxt='n', col="black", bg="darkgoldenrod", cex=0.8)
recline <- lm(invsimpson.x~sample_n, data=recur)
nonline <- lm(invsimpson.x~sample_n, data=non)
reinfline <- lm(invsimpson.x~sample_n, data=reinf)
abline(nonline, col="chartreuse3")
abline(recline, col="orange")
abline(reinfline, col="darkgoldenrod")
legend("topright", c("recurrent", "nonrecurrent", "reinfected"), col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")

# S2B: CHANGE in simpson over sampling, excluding first sampling:
recur2<-subset(subsum4, (subset=group_reinfection.x %in% c("recurrent")))				#subset data by group
non2<-subset(subsum4, (subset=group_reinfection.x %in% c("nonrecurrent")))
reinf2<-subset(subsum4, (subset=group_reinfection.x %in% c("reinfection")))
plot(jitter(recur2$sample_n, amount=0.05), recur2$change.simpson, cex.lab=0.8, cex.axis=0.8,
	pch=21, col="black", bg="orange", cex=0.8, ylim=c(-20,25), mgp=c(2,1,0),
	main="Change in diversity over time", ylab=expression(paste("", Delta, " invsimpson")), xlab="relative sampling")
points(0.25+jitter(non2$sample_n, amount=0.05), non2$change.simpson, pch=21, xaxt='n', col="black", bg="chartreuse3", cex=0.8)
points(0.5+jitter(reinf2$sample_n, amount=0.05), reinf2$change.simpson, pch=21, xaxt='n', col="black", bg="darkgoldenrod", cex=0.8)
recline2 <- lm(change.simpson~sample_n, data=recur2)
nonline2 <- lm(change.simpson~sample_n, data=non2)
reinfline2 <- lm(change.simpson~sample_n, data=reinf2)
abline(nonline2, col="chartreuse3")
abline(recline2, col="orange")
abline(reinfline2, col="darkgoldenrod")
legend("bottomright", c("recurrent", "nonrecurrent", "reinfected"), col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")

### Fig. 2C: change in diversity, over actual time from sampling

# change in diversity over time:
# CHANGE in simpson, excluding first sampling, by relative time:
recur2<-subset(subsum4, (subset=group_reinfection.x %in% c("recurrent")))				#subset data by group
non2<-subset(subsum4, (subset=group_reinfection.x %in% c("nonrecurrent")))
reinf2<-subset(subsum4, (subset=group_reinfection.x %in% c("reinfection")))
plot(recur2$rel_day, recur2$change.simpson, 
	pch=21, col="black", bg="orange", cex=0.8, ylim=c(-20,20),
	main="Change in diversity over time, by day", ylab="delta invsimpson", xlab="relative day")
points(non2$rel_day, non2$change.simpson, pch=21, xaxt='n', col="black", bg="chartreuse3", cex=0.8)
points(reinf2$rel_day, reinf2$change.simpson, pch=21, xaxt='n', col="black", bg="darkgoldenrod", cex=0.8)
recline2 <- lm(change.simpson~rel_day, data=recur2)
nonline2 <- lm(change.simpson~rel_day, data=non2)
reinfline2 <- lm(change.simpson~rel_day, data=reinf2)
abline(nonline2, col="chartreuse3")
abline(recline2, col="orange")
abline(reinfline2, col="darkgoldenrod")
legend("bottomright", c("recurrent", "nonrecurrent", "darkgoldenrod"), col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")


# however, not very many points are present in each group (particularly in the nonrecurrent) after 100 days or so
# let's get the mean relative day after sampling:
summary(subsum4$rel_day)
	   #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    	#1.0    13.0    28.0   101.5    97.5   799.0
# let's subset ONLY days before 120:
subsum5<-subset(subsum4, rel_day < 120)

# now, plot that
# 2B:
par(mai=c(0.75,0.82,0.42,0.22))
recur2<-subset(subsum5, (subset=group_reinfection.x %in% c("recurrent")))				
non2<-subset(subsum5, (subset=group_reinfection.x %in% c("nonrecurrent")))
reinf2<-subset(subsum5, (subset=group_reinfection.x %in% c("reinfection")))
plot(recur2$rel_day, recur2$change.simpson, xlim=c(0,120), cex.lab=0.8, cex.axis=0.8,
	pch=21, col="black", bg="orange", cex=0.8, ylim=c(-20,25), mgp=c(2,1,0),
	main="Change in diversity over time, by day", ylab=expression(paste("", Delta, " invsimpson")), xlab="relative day")
points(non2$rel_day, non2$change.simpson, pch=21, xaxt='n', col="black", bg="chartreuse3", cex=0.8)
points(reinf2$rel_day, reinf2$change.simpson, pch=21, xaxt='n', col="black", bg="darkgoldenrod", cex=0.8)
recline2 <- lm(change.simpson~rel_day, data=recur2)
nonline2 <- lm(change.simpson~rel_day, data=non2)
reinfline2 <- lm(change.simpson~rel_day, data=reinf2)
abline(nonline2, col="chartreuse3")
abline(recline2, col="orange")
abline(reinfline2, col="darkgoldenrod")
legend("topleft", c("recurrent", "nonrecurrent", "darkgoldenrod"), col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")

