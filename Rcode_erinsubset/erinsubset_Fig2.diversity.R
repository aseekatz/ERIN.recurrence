## Making diversity plots (index single-point comparison AND over time)
# 1.4.16
# Anna M. Seekatz

# files used:
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
boxplot(invsimpson~group_reinfection, data=firstdiv2)
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
	if ( n[i] == "1" ) {
	colorvec[i] = "hotpink"
	}
	if ( n[i] == "0" ) {
	colorvec[i] = "thistle1"
	}
	}
	c(colorvec)
	}

# plot it:
par(mfrow=c(1,2))
par(mai=c(0.75,0.82,0.42,0.22))
# 2A: diversity between different patient groups (first sample only):
plot<-plot(invsimpson ~ as.factor(group_reinfection), data = firstdiv2, ylab=expression(paste("", gamma, " (diversity)")), 
	xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,20), cex=0.8, cex.lab=0.8, cex.axis=0.8)
points(invsimpson ~ jitter(as.numeric(group_reinfection, factor=0)), data = firstdiv2, bg=group.col(firstdiv2$group_reinfection), col="black", pch=21, cex=0.8)
names<-as.character(unique(firstdiv2$group_reinfection))
text(x =  seq(1,3,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=0.8)
kruskal.test(invsimpson~group_reinfection, data=firstdiv2)
text(1.5,20, labels="Kruskal-Wallis test \np=0.070", cex=0.6)
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
plot<-plot(invsimpson ~ as.factor(ycat), data = data, ylab=expression(paste("", gamma, " (diversity)")), 
	xlab="", outline=FALSE, xaxt="n", ylim=c(0,20), cex=0.8, cex.lab=0.8, cex.axis=0.8)
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
par(mai=c(0.75,0.52,0.2,0.1))

# S1A: by clinical result (not significant):
# using all samples (subsum2 file)	
ycat<-subsum2$clinical_result
plot<-plot(invsimpson ~ as.factor(ycat), data = subsum2, ylab=expression(paste("", gamma, " (diversity)")), xlab="", 
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
plot<-plot(invsimpson ~ as.factor(ycat), data = data, ylab=expression(paste("", gamma, " (diversity)")), xlab="", 
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
plot<-plot(invsimpson ~ as.factor(ycat), data = data, ylab=expression(paste("", gamma, " (diversity)")), xlab="", outline=FALSE, xaxt="n", 
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

# take first day for each patient, then calculate the change from that...
firstdiv<-subsum2[which(subsum2$sample_n==1), c("patientID","invsimpson","group2", "group_reinfection")]
duplicated(firstdiv$patientID)		# patient EIN191 is duplicated?
	#fix this:
	#subsum2[subsum2$patientID==c("EIN00191"),]
	#subsum2$sample_n[subsum2$patientID==c("EIN00191") & subsum2$rel_day==7]<-2
# now, create 'change in diversity' column:
	# combined to get a column with the diversity upon first sampling
subsum3<-merge(subsum2, firstdiv, by.x="patientID", by.y="patientID")
	# get difference in diversity:
subsum3$change.simpson<-subsum3$invsimpson.y-subsum3$invsimpson.x

## Fig. 2C: diversity over time

# for some of the graphs, we want to eliminate the first samples, since that is 0 for change.simpson
subsum4<-subsum3[-which(subsum3$sample_n==1), ]

### Fig. S2: 
	#A: diversity over time, by sampling period, all samples:
	#B: Change in diversity per patient, by sampling period, ignore first change.simpson point:

# S2A: simpson over time:
recur<-subset(subsum3, (subset=group_reinfection.x %in% c("recurrent")))				#subset data by group
non<-subset(subsum3, (subset=group_reinfection.x %in% c("nonrecurrent")))
reinf<-subset(subsum3, (subset=group_reinfection.x %in% c("reinfection")))
par(mfrow=c(2,1))
par(mai=c(0.75,0.82,0.42,0.22))
# S2A: flat diversity over time
plot(jitter(recur$sample_n, amount=0.05), recur$invsimpson.x, mgp=c(2,1,0),
	pch=21, col="black", bg="orange", cex=0.8, ylim=c(0,50), cex.lab=0.8, cex.axis=0.8,
	main="", ylab=expression(paste("", gamma, " (diversity)")), xlab="relative sampling")
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
	main="", ylab=expression(paste("", Delta, gamma, " (", Delta, "invsimpson)")), xlab="relative sampling")
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
plot(recur2$rel_day, recur2$change.simpson, cex.lab=0.8, cex.axis=0.8,
	pch=21, col="black", bg="orange", cex=0.8, ylim=c(-15,15),
	main="", ylab=expression(paste("", Delta, ,gamma, " (", Delta, "diversity)")), xlab="relative day", mgp=c(2,1,0))
points(non2$rel_day, non2$change.simpson, pch=21, xaxt='n', col="black", bg="chartreuse3", cex=0.8)
points(reinf2$rel_day, reinf2$change.simpson, pch=21, xaxt='n', col="black", bg="darkgoldenrod", cex=0.8)
recline2 <- lm(change.simpson~rel_day, data=recur2)
nonline2 <- lm(change.simpson~rel_day, data=non2)
reinfline2 <- lm(change.simpson~rel_day, data=reinf2)
abline(nonline2, col="chartreuse3")
abline(recline2, col="orange")
abline(reinfline2, col="darkgoldenrod")
legend("bottomright", c("recurrent", "nonrecurrent", "reinfected"), col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")


# however, not very many points are present in each group (particularly in the nonrecurrent) after 100 days or so
# let's get the mean relative day after sampling:
summary(subsum4$rel_day)
	   #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    	#1.0    13.0    28.0   101.5    97.5   799.0
# let's subset ONLY days before 120:
subsum5<-subset(subsum4, rel_day < 120)

# now, plot that
# 2C, replotted:
par(mai=c(0.75,0.82,0.42,0.22))
recur2<-subset(subsum5, (subset=group_reinfection.x %in% c("recurrent")))				#subset data by group
non2<-subset(subsum5, (subset=group_reinfection.x %in% c("nonrecurrent")))
reinf2<-subset(subsum5, (subset=group_reinfection.x %in% c("reinfection")))
plot(recur2$rel_day, recur2$change.simpson, cex.lab=0.8, cex.axis=0.8,
	pch=21, col="black", bg="orange", cex=0.8, ylim=c(-15,15), xlim=c(0,120),
	main="", ylab=expression(paste("", Delta, ,gamma, " (", Delta, "diversity)")), xlab="relative day", mgp=c(2,1,0))
points(non2$rel_day, non2$change.simpson, pch=21, xaxt='n', col="black", bg="chartreuse3", cex=0.8)
points(reinf2$rel_day, reinf2$change.simpson, pch=21, xaxt='n', col="black", bg="darkgoldenrod", cex=0.8)
recline2 <- lm(change.simpson~rel_day, data=recur2)
nonline2 <- lm(change.simpson~rel_day, data=non2)
reinfline2 <- lm(change.simpson~rel_day, data=reinf2)
abline(nonline2, col="chartreuse3")
abline(recline2, col="orange")
abline(reinfline2, col="darkgoldenrod")
legend("bottomright", c("recurrent", "nonrecurrent", "reinfected"), col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")

#### stats on change in diversity:

# set up a new variable for the nonrecurrent group ('0):
subsum5$variable[subsum5$group_reinfection.x==c("nonrecurrent")]<-'0'
subsum5$variable[subsum5$group_reinfection.x==c("recurrent")]<-'1'
subsum5$variable[subsum5$group_reinfection.x==c("reinfection")]<-'2'
lm(formula = change.simpson ~ variable, data = subsum5)
lm(formula = change.simpson ~ variable, data = subsum5)
	#Residuals:
     #Min       1Q   Median       3Q      Max 
	#-197.181    2.168    6.608   10.774   20.235 

	#Coefficients:
    #       Estimate Std. Error t value Pr(>|t|)  
#(Intercept)   -7.628      4.211  -1.811   0.0731 .
#variable1      3.209      7.331   0.438   0.6626  
#variable2      9.418     15.755   0.598   0.5514  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 33.95 on 99 degrees of freedom
#Multiple R-squared:  0.004824,	Adjusted R-squared:  -0.01528 
#F-statistic: 0.2399 on 2 and 99 DF,  p-value: 0.7871

### sent to krishna for mixed model in SAS:

# all samples:
all_samples.div<-subsum3[, c("patientID", "group_reinfection.x", "rel_day", "sample_n", "invsimpson.x", "change.simpson")]
write.table(all_samples.div, file="all_samples.div.txt", sep="\t", quote=FALSE, col.names=NA)

# samples within 120 days:
d120<-subset(subsum3, rel_day < 120)
d120_samples.div<-subsum3[, c("patientID", "group_reinfection.x", "rel_day", "sample_n", "invsimpson.x", "change.simpson")]
write.table(d120_samples.div, file="d120_samples.div.txt", sep="\t", quote=FALSE, col.names=NA)

# samples within 120 days, no timepoint 0 (for change in simpson):
d120_samples_noindex.div<-subsum5[, c("patientID", "group_reinfection.x", "rel_day", "sample_n", "change.simpson")]
write.table(d120_samples_noindex.div, file="d120_samples_noindex.div.txt", sep="\t", quote=FALSE, col.names=NA)

### GEE model:
library(geepack)

recdata<-read.table(file="all_samples.div.txt", header=TRUE)
	# assign variable to sample based on group (y must be a number category, apparently):
recdata$variable[recdata$group_reinfection.x==c("recurrent")]<-1
recdata$variable[recdata$group_reinfection.x==c("nonrecurrent")]<-0
	# take out reinfection group (for ease):
sub.rec<-recdata[recdata$group_reinfection.x %in% c("nonrecurrent", "recurrent"),]
sub.rec<-droplevels(sub.rec)

# run geeglm
# options for corstr= exchangeable, unstructured, ar1, independence

### ar1 is probably the most relevant sorrelation structure to use, since you would expect a correlation over time, 
#but less as time goes by
#ar1:
# inverse simpson:
inv.fit.ar1 <- geeglm(variable~invsimpson.x+sample_n+invsimpson.x*sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "ar1", std.err="san.se")
summary(inv.fit.ar1)			#significant with Wald, now run anova:

# create a null model (using only time):
inv.fit.ar1 <- geeglm(variable~invsimpson.x+sample_n+invsimpson.x*sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "ar1", std.err="san.se")
inv.null <- geeglm(variable~invsimpson.x+sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "ar1", std.err="san.se")
anova(inv.null, inv.fit.ar1, test="Chisq")	
	#Analysis of 'Wald statistic' Table
		#Model 1 variable ~ invsimpson.x + sample_n + invsimpson.x * sample_n 
		#Model 2 variable ~ invsimpson.x + sample_n
 		# Df   X2 P(>|Chi|)    
		# 1  1 10.9   0.00098 ***

# for change in simpson:
chn.fit.ar1 <- geeglm(variable~change.simpson+sample_n+change.simpson*sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "ar1", std.err="san.se")
summary(chn.fit.ar1)			#not significant
chn.null <- geeglm(variable~change.simpson+sample_n+change.simpson, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "ar1", std.err="san.se")
anova(chn.null, chn.fit.ar1, test="Chisq")	
	#Analysis of 'Wald statistic' Table
		#Model 1 variable ~ change.simpson + sample_n + change.simpson * sample_n 
		#Model 2 variable ~ change.simpson + sample_n
 		# Df   X2 P(>|Chi|)    
		# 1  1 10.9   0.0083 ***

# testing model on 'reinfection' group vs 'recurrent':
# not significant, so that is good
recdata$variable[recdata$group_reinfection.x==c("nonrecurrent")]<-0
recdata$variable[recdata$group_reinfection.x==c("recurrent")]<-1
reinf<-recdata[recdata$group_reinfection.x %in% c("reinfection", "recurrent"),]
reinf$variable[reinf$group_reinfection.x==c("reinfection")]<-0
reinf<-droplevels(reinf)

reinf.ar1 <- geeglm(variable~invsimpson.x+sample_n+invsimpson.x*sample_n, family=binomial(link="logit"),
	data=reinf, id=patientID, corstr = "ar1", std.err="san.se")
summary(reinf.ar1)
null.reinf <- geeglm(variable~invsimpson.x+sample_n, family=binomial(link="logit"),
	data=reinf, id=patientID, corstr = "ar1", std.err="san.se")
anova(reinf.ar1, null.reinf, test="Chisq")
	#Analysis of 'Wald statistic' Table
		#Model 1 variable ~ invsimpson.x + sample_n + invsimpson.x * sample_n 
		#Model 2 variable ~ invsimpson.x + sample_n
  	#Df     X2 P(>|Chi|)
	#1  1 0.0264      0.87
	# not significant

reinf2.ar1 <- geeglm(variable~change.simpson+sample_n+change.simpson*sample_n, family=binomial(link="logit"),
	data=reinf, id=patientID, corstr = "ar1", std.err="san.se")
summary(reinf2.ar1)
null.reinf2 <- geeglm(variable~change.simpson+sample_n, family=binomial(link="logit"),
	data=reinf, id=patientID, corstr = "ar1", std.err="san.se")
anova(reinf2.ar1, null.reinf2, test="Chisq")

	# not significant
# if you test 'reinfection' vs. 'nonrecurrent', you also get significance
	#Analysis of 'Wald statistic' Table
		#Model 1 variable ~ change.simpson + sample_n + change.simpson * sample_n 
		#Model 2 variable ~ change.simpson + sample_n
  	#Df     X2 P(>|Chi|)
	#1  1 0.0219      0.88
	# not significant

###---
## other correlation structure options (but probably not as informative)
# using 'exchangeable':
## use this model if you think there is a relationship over time
inv.fit.exch <- geeglm(variable~invsimpson.x+sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "exchangeable", std.err="san.se")
summary(inv.fit.exch)
		# this is significant for flat invsimpson
chn.fit.exch <- geeglm(variable~change.simpson+sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "exchangeable", std.err="san.se")
summary(chn.fit.exch)
		# exchangeable is not significant for change.simpson

# using 'unstructured':
## use this model if you think there is not an effect over time:
inv.fit.unstr <- geeglm(variable~invsimpson.x+sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "unstructured", std.err="san.se")
summary(inv.fit.unstr)
		# this is not significant for flat invsimpson
chn.fit.unstr <- geeglm(variable~change.simpson+sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "unstructured", std.err="san.se")
summary(chn.fit.unstr)
		# this IS significant for change.simpson

# using independence:
## not independent over time
inv.fit.ind <- geeglm(variable~invsimpson.x+sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "independence", std.err="san.se")
summary(inv.fit.ind)			# not significant
chn.fit.ind <- geeglm(variable~change.simpson+sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "independence", std.err="san.se")
summary(chn.fit.ind)			#not significant

