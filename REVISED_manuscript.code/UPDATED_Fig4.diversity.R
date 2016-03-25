## updating diversity figures: Fig. 6 for manuscript
# Anna M. Seekatz
# 3.16.16


# figures used:
	# updated_erinsubset_summary.txt
	# ../erinfmt_extra.meta.txt
	# ../revisionvars_metadata.txt
	# updated_erinsubset_extra.summary.txt # created in this analysis

### Fig. 6A, 6C: diversity between patient groups, severity index

extra<-read.table(file="../erinfmt_extra.meta.txt", header=TRUE)
	# from this, we want to include severeidsa, abx_prior, abx_allprior, ppi, final_plating
filtered.extra<-extra[, c("seqID", "abx_prior", "abx_typeprior", "ppi", "LAB48_WBC_PEAK", "LAB48_CREAT_PEAK", "Abdominal.pain.y.n", "imaging.abnormal.y.n")]

extra2<-read.table(file='../revisionvars_metadata.txt', header=TRUE)
filtered.extra2<-extra2[, c('seqID', 'IBD')]
filtered.extra2$IBD[filtered.extra2$IBD=='TRUE']<-'1'
filtered.extra2$IBD[filtered.extra2$IBD=='FALSE']<-'0'

sum<-read.table(file='updated_erinsubset_summary.txt', header=TRUE)
subsum<-sum[sum$group_reinfection %in% c("recurrent", "nonrecurrent", "reinfection"), ]
subsum<-droplevels(subsum)
subsum2<-merge(subsum, filtered.extra, by.x="seqID", by.y="seqID", all.x=TRUE)
subsum3<-merge(subsum2, filtered.extra2, by.x="seqID", by.y="seqID", all.x=TRUE)
subsum<-subsum3
# subsum = all samples
# testing if removing these samples makes a difference:
#subsum2<-subsum[-which(subsum$seqID %in% c("DA3240", "DA3260_P", "DA3298", "DA3299", "DA3376")), ]
#subsum<-subsum2

# let's add a change in simpson diversity:
firstdiv<-subsum2[which(subsum2$index_sample_n==1), c("patientID","invsimpson")]
# now, create 'change in diversity' column:
	# combined to get a column with the diversity upon first sampling
subsum3<-merge(subsum, firstdiv, by.x="patientID", by.y="patientID")
	# get difference in diversity:
subsum3$change.simpson<-subsum3$invsimpson.y-subsum3$invsimpson.x
#write.table(subsum3, file='updated_erinsubset_extra.summary.txt', sep="\t", quote=FALSE, col.names=NA)
subsum<-subsum3

## going forward, can use the above file for all:
#subsum<-read.table(file="updated_erinsubset_extra.summary.txt", header=TRUE)

# if you want to eliminate some of the data:
subsum2<-subsum[which(subsum$index_sample_n==1), ]	# only has 'index' samples (first POSITIVE sample
subsum2.longi<-subsum2[subsum2$total_sample_n>1, ]		# only has longitudinally sampled individuals
#firstdiv2<-subsum2.longi[which(subsum2.longi$sample_n==1), c("patientID","invsimpson","group_reinfection", "clinical_result", "final_plating", "severeidsa", "abx_prior", "abx_typeprior", "ppi")]
subsum.later<-subsum[-which(subsum$index_sample_n==1), ]

# first looks:
boxplot(invsimpson.x~group_reinfection, data=subsum, ylim=c(0,20))		# all samples
boxplot(invsimpson.x~group_reinfection, data=subsum2, ylim=c(0,20))		# all index samples
boxplot(invsimpson.x~group_reinfection, data=subsum2.longi, ylim=c(0,20))	# only longitudinally sampled

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
prior.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "yes" ) {
	colorvec[i] = "red"
	}
	if ( n[i] == "none" ) {
	colorvec[i] = "blue"
	}
	}
	c(colorvec)
	}
abx.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "yes" ) {
	colorvec[i] = "magenta"
	}
	if ( n[i] == "none" ) {
	colorvec[i] = "thistle1"
	}
	}
	c(colorvec)
	}
prior.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "yes" ) {
	colorvec[i] = "red"
	}
	if ( n[i] == "none" ) {
	colorvec[i] = "blue"
	}
	}
	c(colorvec)
	}
clin.col <- function(n) {
	colorvec <- vector(mode="character", length=length(n))
	for (i in 1:length(n)) {
	colorvec[i] = "light grey"
	if ( n[i] == "negative" ) {
	colorvec[i] = "grey47"
	}
	if ( n[i] == "positive" ) {
	colorvec[i] = "magenta"
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

# plot it:
## Fig. 5A:
par(mfrow=c(1,2))
par(mai=c(0.75,0.82,0.42,0.22))
# 2A: diversity between different patient groups (first sample only):
data<-subsum[subsum$index_sample_n==1, ]		# ALL index
#data<-subsum2.longi	# only longi sampled patients (index)
plot<-plot(invsimpson.x ~ as.factor(group_reinfection), data = data, ylab=expression(paste("", gamma, " (diversity)")), 
	xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,20), cex=0.8, cex.lab=0.8, cex.axis=0.8)
points(invsimpson.x ~ jitter(as.numeric(group_reinfection, factor=0)), data = data, bg=group.col(data$group_reinfection), col="black", pch=21, cex=0.8)
names<-c("nonrecurrent", "recurrent", "reinfection")
text(x =  seq(1,3,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=0.8)
kruskal.test(invsimpson.x~group_reinfection, data=data)
text(1.5,19, labels="Kruskal-Wallis test \np=ns", cex=0.6)
	#Kruskal-Wallis rank sum test
	#data:  invsimpson by group_reinfection
	#Kruskal-Wallis chi-squared = 1.42, df = 2, p-value = 0.4914
	summary(as.factor(firstdiv$group_reinfection))		# uneven proportions
	
## Fig. 5B:
subsum$severeidsa<-as.factor(subsum$severeidsa)	#note: these are only index samples
subsum3<-subsum[subsum$severeidsa %in% c("0", "1") & subsum$index_sample_n==1, ]	# if index only
#subsum3<-subsum[subsum$severeidsa %in% c("0", "1"), ]	# if all samples
ycat<-subsum3$severeidsa
data<-subsum3
subsum3<-droplevels(subsum3)
plot<-plot(invsimpson.x ~ as.factor(ycat), data = data, ylab=expression(paste("", gamma, " (diversity)")), 
	xlab="", outline=FALSE, xaxt="n", ylim=c(0,20), cex=0.8, cex.lab=0.8, cex.axis=0.8)
points(invsimpson.x ~ jitter(as.numeric(ycat, factor=0)), data = data, bg=sev.col(data$severeidsa), cex=0.8, col="black", pch=21, mgp=c(1,2,2))
names<-c("no", "severe")
text(x =  seq(1,2,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=0.8)
wilcox.test(invsimpson.x~ycat, data=data)
	## significance!!!!!
	# all samples: W = 3325, p-value = 0.03658
	# index only: W = 1162, p-value = 0.02207
text(1,19, labels="*Wilcoxon test \np=0.022", cex=0.6)

###
#---
###

# Fig. S2A, B, C, D: Other diversity index comparisons

# Fig. S4A, severity within different groups:
par(mfrow=c(1,4))
par(mai=c(0.75,0.82,0.42,0.22))
subsum$severeidsa<-as.factor(subsum$severeidsa)
data<-subsum[subsum$severeidsa %in% c("0", "1"), ]
data<-droplevels(data)
ycat<-data$severeidsa
ycat2<-data$group_reinfection
# for now, let's do it this way:
group<-as.character(data$group_reinfection)
sev<-as.character(data$severeidsa)
data$combo<-paste(group, sev, sep="_")
data$combo<-as.factor(data$combo)
ycat<-data$combo
plot<-plot(invsimpson.x ~ as.factor(ycat), data = data, ylab=expression(paste("", gamma, " (diversity)")), xlab="", outline=FALSE, xaxt="n", 
	ylim=c(0,20), cex.lab=1, cex.axis=1)
		# although the at function can be used here to separate the groups, I can't figure out how to add the points when the 'at' is specified...
points(invsimpson.x ~ jitter(as.numeric(ycat, factor=0)), data = data, 
	bg=group.col(group_reinfection), col="black", pch=sev.pch(severeidsa), cex=0.9)
names<-as.character(unique(subsum3$group_reinfection))
text(x =  c(1.5, 3.5, 5.5), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
legend("topright", c("not severe", "severe"), col=c("black"), pch=c(21,24), cex=0.8)
kruskal.test(invsimpson.x~combo, data=data)
	#Kruskal-Wallis chi-squared = 7.8699, df = 5, p-value = 0.1636
	#let's try within groups:
data2<-data[data$combo %in% c("recurrent_0", "recurrent_1"),]
wilcox.test(invsimpson.x~combo, data=data2)		#not significant
data3<-data[data$combo %in% c("nonrecurrent_0", "nonrecurrent_1"),]
wilcox.test(invsimpson.x~combo, data=data3)
data4<-data[data$combo %in% c("reinfection_0", "reinfection_1"),]
wilcox.test(invsimpson.x~combo, data=data4)		#not significant
text(2,19.7, labels="*Kruskal-Wallis test \np=0.04", cex=0.8)
	## trending towards significant

# Fig. S4B, comparison of pos/neg samples (clinical lab):
data<-subsum[subsum$POS_NEG %in% c("negative", "positive"),]		# ALL index
plot<-plot(invsimpson.x ~ as.factor(POS_NEG), data = data, ylab=expression(paste("", gamma, " (diversity)")), 
	xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,20), cex=1, cex.lab=1, cex.axis=1)
points(invsimpson.x ~ jitter(as.numeric(POS_NEG, factor=0)), data = data, bg=clin.col(data$POS_NEG), col="black", pch=21, cex=1)
names<-c("negative", "positive")
text(x =  seq(1,2,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
kruskal.test(invsimpson.x~POS_NEG, data=data)
text(1,19.7, labels="Kruskal-Wallis test \np=ns", cex=0.8)

# Fig. S4C, comparison of abx:
data<-subsum[subsum$abx_prior %in% c("no", "yes") & subsum$index_sample_n==1, ]		# ALL index
plot<-plot(invsimpson.x ~ as.factor(abx_prior), data = data, ylab=expression(paste("", gamma, " (diversity)")), 
	xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,20), cex=1, cex.lab=1, cex.axis=1)
points(invsimpson.x ~ jitter(as.numeric(abx_prior, factor=0)), data = data, bg=abx.col(data$abx_prior), col="black", pch=21, cex=1)
names<-c("none", "prior abx")
text(x =  seq(1,2,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
kruskal.test(invsimpson.x~abx_prior, data=data)
text(1,19.7, labels="Kruskal-Wallis test \np=ns", cex=0.8)

# Fig. S4D, comparison of prior cdi:
data<-subsum[subsum$priorcdi %in% c("none", "yes") & subsum$index_sample_n==1, ]		# ALL index
plot<-plot(invsimpson.x ~ as.factor(priorcdi), data = data, ylab=expression(paste("", gamma, " (diversity)")), 
	xlab="", xaxt="n", outline=FALSE, mpg=c(2,1,0), ylim=c(0,20), cex=1, cex.lab=1, cex.axis=1)
points(invsimpson.x ~ jitter(as.numeric(priorcdi, factor=0)), data = data, bg=prior.col(data$priorcdi), col="black", pch=21, cex=1)
names<-c("none", "prior CDI")
text(x =  seq(1,2,by=1), y = par("usr")[3]-1, srt = 45, adj = 1, labels = names, xpd = TRUE, cex=1)
kruskal.test(invsimpson.x~priorcdi, data=data)
text(1,19.7, labels="Kruskal-Wallis test \np=ns", cex=0.8)

###
#---
###

### For over time graphing, need to subset groups:

# let's test out the longi first:
sum<-read.table(file='updated_erinsubset_summary.txt', header=TRUE)
subsum<-sum[sum$group_reinfection %in% c("recurrent", "nonrecurrent", "reinfection"), ]
levels(as.factor(subsum$total_sample_n))
subsum2<-subsum[subsum$total_sample_n > 1, ]
subsum2<-droplevels(subsum2)

# take first day for each patient, then calculate the change from that...
firstdiv<-subsum2[which(subsum2$index_sample_n==1), c("patientID","invsimpson","group_reinfection")]
duplicated(firstdiv$patientID)		# check for duplications

# now, create 'change in diversity' column:
	# combined to get a column with the diversity upon first sampling
subsum3<-merge(subsum2, firstdiv, by.x="patientID", by.y="patientID")
	# get difference in diversity:
subsum3$change.simpson<-subsum3$invsimpson.y-subsum3$invsimpson.x

# for some of the graphs, we want to eliminate the first samples, since that is 0 for change.simpson
subsum4<-subsum3[-which(subsum3$index_sample_n==1), ]

### Fig. S3: 
	#A: diversity over time, by sampling period, all samples:
	#B: Change in diversity per patient, by sampling period, ignore first change.simpson point:

# S3A: simpson over time:
recur<-subset(subsum3, (subset=group_reinfection.x %in% c("recurrent")))				#subset data by group
non<-subset(subsum3, (subset=group_reinfection.x %in% c("nonrecurrent")))
reinf<-subset(subsum3, (subset=group_reinfection.x %in% c("reinfection")))
par(mfrow=c(2,1))
par(mai=c(0.75,0.82,0.42,0.22))
# S3A: flat diversity over time
plot(jitter(recur$index_sample_n, amount=0.05), recur$invsimpson.x, mgp=c(2,1,0),
	pch=21, col="black", bg="orange", cex=0.8, ylim=c(0,50), cex.lab=0.8, cex.axis=0.8,
	main="", ylab=expression(paste("", gamma, " (diversity)")), xlab="relative sampling")
points(0.25+jitter(non$index_sample_n, amount=0.05), non$invsimpson.x, pch=21, xaxt='n', col="black", bg="chartreuse3", cex=0.8)
points(0.5+jitter(reinf$index_sample_n, amount=0.05), reinf$invsimpson.x, pch=21, xaxt='n', col="black", bg="darkgoldenrod", cex=0.8)
recline <- lm(invsimpson.x~index_sample_n, data=recur)
nonline <- lm(invsimpson.x~index_sample_n, data=non)
reinfline <- lm(invsimpson.x~index_sample_n, data=reinf)
abline(nonline, col="chartreuse3")
abline(recline, col="orange")
abline(reinfline, col="darkgoldenrod")
legend("topright", c("recurrent", "nonrecurrent", "reinfected"), col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")
	### this still looks good for recurrent/nonrecurrent patients!
	
# S3B: CHANGE in simpson over sampling, excluding first sampling:
recur2<-subset(subsum4, (subset=group_reinfection.x %in% c("recurrent")))				#subset data by group
non2<-subset(subsum4, (subset=group_reinfection.x %in% c("nonrecurrent")))
reinf2<-subset(subsum4, (subset=group_reinfection.x %in% c("reinfection")))
plot(jitter(recur2$index_sample_n, amount=0.05), recur2$change.simpson, cex.lab=0.8, cex.axis=0.8,
	pch=21, col="black", bg="orange", cex=0.8, ylim=c(-20,25), mgp=c(2,1,0),
	main="", ylab=expression(paste("", Delta, gamma, " (", Delta, "invsimpson)")), xlab="relative sampling")
points(0.25+jitter(non2$index_sample_n, amount=0.05), non2$change.simpson, pch=21, xaxt='n', col="black", bg="chartreuse3", cex=0.8)
points(0.5+jitter(reinf2$index_sample_n, amount=0.05), reinf2$change.simpson, pch=21, xaxt='n', col="black", bg="darkgoldenrod", cex=0.8)
recline2 <- lm(change.simpson~sample_n, data=recur2)
nonline2 <- lm(change.simpson~sample_n, data=non2)
reinfline2 <- lm(change.simpson~sample_n, data=reinf2)
abline(nonline2, col="chartreuse3")
abline(recline2, col="orange")
abline(reinfline2, col="darkgoldenrod")
legend("bottomright", c("recurrent", "nonrecurrent", "reinfected"), col=c("chartreuse3", "orange", "darkgoldenrod"), pch=19, cex=0.6, bg="white")
	### looks OK...fairly similar, line is below, however
	
#### Note: did not include an over time graph of the data in the manuscript: it was too hard to visualize the trend lines
### however, we did include the GEE modeling results below:
	
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


### GEE model:
library(geepack)

recdata<-read.table(file="updated_all_samples.div.txt", header=TRUE)
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
inv.fit.ar1 <- geeglm(variable~invsimpson.x+index_sample_n+invsimpson.x*index_sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "ar1", std.err="san.se")
summary(inv.fit.ar1)			#significant with Wald, now run anova:

# create a null model (using only time):
inv.fit.ar1 <- geeglm(variable~invsimpson.x+index_sample_n+invsimpson.x*index_sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "ar1", std.err="san.se")
inv.null <- geeglm(variable~invsimpson.x+index_sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "ar1", std.err="san.se")
anova(inv.null, inv.fit.ar1, test="Chisq")	
	#Analysis of 'Wald statistic' Table
		#Model 1 variable ~ invsimpson.x + index_sample_n + invsimpson.x * index_sample_n 
		#Model 2 variable ~ invsimpson.x + index_sample_n
 		# Df   X2 P(>|Chi|)    
		#1  1 10.3    0.0013 **

# for change in simpson:
# take out index samples (to look at just change):
#sub.rec2<-sub.rec[-which(sub.rec$index_sample_n==1), ]
chn.fit.ar1 <- geeglm(variable~change.simpson+index_sample_n+change.simpson*index_sample_n, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "ar1", std.err="san.se")
summary(chn.fit.ar1)			#not significant
chn.null <- geeglm(variable~change.simpson+index_sample_n+change.simpson, family=binomial(link="logit"),
	data=sub.rec, id=patientID, corstr = "ar1", std.err="san.se")
anova(chn.null, chn.fit.ar1, test="Chisq")	
	#Analysis of 'Wald statistic' Table
		#Model 1 variable ~ change.simpson + sample_n + change.simpson * sample_n 
		#Model 2 variable ~ change.simpson + sample_n
 		# Df   X2 P(>|Chi|)    
		# 1  1 8.21    0.16
		# not significant with change anymore....

# testing model on 'reinfection' group vs 'recurrent':
# not significant, so that is good
recdata$variable[recdata$group_reinfection.x==c("nonrecurrent")]<-0
recdata$variable[recdata$group_reinfection.x==c("recurrent")]<-1
reinf<-recdata[recdata$group_reinfection.x %in% c("reinfection", "recurrent"),]
reinf$variable[reinf$group_reinfection.x==c("reinfection")]<-0
reinf<-droplevels(reinf)

reinf.ar1 <- geeglm(variable~invsimpson.x+index_sample_n+invsimpson.x*index_sample_n, family=binomial(link="logit"),
	data=reinf, id=patientID, corstr = "ar1", std.err="san.se")
summary(reinf.ar1)
null.reinf <- geeglm(variable~invsimpson.x+index_sample_n, family=binomial(link="logit"),
	data=reinf, id=patientID, corstr = "ar1", std.err="san.se")
anova(reinf.ar1, null.reinf, test="Chisq")
	#Analysis of 'Wald statistic' Table
		#Model 1 variable ~ invsimpson.x + index_sample_n + invsimpson.x * index_sample_n 
		#Model 2 variable ~ invsimpson.x + index_sample_n
  	#Df     X2 P(>|Chi|)
	#1  1 0.378      0.54
	# not significant

reinf2.ar1 <- geeglm(variable~change.simpson+index_sample_n+change.simpson*index_sample_n, family=binomial(link="logit"),
	data=reinf, id=patientID, corstr = "ar1", std.err="san.se")
summary(reinf2.ar1)
null.reinf2 <- geeglm(variable~change.simpson+index_sample_n, family=binomial(link="logit"),
	data=reinf, id=patientID, corstr = "ar1", std.err="san.se")
anova(reinf2.ar1, null.reinf2, test="Chisq")

	# not significant
# if you test 'reinfection' vs. 'nonrecurrent', you also get significance
	#Analysis of 'Wald statistic' Table
		#Model 1 variable ~ change.simpson + sample_n + change.simpson * sample_n 
		#Model 2 variable ~ change.simpson + sample_n
  	#Df     X2 P(>|Chi|)
	#1  1 0.0795      0.78
	# not significant

