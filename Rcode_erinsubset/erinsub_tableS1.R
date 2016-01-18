### Make metadata file for Table S1:
# 1.12.16
# Anna M. Seekatz

## Read in metadata file:
sum<-read.table(file="suberin_summary.phylotyped.txt", header=TRUE)
## extra metadata:
extra<-read.table(file="erinfmt_extra.meta.txt", header=TRUE)
	# combine:
sum2<-merge(sum, extra, by.x="sampleID", by.y="seqID", all.x=TRUE)
	# get relevant metadata:
sum3<-sum2[, c("sampleID", "patientID.x", "index.x", "recur.x", "group_reinfection", 
	"sample_n.x", "clinical_result.x", "age.x", "sex.x", "rel_day.x", "phylopam4_cluster", 
	"invsimpson", "severeidsa", "abx_prior")]

# read in OTUs (most abundant):
otus<-read.table(file="suberin.relOTUs.txt", header=TRUE)
	# this is quite large-lets just get the top 100:
topotus<- otus[, order(-colSums(otus))]	
top100<-topotus[, 1:100]
top100.sorted<-top100[, order(colnames(top100))]

# let's also add the taxonomic info:
tax<-read.table(file="erinfmt.new.taxonomy.names.txt", header=TRUE)
tax$OTU<-as.character(tax$OTU)
	# get list of top100 otus:
topotu.list<-colnames(top100.sorted)												#define the OTUs you want to keep (those in the top 98%)
filtered.tax<-tax[which(tax$OTU %in% topotu.list) , ]							#filter out the taxonomy metadata file
filtered.tax<-droplevels(filtered.tax)
colnames(top100.sorted)<-filtered.tax$taxname
top100.sorted$sampleID<-rownames(top100.sorted)

	# merge with the selected summary:
sum4<-merge(sum3, top100.sorted, by="sampleID")
#write.table(sum4, file="TableS1_sample.descr.txt", sep="\t", quote=FALSE, col.names=NA)

## Creating Table 1 (patient metadata):

## all samples together:
# n total samples in each group:
summary(sum$total_sample_n)
   #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #1.00    2.00    3.00    2.99    4.00    7.00
summary(sum$group_reinfection)
	#nonrecurrent    recurrent  reinfection 
    #     138           66           26 

prop.table(table(sum$group_reinfection, sum$total_sample_n))
	# gives prop
	
# let's split the groups, then explore within:
groups<-split(sum, sum$group_reinfection)
rec<-groups$'recurrent'
non<-groups$'nonrecurrent'
reinf<-groups$'reinfection'

##recurrent group:
length(unique(rec$patientID))
	#[1] 19
summary(rec$sex)
 	#F  M 
	#49 17 
summary(rec$total_sample_n)
   #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   #2.00    3.00    4.00    3.97    4.00    7.00
   
##nonrecurrent group:
length(unique(non$patientID))
	#[1] 70
summary(non$sex)
 	#F  M 
	#70 68
summary(non$total_sample_n)
   #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 1.00    2.00    2.00    2.47    3.00    5.00 
   
##reinfected group:
length(unique(reinf$patientID))
	#[1] 9
summary(reinf$sex)
 	#F  M 
	#17  9
summary(reinf$total_sample_n)
   #Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   # 2.00    3.00    3.00    3.27    4.00    4.00 

# to get no. patients within each group that are male/female:
data <- ddply(sum, c("group_reinfection", "sex"), summarise,
               N    = length(unique(patientID))
               )
               
data <- ddply(sum, c("group_reinfection"), summarise,
               N    = length(rel_day),
               min = min(rel_day),
               mean = mean(rel_day),
               max = max(rel_day),
               median = median(rel_day) )
# didn't really work

# do this:
sum0<-subset(sum, rel_day>0)		# to get real minimum
groups<-split(sum0, sum0$group_reinfection)
rec0<-groups$'recurrent'
non0<-groups$'nonrecurrent'
reinf0<-groups$'reinfection'

summary(sum0$rel_day)
summary(non0$rel_day)
summary(rec0$rel_day)
summary(reinf0$rel_day)
