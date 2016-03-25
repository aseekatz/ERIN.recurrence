## Updating: PAM-clustering of samples and making a heatmap
# 3.16.16
# Anna M. Seekatz

##-----------------------------------------
#### Part 1: creating a phylotype file:
	# note: this was done previously with a larger data set
	# please contact me at aseekatz@umich.edu for more info on this step

# phylotype analysis

# files used:
	# erinfmt.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary
	# erinfmt_genbar.txt					# only using the data
	# erinfmt_genfrac2p_ordered.txt			# only using the data
	# updated_erinsubset_extra.summary.txt	# new file for metadata

#used mothur file:
	# erinfmt.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary
	# as erinfmt_phylotype.xlsx
	# in 'genus_only' sheet:
		# removed all phylotypes with < 100 seqs
		# added unclassified info (to whatever level it classified to)
		# added phylum info
		# ordered by phylum, total
		# saved as 'erinfmt_genbar.txt'

# in R:
genus<-read.table(file="../../erinfmt_genbar.txt", header=TRUE, row.names=3)
genus.df<-subset(genus, select=-c(taxlevel, phylum, rankID, total))
genus.df<-as.data.frame(t(genus.df))
genus.fr<-genus.df/rowSums(genus.df)*100
genus.fr<-t(genus.fr)
#write.table(genus.fr, file="erinfmt_genfrac.txt", quote=FALSE, sep="\t", col.names=NA)
genus1<- genus.fr[rowSums(genus.fr>=1)>=1,] #will list all genera above 1%
genus2<- genus.fr[rowSums(genus.fr>=2)>=2,] #will list all genera above 2%
	# note: still about 100 phylotypes! yikes
#write.table(genus2, 'erinfmt_genfrac2p.txt',quote=FALSE,sep="\t", col.names=NA)

# combine with metadata:
meta<-read.table(file="updated_erinsubset_extra.summary.txt", header=TRUE)
genbar<-read.table(file="../../erinfmt_genfrac.txt", header=TRUE, row.names=1)
	#genbar<-subset(rm_g, select =-c(color) )
	barg<-as.data.frame(t(genbar))
	barg$sampleID<-rownames(barg)
	#col.gen<-as.character(rm_g$color)
# test that samples match...
	#test.barg<-barg[order(barg$sampleID),]
	#test.meta<-meta[order(meta$seqID),]
	#unique(test.meta$sampleID)						#these match, but why doesn't seqID match?
	#unique(test.barg$sampleID)			

bar<-merge(barg, meta, by.y=c("seqID"), by.x=c("sampleID"))
	#note: you MUST put barg first! otherwise, it will merge incorrectly
bar_ordered<- bar[order(bar$total_sample_n, bar$group_reinfection, bar$patientID, bar$sampleID, bar$sample_n),]
#write.table(bar_ordered, 'updated_erinfmt_genfrac.all_w.meta.txt',quote=FALSE,sep="\t", col.names=NA)



##-----------------------------------------
#### Part 2: PAM-clustering of samples
# based on Jensen-Shannon index

# files used:
	# erinsubset.0.03.shared
	# updated_erinfmt_genfrac.all_w.meta.txt
	# updated_suberin_summary.phylotyped.txt 	#created from this section
		# if using OTU file, not phylotype file: your .shared file
	
library(cluster)
library(vegan)
library(labdsv)

names.subset<-read.table(file="../erinsubset_mothurfiles/erinsubset.0.03.shared", header=TRUE, row.names=2)
names<-as.character(rownames(names.subset))		# these are all of the subset samples
phylos<-read.table(file="updated_erinfmt_genfrac.all_w.meta.txt", header=TRUE)
	# note: this was a file recreated from your 'phylotype' file, or 
rownames(phylos)<-phylos$sampleID
sub.phylos.df<-phylos[which(phylos$sampleID %in% names), 2:277]
	# if you want to do this by OTU, use your .shared file directly:
#dd.all<- subset(names.subset, select = -c(label, numOtus ))					#removing the columns you do not need
#dd.all<-dd.all[,colSums(dd.all) > 100]										#this removes all OTUs that are present at < 100 seqs
#sub.dd.frac<-dd.all/rowSums(dd.all)

# code for jensen-shannon distance:
all.pairs.jensen.shannon <- function(d) {
  kl.divergence <- function(pq, midpoint) {
    sum(ifelse(pq == 0.0, 0.0, pq * log(pq / midpoint)))
  }
 
  jensen.shannon <- function(row1, row2) {
    midpoint <- 0.5 * (row1 + row2)
    kl.divergence(row1, midpoint) + kl.divergence(row2, midpoint)
  }
 
  d <- as.matrix(d)
  num.rows <- nrow(d)
  result <- matrix(0, nrow=num.rows, ncol=num.rows, dimnames=list(rownames(d), rownames(d)))
  for (j in 1:(num.rows - 1)) {
    for (i in (j + 1):num.rows) {
      distance <- jensen.shannon(c(d[i, ]), c(d[j, ]))
      result[i, j] <- distance
      result[j, i] <- distance
    }
  }
  as.dist(result)
 }

# run your data through the code:
sub.js.phylo <- all.pairs.jensen.shannon(sub.phylos.df)
#sub.js.otu <- all.pairs.jensen.shannon(sub.dd.frac)		#if using shared file

#for phylos:
postscript('updated_ALLphylos.sub.erin.phylo_pamtest.ps',width=8,height=8,horizontal=TRUE)
for (nc in 2:10) {
plot(pam(sub.js.phylo, k=nc))
}
dev.off()
	#for OTUs:
#postscript('test.sub.erin.otu_pamtest.ps',width=8,height=8,horizontal=TRUE)
#for (nc in 2:10) {
#plot(pam(js.otu, k=nc))
#}
#dev.off()
	# looking at the data, an n=4 has the highest silhouette score (mean score=0.49)--not bad!

# create a phylotyped file
# if you wanted to add pam=2 cluster just in case (silhouette score=0.37)
phylotype2 <- pam(sub.js.phylo, k=2)										#the k= , will give you however many clusters you are defining
phylo2<-as.data.frame(phylotype2$silinfo$widths[,1])   											#this will list the samples, their assigned cluster, and the silhouette info for each sample
phylo2$SEQ_NAME<-rownames(phylo2)
colnames(phylo2)[1] <- c("phylopam2_cluster")
suberin.summary<-merge(phylos, phylo2, by.x=c("sampleID"), by.y=c("SEQ_NAME"))
phylotype4 <- pam(sub.js.phylo, k=4)										#the k= , will give you however many clusters you are defining
phylo4<-as.data.frame(phylotype4$silinfo$widths[,1])   											#this will list the samples, their assigned cluster, and the silhouette info for each sample
phylo4$SEQ_NAME<-rownames(phylo4)
colnames(phylo4)[1] <- c("phylopam4_cluster")
suberin.summary2<-merge(suberin.summary, phylo4, by.x=c("sampleID"), by.y=c("SEQ_NAME"))
#write.table(suberin.summary2, file="updated_suberin_summary.phylotyped.txt", sep="\t", quote=FALSE, col.names=NA)


##-----------------------------------------
#### Part 3: Heatmap of data (Figure 2)
	# note: some of the sample names in the .taxonomy file and .shared files were processed together originally
	# thus, there is a subset step at some point to eliminate these samples
# files used:
	# updated_suberin_summary.phylotyped.txt
	# erinsubset.0.03.shared					# data
	# updated_erinfmt_genfrac.all_w.meta.txt	# more meta, if needed
	# erinfmt_otu.taxonomy.txt (since I subsetted data from the shared file, they should all be the same)
	# erinfmt.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy
		# later produced the erinfmt.new.taxonomy.names.txt that is used in other analyses
		# suberin.relOTUs.txt: created from this file (relative abundance of most common OTUs)

library(RColorBrewer)
library(gplots)
library(vegan)
library(plyr)

sum<-read.table(file='updated_suberin_summary.phylotyped.txt', header=TRUE)
shared<-read.table(file="../erinsubset_mothurfiles/erinsubset.0.03.shared", header=TRUE, row.names=2)
shared2<-shared[-which(rownames(shared) %in% c("DA3240", "DA3260_P", "DA3298", "DA3299", "DA3376")), ]
	# eliminated samples that are unconfirmed...
shared<-shared2

otu<-subset(shared, select =-c(label, numOtus) )								#gets rid of extra columns in R
otu.filtered<-otu[,which(colSums(otu)>=100)]									#eliminates any OTUs present at less than a count of 100
otu.filtered<- otu.filtered[ order(row.names(otu.filtered)), ]					#order by the rownames (your samples)
otu.matrix<-as.matrix(otu.filtered)
	
	# get relative abundance of file instead of raw counts:
#otu.rel<-otu.matrix/rowSums(otu.matrix)
otu.rel<-(otu.matrix/rowSums(otu.matrix))*100
#write.table(otu.rel, file="suberin.relOTUs.txt", sep="\t", quote=FALSE, col.names=NA)

	# to filter out OTUs that are present in less than x%:
otu.rel.max<-apply(otu.rel,2,max)
otu.rel.filtered<-otu.rel[,otu.rel.max>0.020]
	#median, if you wanted to do that...
otu.rel.med<-apply(otu.rel,2,median)
otu.rel.filtered.med<-otu.rel[,otu.rel.med>0.0001]

# graph it:			
my.col <- colorRampPalette(c("aliceblue","blue"))(n = 249)			#you can put in as many colors as you want here
my.breaks = c(seq(0,0.001,length=50),  								#just remember: your breaks must always = 1 more than the number of colors you are defining!
               seq(0.001,0.01,length=50),
               seq(0.01,0.10,length=50),
               seq(0.10,0.50,length=50),
               seq(0.50,1,length=50))
#heatmap.2(otu.rel.filtered, 
			col=my.col, 
			breaks=my.breaks, 
			cexRow=0.5, 
			cexCol=0.5, 
			trace="none"
			)
# there are some interesting little patches in here...			
			
	# heatmap of top 30 (?) OTUs:
topotus<- otu.rel.filtered[, order(-colSums(otu.rel.filtered))]	
top30<-topotus[, 1:30]
top40<-topotus[, 1:40]

my.col <- colorRampPalette(c("aliceblue","blue"))(n = 249)			#you can put in as many colors as you want here
my.breaks = c(seq(0,0.1,length=50),  								#just remember: your breaks must always = 1 more than the number of colors you are defining!
               seq(0.1,1,length=50),
               seq(1,10,length=50),
               seq(10,50,length=50),
               seq(50,100,length=50))
#heatmap.2(top30, 
			col=my.col, 
			breaks=my.breaks, 
			trace="none", 
			density.info="none", 
			Colv=F, 			
			Rowv=T,
			dendrogram="none",
			#RowSideColors = as.character(sum2$color1),
			#ColSideColors = as.character(phyla.col)
			)
#legend("left", c("control","patient"), col=c("red","orange"), pch=15, cex=0.8)

# add row metadata:
sum<-read.table(file="updated_suberin_summary.phylotyped.txt", header=TRUE)
# need to replace some NA's with 'not_tested':
levels(sum$POS_NEG) <- c(levels(sum$POS_NEG), "not_tested")
sum$POS_NEG[is.na(sum$POS_NEG)]<-'not_tested'
# let's also replace some others in the final_plating
levels(sum$final_plating) <- c(levels(sum$final_plating), "unknown")
sum$final_plating[is.na(sum$final_plating)]<-'unknown'
sum$final_plating[sum$final_plating==c("positive_maybe")]<-'unknown'
sum$final_plating[sum$final_plating==c("No_sample")]<-'unknown'

sum<- sum[order(sum$sampleID),]									#remember that the order of this must equal the order of the matrix read into the heatmap command
samples<-rownames(top40)										#this is the list of sample names in the heatmap matrix
sum2<-sum[sum$sampleID %in% samples , ]							#eliminate any samples from the original meta file to get the exact list (since not doing a merge)
top40<-top40[order(rownames(top40)), ]							#must be in same order
cbind(rownames(top40), sum2$sampleID)									#check that this is right

# giving some options to color by:
myCol4 <- colorRampPalette(brewer.pal(4,"BrBG"))(4)			# only for more than 3 colors
myCol3<-c("yellow2", "blue2")
sum2$color<-mapvalues(sum2$group_reinfection, from = c("nonrecurrent", "recurrent", "reinfection"), to = c("chartreuse3", "orange", "darkgoldenrod"))
sum2$phylopam2_cluster<-as.factor(sum2$phylopam2_cluster)
sum2$color2<-mapvalues(sum2$phylopam2_cluster, from = c("1", "2"), to = myCol3)
sum2$phylopam4_cluster<-as.factor(sum2$phylopam4_cluster)
sum2$color4<-mapvalues(sum2$phylopam4_cluster, from = c("1", "2", "3", "4"), to = myCol4)
#sum2$index<-as.factor(sum2$index)
#sum2$color.index<-mapvalues(sum2$index, from = c("0", "1"), to = c("green", "red"))
#sum2$recur<-as.factor(sum2$recur)
#sum2$color.recur<-mapvalues(sum2$recur, from = c("0", "1"), to = c("green", "red"))
#sum2$patient_status<-paste(sum2$group2, sum2$index, sep="_")
#sum2$color.status<-mapvalues(sum2$patient_status, from = c("nonrecurrent_0", "nonrecurrent_1", "recurrent_0", "recurrent_1"), to = c("lightgreen", "lightpink", "green3", "red"))
#sum2$patient_status<-paste(sum2$group2, sum2$index, sep="_")
sum2$color.clinical<-mapvalues(sum2$POS_NEG, from = c("negative", "positive", "not_tested"), to = c("grey47","magenta", "black"))
sum2$color.plating<-mapvalues(sum2$final_plating, from = c("negative", "positive", "unknown"), to = c("grey47","magenta", "black"))
sum2$severeidsa[is.na(sum2$severeidsa)]<-"unknown"
sum2$severeidsa<-as.factor(sum2$severeidsa)
sum2$color.severe<-mapvalues(sum2$severeidsa, from = c("0", "1", "unknown"), to = c("blue","red", "black"))
sum2$col.status<-mapvalues(sum2$sample_status, from = c("index", "recovery", "recurrence", "reinfection", "treatment"), to = c("red", "chartreuse3", "orange", "darkgoldenbrown", "pink"))

# add column metadata (taxonomy info):
taxonomy_file<-read.table(file="../../mothur.files_3.18.15/erinfmt.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy", header=TRUE)
tax <- taxonomy_file$Taxonomy
tax <- gsub("\\(\\d*\\)", "", tax)
tax <- gsub(";unclassified", "", tax)
tax <- gsub("_1", "", tax)
tax <- gsub(";$", "", tax)
tax <- gsub("/.*", "", tax)
tax <- gsub(".*;", "", tax)
tax.names <-paste(taxonomy_file$OTU, tax)
tax.names <-gsub("000", "", tax.names)
taxonomy_file$taxname<-tax.names
phylum <- taxonomy_file$Taxonomy
phylum <- gsub("\\(\\d*\\)", "", phylum)
phylum <- gsub("Bacteria;", "", phylum)
phylum <- gsub(";$", "", phylum)
phylum <- gsub(";.*", "", phylum)
taxonomy_file$phylum<-phylum
tax<-taxonomy_file
tax<-tax[order(tax$OTU),]
#write.table(tax, file="erinfmt.new.taxonomy.names.txt", sep="\t", quote=FALSE, col.names=NA)
#tax<-read.table(file="erinfmt.new.taxonomy.names.txt", header=TRUE)
tax$taxname <- gsub("000", "", tax$taxname)
tax$taxname <- gsub("_", " ", tax$taxname)

topotu.list<-colnames(top40)												#define the OTUs you want to keep (those in the top 98%)
filtered.tax<-tax[tax$OTU %in% topotu.list , ]							#filter out the taxonomy metadata file
filtered.tax<-droplevels(filtered.tax)										#sometimes R keeps levels that were discarded--this gets rid of them completely

filtered.tax$phylum<-factor(filtered.tax$phylum, levels=c("Bacteroidetes", "Firmicutes", "Proteobacteria", "Actinobacteria", "Verrucomicrobia", "Fusobacteria", "unclassified"))
tax.col<- filtered.tax[order(filtered.tax$phylum, -filtered.tax$Size) , ]
tax.col$color<-mapvalues(tax.col$phylum, from = c("Bacteroidetes", "Firmicutes", "Proteobacteria", "Actinobacteria", "Verrucomicrobia", "Fusobacteria", "unclassified"), to = c("green4", "dodgerblue1","gold", "firebrick1", "firebrick4", "hotpink", "grey47"))
tax.col2<-t(tax.col)
phyla.col<-tax.col2[6, ]

col.order<-as.character(tax.col2[1, ])								#convert the order of the OTUs into a list
top40<-top40[,col.order]											#order your matrix by the ordered list
rbind(colnames(top40), tax.col2)									#check order of the matrix and list to ensure that colors will be correct (they should match!)

# by group:
top40<-top40[order(sum2$group_reinfection), ]
sum2<-sum2[order(sum2$group_reinfection), ]
# OR by sample status:
#top40<-top40[order(sum2$sample_status), ]
#sum2<-sum2[order(sum2$sample_status), ]
# OR by phylotype:
top40<-top40[order(sum2$phylopam4_cluster), ]
sum2<-sum2[order(sum2$phylopam4_cluster), ]

# graph--subset samples:
heatmap.2(top40,
                    notecol="black",     
                    density.info="none",
                    key.xlab="",
                    key.title="",  
                    trace="none",         
                    margins =c(2,4),    
                    col=my.col,        
                    breaks=my.breaks,    
                    dendrogram="none",    
                    Colv=F,
                    Rowv=F,
                    srtCol=45,
                    cexRow= 0.7,
                    cexCol = 0.8,
                    keysize=1,
                    lwid = c(3,5),
                    lhei = c(2,5),
                    labCol=tax.col$taxname,
                    labRow="",
					RowSideColors = as.character(sum2$color.severe),
					ColSideColors = as.character(phyla.col)
					)
legend("topleft",legend=c("recurrent", "nonrecurrent", "reinfected"), col=c("orange", "chartreuse3", "darkgoldenrod"), cex=0.6, pch=19)
#legend("topleft",legend=rep(1:4, 1), col=myCol4, cex=0.8, pch=19)
legend("bottomleft", c("Bacteroidetes", "Firmicutes", "Proteobacteria", "Actinobacteria", "Verrucomicrobia", "Fusobacteria", "unclassified"), 
	col=c("green4", "dodgerblue1","gold", "firebrick1", "firebrick4", "hotpink", "grey47"), pch=15, cex=0.6)
legend("topright",legend=c("negative", "positive", "unknown"), col=c("grey47","magenta", "grey87"), cex=0.6, pch=19)
legend("topmiddle",legend=c("not severe", "severe", "unknown"), col=c("blue","red", "black"), cex=0.6, pch=19)
