## Collating files from mothur results
# 1.4.16
# Anna M. Seekatz

# the following files and directories were used to generate the figures in the manuscript '':

# description of directories:
	# Data.files_erinsubset: contains metafiles, filtered/edited mothur-generated files
	# Mothurfiles_erinsubset: contains raw mothur results produced in mothur
		#including the 'mbatch' file, which contains the mothur commands
	# Rcode_erinsubset: contains code for each figure
	# Figures_erinsubset: contains several reiterations of 

# Fig. 2 (diversity):	
	# erinfmt_summary.txt
	# erinfmt_extra.meta.txt
	
# Fig. 3 (PAM-clustering and heatmap):
# part 1 (creating a phylotype file):
	# erinfmt.trim.contigs.good.unique.good.filter.unique.precluster.pick.rdp.wang.tax.summary.txt: mothur-generated file
		# erinfmt_genbar.txt: created from the mothur-generated file to only include genus-level assignment
		# erinfmt_genfrac.all_ordered.txt: contains top genera (98% of the total sequences)
	# ERIN.UMFMT_metadata.txt: metadata file
# part 2: Clustering:
	# erinsubset.0.03.shared: use only to get list of names; OR, directly use this file if doing OTU-based community type clustering
	# erinfmt_genfrac.all_w.meta.txt: your phylotype relative abundance
# part 3: Heatmap:
	# suberin_summary.phylotyped.txt: contains phylotype information, clustering info, some metadata
	# erinsubset.0.03.shared: mothur-generated .shared file of subsetted samples
	# ERIN.UMFMT_metadata_filtered.txt: an older metafile (has other samples in it, as well)
	# erinfmt.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.0.03.cons.taxonomy: original taxonomy identification file
		# later produced the erinfmt.new.taxonomy.names.txt that is used in other analyses
		# suberin.relOTUs.txt: created from this file (relative abundance of most common OTUs)
	
# Fig. 4 (PCOA biplot):
	# erinsubset.thetayc0.03.pcoa.axes: mothur-generated
	# erinsubset.thetayc0.03.nmds.axes: mothur-generated
	# erinfmt.unique_list.groups.summary.txt: mothur-generated (from a larger batch)
	# ERIN.UMFMT_metadata_filtered.txt: an older metafile (has other samples in it, as well)
	# erinfmt_extra.meta.txt: an extra metadata file, queried from available samples
	# erinsubset_summary.txt #created from the mothur files
	
# Fig. 5 (dotplot of lefse results):
	# suberin_all.alpha.phylo.txt: contains summary and phylotype information collated previously
	# suberin.relOTUs.txt: relative abundance of OTUs, filtered
	# erinfmt.new.taxonomy.names.txt: edited taxonomy files with OTU info
	# erinsubset_all.lefse.results.txt: lefse results, edited from mothur
	
# Fig. 6 (thetayc over time and between clinical conversion):
	# erinsubset.0.03.summary: mothur-generated file from summary.shared (note: fix columns before reading in)
	# erinsub_all.extrameta.summary.txt: generated previously 
		# erinsub_ind.dist.byindex.txt: (created)
		# erinsub_shared.meta.dist.txt: (created)
		# erinsub_ind.dist.overtime.txt: (created)
