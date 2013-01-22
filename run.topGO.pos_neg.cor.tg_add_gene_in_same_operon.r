# for positive and negative correlation: need to change the background genelist (geneNames)
# load("/users/slo/R/work/2012-08-03-DataforCorrelation.RData") # load EO term, Lrpnames, ratios
load("2012-12-05-Result.for.Correlation_68.environmental.EOs_Lrp.fold.change.1.25_p.val.0.05_median.r.NA_117EOs_up.250.bp_in.coding.50.bp.RData")
load("operons.RData")
source("Modified_MeDiChI_functions.R")
source("medichi.analysis.misc.r")

target.type = "correlated.ChIP.targets.overlapwith.KO"
# target.type = "correlated.ChIP.targets"


median.r.cutoff = NA
fold.change = 1.25
p.adjusted.cutoff = 0.05
dir.name = paste(substr(Sys.time(),1,10),"-topGO.68.Env.EO_Lrp.fold.change.",fold.change,"_p.cutoff.0.05",sep="")

# load(gene.list.fl) # load ChIP-chip datasets for geneList

library(topGO)
vng2GO <- readMappings(file="/users/slo/R/work/halo.go.list.txt") #load halo GO file
# vng2GO <- readMappings(file="U:/R/work/halo.go.list.txt")

##################
# set directory
##################
# dir.name = paste(substr(Sys.time(),1,10),"-topGO.117EO_Lrp.fold.change.2_p.cutoff.0.05_median.r.0.5",sep="")
dir.create(dir.name)

setwd(paste("/users/slo/R/work/",dir.name,sep=""))
# setwd(paste("U:/R/work/",dir.name,sep=""))
set_dir1 = getwd()

dir.create(paste(set_dir1,"/",substr(Sys.time(),1,10),"_",target.type,sep=""))
set_dir = paste(set_dir1,"/",substr(Sys.time(),1,10),"_",target.type,sep="")
setwd(set_dir)
###################
# genelists
###################
pos.genelist = get(paste("pos.",target.type,sep=""))
neg.genelist = get(paste("neg.",target.type,sep=""))

###################
# include operons
###################
corr.coeff.matrix = sapply(Lrpnames, function(lrp)
						{
							sapply(rownames(ratios.norm), function(gene)
							{
								cor(ratios.norm[lrp,],ratios.norm[gene,] )
							})
						})

#---positive ----
pos.genelist1 = list()
for (i in 1: length(pos.genelist))
{
	if (length(pos.genelist[[i]])>0)
	{
		each.lrp <- list()
		count <- 1
		for (j in 1: length(pos.genelist[[i]]))
		{
			if (length(pos.genelist[[i]][[j]])>1 | (length(pos.genelist[[i]][[j]])==1 && !is.na(pos.genelist[[i]][[j]])) )
			{
				each.lrp[[count]] <- unique(unlist(get.all.operons.list_Lrp.tg(pos.genelist[[i]][[j]],Lrp = names(pos.genelist)[i],
										Lrp.ChIP.targets,corr.coeff.matrix,cor.cutoff = 0.2,DirectionOfCorrelation="pos")))	
				names(each.lrp)[count] <- names(pos.genelist[[i]])[j]
				count <- count +1 ;
			}
		}
		pos.genelist1[[i]] <- each.lrp	
	} else 
	{
		pos.genelist1[[i]] <- pos.genelist[[i]]
	}
}
names(pos.genelist1) <- names(pos.genelist)

#---negative ----
neg.genelist1 = list()
for (i in 1: length(neg.genelist))
{
	if (length(neg.genelist[[i]])>0)
	{
		each.lrp <- list()
		count <- 1
		for (j in 1: length(neg.genelist[[i]]))
		{
			if (length(neg.genelist[[i]][[j]])>1 | (length(neg.genelist[[i]][[j]])==1 && !is.na(neg.genelist[[i]][[j]])) )
			{
				each.lrp[[count]] <- unique(unlist(get.all.operons.list_Lrp.tg(neg.genelist[[i]][[j]],Lrp = names(neg.genelist)[i],
										Lrp.ChIP.targets,corr.coeff.matrix,cor.cutoff = 0.2,DirectionOfCorrelation="neg")))	
				names(each.lrp)[count] <- names(neg.genelist[[i]])[j]
				count <- count +1;
			}
		}
		neg.genelist1[[i]] <- each.lrp	
	} else 
	{
		neg.genelist1[[i]] <- neg.genelist[[i]]
	}
}
names(neg.genelist1) <- names(neg.genelist)




pos.genelist <- pos.genelist1 
neg.genelist <- neg.genelist1

#################
# positive
#################

for (i in 1: length(pos.genelist)) {
	if (length(pos.genelist[[i]])!=0) {
		geneSets <- pos.genelist[[i]]
		Lrp = names(pos.genelist)[i]
		print (Lrp)
		
		################
		# Run analysis #
		################
		# # Make Gene list
		# geneNames <- rownames(ratios.norm) # Background geneset. will be different for positive and negative correlation
		#----------------------------------------------------------------------------------------------------------------
		# m2.BP is the output:column 1:the top 10 GOs, column 2: ids of GOs that passed BH.adjusted p value<=0.05;
		# column 3: same as column2, use the names og GOs instead of ids; column 4: names of of GOs that passed BH.adjusted p value<=0.1;
		# column 5:names of GOs that passed (non-corrected) p values <=0.01;
		m2.BP <- matrix(nrow = length(geneSets), ncol= 13, dimnames = list(1:length(geneSets), c("#of.conditions","Top10.Terms.BP", "BH.0.05.GO.Ids.BP","BH.0.05.GO.names.BP","BH.0.1.GO.Ids.BP","BH.0.1.GO.names.BP",
					"p-val.0.01.GO.Ids.BP","p-val.0.01.GO.names.BP","p-val.0.001.GO.Ids.BP","p-val.0.001.GO.names.BP","p-val.0.0001.GO.Ids.BP","p-val.0.0001.GO.names.BP",
					"#of.genes.in.this.list")))
		#---------------------------------------------------------------------------------------------------------------
		for( geneSet in (1:length(geneSets)) ) {
			print(paste(geneSet," in ",length(geneSets),sep=""))
		#----------------------------------------------------------------
		#    add this part for setting gene background
		#----------------------------------------------------------------
			# set geneNames ( background genelist. if using absolute corr.is all 2400 genes. but changed when switch to pos/neg cor)
			set.cor.coeff = sapply(rownames(ratios.norm), function(gene){cor(ratios.norm[Lrp,term[[names(geneSets)[geneSet]]]], ratios.norm[gene,term[[names(geneSets)[geneSet]]]])})
			geneNames =  names(set.cor.coeff[set.cor.coeff> 0])

			# initiate Gene List
			tmp1 <- c(geneNames[1],geneNames[2]) # put in a couple of gene ids to build structure for automation
			geneList <- factor(as.integer(geneNames %in% tmp1))
			names(geneList) <- geneNames	
			# Make Biological Process GOData object
			GOdata.BP <- new("topGOdata", ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = vng2GO)
			#m1.BP <- matrix(nrow = length(GOdata.BP@graph@nodes), ncol = length(geneSets), dimnames = list(GOdata.BP@graph@nodes, 1:length(geneSets)))
		#-----------------------------------------------------------------	
			# Expand gene list and change factor in GOdata
			genes <- geneSets[[geneSet]]
			GOdata.BP@allScores <- factor(as.integer(geneNames %in% genes))
			#GOdata.MF@allScores <- factor(as.integer(geneNames %in% genes))
			#GOdata.CC@allScores <- factor(as.integer(geneNames %in% genes))
			# Biological process
			r1.BP = runTest(GOdata.BP, algorithm = 'classic', statistic = 'fisher')
			# m1.BP[,geneSet] = r1.BP@score
			m2.BP[geneSet,1] = length(term[[names(geneSets)[geneSet]]]) # number of conditions in each EO
			m2.BP[geneSet,2] = paste(GenTable(GOdata.BP, r1.BP)[,2],collapse=', ')  # top 10 GO terms
			m2.BP[geneSet,3] = paste(names(which(p.adjust(r1.BP@score,method='BH')<=0.05)),collapse=',') # GO.id passed BH.p-val<=0.05
			# create a conversion table for GO terms that have BH.p.value<=0.05
			temp.dataframe = GenTable(GOdata.BP, r1.BP, topNodes = length(GOdata.BP@graph@nodes))
			rownames(temp.dataframe)=temp.dataframe[,1]
			temp.dataframe[,6] <- as.numeric(temp.dataframe[,6])
			# temp.dataframe1 = temp.dataframe[names(which(p.adjust(r1.BP@score,method='BH')<=0.05)),]
			m2.BP[geneSet,4] = paste(temp.dataframe[names(which(p.adjust(r1.BP@score,method='BH')<=0.05)),][,2],collapse=',') # GO.name passed BH.p-val<=0.05
			m2.BP[geneSet,5] = paste(names(which(p.adjust(r1.BP@score,method='BH')<=0.1)),collapse=',') # GO.id passed BH.p-val<=0.1
			m2.BP[geneSet,6] = paste(temp.dataframe[names(which(p.adjust(r1.BP@score,method='BH')<=0.1)),][,2],collapse=',') # GO.name passed BH.p-val<=0.1
			m2.BP[geneSet,7] = paste(names(which(r1.BP@score<=0.01)),collapse=',')
			m2.BP[geneSet,8] = paste(temp.dataframe[names(which(r1.BP@score<=0.01)),][,2],collapse=',')
			m2.BP[geneSet,9] = paste(names(which(r1.BP@score<=0.001)),collapse=',')
			m2.BP[geneSet,10] = paste(temp.dataframe[names(which(r1.BP@score<=0.001)),][,2],collapse=',')
			m2.BP[geneSet,11] = paste(names(which(r1.BP@score<=0.0001)),collapse=',')
			m2.BP[geneSet,12] = paste(temp.dataframe[names(which(r1.BP@score<=0.0001)),][,2],collapse=',')
			m2.BP[geneSet,13] = length(genes)
			
			# Plot the structure of GO enrichment if there is any GO term that is enriched with p-value < 0.01
			if (sum(temp.dataframe[,6]<0.001)) {
				dir.create(paste(Lrp,"_pos.",target.type,"_pv_0.001_",substr(Sys.time(),1,10),sep=""))
				setwd(paste(set_dir,"/",paste(Lrp,"_pos.",target.type,"_pv_0.001_",substr(Sys.time(),1,10),sep=""),sep=""))
				printGraph(GOdata.BP, r1.BP, firstSigNodes = 5, useInfo = 'all', fn.prefix = paste(Lrp,"_pos.",target.type,"_pv_0.001_",names(geneSets)[geneSet],sep=""),pdfSW=TRUE)
				setwd(set_dir)
			}
			if (length(which(p.adjust(r1.BP@score,method='BH')<=0.1))> 0) {
				dir.create(paste(Lrp,"_pos.",target.type,"_BH_0.1_",substr(Sys.time(),1,10),sep=""))
				setwd(paste(set_dir,"/",paste(Lrp,"_pos.",target.type,"_BH_0.1_",substr(Sys.time(),1,10),sep=""),sep=""))
				printGraph(GOdata.BP, r1.BP, firstSigNodes = 5, useInfo = 'all', fn.prefix = paste(Lrp,"_pos.",target.type,"_BH_0.1_",names(geneSets)[geneSet],sep=""),pdfSW=TRUE)
				setwd(set_dir)
			}
		}
		rownames(m2.BP) = names(geneSets)
		# Write out the results, can add the others by doing cbind to add the other matrices as columns
		write.csv(m2.BP, file = paste(names(pos.genelist)[i],"_pos.",target.type,"_",substr(Sys.time(),1,10),".csv",sep="")) 
	}
}

################
# Negative
################
for (i in 1: length(neg.genelist)) {
	if (length(neg.genelist[[i]])!=0) {
		geneSets <- neg.genelist[[i]]
		Lrp = names(neg.genelist)[i]
		print (Lrp)

		################
		# Run analysis #
		################
		# # Make Gene list
		#----------------------------------------------------------------------------------------------------------------
		# m2.BP is the output:column 1:the top 10 GOs, column 2: ids of GOs that passed BH.adjusted p value<=0.05;
		# column 3: same as column2, use the names og GOs instead of ids; column 4: names of of GOs that passed BH.adjusted p value<=0.1;
		# column 5:names of GOs that passed (non-corrected) p values <=0.01;
		m2.BP <- matrix(nrow = length(geneSets), ncol= 13, dimnames = list(1:length(geneSets), c("#of.conditions","Top10.Terms.BP", "BH.0.05.GO.Ids.BP","BH.0.05.GO.names.BP","BH.0.1.GO.Ids.BP","BH.0.1.GO.names.BP",
					"p-val.0.01.GO.Ids.BP","p-val.0.01.GO.names.BP","p-val.0.001.GO.Ids.BP","p-val.0.001.GO.names.BP","p-val.0.0001.GO.Ids.BP","p-val.0.0001.GO.names.BP",
					"#of.genes.in.this.list")))
		#---------------------------------------------------------------------------------------------------------------

		for( geneSet in (1:length(geneSets)) ) {
			print(paste(geneSet," in ",length(geneSets),sep=""))
		#----------------------------------------------------------------
		#    add this part for setting gene background
		#----------------------------------------------------------------
			# set geneNames ( background genelist. if using absolute corr.is all 2400 genes. but changed when switch to pos/neg cor)
			set.cor.coeff = sapply(rownames(ratios.norm), function(gene){cor(ratios.norm[Lrp,term[[names(geneSets)[geneSet]]]], ratios.norm[gene,term[[names(geneSets)[geneSet]]]])})
			geneNames =  names(set.cor.coeff[set.cor.coeff< 0])
			# initiate Gene List
			tmp1 <- c(geneNames[1],geneNames[2]) # put in a couple of gene ids to build structure for automation
			geneList <- factor(as.integer(geneNames %in% tmp1))
			names(geneList) <- geneNames	
			# Make Biological Process GOData object
			GOdata.BP <- new("topGOdata", ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = vng2GO)
			#m1.BP <- matrix(nrow = length(GOdata.BP@graph@nodes), ncol = length(geneSets), dimnames = list(GOdata.BP@graph@nodes, 1:length(geneSets)))
		#-----------------------------------------------------------------	
			# Expand gene list and change factor in GOdata
			genes <- geneSets[[geneSet]]
			GOdata.BP@allScores <- factor(as.integer(geneNames %in% genes))
			#GOdata.MF@allScores <- factor(as.integer(geneNames %in% genes))
			#GOdata.CC@allScores <- factor(as.integer(geneNames %in% genes))
			# Biological process
			r1.BP = runTest(GOdata.BP, algorithm = 'classic', statistic = 'fisher')
			# m1.BP[,geneSet] = r1.BP@score
			m2.BP[geneSet,1] = length(term[[names(geneSets)[geneSet]]]) # number of conditions in each EO
			m2.BP[geneSet,2] = paste(GenTable(GOdata.BP, r1.BP)[,2],collapse=', ')  # top 10 GO terms
			m2.BP[geneSet,3] = paste(names(which(p.adjust(r1.BP@score,method='BH')<=0.05)),collapse=',') # GO.id passed BH.p-val<=0.05
			# create a conversion table for GO terms that have BH.p.value<=0.05
			temp.dataframe = GenTable(GOdata.BP, r1.BP, topNodes = length(GOdata.BP@graph@nodes))
			rownames(temp.dataframe)=temp.dataframe[,1]
			temp.dataframe[,6] <- as.numeric(temp.dataframe[,6])
			# temp.dataframe1 = temp.dataframe[names(which(p.adjust(r1.BP@score,method='BH')<=0.05)),]
			m2.BP[geneSet,4] = paste(temp.dataframe[names(which(p.adjust(r1.BP@score,method='BH')<=0.05)),][,2],collapse=',') # GO.name passed BH.p-val<=0.05
			m2.BP[geneSet,5] = paste(names(which(p.adjust(r1.BP@score,method='BH')<=0.1)),collapse=',') # GO.id passed BH.p-val<=0.1
			m2.BP[geneSet,6] = paste(temp.dataframe[names(which(p.adjust(r1.BP@score,method='BH')<=0.1)),][,2],collapse=',') # GO.name passed BH.p-val<=0.1
			m2.BP[geneSet,7] = paste(names(which(r1.BP@score<=0.01)),collapse=',')
			m2.BP[geneSet,8] = paste(temp.dataframe[names(which(r1.BP@score<=0.01)),][,2],collapse=',')
			m2.BP[geneSet,9] = paste(names(which(r1.BP@score<=0.001)),collapse=',')
			m2.BP[geneSet,10] = paste(temp.dataframe[names(which(r1.BP@score<=0.001)),][,2],collapse=',')
			m2.BP[geneSet,11] = paste(names(which(r1.BP@score<=0.0001)),collapse=',')
			m2.BP[geneSet,12] = paste(temp.dataframe[names(which(r1.BP@score<=0.0001)),][,2],collapse=',')
			m2.BP[geneSet,13] = length(genes)	
			# Plot the structure of GO enrichment if there is any GO term that is enriched with p-value < 0.01
			if (sum(temp.dataframe[,6]<0.001)) {
				dir.create(paste(Lrp,"_neg.",target.type,"_pv_0.001_",substr(Sys.time(),1,10),sep=""))
				setwd(paste(set_dir,"/",paste(Lrp,"_neg.",target.type,"_pv_0.001_",substr(Sys.time(),1,10),sep=""),sep=""))
				printGraph(GOdata.BP, r1.BP, firstSigNodes = 5, useInfo = 'all', fn.prefix = paste(Lrp,"_neg.",target.type,"_pv_0.001_",names(geneSets)[geneSet],sep=""),pdfSW=TRUE)
				setwd(set_dir)
			}
			if (length(which(p.adjust(r1.BP@score,method='BH')<=0.1))> 0) {
				dir.create(paste(Lrp,"_neg.",target.type,"_BH_0.1_",substr(Sys.time(),1,10),sep=""))
				setwd(paste(set_dir,"/",paste(Lrp,"_neg.",target.type,"_BH_0.1_",substr(Sys.time(),1,10),sep=""),sep=""))
				printGraph(GOdata.BP, r1.BP, firstSigNodes = 5, useInfo = 'all', fn.prefix = paste(Lrp,"_neg.",target.type,"_BH_0.1_",names(geneSets)[geneSet],sep=""),pdfSW=TRUE)
				setwd(set_dir)
			}
		}
		rownames(m2.BP) = names(geneSets)
		# Write out the results, can add the others by doing cbind to add the other matrices as columns
		write.csv(m2.BP, file = paste(names(neg.genelist)[i],"_neg.",target.type,"_",substr(Sys.time(),1,10),".csv",sep="")) 
	}
}

########################################### Functions #############################################
##                                                                                               ##
##                            Get all genes  in the same operon                                  ##
##                                                                                               ##
###################################################################################################

corr.coeff.matrix = sapply(Lrpnames, function(lrp)
						{
							sapply(rownames(ratios.norm), function(gene)
							{
								cor(ratios.norm[lrp,],ratios.norm[gene,] )
							})
						})

`get.all.operons.list_Lrp.tg` <- function(genelist,Lrp,Lrp.ChIP.targets,corr.coeff.matrix,cor.cutoff = 0.2,DirectionOfCorrelation="pos")
{
# Include genes that are (1) on the same operon & downstream of the input gene(s); 
#                        (2) on the same operon & upstream of the input gene(s) & ChIP target;
#                        (3) on the same operon & upstream of the input gene(s) & correlated with lrp ( same direction and pass a cutoff); 
# genelist: a character vector of genes canonical name
# operons: operon list
# DirectionOfCorrelation = "pos" or "neg"
    
	operons <- get.operons.halo()
	coords <- get.gene.coords.halo()
	gene.matrix <- coords[coords$canonical_Name %in% genelist,]
	is.for <- gene.matrix$Orientation == "For"
	targets <- Lrp.ChIP.targets[[Lrp]]
	
	
	operon.list=list()
	for (i in 1:dim(gene.matrix)[1])
	{
		if (!(gene.matrix[i,1] %in% unlist(operons))) operon.list[[i]]<-as.character(gene.matrix[i,1])
		else{	
				for (j in 1:length(operons))
				{
					if (gene.matrix[i,1] %in% operons[[j]])
					{
						indx = which (operons[[j]] %in% gene.matrix[i,1])
						# operon.list[[i]] = operons[[j]]
						if (is.for[i])
						{
							downsreamhits <- operons[[j]][indx:length(operons[[j]])]
							upstreamhits <- character()
							
							if (indx > 1)
							{
								count <- 1
								for (k in 1: (indx-1))
								{
									# if it's a target:included. correlation coeff.
									if (operons[[j]][k] %in% rownames(corr.coeff.matrix))
									{
										if ((operons[[j]][k] %in% targets)|
										((corr.coeff.matrix[operons[[j]][k],Lrp] >= cor.cutoff)&(DirectionOfCorrelation =="pos"))|
										((corr.coeff.matrix[operons[[j]][k],Lrp] <= -cor.cutoff)&(DirectionOfCorrelation =="neg")))
										{
											upstreamhits[count] <- operons[[j]][k];
											count <- count+1;
										}
									} else if (!operons[[j]][k] %in% rownames(corr.coeff.matrix)) 
									{
										if (operons[[j]][k] %in% targets)
										{
											upstreamhits[count] <- operons[[j]][k];
											count <- count+1;
										}
									} 
								}
							}
							operon.list[[i]] <- c(upstreamhits,downsreamhits)
						  
						} else 
						{
							downsreamhits <- operons[[j]][1:indx]
							upstreamhits <- character()
							if (indx <length(operons[[j]]) )
							{
								count <- 1
								for (k in (indx+1): length(operons[[j]]))
								{
									if (operons[[j]][k] %in% rownames(corr.coeff.matrix))
									{
										if ((operons[[j]][k] %in% targets)|
										((corr.coeff.matrix[operons[[j]][k],Lrp] >= cor.cutoff)&(DirectionOfCorrelation =="pos"))|
										((corr.coeff.matrix[operons[[j]][k],Lrp] <= -cor.cutoff)&(DirectionOfCorrelation =="neg")))
										{
											upstreamhits[count] <- operons[[j]][k];
											count <- count+1;
										}
									} else if (!operons[[j]][k] %in% rownames(corr.coeff.matrix)) 
									{
										if (operons[[j]][k] %in% targets)
										{
											upstreamhits[count] <- operons[[j]][k];
											count <- count+1;
										}
									}
								}
							}
							operon.list[[i]] <- c(upstreamhits,downsreamhits)
						} 
					} 
				}
			}
	}
	return(operon.list)
}
