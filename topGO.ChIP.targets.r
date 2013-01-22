library(topGO)

# Read in GO mappings to gene ids
# <gene>\t<GOID>,<GOID>,...
vng2GO <- readMappings(file="U:/R/work/halo.go.list.txt") #load halo GO file

# Load up cluster or gene set identifiers
geneSets <- Lrp.ChIP.targets

################
# Run analysis #
################
# Make Gene list
geneNames <- rownames(ratios)
tmp1 <- geneSets[[1]][c(1:3)] # put in a couple of gene ids to build structure for automation
geneList <- factor(as.integer(geneNames %in% tmp1))
names(geneList) <- geneNames
# Make Biological Process GOData object
GOdata.BP <- new("topGOdata", ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = vng2GO)
m1.BP <- matrix(nrow = length(GOdata.BP@graph@nodes), ncol = length(geneSets), dimnames = list(GOdata.BP@graph@nodes, 1:length(geneSets)))
m2.BP <- matrix(nrow = length(geneSets), ncol= 11, dimnames = list(1:length(geneSets), c("Top10.Terms.BP", "BH.0.05.GO.Ids.BP","BH.0.05.GO.names.BP","BH.0.1.GO.Ids.BP","BH.0.1.GO.names.BP",
				"p-val.0.01.GO.Ids.BP","p-val.0.01.GO.names.BP","p-val.0.001.GO.Ids.BP","p-val.0.001.GO.names.BP","p-val.0.0001.GO.Ids.BP","p-val.0.0001.GO.names.BP")))

for( geneSet in (1:length(geneSets)) ) {
    # Expand gene list and change factor in GOdata
    genes <- geneSets[[geneSet]]
    GOdata.BP@allScores <- factor(as.integer(geneNames %in% genes))
    #GOdata.MF@allScores <- factor(as.integer(geneNames %in% genes))
    #GOdata.CC@allScores <- factor(as.integer(geneNames %in% genes))
    # Biological process
	r1.BP = runTest(GOdata.BP, algorithm = 'elim', statistic = 'fisher')
	m1.BP[,geneSet] = r1.BP@score
	m2.BP[geneSet,1] = paste(GenTable(GOdata.BP, r1.BP)[,2],collapse=', ')  # top 10 GO terms
	m2.BP[geneSet,2] = paste(names(which(p.adjust(r1.BP@score,method='BH')<=0.05)),collapse=',') # GO.id passed BH.p-val<=0.05
	# create a conversion table for GO terms that have BH.p.value<=0.05
	temp.dataframe = GenTable(GOdata.BP, r1.BP, topNodes = length(GOdata.BP@graph@nodes))
	rownames(temp.dataframe)=temp.dataframe[,1]
	temp.dataframe[,6] <- as.numeric(temp.dataframe[,6])
	# temp.dataframe1 = temp.dataframe[names(which(p.adjust(r1.BP@score,method='BH')<=0.05)),]
	m2.BP[geneSet,3] = paste(temp.dataframe[names(which(p.adjust(r1.BP@score,method='BH')<=0.05)),][,2],collapse=',') # GO.name passed BH.p-val<=0.05
	m2.BP[geneSet,4] = paste(names(which(p.adjust(r1.BP@score,method='BH')<=0.1)),collapse=',') # GO.id passed BH.p-val<=0.1
	m2.BP[geneSet,5] = paste(temp.dataframe[names(which(p.adjust(r1.BP@score,method='BH')<=0.1)),][,2],collapse=',') # GO.name passed BH.p-val<=0.1
	m2.BP[geneSet,6] = paste(names(which(r1.BP@score<=0.01)),collapse=',')
	m2.BP[geneSet,7] = paste(temp.dataframe[names(which(r1.BP@score<=0.01)),][,2],collapse=',')
	m2.BP[geneSet,8] = paste(names(which(r1.BP@score<=0.001)),collapse=',')
	m2.BP[geneSet,9] = paste(temp.dataframe[names(which(r1.BP@score<=0.001)),][,2],collapse=',')
	m2.BP[geneSet,10] = paste(names(which(r1.BP@score<=0.0001)),collapse=',')
	m2.BP[geneSet,11] = paste(temp.dataframe[names(which(r1.BP@score<=0.0001)),][,2],collapse=',')
}

# Write out the results, can add the others by doing cbind to add the other matrices as columns
write.csv(m2.BP,'2012-11-13-ChIP.targets_GOBP_up.250bp_in.50bp_elim.csv') 