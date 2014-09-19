load("U:/R/work/2012-07-20-median correlation and p-values of 8Lrps in 148EOs.RData")
load("U:/R/work/2012-06-27-DataforCorrelation.RData")
load("U:/R/work/2012-07-17-diff.genes.BH.0.05.KO.RData")
load("U:/R/work/2012-06-30-GO2vng.RData") # load from GO2VNG gene names
load("U:/R/work/2012-07-19-CCgenelist_adjusted.p.value.cutoff.0.05_correlated.CC.targets_correlated.genes.overlap.w.KO.RData")

library(topGO)
vng2GO <- readMappings(file="U:/R/work/halo.go.list.txt") #load halo GO file
geneNames <- rownames(ratios.norm) # Background geneset. will be different for positive and negative correlation
tmp1 <- c('VNG1816G','VNG1814G') # put in a couple of gene ids to build structure for automation
geneList <- factor(as.integer(geneNames %in% tmp1))
names(geneList) <- geneNames
# Make Biological Process GOData object
GOdata.BP.all <- new("topGOdata", ontology='BP', allGenes = geneList, annot = annFUN.gene2GO, gene2GO = vng2GO)


Lrpnames <- c("VNG1377G","VNG1179C", "VNG1237C","VNG1285G", "VNG1816G","VNG2094G","VNG1351G", "VNG1123G")
##############################################################################################################################
# Analyze topGO data 
# Positive correlated genes overlap with KO
# topGO enrichment try using BH corrected p-value 0.1 as cut off
###################################################################################################################
set_dir = "U:/R/work/2012-07-17-topGO-functional enrichment.148EO/2012-07-17-corrected.permuted.BH.pv.0.05/2012-07-19_pos.correlated.genes.overlap.w.KO"
genels = pos.correlated.genes.overlapwith.KO

temp.list = list.files(set_dir)
temp.list = temp.list[grep(".csv", temp.list)]
temp.list = temp.list[grep("VNG",temp.list)]

setwd(set_dir)
dirpath = "for.cytoscape_topGO"
# dir.create(dirpath)
genes.in.GO = list()

for (j in 1: length(temp.list)) {
	# fl.path.0.001 = paste(set_dir,"/",dirpath,"/",substr(temp.list[j],1,16),".genes.KO.GO.0.001.csv",sep="")
	fl.path.BH.0.1 = paste(set_dir,"/",dirpath,"/",substr(temp.list[j],1,16),".genes.KO.GO.BH.0.1.csv",sep="")
	
	Lrp =  substr(temp.list[j],1,8)
	
	m1 = read.csv(temp.list[j])
	m1 = m1[!duplicated(m1[,1]),] # getting rid of duplicated rows
	rownames(m1) = m1[,1]
	m1 = m1[,-1]
	m2 = m1[,c("X.of.conditions" ,"BH.0.1.GO.Ids.BP","BH.0.1.GO.names.BP","p.val.0.001.GO.Ids.BP","p.val.0.001.GO.names.BP","X.of.genes.in.this.list")]
	
	# Add permuted p-values and BH corrected p-values to the output csv files:
	if (sum(grep('pos',temp.list[j]))>0) {
		permuted.p.val = all.median.r.p.values[[Lrp]][rownames(m2),c("BH.pos.p.p-values","pos.p.p-values")]
		} else if (sum(grep('neg',temp.list[j]))>0){
			permuted.p.val = all.median.r.p.values[[Lrp]][rownames(m2),c("BH.neg.p.p-values","neg.p.p-values")]}
	if (length(rownames(m2))==1) permuted.p.val = t(permuted.p.val)
	m2 = cbind(m2, permuted.p.val)
	#-------------------------------------------------------------------------------------------------
	# # if ('thr4.ko' %in% rownames(m2)){rownames(m2)[which(rownames(m2)=="thr4.ko")]<- "trh4.ko"}
	# cond.num = sapply(rownames(m2), function(x)length(term[[x]]))
	#----- out put files ---------------------------------------------------------------------------------------
	GoTermList =list()
	genes.in.GO.each.eo = list()
	for (eo in rownames(m2)) {
		if (sum(is.na(m2))>0) {
			if(!is.na(m2[eo,2])) {
				go.names = unlist(strsplit(as.character(m2[eo,"BH.0.1.GO.names.BP"]),","));
				go.ids = unlist(strsplit(as.character(m2[eo,"BH.0.1.GO.Ids.BP"]),","));
				m3 = array(dim = c(length(go.names),2));
				m3[c(1:length(go.names)),1] = paste(eo,"(",m2[eo,"X.of.conditions"],")(p-value ",0.00001*round(100000*m2[eo,7]),")",sep="");
				m3[,2] = go.names;
				GoTermList[[eo]] = m3
				# genes in GO terms:
				GO.genes =list()
				for (id in 1:length(go.ids)) {
					GO.genes[[id]] = intersect(genels[[Lrp]][[eo]], genesInTerm(GOdata.BP.all,go.ids[id])[[1]])
					names(GO.genes)[id] = go.ids[id]			
				} 
				genes.in.GO.each.eo[[eo]] = GO.genes
			}
		} else if(!as.character(m2[eo,2]) =="") {
			go.names = unlist(strsplit(as.character(m2[eo,"BH.0.1.GO.names.BP"]),","));
			go.ids = unlist(strsplit(as.character(m2[eo,"BH.0.1.GO.Ids.BP"]),","));
			m3 = array(dim = c(length(go.names),2));
			m3[c(1:length(go.names)),1] = paste(eo,"(",m2[eo,"X.of.conditions"],")(p-value ", 0.00001*round(100000*m2[eo,7]),")",sep="");
			# intersect(genels[[Lrp]][[eo]], unlist(genesInTerm(GOdata.BP.all,go.ids[1])))
			m3[,2] = go.names;
			GoTermList[[eo]] = m3 # for writing to csv files
			# genes in GO terms:
			GO.genes =list()
			for (id in 1:length(go.ids)) {
				GO.genes[[id]] = intersect(genels[[Lrp]][[eo]], genesInTerm(GOdata.BP.all,go.ids[id])[[1]])
				names(GO.genes)[id] = go.ids[id]			
			} 
			genes.in.GO.each.eo[[eo]] = GO.genes	
		}
		
	}
	genes.in.GO [[Lrp]] = genes.in.GO.each.eo
	#----
	GO1 = sapply(rownames(m2),function(x){paste(x,"(",m2[x,"X.of.conditions"],")(p-value ",0.00001*round(100000*m2[eo,7]),")",sep="")})
	
	out_file <- file(fl.path.BH.0.1, open="a") 
	if (length(GoTermList)>0) {
		for (i in 1: length(GoTermList)){
			if (length(GoTermList[[i]])>0) {
				write.table(GoTermList[[i]], file=out_file, sep=",",quote=FALSE, col.names=FALSE, row.names=FALSE) 
			}
		}
	}
	m4 = array(dim= c(length(rownames(m2)),2))
	m4[c(1:length(rownames(m2))),1] = Lrp
	m4[c(1:length(rownames(m2))),2] = GO1
	write.table(m4, file=out_file, sep=",",quote=FALSE, col.names=FALSE, row.names=FALSE)
	close(out_file)
}
save(genes.in.GO, file = paste(substr(Sys.time(),1,10),"_genes.in.enrichedGOTerms_","pos.cor.KO",".RData",sep = ""))
