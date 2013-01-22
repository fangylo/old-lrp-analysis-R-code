# load("2012-11-09-run.corr.pos.mcl.mean_median_x1000_w.modified.resampling.methods.RData")
# load("/users/slo/R/work/2012-11-09-Re-calculation of Lrp targets/2012-11-09-DataforCorrelation.RData")
# load("2012-11-10-rerun.smaller.pos.p.valsx100000_slice to every 1000.RData")

load("/users/slo/R/work/2012-11-09-Re-calculation of Lrp targets/2012-11-09-DataforCorrelation.RData")
# load("2012-11-27-FixedTerms_with_controls_included.RData")
load("2012-12-01-FixedTerms_with_controls_included.RData")
# ##################################################################################################################
# #        Does not need this part for ChIP targets that were defined by up 250bp & down 50bp(default, 
# #        saved in the above RData file) of a gene start site
# # to overide the original Lrp.ChIP.targets based on up250bp,incoding50bp:
# load("2012-11-12-Lrp.ChIP.targets_up_225.bp_incoding_25.bp.RData")
# ##################################################################################################################

require(multicore)
# Lrpnames <- c("VNG1377G","VNG1179C", "VNG1237C","VNG1285G", "VNG1816G","VNG2094G","VNG1351G", "VNG1123G")

term <- FixedTerm

sumConditions <- sapply(names(term), function(x)length(term[[x]]))
names.select.cond <- names(sumConditions[sumConditions >=4])# select a subset of EO that pass the above limitations. select.term is a list
select.term <- term [names.select.cond] 


all.genes.cor.coeff_ls = list()
targets.cor.coeff_ls = list()
neg.med.cor.coeff_ls = list()
pos.med.cor.coeff_ls = list()

for (lrp.x in Lrpnames){
	cat(paste(lrp.x,"\n",sep=""))
	
	Lrp = lrp.x
	
	targets = Lrp.ChIP.targets[[Lrp]]
	targets = targets[targets %in% rownames(ratios.norm)]
	l1 = length(targets)
	
	all.genes.cor.coeff_ls[[lrp.x]] <- mclapply(names(select.term), function(eo.term) 
	{
		cat(paste(eo.term,'\n',sep=''));return(
		sapply(rownames(ratios.norm), function(gene) { 
		cor(ratios.norm[Lrp,select.term[[eo.term]]], ratios.norm[gene,select.term[[eo.term]]]) }) )   
	})
	all.genes.cor.coeff_ls[[lrp.x]] <- do.call(cbind, all.genes.cor.coeff_ls[[lrp.x]])
	colnames(all.genes.cor.coeff_ls[[lrp.x]]) <- names(select.term)
	
	targets.cor.coeff_ls[[lrp.x]] <- all.genes.cor.coeff_ls[[lrp.x]][targets,]

	neg.med.cor.coeff_ls[[lrp.x]] <- sapply(names(select.term), function(eo.term){median(targets.cor.coeff_ls[[lrp.x]][,eo.term][targets.cor.coeff_ls[[lrp.x]][,eo.term]<0])})
	pos.med.cor.coeff_ls[[lrp.x]] <- sapply(names(select.term), function(eo.term){median(targets.cor.coeff_ls[[lrp.x]][,eo.term][targets.cor.coeff_ls[[lrp.x]][,eo.term]>0])})
	# abs.med.cor.coeff <- sapply(names(select.term), function(eo.term){median(abs(targets.cor.coeff[,eo.term]))})
}
save(all.genes.cor.coeff_ls,targets.cor.coeff_ls,neg.med.cor.coeff_ls,pos.med.cor.coeff_ls, file = "2012-12-01-all_median_correlation.coeff_128EO_up_250.bp_incoding_50.bp.RData")