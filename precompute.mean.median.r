# re-run at deeper permutations for positive and negative correlation
# store both mean and median
# try doing them at separate level so it runs faster
# reduce some environmental ontology

load("/users/slo/R/work/2012-11-09-Re-calculation of Lrp targets/2012-11-09-DataforCorrelation.RData")
##################################################################################################################
# load Fixed term (128EOs, with controls) 
# load("2012-11-27-FixedTerms_with_controls_included.RData") # Use the EO list whose genetic EOs also does not include the parental strain controls
load("2012-12-01-FixedTerms_with_controls_included.RData") # Use the EO list whose genetic EOs also include the parental strain controls
term <- FixedTerm
##################################################################################################################
# #        Does not need this part for ChIP targets that were defined by up 250bp & down 50bp(default, 
# #        saved in the above RData file) of a gene start site
# # to overide the original Lrp.ChIP.targets based on up250bp,incoding50bp:
# load("2012-11-13-Lrp.ChIP.targets_up_250.bp_incoding_25.bp.RData")
##################################################################################################################
require(multicore)
# Lrpnames <- c("VNG1377G","VNG1179C", "VNG1237C","VNG1285G", "VNG1816G","VNG2094G","VNG1351G", "VNG1123G")

sumConditions <- sapply(names(term), function(x)length(term[[x]]))
names.select.cond <- names(sumConditions[sumConditions >=4])# select a subset of EO that pass the above limitations. select.term is a list
select.term <- term [names.select.cond]

permuted.time = 1000

# abs.mean.perm.p.val = list()
# neg.mean.perm.p.val = list()
pos.mean.perm.p.val = list()
# abs.median.perm.p.val = list()
# neg.median.perm.p.val = list()
pos.median.perm.p.val = list()

for ( lrp.x in Lrpnames ){

	cat(paste(lrp.x,"\n",sep=""))
	Lrp = lrp.x
	targets = Lrp.ChIP.targets[[Lrp]]
	targets = targets[targets %in% rownames(ratios.norm)]

	# the correlation coefficients of all 2400 genes relative to Lrp in each selected environmental ontology conditions.
	# This way we can do re-sampling on only the genes that have either positive or negative correlations within a certain EO
	# produce a matrix: dim = c(2400, length(select.term))
	all.genes.cor.coeff<- mclapply(names(select.term), function(eo.term) {cat(paste(eo.term,'\n',sep=''));return(
	sapply(rownames(ratios.norm), function(gene) { 
	cor(ratios.norm[Lrp,select.term[[eo.term]]], ratios.norm[gene,select.term[[eo.term]]]) }) )   })
	all.genes.cor.coeff <- do.call(cbind, all.genes.cor.coeff)
	colnames(all.genes.cor.coeff) <- names(select.term)

	# the correlation coefficient of each target relative to Lrp in each selected environmental ontology conditions.
	# produce a matrix: dim = c(length(targets), length(select.term))
	targets.cor.coeff <- all.genes.cor.coeff[targets,]
	
	# neg.mean.cor.coeff <- sapply(names(select.term), function(eo.term){mean(targets.cor.coeff[,eo.term][targets.cor.coeff[,eo.term]<0])})
	pos.mean.cor.coeff <- sapply(names(select.term), function(eo.term){mean(targets.cor.coeff[,eo.term][targets.cor.coeff[,eo.term]>0])})
	# abs.mean.cor.coeff <- sapply(names(select.term), function(eo.term){mean(abs(targets.cor.coeff[,eo.term]))})

	# neg.med.cor.coeff <- sapply(names(select.term), function(eo.term){median(targets.cor.coeff[,eo.term][targets.cor.coeff[,eo.term]<0])})
	pos.med.cor.coeff <- sapply(names(select.term), function(eo.term){median(targets.cor.coeff[,eo.term][targets.cor.coeff[,eo.term]>0])})
	# abs.med.cor.coeff <- sapply(names(select.term), function(eo.term){median(abs(targets.cor.coeff[,eo.term]))})


	# Permuted vector
	l1 = length(targets[targets %in% rownames(ratios.norm)])

	# perm.abs.mean.med <-  sapply(1:permuted.time, function(x) { 
		# cat(paste(lrp.x,':abs:',x,' ',sep='')); 
		# unlist(mclapply(names(select.term), function(eo.term) {
		# sampled.genes = sample(rownames(ratios.norm),l1);
		# all.cor = abs(sapply(sampled.genes, function(gene) {cor(ratios.norm[Lrp,select.term[[eo.term]]], ratios.norm[gene,select.term[[eo.term]]])}));
		# tmp.mean = mean(all.cor);
		# tmp.med = median(all.cor);
		# return(list("mean" = tmp.mean,"median" = tmp.med))
	# }))	})
	# perm.abs.mean = perm.abs.mean.med [rownames(perm.abs.mean.med )=="mean",]
	# perm.abs.med = perm.abs.mean.med [rownames(perm.abs.mean.med )=="median",]

	# perm.neg.mean.med <- sapply(1:permuted.time, function (x){
		# cat(paste(lrp.x,':neg:',x,' ',sep='')); 
		# unlist(mclapply(names(select.term) ,function(eo.term){
		# sampled.genes = sample(names(all.genes.cor.coeff[,eo.term][all.genes.cor.coeff[,eo.term]<0]),sum(targets.cor.coeff[,eo.term]<0)) ;
		# all.cor = sapply(sampled.genes, function(gene) {cor(ratios.norm[Lrp,select.term[[eo.term]]], ratios.norm[gene,select.term[[eo.term]]])});
		# tmp.mean = mean(all.cor);
		# tmp.med = median(all.cor);
		# return(list("mean" = tmp.mean,"median" = tmp.med))
	# }))})
	# perm.neg.mean = perm.neg.mean.med [rownames(perm.neg.mean.med )=="mean",]
	# perm.neg.med = perm.neg.mean.med [rownames(perm.neg.mean.med )=="median",]

	perm.pos.mean.med <- sapply(1:permuted.time, function (x){
		cat(paste(lrp.x,':pos:',x,' ',sep='')); 
		unlist(mclapply(names(select.term) ,function(eo.term){
		sampled.genes = sample(names(all.genes.cor.coeff[,eo.term][all.genes.cor.coeff[,eo.term]>0]),sum(targets.cor.coeff[,eo.term]>0)) ;
		all.cor = sapply(sampled.genes, function(gene) {cor(ratios.norm[Lrp,select.term[[eo.term]]], ratios.norm[gene,select.term[[eo.term]]])});
		tmp.mean = mean(all.cor);
		tmp.med = median(all.cor);
		return(list("mean" = tmp.mean,"median" = tmp.med))
	}))})
	perm.pos.mean = perm.pos.mean.med [rownames(perm.pos.mean.med )=="mean",]
	perm.pos.med = perm.pos.mean.med [rownames(perm.pos.mean.med )=="median",]

	rownames(perm.pos.mean) <- rownames(perm.pos.med) <- names(select.term)

	# Calculate perumted p-values
	# abs.mean.perm.p.val[[lrp.x]] = unlist(sapply(names(select.term), function(eo.term) { sum(perm.abs.mean[eo.term,]>=abs.mean.cor.coeff[eo.term])/permuted.time } ))
	# neg.mean.perm.p.val[[lrp.x]] = unlist(sapply(names(select.term), function(eo.term) { sum(perm.neg.mean[eo.term,]<=neg.mean.cor.coeff[eo.term])/permuted.time } ))
	pos.mean.perm.p.val[[lrp.x]] = unlist(sapply(names(select.term), function(eo.term) { sum(perm.pos.mean[eo.term,]>=pos.mean.cor.coeff[eo.term])/permuted.time } ))
	
	# abs.median.perm.p.val[[lrp.x]] = unlist(sapply(names(select.term), function(eo.term) { sum(perm.abs.med[eo.term,]>=abs.med.cor.coeff[eo.term])/permuted.time } ))
	# neg.median.perm.p.val[[lrp.x]] = unlist(sapply(names(select.term), function(eo.term) { sum(perm.neg.med[eo.term,]<=neg.med.cor.coeff[eo.term])/permuted.time } ))
	pos.median.perm.p.val[[lrp.x]] = unlist(sapply(names(select.term), function(eo.term) { sum(perm.pos.med[eo.term,]>=pos.med.cor.coeff[eo.term])/permuted.time } ))

}
save(pos.mean.perm.p.val,pos.median.perm.p.val,file = "2012-12-01-run.corr.pos.mcl.mean_median_x1000_w.modified.resampling.methods_up_250.bp_incoding_50.bp.RData")










