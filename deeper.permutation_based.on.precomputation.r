# load("/users/slo/R/work/2012-11-09-run.corr.pos.mcl.mean_median_x1000_w.modified.resampling.methods.RData")
# load("/users/slo/R/work/2012-11-09-Re-calculation of Lrp targets/2012-11-09-DataforCorrelation.RData")

load("/users/slo/R/work/2012-11-27-run.corr.pos.mcl.mean_median_x1000_w.modified.resampling.methods_up_250.bp_incoding_50.bp.RData")
load("/users/slo/R/work/2012-11-09-Re-calculation of Lrp targets/2012-11-09-DataforCorrelation.RData")
##################################################################################################################
#        Does not need this part for ChIP targets that were defined by up 250bp & down 50bp(default, 
#        saved in the above RData file) of a gene start site
# to overide the original Lrp.ChIP.targets based on up250bp,incoding50bp:
# load("2012-11-12-Lrp.ChIP.targets_up_225.bp_incoding_25.bp.RData")
##################################################################################################################
# load Fixed term (128EOs, with controls) 
load("2012-11-27-FixedTerms_with_controls_included.RData")
term <- FixedTerm

# term = term2
require(multicore)

# Lrpnames <- c("VNG1377G","VNG1179C", "VNG1237C","VNG1285G", "VNG1816G","VNG2094G","VNG1351G", "VNG1123G")

sumConditions <- sapply(names(term), function(x)length(term[[x]]))
names.select.cond <- names(sumConditions[sumConditions >=4])# select a subset of EO that pass the above limitations. select.term is a list
select.term <- term [names.select.cond] 

permuted.time <- 100000
permuted.slice <- 1000
cutoff.p.val = 0.0005

# abs.smaller.pv = sapply(names(abs.mean.perm.p.val), function(lrpx)which(abs.mean.perm.p.val[[lrpx]]<=0.1))
# pos.smaller.mean.pv = sapply(names(pos.mean.perm.p.val), function(lrpx)which(pos.mean.perm.p.val[[lrpx]]== 0))
pos.smaller.median.pv = sapply(names(pos.median.perm.p.val), function(lrpx)which(pos.median.perm.p.val[[lrpx]]== 0))
# neg.smaller.mean.pv = sapply(names(neg.mean.perm.p.val), function(lrpx)which(neg.mean.perm.p.val[[lrpx]] == 0))
# neg.smaller.median.pv = sapply(names(neg.mean.perm.p.val), function(lrpx)which(neg.median.perm.p.val[[lrpx]]==0))

# abs.perm.p.values1 = list()
# pos.mean.perm.p.val1 = list()
pos.median.perm.p.val1 =list()
# neg.mean.perm.p.val1 = list()
# neg.median.perm.p.val1 =list()
final.permuted.times_med.ls =list()
perm.pos.median_part_med.ls =list()
sum_permuted.over.target_med.ls =list()
collapsed.sum_med.ls =list()
# final.permuted.times_mean.ls =list()
# perm.pos.median_part_mean.ls =list()
# sum_permuted.over.target_mean.ls =list()
# collapsed.sum_mean.ls =list()

for (lrp.x in Lrpnames){
	cat(paste(lrp.x,"\n",sep=""))
	Lrp = lrp.x
	targets = Lrp.ChIP.targets[[Lrp]]
	targets = targets[targets %in% rownames(ratios.norm)]
	l1 = length(targets)
	
	all.genes.cor.coeff<- mclapply(names(select.term), function(eo.term) {cat(paste(eo.term,'\n',sep=''));return(
	sapply(rownames(ratios.norm), function(gene) { 
	cor(ratios.norm[Lrp,select.term[[eo.term]]], ratios.norm[gene,select.term[[eo.term]]]) }) )   })
	all.genes.cor.coeff <- do.call(cbind, all.genes.cor.coeff)
	colnames(all.genes.cor.coeff) <- names(select.term)
	
	targets.cor.coeff <- all.genes.cor.coeff[targets,]
	# neg.mean.cor.coeff <- sapply(names(select.term), function(eo.term){mean(targets.cor.coeff[,eo.term][targets.cor.coeff[,eo.term]<0])})
	pos.mean.cor.coeff <- sapply(names(select.term), function(eo.term){mean(targets.cor.coeff[,eo.term][targets.cor.coeff[,eo.term]>0])})
	# abs.mean.cor.coeff <- sapply(names(select.term), function(eo.term){mean(abs(targets.cor.coeff[,eo.term]))})

	# neg.med.cor.coeff <- sapply(names(select.term), function(eo.term){median(targets.cor.coeff[,eo.term][targets.cor.coeff[,eo.term]<0])})
	pos.med.cor.coeff <- sapply(names(select.term), function(eo.term){median(targets.cor.coeff[,eo.term][targets.cor.coeff[,eo.term]>0])})
	# abs.med.cor.coeff <- sapply(names(select.term), function(eo.term){median(abs(targets.cor.coeff[,eo.term]))})

	#----   positive correlation / slice the data and precompute first for faster performance ------------------------------------ #
	
	# median #
	if(length(pos.smaller.median.pv[[lrp.x]])>0) {
		smaller.eo.term <- names(pos.smaller.median.pv[[lrp.x]])
		reduced.smaller.eo.term = smaller.eo.term
		select.term1 <- term [reduced.smaller.eo.term] # only pick the EO that has p value =0
				
		perm.pos.median_part =list()
		sum_permuted.over.target =list()
		collapsed.sum = array(0, dim = c(length(smaller.eo.term),1))
		rownames(collapsed.sum) = smaller.eo.term
		final.permuted.times = array(0, dim = c(length(smaller.eo.term),1))
		rownames(final.permuted.times) = smaller.eo.term
		
		for (k in 1:(permuted.time/permuted.slice)) {
			
			perm.pos.median_part[[k]] <- sapply(1:permuted.slice, function (x){
			cat(paste(lrp.x,":pos_median:slice #",k,":",x,' ',sep='')); 
			unlist(mclapply(names(select.term1) ,function(eo.term){
			sampled.genes = sample(names(all.genes.cor.coeff[,eo.term][all.genes.cor.coeff[,eo.term]>0]),sum(targets.cor.coeff[,eo.term]>0)) ;
			permuted.cor = sapply(sampled.genes, function(gene) {cor(ratios.norm[Lrp,select.term[[eo.term]]], ratios.norm[gene,select.term[[eo.term]]])});
			# tmp.mean = mean(permuted.cor);
			tmp.med = median(permuted.cor);
			return(tmp.med)
			}))})
			if(!is.matrix(perm.pos.median_part[[k]]))   perm.pos.median_part[[k]] <- t(as.matrix(perm.pos.median_part[[k]]))
			rownames(perm.pos.median_part[[k]]) <- reduced.smaller.eo.term
			
			sum_permuted.over.target[[k]] = sapply(reduced.smaller.eo.term, function(eo.term)sum(perm.pos.median_part[[k]][eo.term,]>= pos.med.cor.coeff[eo.term]))
			
			for (eo in names(sum_permuted.over.target[[k]])) {
				collapsed.sum[eo,1] = collapsed.sum[eo,1] + sum_permuted.over.target[[k]][[eo]]
			}
			for (eo in names(sum_permuted.over.target[[k]])) {
				final.permuted.times[eo,1] = final.permuted.times[eo,1] + (permuted.slice)
			}
			reduced.smaller.eo.term = setdiff(rownames(collapsed.sum),  rownames(collapsed.sum)[which(collapsed.sum > (cutoff.p.val* final.permuted.times))])
			select.term1 <- term [reduced.smaller.eo.term]
		}
		pos.median.perm.p.val1[[lrp.x]] = collapsed.sum/final.permuted.times
		final.permuted.times_med.ls[[lrp.x]] =final.permuted.times
		perm.pos.median_part_med.ls[[lrp.x]] = perm.pos.median_part
		sum_permuted.over.target_med.ls[[lrp.x]] = sum_permuted.over.target
		collapsed.sum_med.ls[[lrp.x]] = collapsed.sum    
	}
	
	# # mean #
	# if(length(pos.smaller.mean.pv[[lrp.x]])>0) {
		# smaller.eo.term <- names(pos.smaller.mean.pv[[lrp.x]])
		# reduced.smaller.eo.term = smaller.eo.term
		# select.term1 <- term [reduced.smaller.eo.term] # only pick the EO that has p value =0
				
		# perm.pos.mean_part =list()
		# sum_permuted.over.target =list()
		# collapsed.sum = array(0, dim = c(length(smaller.eo.term),1))
		# rownames(collapsed.sum) = smaller.eo.term
		# final.permuted.times = array(0, dim = c(length(smaller.eo.term),1))
		# rownames(final.permuted.times) = smaller.eo.term
		
		# for (k in 1:(permuted.time/permuted.slice)) {
		
			# perm.pos.mean_part[[k]] <- sapply(1:permuted.slice, function (x){
			# cat(paste(lrp.x,":pos_mean:slice #",k,":",x,' ',sep='')); 
			# unlist(mclapply(names(select.term1) ,function(eo.term){
			# sampled.genes = sample(names(all.genes.cor.coeff[,eo.term][all.genes.cor.coeff[,eo.term]>0]),sum(targets.cor.coeff[,eo.term]>0)) ;
			# permuted.cor = sapply(sampled.genes, function(gene) {cor(ratios.norm[Lrp,select.term[[eo.term]]], ratios.norm[gene,select.term[[eo.term]]])});
			# # tmp.mean = mean(permuted.cor);
			# tmp.mean = mean(permuted.cor);
			# return(tmp.mean)
			# }))})
			# if(!is.matrix(perm.pos.mean_part[[k]]))   perm.pos.mean_part[[k]] <- t(as.matrix(perm.pos.mean_part[[k]]))
			# rownames(perm.pos.mean_part[[k]]) <- reduced.smaller.eo.term
			
			# sum_permuted.over.target[[k]] = sapply(reduced.smaller.eo.term, function(eo.term)sum(perm.pos.mean_part[[k]][eo.term,]>= pos.mean.cor.coeff[eo.term]))
			
			# for (eo in names(sum_permuted.over.target[[k]])) {
				# collapsed.sum[eo,1] = collapsed.sum[eo,1] + sum_permuted.over.target[[k]][[eo]]
			# }
			# for (eo in names(sum_permuted.over.target[[k]])) {
				# final.permuted.times[eo,1] = final.permuted.times[eo,1] + (permuted.slice)
			# }
			# reduced.smaller.eo.term = setdiff(rownames(collapsed.sum),  rownames(collapsed.sum)[which(collapsed.sum > (cutoff.p.val* final.permuted.times))])
			# select.term1 <- term [reduced.smaller.eo.term]
		# }
		# pos.mean.perm.p.val1[[lrp.x]] = collapsed.sum/final.permuted.times
		# final.permuted.times_mean.ls[[lrp.x]] =final.permuted.times
		# perm.pos.median_part_mean.ls[[lrp.x]] = perm.pos.median_part
		# sum_permuted.over.target_mean.ls[[lrp.x]] = sum_permuted.over.target
		# collapsed.sum_mean.ls[[lrp.x]] = collapsed.sum  
		    
	# }
}
save(pos.median.perm.p.val1,
	final.permuted.times_med.ls,perm.pos.median_part_med.ls,
	sum_permuted.over.target_med.ls,collapsed.sum_med.ls,
	file = "2012-11-28-rerun.smaller.pos.p.valsx100000_slice to every 1000_up_250.bp_incoding_50.bp.RData")



