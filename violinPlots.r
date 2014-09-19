load("2013-03-15-Result.for.Correlation_35.environmental.EOs_Lrp.fold.change.1.75_p.val.0.05_median.r.NA_35EOs_up.250.bp_in.coding.50.bp.RData")

library(vioplot)
source("my.vioplot.r")

pos.select.pv <- sapply(pos.correlated.ChIP.targets, names)
neg.select.pv <- sapply(neg.correlated.ChIP.targets, names)

#########################
direction = "neg"
#########################
#                       #
#  Positive / negative  # 
#                       #
#########################
if (direction == "pos"){select.pv <- pos.select.pv} else {select.pv <- neg.select.pv}

for (i in 1: length(select.pv))
{
	if(!is.null(select.pv[[i]]))
	{
		lrpx <- names(select.pv)[i]
		targets = Lrp.ChIP.targets[[lrpx]]
		targets = targets[targets %in% rownames(ratios.norm)]
		select.eo.set <- select.pv[[lrpx]]
		
		all.corr.coeff = sapply(rownames(ratios.norm), function(gene){cor( ratios.norm[lrpx,],ratios.norm[gene,] )})

		hyper.pv.coex.targets = array()
		hyper.pv.coex.allgenes = array()
		tg.qmnk_check = array(dim = c(length(select.eo.set),5))
		rownames(tg.qmnk_check)= select.eo.set
		colnames(tg.qmnk_check)= c("q","m","n","k","p.value")
		cor.qmnk_check = array(dim = c(length(select.eo.set),5))
		rownames(cor.qmnk_check)= select.eo.set
		colnames(cor.qmnk_check)= c("q","m","n","k","p.value")
		
		for (j in 1: length(select.pv[[i]]))
		{
			select.eo <- select.pv[[i]][j]
			set.cor.coeff = sapply(rownames(ratios.norm), function(gene){cor(ratios.norm[lrpx,term[[select.eo]]], ratios.norm[gene,term[[select.eo]]])})
			targets.cor.coeff = set.cor.coeff[targets]

			
			
			if (direction=="pos")
			{
				##############
				# positive
				##############
				stdv <- sd(targets.cor.coeff[targets.cor.coeff>0])
				median.r.targets <- median(targets.cor.coeff[targets.cor.coeff>0])

				targets_coexpressed = names(targets.cor.coeff[targets.cor.coeff>(median.r.targets-stdv)])
				allgenes_coexpressed = names(set.cor.coeff[set.cor.coeff>(median.r.targets-stdv)])

				allgenes_background = names(set.cor.coeff[set.cor.coeff> 0])
				diff.KO = intersect(allgenes_background, lrp.diff.in.KO[[lrpx]])

			} else if (direction=="neg"){
				##############
				# negative
				##############
				stdv <- sd(targets.cor.coeff[targets.cor.coeff<0])
				median.r.targets <- median(targets.cor.coeff[targets.cor.coeff<0])

				targets_coexpressed = names(targets.cor.coeff[targets.cor.coeff<(median.r.targets+stdv)])
				allgenes_coexpressed = names(set.cor.coeff[set.cor.coeff<(median.r.targets+stdv)])

				allgenes_background = names(set.cor.coeff[set.cor.coeff< 0])
			}

			#########################
			#  plot violin plot     #
			#########################
			pdf(paste(lrpx,"_vioplot_",select.eo,"_",substr(Sys.time(),1,10),".pdf",sep=""))

			x1 <- all.corr.coeff
			x2 <- set.cor.coeff
			x3 <- set.cor.coeff[targets]
			x4 <- set.cor.coeff[allgenes_coexpressed]
			
			par(mar = c(5, 15, 4, 2))
			my.vioplot(x1,x2,x3,x4,col = c("gray","dodgerblue","red","blue"),horizontal=T, names= c("","","",""))
			labels = c("All Genes and Conditions",paste("All Genes ",select.eo,sep=""),paste(lrpx," Target Genes ",select.eo,sep=""), paste("Correlated ",lrpx," Target Genes ",select.eo,sep=""))
			text(y=c(1:4), par("usr")[1]-0.1, adj = 1, srt = 45, labels = labels, xpd = TRUE)
			title(paste(lrpx,"_",select.eo,sep=""))

			dev.off()	
		}
	}
}
