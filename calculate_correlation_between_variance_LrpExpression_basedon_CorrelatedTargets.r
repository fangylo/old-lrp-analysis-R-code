load("U:/R/work/2012-11-29-DataforCorrelation.RData")
load("2012-12-17-Env.EOTermlist_Published.RData")
term <- FixedTerm

### Load data from permutation run  ###
load("2012-12-18-un_corrected_permuted.pvalue_100000permutations.RData")

### Load data of median correlation coeff between each gene and Lrp ###
load("2012-12-17-all_median_correlation.coeff_48EO_up_250.bp_incoding_50.bp.RData")

###########################################################################
# (1) Construct Matrixes with correlation between the expression          #
#     of the Lrp and the variance among target genes whose expression     # 
#     correlated with the Lrp in each EO.                                 #
# (2) Matrix: nrow = 128EOs; ncol = 3. First column is the correlation    #
#     coefficient between the expression of the Lrp and the variance      #
#     among the targets whose expression correlated with the Lrp (defined #
#     as the targets that passed one stdv away from the median); second   #
#     column is the (uncorrected) p value of correlation (from cor.test); # 
#     third column is the list of the above correlated targets.           #
#     Also, here only EOs whose uncorrelcted permuted p-values from the   #
#     correlation analysis less than 0.05 run are calculated.             #
# (3) We want to see a negative correlation between the variances         #
#     among targets and the expression of the TF                          #
###########################################################################

p.cutoff = 0.2
Pos.CorVarMatrix = list()
Neg.CorVarMatrix = list()

for (lrpx in Lrpnames)
	{
		targets = Lrp.ChIP.targets[[lrpx]]
		targets = targets[targets %in% rownames(ratios.norm)]
		
		# Calculate variance based on positive correlated targets
		pos.sig.EOs <- names(pos.median.perm.p.val[[lrpx]])[pos.median.perm.p.val[[lrpx]]<=p.cutoff];
		Pos.CorVarMatrix[[lrpx]] <- matrix(nrow = length(term), ncol = 3, dimnames=list(names(term), c("Pos.CorrCoeff_corTargets","p.value","correlatedTargets")))
		# Only look at the EOs in which the expression of Lrp (significantly? uncorrected p value<=0.05) correlated w/ the expression of targets
		for (i in 1: length(pos.sig.EOs))
		{
			EOTerm <- pos.sig.EOs[i];
			
			TargetsCorrCoeffs <- all.genes.cor.coeff_ls[[lrpx]][targets,EOTerm]
			median.TargetsCorrCoeffs <- median(TargetsCorrCoeffs[TargetsCorrCoeffs>0])
			stdv.TargetsCorrCoeffs <- sd(TargetsCorrCoeffs[TargetsCorrCoeffs>0])
			
			CorrelatedTargets <- names(TargetsCorrCoeffs)[TargetsCorrCoeffs>(median.TargetsCorrCoeffs-stdv.TargetsCorrCoeffs)]
			variance <- sapply(1: length(term[[EOTerm]]), function(x){var(ratios[CorrelatedTargets, term[[EOTerm]][x]])})
			cor.result <- cor.test(variance, ratios[lrpx, term[[EOTerm]]])
			Pos.CorVarMatrix[[lrpx]][EOTerm,1] = cor.result$estimate
			Pos.CorVarMatrix[[lrpx]][EOTerm,2] = cor.result$p.value
			Pos.CorVarMatrix[[lrpx]][EOTerm,3] = paste(CorrelatedTargets,collapse=',')
		}
		
		# Calculate variance based on negative correlated targets
		neg.sig.EOs <- names(neg.median.perm.p.val[[lrpx]])[neg.median.perm.p.val[[lrpx]]<=p.cutoff];
		Neg.CorVarMatrix[[lrpx]] <- matrix(nrow = length(term), ncol = 3, dimnames=list(names(term), c("Neg.CorrCoeff_corTargets","p.value","correlatedTargets")))
		# Only look at the EOs in which the expression of Lrp (significantly? uncorrected p value<=0.05) correlated w/ the expression of targets
		for (i in 1: length(neg.sig.EOs))
		{
			EOTerm <- neg.sig.EOs[i];
			
			TargetsCorrCoeffs <- all.genes.cor.coeff_ls[[lrpx]][targets,EOTerm]
			median.TargetsCorrCoeffs <- median(TargetsCorrCoeffs[TargetsCorrCoeffs<0])
			stdv.TargetsCorrCoeffs <- sd(TargetsCorrCoeffs[TargetsCorrCoeffs<0])
			
			CorrelatedTargets <- names(TargetsCorrCoeffs)[TargetsCorrCoeffs<(median.TargetsCorrCoeffs+stdv.TargetsCorrCoeffs)]
			variance <- sapply(1: length(term[[EOTerm]]), function(x){var(ratios[CorrelatedTargets, term[[EOTerm]][x]])})
			cor.result <- cor.test(variance, ratios[lrpx, term[[EOTerm]]])
			Neg.CorVarMatrix[[lrpx]][EOTerm,1] = cor.result$estimate
			Neg.CorVarMatrix[[lrpx]][EOTerm,2] = cor.result$p.value
			Neg.CorVarMatrix[[lrpx]][EOTerm,3] = paste(CorrelatedTargets,collapse=',')
		}
	}

save(Pos.CorVarMatrix, Neg.CorVarMatrix, file = "2012-12-18-Correlation_CoeffVar_LrpExpression_basedon_CorrelatedTargets.RData")