load("U:/R/work/2012-11-29-DataforCorrelation.RData")
## Only consider published data:
load("2012-12-17-Env.EOTermlist_Published.RData")
term <- FixedTerm

##########################################################################
# (1) Construct Matrixes with correlation between the expression         #
#     of the Lrp and the variance among target genes.                    #
# (2) We want to see a negative correlation between the variances        #
#     among targets and the expression of the TF                         #
##########################################################################

Correlation_Variances_LrpExpression = list()

for (lrpx in Lrpnames)
{
	Lrp <- lrpx
	targets <- Lrp.ChIP.targets[[Lrp]] %in% rownames(ratios)
	
	CorVarMatrix <- matrix(nrow = length(term), ncol = 2, dimnames=list(names(term), c("cor","p.value")))
	for (i in 1:length(term))
	{
		EOTerm <- names(term)[i]
		variance = sapply(1: length(term[[EOTerm]]), function(x){var(ratios[targets, term[[EOTerm]][x]])})
		cor.result <- cor.test(variance, ratios[Lrp, term[[EOTerm]]])
		CorVarMatrix[i,1] = cor.result$estimate
		CorVarMatrix[i,2] = cor.result$p.value
	}

	# The EO terms that have negative correlation between the expression of Lrp and the variance among the target genes
	EO_negative <- rownames(CorVarMatrix)[CorVarMatrix[,"cor"]<=0]
	# These are the EOs that we are not considering. Set the p values to 1
	EO_positive <- rownames(CorVarMatrix)[CorVarMatrix[,"cor"]>0]

	adjusted.p <-  matrix(nrow = length(term), ncol = 1, dimnames = list(names(term), "BH.p.value"))
	adjusted.p[EO_negative,] <- p.adjust(CorVarMatrix[EO_negative,"p.value"], method = "BH")
	adjusted.p[EO_positive,] <- 1

	Correlation_Variances_LrpExpression[[lrpx]] <- cbind(CorVarMatrix,adjusted.p)
}


# sapply(Lrpnames, function(x){Correlation_Variances_LrpExpression[[x]][Correlation_Variances_LrpExpression[[x]][,3]<=0.05,]})
save(Correlation_Variances_LrpExpression, file ="2012-12-18-Correlation_Variances_LrpExpression.RData")

