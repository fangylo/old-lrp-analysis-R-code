load("U:/R/work/2012-11-29-DataforCorrelation.RData")
load("2012-12-17-Env.EOTermlist_Published.RData")
term <- FixedTerm
allexp = unique(unlist(term))

FixedEO_withGroupID = read.csv2("2012-11-29-fixedEOs_withGroupID.csv",check.names=FALSE, sep=",")
rownames(FixedEO_withGroupID) <- as.character(FixedEO_withGroupID[,1])
FixedEO_withGroupID <- FixedEO_withGroupID[,-1]
FixedEO_withGroupID <- FixedEO_withGroupID[allexp,]


GroupIDIndex = as.numeric(FixedEO_withGroupID[,1])
names(GroupIDIndex) = rownames(FixedEO_withGroupID)


####################################################################################
##  Goal: To calculate Fold Change of the expression of each Lrp under each EO    ##     
##  (1) Based on the 'experiments'(recorded as GroupID in EnvMap). In each        ##
##      EO, calculate the fold change from different experiments seperately.      ##
##  (2) For each EO, record the maximum FC from one of the experiment.            ##
##  (3) Check if the control experiment is included in the above selected         ##
##      experiment in this EO. If so, record the highest and lowest expression    ##
##      level relative to the expression of Lrp at the control condition          ##
####################################################################################

lrpx <- "VNG1377G"

###############################################################
#  For each Lrp, contruct a 128(EO) x 3 matrix:               #
#  1st column: the value of maximun FC of Lrp expression      # 
#  2nd column and 3rd column are meaningful only when a       # 
#  control condition was included in the experiment that max  #
#  FC was observed within this EO. 2nd column is the max      #
#  expression of Lrp relative to its expression in the        #
#  control condition, and the 3rd col is the min.             #
###############################################################

Lrp.FoldChange = list()
for (lrpx in Lrpnames)
{
	FoldChange = matrix(nrow = length(term), ncol = 3, dimnames = list(names(term),c("maxFC","maxExp_RelativeToCtrl","minExp_RelativeToCtrl")))

	for (EOterm in names(term))
	{
		### Separate the experiments based on GroupIDs ###
		Conditions_thisEO = term[[EOterm]]
		ExperimentsID = unique(GroupIDIndex[Conditions_thisEO])

		### Calculate fold change of the TF expresison (within each experiment) within this EO ###
		FC_thisEO <- sapply(ExperimentsID, function(ExpID){abs(max(ratios[lrpx,Conditions_thisEO[GroupIDIndex[Conditions_thisEO]==ExpID]])-min(ratios[lrpx,Conditions_thisEO[GroupIDIndex[Conditions_thisEO]==ExpID]]))})

		### Maximun fold change of this EO ###
		FoldChange[EOterm,"maxFC"] = 10^(max(FC_thisEO))
		
		### Check if there is control condition included in this experiment in this EO ###
		
		# GroupID of this experiment
		ExperimentsID_maxFC <- ExperimentsID[which(FC_thisEO == max(FC_thisEO))]

		# Select the conditions from the same experiment:
		Conditions_thisEO_thisExperiment <- Conditions_thisEO[GroupIDIndex[Conditions_thisEO] == ExperimentsID_maxFC]

		# If there is Control, fill in the 2nd and the 3rd columns of the matrix:
		if ( 2 %in% FixedEO_withGroupID[Conditions_thisEO_thisExperiment,EOterm])
		{
			CtrlCondition_thisEO_thisExperiment <- Conditions_thisEO_thisExperiment[which(FixedEO_withGroupID[Conditions_thisEO_thisExperiment,EOterm]==2)]
			CtrlExpressionLevel <- mean(ratios[lrpx, CtrlCondition_thisEO_thisExperiment]) # in case there are more than one controls
			FoldChange[EOterm,"maxExp_RelativeToCtrl"] <- 10^(max(ratios[lrpx,Conditions_thisEO_thisExperiment]) - CtrlExpressionLevel)
			FoldChange[EOterm,"minExp_RelativeToCtrl"] <- 10^(min(ratios[lrpx,Conditions_thisEO_thisExperiment]) - CtrlExpressionLevel)
		}
	}
	
	Lrp.FoldChange[[lrpx]] <- FoldChange
}

save(Lrp.FoldChange, file = "2012-12-18-FoldChange_relativetoCtrls.RData")




fold.change = 1.25
EO2test = sapply(names(Lrp.FoldChange), function(x){names(Lrp.FoldChange[[x]][,1])[Lrp.FoldChange[[x]][,1]>= fold.change ]})
