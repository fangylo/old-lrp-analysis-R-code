#############################################################################################################
##  Calculate the length of intergenic region based on each transcription unit                             ##
##  (1) operon.coords1 is the modified genome coordinates that take the operons into consideration         ## 
##  so that the gaps between genes in the same operon are ignored; (2) genome.length and operons.coords    ## 
##  can be calculated from  U:/R/work/2012-11-08-calculate.average.intergenic.based.on.genes.r;            ##
##  (3) chr is "HALCHR", "pNRC200" or "pNRC100" ; (4) output a vector with the length of each intergenic   ##
##  fragments;                                                                                             ##
#############################################################################################################

# chr = "HALCHR"
# coords = operon.coords1[[chr]]

`get.operon.coords`<- function() 
{
  if ( ! exists( "operon.coords1" ) ) try( load("2012-11-08-for.calculating.intergenic.length.RData")  )
  if ( exists( "operon.coords1" ) ) return( operon.coords1 )
}

`get.genome.length`<- function()
{
  if ( ! exists( "genome.length" ) ) try( load("2012-11-08-for.calculating.intergenic.length.RData")  )
  if ( exists( "genome.length" ) ) return( genome.length )
}


`get.all.intergenic.fragment.length` <- function(chr){

	####################################################################
	#  Create a vector composed of 1 or 0:[ 0 0 1 1 1 0 0 1 1 ....],   # 
	#  where 1 indicates intergenic region; 0 indicates this position  #
	#  lies within a gene on either forward or reverse strand          #
	####################################################################
	## load data
	operon.coords1 <- get.operon.coords()
	genome.length <- get.genome.length()
	
	coords <- operon.coords1[[chr]]
	## First modify the genome coordinates so that they include stop condons as part of the genes:
	coords[coords[,"Orientation"]=="For","Stop"] <- coords[coords[,"Orientation"]=="For","Stop"]+3
	coords[coords[,"Orientation"]=="Rev","Stop"] <- coords[coords[,"Orientation"]=="Rev","Stop"]-3


	intergenic.region <- array(1,genome.length[[chr]])
	for (i in 1: dim(coords)[1])
	{	
		if (coords[i,"Orientation"]=="For")
		{
			intergenic.region[coords[i,"Start"]:coords[i,"Stop"]] <- 0
		} else if (coords[i,"Orientation"]=="Rev")
		{
			intergenic.region[coords[i,"Stop"]:coords[i,"Start"]] <- 0
		}
	}

	###############################################
	#  Sort the coordinates on FORWARD and REVERSE 
	#  strands seperately 
	###############################################
	## Sort the genes based on Start sites
	for.coords = coords[coords[,"Orientation"]=="For",]
	index = order(for.coords[,"Start"])
	for.ordered.coords = for.coords[index,]

	# Reverse strand:
	# sort the genes based on Stop sites
	rev.coords = coords[coords[,"Orientation"]=="Rev",]
	index = order(rev.coords[,"Stop"])
	rev.ordered.coords = rev.coords[index,]

	###############################################################
	#  Go though each gene and search for upstream region for     #
	#  intergenic length. Calculate genes on FORWARD and REVERSE  #
	#  strand seperately                                          #
	###############################################################

	##################
	# FORWARD STRAND #
	##################

	FORWARD.intergenic.length <- array()

	## First gene: consider circular genome structure ##
	i <- 1
	in.start <- for.ordered.coords[i,"Start"]-1 # Start of the intergenic region upstream of the current gene
	if (intergenic.region[in.start]==1) # Check if this is within any intergenic region
	{
		## Calculate the shortest intergenic distance to the first gene on FORWARD strand
		for.upstream.gap <- genome.length[[chr]]-for.ordered.coords[dim(for.ordered.coords)[1],'Stop']+1+in.start
		
		## Calculate the shortest intergenic distance to the first gene on REVERSE strand
		if (sum(in.start > rev.ordered.coords[,"Start"])>0) 
		{
			index.r <- which (in.start > rev.ordered.coords[,"Start"] )
			rev.upstream.gap <- min(in.start - rev.ordered.coords[index.r,"Start"] + 1)
			
		} else # go back to the last gene on the reverse strand (consider circular genome)
		{
			rev.upstream.gap <- genome.length[[chr]] - rev.ordered.coords[dim(rev.ordered.coords)[1],"Start"] +1 + in.start
		}
		
		## Determine which one ( FORWARD or REVERSE ) is closer
		FORWARD.intergenic.length[i] <- min(for.upstream.gap,rev.upstream.gap)
		
	} else {
		FORWARD.intergenic.length[i] <- 0
	}

	## From the second gene to the last ##
	for (i in 2: dim(for.ordered.coords)[1])
	{	
		## Start of the intergenic region
		in.start <- for.ordered.coords[i,"Start"]-1
		
		## First check if this is within any intergenic region
		if (intergenic.region[in.start]==1)
		{
			## Calculate the shortest intergenic distance to the first gene on FORWARD strand
			for.upstream.gap <- in.start - for.ordered.coords[i-1,"Stop"] + 1
			
			## Calculate the shortest intergenic distance to the first gene on REVERSE strand
			# Only consider genes that have Start site upstream of for.ordered.coords[i,"Start"]
			if (sum(in.start > rev.ordered.coords[,"Start"])>0) 
			{
				index.r <- which (in.start > rev.ordered.coords[,"Start"] )
				rev.upstream.gap <- min(in.start - rev.ordered.coords[index.r,"Start"] + 1)
				
			} else # go back to the last gene on the reverse strand (consider circular genome)
			{
				rev.upstream.gap <- Inf
			}
					
			## Determine which one ( FORWARD or REVERSE ) is closer
			FORWARD.intergenic.length[i] <- min(for.upstream.gap,rev.upstream.gap)
				
			} else {
			
				FORWARD.intergenic.length[i] <- 0
			}
	}

	##################
	# REVERSE STRAND #
	##################
	REVERSE.intergenic.length <- array()

	## Last gene: consider circular genome structure
	i <- dim(rev.ordered.coords)[1]
	in.start <- rev.ordered.coords[i,"Start"]+1 # Start of the intergenic region upstream of the current gene
	if (intergenic.region[in.start]==1) # Check if this is within any intergenic region
	{
		## Calculate the shortest intergenic distance to the first gene on REVERSE strand
		rev.upstream.gap <- genome.length[[chr]] - in.start +1 + rev.ordered.coords[1,"Stop"]
		
		## Calculate the shortest intergenic distance to the first gene on FORWARD strand
		if (sum(in.start < for.ordered.coords[,"Start"])>0) 
		{
			index.f <- which (in.start < for.ordered.coords[,"Start"])
			for.upstream.gap <- min(for.ordered.coords[index.f,"Start"] - in.start + 1)
			
		} else # go back to the last gene on the reverse strand (consider circular genome)
		{
			for.upstream.gap <- genome.length[[chr]] - in.start +1 + for.ordered.coords[1,"Start"]
		}
		
		## Determine which one ( FORWARD or REVERSE ) is closer
		REVERSE.intergenic.length[i] <- min(rev.upstream.gap,for.upstream.gap)
		
	} else {
		REVERSE.intergenic.length[i] <- 0
	}

	## From the first gene to second to the last gene:
	for (i in 1:(dim(rev.ordered.coords)[1]-1))
	{	
		## Start of the intergenic region
		in.start <- rev.ordered.coords[i,"Start"]+1
		
		## First check if this is within any intergenic region
		if (intergenic.region[in.start]==1)
		{
			## Calculate the shortest intergenic distance to the first gene on REVERSE strand
			rev.upstream.gap <- rev.ordered.coords[i+1,"Stop"] - in.start  + 1
			
			## Calculate the shortest intergenic distance to the first gene on FORWARD strand
			# Only consider genes that have Start site upstream of rev.ordered.coords[i,"Start"]
			if (sum(in.start < for.ordered.coords[,"Start"])>0) 
			{
				index.f <- which (in.start < for.ordered.coords[,"Start"] )
				for.upstream.gap <- min(for.ordered.coords[index.f,"Start"] - in.start  + 1)
				
			} else # go back to the last gene on the reverse strand (consider circular genome)
			{
				for.upstream.gap <- Inf
			}
					
			## Determine which one ( FORWARD or REVERSE ) is closer
			REVERSE.intergenic.length[i] <- min(rev.upstream.gap, for.upstream.gap)
				
			} else {
			
				REVERSE.intergenic.length[i] <- 0
			}
	}

	######################################
	#  Combine both FORWARD and REVERSE  #
	######################################
	All.Intergenic.Length <- c(FORWARD.intergenic.length,REVERSE.intergenic.length)
	return(All.Intergenic.Length )
}