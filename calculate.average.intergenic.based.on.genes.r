load("U:/R/work/gene.coords.new.RData")
load("U:/R/work/operons.RData")
load("U:/R/work/halo_genome_strings.RData")


genome.length = list()
genome.length$HALCHR <- length(unlist(strsplit(seq.HALCHR,'')))
genome.length$pNRC100 <- length(unlist(strsplit(seq.PNRC100,'')))
genome.length$pNRC200 <- length(unlist(strsplit(seq.PNRC200,'')))

##################################################################
#   Combine genes in Operons into the same Transcription Units:  #
#   Output new coordinates that treat each operon as a whole     #
#   and ignore the gap between genes within operons              #
##################################################################
# First deal with the differences in names ( e.g. VNG1234C, VNG1234a,etc. ): 
# convert all gene names to VNG1234 for easy matching between genome coordinates and operons
all.genes <- rownames(gene.coords1) <- gsub("[1-9]{1}$","", rownames(gene.coords1), perl=T)
operons <- sapply(1:length(operons), function(x){gsub("[a-zA-Z]{1}$","", operons[[x]], perl=T)})

# Genes that are not part of any operon stay the same
GenesNotInOperon <- setdiff(all.genes,unique(unlist(operons)))
GenesNotInOperon.coords <- gene.coords1[GenesNotInOperon,]

# Combine genes that are in the same operon as one unit
# Make sure that only include the genes that are found within all.genes
GenesInOperon <- intersect(unique(unlist(operons)),all.genes)
names(operons) <- sapply(1:length(operons), function(x)operons[[x]][1])
operons <- operons[intersect(names(operons) ,all.genes)]
GenesInOperon.coords <- gene.coords1[names(operons),]

for (i in 1: length(operons))
{	# Genes on the forward strand
	if (as.character(gene.coords1[operons[[i]][1],"Orientation"])=="For")
	{
		operon.Start <- min(gene.coords1[operons[[i]],"Start"])
		operon.Stop <- max(gene.coords1[operons[[i]],"Stop"])
		GenesInOperon.coords[names(operons)[i],"Start"] <- operon.Start
		GenesInOperon.coords[names(operons)[i],"Stop"] <- operon.Stop
		
	} 
	# Genes on the reverse strand
	if (as.character(gene.coords1[operons[[i]][1],"Orientation"])=="Rev")
	{
		operon.Start <- max(gene.coords1[operons[[i]],"Start"])
		operon.Stop <- min(gene.coords1[operons[[i]],"Stop"])
		GenesInOperon.coords[names(operons)[i],"Start"] <- operon.Start
		GenesInOperon.coords[names(operons)[i],"Stop"] <- operon.Stop
	}
}
operon.coords <- rbind(GenesNotInOperon.coords,GenesInOperon.coords)

## Seperate the three chromosomes
operon.coords1 = list("HALCHR" = operon.coords[operon.coords[,"where"]=="HALCHR",], 
					  "pNRC100" = operon.coords[operon.coords[,"where"]=="pNRC100",], 
					  "pNRC200" = operon.coords[operon.coords[,"where"]=="pNRC200",])

save(operon.coords1,genome.length, file = "2012-11-08-for.calculating.intergenic.length.RData")
#-----------------------------------------------------------------------------------------
`get.operon.coords`<-
function() {
  if ( ! exists( "operon.coords1" ) ) try( load("2012-11-08-for.calculating.intergenic.length.RData")  )
  if ( exists( "operon.coords1" ) ) return( operon.coords1 )
}

`get.genome.length`<-
function() {
  if ( ! exists( "genome.length" ) ) try( load("2012-11-08-for.calculating.intergenic.length.RData")  )
  if ( exists( "genome.length" ) ) return( genome.length )
}

`get.all.intergenic.fragment.length` <- function(chr){
	####################################################################
	#  Create a vector composed of 1 or 0:[ 0 0 1 1 1 0 0 1 1 ....],   # 
	#  where 1 indicates intergenic region; 0 indicates this position  #
	#  lies within a gene on either forward or reverse strand          #
	####################################################################
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