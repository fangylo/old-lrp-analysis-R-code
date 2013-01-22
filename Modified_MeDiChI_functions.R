##### Some modified MeDiChI functions and other miscellaneous functions
##### Modified from original David's functions
##### Should be loaded lastly, since this modify some original functions

source("http://pinnacle/~dreiss/medichi.utils.R")

`medichi.clone.files1`<-
## Modified 'medichi.clone.files.r'for new high res 60 mer tiling array data. 
## load new kernel learned from 60mer array data Chr only.
## Do this first: kernel1=kernel$kernel
## For high res array set fit.res=10
function( files, fit.res=10, n.boot=10, verbose=F, kernel=kernel1,... ) {  
  data <- NULL
  for ( ff in files ) {
    cat( "LOADING:", ff, "\n" )
    data.tmp <- try( load.data.halo( "isb", data.file=ff, verbose=T, cor.cutoff=0.5, ... ) )
    if ( class( data.tmp ) == "try-error" ) { cat( "UH OH (1) - cant load data!\n" ); next }
    if ( nrow( data.tmp ) <= 5 ) { cat( "UH OH (2) - cant load data!\n" ); next }
    if ( unique( rownames( data.tmp ) ) > 3 )
      data.tmp <- data.tmp[ rownames( data.tmp ) %in% c( "HALCHR", "pNRC100", "pNRC200","NC"), ]
    if ( nrow( data.tmp ) <= 5 ) { cat( "UH OH (2) - cant load data!\n" ); next }
    data <- rbind( data, data.tmp )
  }
  if ( is.null( data ) || nrow( data ) <= 5 ) { cat( "UH OH (3) - cant load data!\n" ); next }
  # kernel1=get("kernel",envir=environment(medichi.clone.files1))
  # kernel1=kernel$kernel
  fits<-deconv.entire.genome( data,kernel=kernel,verbose=verbose,fit.res=fit.res, n.boot=n.boot,...);
  #try( data( "halo.lowres", package="MeDiChI" ) ) ## Load the deconv. kernel and gene.coords (halo)
  #fits <- deconv.entire.genome( data=data, max.steps=100, fit.res=fit.res, n.boot=n.boot,
                              # boot.sample="residual", kernel=kernel.halo.lowres, verbose=verbose, ... )
  #if ( class( fits ) == "try-error" ) { cat( "UH OH - cant load data!\n" ); print( fits ); next }
  #cat( "Saving fits to", rdata.file, "\n" )
  #save( data, fits, file=rdata.file, compress=T )
  fits
}

`plot.fractions.in.coding.rgns`<-function(fits,p.seq=NA,gap=0.05)
{
coeffs <- get.strongest.hits( fits, 0.999)
p.min=min(coeffs[,3])
p.max=max(coeffs[,3])
if (is.na(p.seq)) p.seq<-seq(p.min,p.max,by=gap)
hits.sel<-sapply(p.seq,function(i)fraction.in.coding.rgns(fits,p.cutoff=i))
plot(p.seq,hits.sel[6,],type="o",xlab="p.value cut off",ylab="*100=% within coding region")
}
`fraction.in.coding.rgns`<-
function( fits, p.cutoff=0.05, slop=0, chr.in=NA, ... ) {
# if feed in coeff instead of fits
	if (class(fits)=="matrix") coeffs <- fits
	else  coeffs <- get.strongest.hits( fits, p.cutoff )
  get.fraction.in.coding.rgns( coeffs, slop=slop, chr=chr.in, ... )
}

`get.fraction.in.coding.rgns`<-
#Original out and in were inverted
function( coeffs, max.npeaks=NA, slop=0, chr=NA ) {
  ##slop <- 20 ## Allow +/- 20 bp "slop" in gene start/stop sites -- see if this does anything.
  coords <- get.gene.coords.halo()
  if ( is.na( max.npeaks ) || max.npeaks > nrow( coeffs ) ) max.npeaks <- nrow( coeffs )
  if ( is.na( chr ) ) chr <- unique( as.character( coords$where ) )
  out <- NULL
  out.chr <- character()
  hits.in <- hits.out <- bg.in <- bg.out <- total.len <- coe.len <- 0
  for ( w in chr ) {
    if ( "where" %in% colnames( coeffs ) ) coe <- coeffs[ coeffs$where == w, ,drop=F ]
    else coe <- coeffs[ rownames( coeffs ) == w, ,drop=F ]
    if ( nrow( coe ) <= 0 ) next
    if ( nrow( coe ) > max.npeaks ) coe <- coe[ order( coe[ ,2 ], decreasing=TRUE ), ][ 1:max.npeaks, ]
    max.np <- nrow( coe )
    cc <- coords[ as.character( coords$where ) == w, ]
    coo <- rep( FALSE, max( c( cc$Start, cc$Stop ) ) )
    for ( i in 1:nrow( cc ) ) coo[ cc$Start[ i ]:cc$Stop[ i ] ] <- TRUE
    if ( slop > 0 ) {
      for ( i in 1:slop ) {
        c.up <- c( FALSE, coo[ 1:( length( coo ) - 1 ) ] )
        c.down <- c( coo[ 2:length( coo ) ], FALSE )
        coo <- coo & c.up & c.down
      }
    }
    ##if ( verbose ) cat( w, sum( coo ) / length( coo ), sum( ! coo ) / length( coo ), "\t" )
    coe[ ,1 ] <- round( coe[ ,1 ] )
    hits.in <- hits.in + sum( coo[ coe[ ,1 ] ], na.rm=T )
    hits.out <- hits.out + sum( ! coo[ coe[ ,1 ] ], na.rm=T )
    bg.in <- bg.in + sum( coo, na.rm=T )
    bg.out <- bg.out + sum( ! coo, na.rm=T )
    total.len <- total.len + length( coo )
    coe.len <- coe.len + nrow( coe )
##     if ( verbose ) {
##       cat( hits.in, hits.out, "\t" )
##       cat( pbinom( sum( coo[ coe[ ,1 ] ] ), nrow( coe ), sum( coo ) / length( coo ) ), "" )
##       cat( pbinom( sum( ! coo[ coe[ ,1 ] ] ), nrow( coe ), sum( ! coo ) / length( coo ), lower=F ), "\n" )
##     }
##     out.chr <- c( out.chr, w )
##     out <- rbind( out, c( max.np, slop, sum( coo ) / length( coo ), sum( ! coo ) / length( coo ),
##                          hits.in, hits.out,
##                          pbinom( sum( coo[ coe[ ,1 ] ] ), nrow( coe ), sum( coo ) / length( coo ) ),
##                          pbinom( sum( ! coo[ coe[ ,1 ] ] ), nrow( coe ), sum( ! coo ) / length( coo ), lower=F ) ) )
  }

  out <- c( coe.len, slop,  bg.out / total.len,bg.in / total.len,
           hits.out / coe.len,hits.in / coe.len, 
		   log10( pbinom( hits.out, coe.len, bg.out / total.len, lower=F ) ),
		   log10( pbinom( hits.in, coe.len, bg.in / total.len ) )
           )
  
  names( out ) <- c( "N. peaks", "Slop", "Expected Out", "Expected In", "Observed Out", "Observed In",
                       "log10-P-value Out", "log10-P-value In" )
  t( t( out ) )
}


# `get.genes.hit`<-
# # Original function asigned wrong coding/downstream region.
# function( fits, coeffs=NULL, p.cutoff=0.05, dist.cut=2000 ) {
  
  # if ( is.null( coeffs ) ) coeffs <- get.strongest.hits( fits, p.cutoff )
  # else coeffs <- get.strongest.hits( coeffs, p.cutoff )

  # coords <- get.gene.coords.halo()

  # is.for <- coords$Orientation == "For"
  # start <- as.integer( as.vector( coords$Start ) )
  # end <- as.integer( as.vector( coords$Stop ) )
  # wheres <- as.character( coords$where )
  # genes <- as.character( coords$canonical_Name )
  # name <- as.character( coords$Gene_Name )
  # chr <- as.character( coords$where )
  # names( start ) <- names( end ) <- names( is.for ) <- names( wheres ) <- names( name ) <- names( chr ) <- genes
 
  # out <- data.frame()
  # for( i in 1:nrow( coeffs ) ) {
    # dists <- coeffs[ i, 1 ] - start
    # hits <- which( abs( dists ) <= dist.cut & wheres == rownames( coeffs )[ i ] )
    # if ( length( hits ) <= 0 ) next
	
    # where.hit <- 
        # ifelse( dists[ hits ] <= 0 & is.for[ hits ], "upstream",
            # ifelse( dists[ hits ] >= 0 & ! is.for[ hits ], "upstream",
                # ifelse( dists[ hits ] > 0 & is.for[ hits ] & dists[ hits ] <= (end[ hits ]-start[ hits]), "coding",
                    # ifelse( dists[ hits ] < 0 & ! is.for[ hits ] & dists[ hits ] >= (end[ hits ]-start[ hits]), "coding",
                        # ifelse( dists[ hits ] > (end[ hits ]-start[ hits]) & is.for[ hits ], "downstream",
                            # ifelse( dists[ hits ] < (end[ hits ]-start[ hits])& ! is.for[ hits ], "downstream", "" ) ) ) ) ) )

    # out <- rbind( out, data.frame( name[ hits ], round( abs( dists[ hits ] ) ), where.hit, chr[ hits ],
                             # round( rep( coeffs[ i, 1 ], length( hits ) ) ),
                             # sprintf( "%.3f", rep( coeffs[ i, 2 ], length( hits ) ) ),
                             # sprintf( "%.5f", rep( coeffs[ i, 3 ], length( hits ) ) ) ) )
  # }
  # if ( nrow( out ) > 0 ) colnames( out ) <- c( "Gene", "Distance", "Where", "Chr",
                                              # "Pk.coord", "Pk.intens", "Pk.p.val" )
  # out
# }

####################################################################################
#  Modify the original get.genes.hit function so that how far a peak go into 
#  the coding region ( dist.cut.in ) can be specified. 
#  dist.cut.up: how far upstream a peak can be from the gene start site
#  dist.cut.in: how far into the coding region a peak can be from a gene start site
####################################################################################

`get.genes.hit`<-
function( fits, coeffs=NULL, p.cutoff=0.05, dist.cut.in =50, dist.cut.up =250 ) 
{  
  if ( is.null( coeffs ) ) coeffs <- get.strongest.hits( fits, p.cutoff )
  else coeffs <- get.strongest.hits( coeffs, p.cutoff )

  coords <- get.gene.coords.halo()

  is.for <- coords$Orientation == "For"
  start <- as.integer( as.vector( coords$Start ) )
  end <- as.integer( as.vector( coords$Stop ) )
  wheres <- as.character( coords$where )
  genes <- as.character( coords$canonical_Name )
  name <- as.character( coords$Gene_Name )
  chr <- as.character( coords$where )
  names( start ) <- names( end ) <- names( is.for ) <- names( wheres ) <- names( name ) <- names( chr ) <- genes
 
  out <- data.frame()
  for( i in 1:nrow( coeffs ) ) {
    dists <- coeffs[ i, 1 ] - start
	
    ###################################### this part is modified  ######################################
	# modify so that take upstream and downstream into consideration
	
	hits <- which(((abs(dists) <= dist.cut.up) & (dists <=0) &(wheres == rownames( coeffs )[ i ]) & (is.for)) |
	((abs(dists) <= dist.cut.up) & (dists >=0) &(wheres == rownames( coeffs )[ i ]) & (!is.for))|
	((abs(dists) <= dist.cut.in) & (dists >=0) &(wheres == rownames( coeffs )[ i ]) & (is.for))|
	((abs(dists) <= dist.cut.in) & (dists <=0) &(wheres == rownames( coeffs )[ i ]) & (!is.for)) )
	# hits <- which( abs( dists ) <= dist.cut & wheres == rownames( coeffs )[ i ] )
	####################################################################################################
	
    if ( length( hits ) <= 0 ) next

    where.hit <- 
        ifelse( dists[ hits ] <= 0 & is.for[ hits ], "upstream",
            ifelse( dists[ hits ] >= 0 & ! is.for[ hits ], "upstream",
                ifelse( dists[ hits ] > 0 & is.for[ hits ] & dists[ hits ] <= (end[ hits ]-start[ hits]), "coding",
                    ifelse( dists[ hits ] < 0 & ! is.for[ hits ] & dists[ hits ] >= (end[ hits ]-start[ hits]), "coding",
                        ifelse( dists[ hits ] > (end[ hits ]-start[ hits]) & is.for[ hits ], "downstream",
                            ifelse( dists[ hits ] < (end[ hits ]-start[ hits])& ! is.for[ hits ], "downstream", "" ) ) ) ) ) )

    out <- rbind( out, data.frame( name[ hits ], round( abs( dists[ hits ] ) ), where.hit, chr[ hits ],
                             round( rep( coeffs[ i, 1 ], length( hits ) ) ),
                             sprintf( "%.3f", rep( coeffs[ i, 2 ], length( hits ) ) ),
                             sprintf( "%.5f", rep( coeffs[ i, 3 ], length( hits ) ) ) ) )
  }
  if ( nrow( out ) > 0 ) colnames( out ) <- c( "Gene", "Distance", "Where", "Chr",
                                              "Pk.coord", "Pk.intens", "Pk.p.val" )
  out
  
}



