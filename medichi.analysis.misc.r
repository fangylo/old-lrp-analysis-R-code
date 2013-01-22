# setlist=list(asnC_1.1=coef(fitslist.each$asnC_1.1$fits.fin$HALCHR),trh4.1.2=coef(fitslist.each$trh4_1.2$fits.fin$HALCHR),trh6.2.2=coef(fitslist.each$trh6_1.2$fits.fin$HALCHR))
# setlist = list(hitslist_operons$)
# modified from [http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R <-more details see here] for ChIPChip peaks data
# much more scalable than the previous hist.overlap.parts.r function
# Examples: 2 way: setlist2 <- setlist[1:2]; OLlist2 <- overLapper(setlist=setlist2, p.cutoff=0.05,dis=200,complexity=1:length(setlist), sep="_",distance_to_geneStart=150, type="vennsets");
# OLlist2$Venn_List; counts <- sapply(OLlist2$Venn_List, length); vennPlot(counts=counts) 
# 3 way: OLlist_3 <- overLapper(setlist=setlist3, p.cutoff=0.05,dis=200,complexity=1:length(setlist), sep="_", distance_to_geneStart=NA,type="vennsets");
# counts <- sapply(OLlist_3$Venn_List, length); vennPlot(counts=counts, mysub="Top: var1; Bottom: var2", yoffset=c(0.3, -0.2))

overLapper <- function(setlist=setlist, p.cutoff=0.05,dis=200,complexity=1:length(setlist), sep="-",distance_to_geneStart=NA,type="vennsets") 
{
	
	## Create intersect matrix (removes duplicates!)
	# Filter data by p.cutoff (now each element in the setlist.p list is still a 3-column matrix with position, intensity ,p value)
	setlist.p <- sapply(names(setlist), function(x)  setlist[[x]][setlist[[x]][,3]<=p.cutoff,]) 
	# Now only keep the position part
	setlist.p <- sapply(names(setlist.p), function(x)  setlist.p[[x]][,1])
	
	# Filter the peak location based on their distance to gene start sites
	if (!is.na(distance_to_geneStart)){
	setlist.p<- sapply(names(setlist.p), function(x) filter.by.dis.to.geneStarts(setlist.p[[x]],distance=distance_to_geneStart))
	}
		setunion <- sort(unique(unlist(setlist.p)))
	# setmatrix3 <- sapply(names(setlist3), function(x) setunion3 %in% unique(setlist3[[x]])) 
	# setmatrix1 <- sapply(names(setlist1.p), function(x) sum(abs(setunion1[1]-setlist1.p[[1]])<dis)!=0
	setunion <-unique(round(setunion))
	setmatrix=matrix(0,nrow=length(setunion),ncol=length(setlist.p))
	tempA= as.numeric(vector(length=length(setunion)))
	
		for ( i in 1:length(setlist.p)){
			for (j in 1: length(setunion)){
			if (!sum(abs(setunion[j]-setlist.p[[i]])<=dis)) tempA[j]=0
			else tempA[j]=1		
			}
			 setmatrix[,i]=tempA
		}
	colnames(setmatrix)=names(setlist.p)
	rownames(setmatrix)<-setunion
	# rownames(setmatrix) <- setunion
	storage.mode(setmatrix) <- "numeric"

	## Create all possible sample combinations within requested complexity levels
	labels <- names(setlist)
	allcombl <- lapply(complexity, function(x) combn(labels, m=x, simplify=FALSE))
	allcombl <- unlist(allcombl, recursive=FALSE)
	complevels <- sapply(allcombl, length)
	
	## Return intersect list for generated sample combinations 
	if(type=="intersects") {
		OLlist <- sapply(seq(along=allcombl), function(x) setunion[rowSums(setmatrix[, rep(allcombl[[x]], 2)]) == 2 * length(allcombl[[x]])])
		names(OLlist) <- sapply(allcombl, paste, collapse=sep)
		return(list(Set_List=setlist.p, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Intersect_List=OLlist))
	}	

	## Return Venn intersect list for generated sample combinations 
	if(type=="vennsets") {
		vennSets <- function(setmatrix=setmatrix, allcombl=allcombl, index=1) {
			mycol1 <- which(colnames(setmatrix) %in% allcombl[[index]])
			mycol2 <- which(!colnames(setmatrix) %in% allcombl[[index]])
			cond1 <- rowSums(setmatrix[, rep(mycol1, 2)]) == 2 * length(mycol1)
			cond2 <- rowSums(setmatrix[, rep(mycol2, 2)]) == 0
			return(setunion[cond1 & cond2])
		}
		vennOLlist <- sapply(seq(along=allcombl), function(x) vennSets(setmatrix=setmatrix, allcombl=allcombl, index=x))
		names(vennOLlist) <- sapply(allcombl, paste, collapse=sep)
		return(list(Set_List.p=setlist.p, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Venn_List=vennOLlist))
	}
}

###########################################
## Define Venn Diagram Plotting Function ##
###########################################
vennPlot <- function(counts=counts, mymain="Venn Diagram", mysub="default", setlabels="default", yoffset=seq(0,10,by=0.34), ccol=rep(1,31), lcol=c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), lines=c("#FF0000", "#008B00", "#0000FF", "#FF00FF", "#CD8500"), mylwd=3, diacol=1, type="ellipse", ccex=1.0, lcex=1.0, ...) {
	## Enforce list structure to support multiple venn sets 
	if(is.list(counts)==FALSE) {
		counts <- list(counts)
	}
	
	## Check for supported number of Venn counts: 3, 7, 15 and 31
	if(!length(counts[[1]]) %in%  c(3,7,15,31)) stop("Only the counts from 2-5 way venn comparisons are supported.")
	
	## 2-way Venn diagram
	if(length(counts[[1]])==3) {
		## Define subtitle
		if(mysub=="default") {
			sample_counts <- sapply(names(counts[[1]])[1:2], function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), sep="")
		} else { 
			mysub <- mysub 
		}
		
		## Plot venn shapes
		symbols(x=c(4, 6), y = c(6, 6), circles=c(2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=F, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="",  xaxt="n", yaxt="n", bty="n", fg=lines, ...); 
		
		## Add counts
		for(i in seq(along=counts)) {
			olDF <- data.frame(x=c(3.1, 7.0, 5.0), 
                                           y=c(6.0, 6.0, 6.0), 
                                           counts=counts[[i]])
			text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...)
		}
                
		## Add sample labels
		if(length(setlabels)==1 & setlabels[1]=="default") { 
			setlabels <- names(counts[[1]][1:2])
		} else {
			setlabels <- setlabels
		}
		text(c(2.0, 8.0), c(8.8, 8.8), labels=setlabels, col=lcol, cex=lcex, ...)	
	}
 
	## 3-way Venn diagram
	if(length(counts[[1]])==7) { 
		## Define subtitle
		if(mysub=="default") {
			sample_counts <- sapply(names(counts[[1]])[1:3], function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), sep="")
		} else { 
			mysub <- mysub
		}
		
		## Plot venn shapes
		symbols(x=c(4, 6, 5), y=c(6, 6, 4), circles=c(2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=FALSE, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", fg=lines, ...)
		
		## Add counts
		for(i in seq(along=counts)) {
			olDF <- data.frame(x=c(3.0, 7.0, 5.0, 5.0, 3.8, 6.3, 5.0), 
                                           y=c(6.5, 6.5, 3.0, 7.0, 4.6, 4.6, 5.3), 
                                           counts=counts[[i]])
	        	text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...)

		}

                ## Add sample labels
		if(length(setlabels)==1 & setlabels[1]=="default") { 
			setlabels <- names(counts[[1]][1:3])
		} else {
			setlabels <- setlabels
		}
		text(c(2.0, 8.0, 6.0), c(8.8, 8.8, 1.1), labels=setlabels, col=lcol, cex=lcex, ...)	
	}
	
	## 4-way Venn diagram with ellipses
	if(length(counts[[1]])==15 & type=="ellipse") {
		## Define subtitle
		if(mysub=="default") {
			sample_counts <- sapply(names(counts[[1]])[1:4], function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
		} else { 
			mysub <- mysub
		}
		
		## Plot ellipse
		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments  
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])	
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}
		## Plot ellipse as 4-way venn diagram
		ellipseVenn <- function(...) {
			split.screen(c(1,1))
			plotellipse(center=c(3.5,3.6), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=mymain, sub=mysub, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.7,4.4), radius=c(2,4), rotate=-35, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(5.3,4.4), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.5,3.6), radius=c(2,4), rotate=35, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, lwd=mylwd, ...)
			## Add counts
			for(i in seq(along=counts)) {
				olDF <- data.frame(x=c(1.5, 3.5, 6.5, 8.5, 2.9, 3.1, 5.0, 5.0, 6.9, 7.1, 3.6, 5.8, 4.2, 6.4, 5.0), 
                                                   y=c(4.8, 7.2, 7.2, 4.8, 5.9, 2.2, 0.7, 6.0, 2.2, 5.9, 4.0, 1.4, 1.4, 4.0, 2.8), 
                                                   counts=counts[[i]])
				text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...)
			}
			## Add sample labels
			if(length(setlabels)==1 & setlabels[1]=="default") { 
				setlabels <- names(counts[[1]][1:4])
			} else {
				setlabels <- setlabels
			}
			text(c(0.4, 2.8, 7.5, 9.4), c(7.3, 8.3, 8.3, 7.3), labels=setlabels, col=lcol, cex=lcex, ...)
			close.screen(all=TRUE) 
		}
		ellipseVenn(...)
	} 

	## 4-way Venn diagram with circles (pseudo-venn diagram that misses two overlap sectors) 
	if(length(counts[[1]])==15 & type=="circle") {
		## Define subtitle
		if(mysub=="default") {
			sample_counts <- sapply(names(counts[[1]])[1:4], function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), sep="")
		} else { 
			mysub <- mysub
		}
		
		## Plot venn shapes
		symbols(x=c(4, 5.5, 4, 5.5), y = c(6, 6, 4.5, 4.5), circles=c(2, 2, 2, 2), xlim=c(0, 10), ylim=c(0, 10), inches=FALSE, main=mymain, sub=mysub, lwd=mylwd, xlab="", ylab="", xaxt="n", yaxt="n", bty="n", fg=lines, ...)
		
		## Add counts
		for(i in seq(along=counts)) {
		        olDF <- data.frame(x=c(3.0, 6.5, 3.0, 6.5, 4.8, 3.0, 4.8, 4.8, 6.5, 4.8, 3.9, 5.7, 3.9, 5.7, 4.8), 
                                           y=c(7.2, 7.2, 3.2, 3.2, 7.2, 5.2, 0.4, 0.4, 5.2, 3.2, 6.3, 6.3, 4.2, 4.2, 5.2), 
                                           counts=counts[[i]])
			text(olDF$x[-c(7,8)], olDF$y[-c(7,8)] + yoffset[i], olDF$counts[-c(7,8)], col=ccol, cex=ccex, ...) # rows 14-15 of olDF are printed in next step
			text(c(4.8), c(0.8) + yoffset[i], paste("Only in ", names(counts[[1]][1]), " & ", names(counts[[1]][4]), ": ", olDF$counts[7], "; Only in ", names(counts[[1]][2]), " & ", names(counts[[1]][3]), ": ", olDF$counts[8], sep=""), col=diacol, cex=ccex, ...)
		}

                ## Add sample labels
			if(length(setlabels)==1 & setlabels[1]=="default") { 
				setlabels <- names(counts[[1]][1:4])
			} else {
				setlabels <- setlabels
			}
		text(c(2.0, 7.5, 2.0, 7.5), c(8.3, 8.3, 2.0, 2.0), labels=setlabels, col=lcol, cex=lcex, ...)
	} 
	
	## 5-way Venn diagram
	if(length(counts[[1]])==31) {
		## Define subtitle
		if(mysub=="default") {
			sample_counts <- sapply(names(counts[[1]])[1:5], function(x) sum(counts[[1]][grep(x, names(counts[[1]]))]))
			mysub <- paste(paste("Unique objects: All =", sum(counts[[1]])), paste("; S1 =", sample_counts[1]), paste("; S2 =", sample_counts[2]), paste("; S3 =", sample_counts[3]), paste("; S4 =", sample_counts[4]), paste("; S5 =", sample_counts[5]), sep="")
		} else { 
			mysub <- mysub
		}
		
		## Plot ellipse
		plotellipse <- function (center=c(1,1), radius=c(1,2), rotate=1, segments=360, xlab="", ylab="", ...) {
			angles <- (0:segments) * 2 * pi/segments  
			rotate <- rotate*pi/180
			ellipse <- cbind(radius[1] * cos(angles), radius[2] * sin(angles))
			ellipse <- cbind( ellipse[,1]*cos(rotate) + ellipse[,2]*sin(rotate), ellipse[,2]*cos(rotate) - ellipse[,1]*sin(rotate) )
			ellipse <- cbind(center[1]+ellipse[,1], center[2]+ellipse[,2])	
			plot(ellipse, type = "l", xlim = c(0, 10), ylim = c(0, 10), xlab = "", ylab = "", ...)
		}
		## Plot ellipse as 5-way venn diagram
		ellipseVenn <- function(...) {
			split.screen(c(1,1))
			screen(1, new=FALSE)
			plotellipse(center=c(4.83,6.2), radius=c(1.43,4.11), rotate=0, segments=360, xlab="", ylab="", col=lines[1], axes=FALSE, main=mymain, sub=mysub, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.25,5.4), radius=c(1.7,3.6), rotate=66, segments=360, xlab="", ylab="", col=lines[2], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(6.1,3.5), radius=c(1.55,3.9), rotate=150, segments=360, xlab="", ylab="", col=lines[3], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(4.48,3.15), radius=c(1.55,3.92), rotate=210, segments=360, xlab="", ylab="", col=lines[4], axes=FALSE, lwd=mylwd, ...)
			screen(1, new=FALSE)
			plotellipse(center=c(3.7,4.8), radius=c(1.7,3.6), rotate=293.5, segments=360, xlab="", ylab="", col=lines[5], axes=FALSE, lwd=mylwd, ...)

			## Add counts
			for(i in seq(along=counts)) {
				olDF <- data.frame(x=c(4.85, 8.0, 7.1, 3.5, 2.0, 5.90, 4.4, 4.60, 3.60, 7.1, 6.5, 3.2, 5.4, 6.65, 3.40, 5.00, 6.02, 3.60, 5.20, 4.03, 4.20, 6.45, 6.8, 3.39, 6.03, 5.74, 4.15, 3.95, 5.2, 6.40, 5.1), 
                                                   y=c(8.30, 6.2, 1.9, 1.6, 5.4, 6.85, 6.6, 2.45, 6.40, 4.3, 6.0, 4.6, 2.1, 3.40, 3.25, 6.43, 6.38, 5.10, 2.49, 6.25, 3.08, 5.30, 4.0, 3.80, 3.20, 5.95, 5.75, 3.75, 3.0, 4.50, 4.6),
					counts=counts[[i]]) 
				text(olDF$x, olDF$y + yoffset[i], olDF$counts, col=ccol, cex=ccex, ...)
			}
			## Add sample labels
			if(length(setlabels)==1 & setlabels[1]=="default") { 
				setlabels <- names(counts[[1]][1:5])
			} else {
				setlabels <- setlabels
			}
			text(c(5.7, 7.9, 8.5, 4.2, 0.8), c(9.9, 7.9, 1.9, 0.0, 7.3), adj=c(0, 0.5), labels=setlabels, col=lcol, cex=lcex, ...)
			close.screen(all=TRUE) 
		}
		ellipseVenn(...)
	} 
}
# setlist = list(VNG1179C_early = hitslist.lrp2$VNG1179C_early,VNG1179C_late = hitslist.lrp2$VNG1179C_late,trh7_early=hitslist.lrp2$trh7_early,trh7_late = hitslist.lrp2$trh7_late)
#######################################
## overLapper1 function ## overLapper function starts from peak locations. overLapper1 starts from TGs.
#######################################
## Computation of (1) Venn Intersects and (2) Regular Intersects
# setlist = list(trh3_early = hitslist.lrp2$trh3_early, trh6_early = hitslist.lrp2$trh6_early,
				# trh6_late = hitslist.lrp2$trh6_late,asnC_late = hitslist.lrp2$asnC_late,trh4_late = hitslist.lrp2$trh4_late)
# setlist = list(asnC_early = hitslist_operons$asnC_1.1_2.1, asnC_late = hitslist_operons$asnC_2.2,VNG1237C_late = hitslist_operons$p1237_1.2_2.2)
# OLlist1 <- overLapper1(setlist = setlist)
# counts <- sapply(OLlist1$Venn_List, length); vennPlot(counts=counts, mysub="Top: var1; Bottom: var2", yoffset=c(0.3, -0.2))
overLapper1 <- function(setlist=setlist, complexity=1:length(setlist), sep="-", cleanup=FALSE, keepdups=FALSE, type = "vennsets") {
	## Clean up of sample sets to minimize formatting issues 
	if(cleanup==TRUE) {
		## Set all characters to upper case 
		setlist <- sapply(setlist, function(x) gsub("([A-Z])", "\\U\\1", x, perl=T, ignore.case=T))
		## Remove leading and trailing spaces
		setlist <- sapply(setlist, function(x) gsub("^ {1,}| {1,}$", "", x, perl=T, ignore.case=T))
	}
	
	## Append object counter to retain duplicates 
	if(keepdups==TRUE) {
		dupCount <- function(setlist=setlist) {
			count <- table(setlist)
			paste(rep(names(count), count), unlist(sapply(count, function(x) seq(1, x))), sep=".")
		}
		mynames <- names(setlist)
		setlist <- lapply(setlist, function(x) dupCount(x)) # lapply necessary for numeric data!
		names(setlist) <- mynames
	}	

	## Create intersect matrix (removes duplicates!)
	setunion <- sort(unique(unlist(setlist)))
	setmatrix <- sapply(names(setlist), function(x) setunion %in% unique(setlist[[x]])) 
	rownames(setmatrix) <- setunion
	storage.mode(setmatrix) <- "numeric"

	## Create all possible sample combinations within requested complexity levels
	labels <- names(setlist)
	allcombl <- lapply(complexity, function(x) combn(labels, m=x, simplify=FALSE))
	allcombl <- unlist(allcombl, recursive=FALSE)
	complevels <- sapply(allcombl, length)
	
	## Return intersect list for generated sample combinations 
	if(type=="intersects") {
		OLlist <- sapply(seq(along=allcombl), function(x) setunion[rowSums(setmatrix[, rep(allcombl[[x]], 2)]) == 2 * length(allcombl[[x]])])
		names(OLlist) <- sapply(allcombl, paste, collapse=sep)
		return(list(Set_List=setlist, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Intersect_List=OLlist))
	}	

	## Return Venn intersect list for generated sample combinations 
	if(type=="vennsets") {
		vennSets <- function(setmatrix=setmatrix, allcombl=allcombl, index=1) {
			mycol1 <- which(colnames(setmatrix) %in% allcombl[[index]])
			mycol2 <- which(!colnames(setmatrix) %in% allcombl[[index]])
			cond1 <- rowSums(setmatrix[, rep(mycol1, 2)]) == 2 * length(mycol1)
			cond2 <- rowSums(setmatrix[, rep(mycol2, 2)]) == 0
			return(setunion[cond1 & cond2])
		}
		vennOLlist <- sapply(seq(along=allcombl), function(x) vennSets(setmatrix=setmatrix, allcombl=allcombl, index=x))
		names(vennOLlist) <- sapply(allcombl, paste, collapse=sep)
		return(list(Set_List=setlist, Intersect_Matrix=setmatrix, Complexity_Levels=complevels, Venn_List=vennOLlist))
	}
}


######################
##  nearby genes    ## 
######################
# load "genc.RData"
# location_v: length n. a 1xn vector. 
# distance: within this distance find the start site of genes. 
# The last column include the distance from the gene Start_site to the query location

`genes.nearby`<-function(location_v,distance=200){
		load('U://R/work/gene.coords.new.RData')
		gene.coords.HALCHR<-gene.coords[gene.coords$where=="HALCHR",]
			row_id<-function(location,distance=distance,coords=gene.coords.HALCHR){
					id=which(abs(location-coords$Start)<=distance)
					return(id) # represent the corresponding row number on gene,coords.HALCHR matrix. can be matched to get the geneNames.
					}
			ind=list();
			for (i in 1: length(location_v)){
			ind[[i]]=row_id(location_v[i],distance=distance,coords=gene.coords.HALCHR)
	}
index=unique(unlist(ind));
gene_matrix=gene.coords.HALCHR[index,]	
gene_List=as.character(gene_matrix[,1])
return(list(gene_Matrix=gene_matrix,gene_List=gene_List))
}
########################################### 
#Chr = "HALCHR" or "pNRC100" or "pNRC200"

`genes.nearby1`<-function(location_v,distance=200,Chr = "HALCHR"){
		load('U://R/work/gene.coords.new.RData')
		gene.coords.ref <-gene.coords[gene.coords$where == Chr,]
			row_id<-function(location,distance=distance,coords=gene.coords.ref){
					id=which(abs(location-coords$Start)<=distance)
					return(id) # represent the corresponding row number on gene,coords.HALCHR matrix. can be matched to get the geneNames.
					}
			ind=list();
			for (i in 1: length(location_v)){
			ind[[i]]=row_id(location_v[i],distance=distance,coords=gene.coords.ref)
	}
index=unique(unlist(ind));
gene_matrix=gene.coords.ref[index,]	
gene_List=as.character(gene_matrix[,1])
return(list(gene_Matrix=gene_matrix,gene_List=gene_List))
}
####################################################################
##  filter peaks by their location relative to gene start sites   ## 
####################################################################

`filter.by.dis.to.geneStarts`<-function(location_v,distance=distance_to_geneStart){
	gene.coords.HALCHR<-gene.coords[gene.coords$where=="HALCHR",]
	# coords=temp.coords
	location_id=list()
		for (i in 1:dim(gene.coords.HALCHR)[1]){
		location_id[[i]]=which(abs(location_v-gene.coords.HALCHR$Start[i])<=distance)
		}
	location_id=unique(unlist(location_id))
	return(location_v[location_id])
}

###################################################################################################
#### From the folder of all gene.hits.csv files, generate hclust and and do pairwise comparison  ##
###################################################################################################
# 20110218 .have to first setwd to the current folder that have all the hits file generated
# pattern1 : pattern to match the hits. files in the folder.
# pattern2 : pattern to match to the original long file names and to be trimmed


`hits.overlap.hclust`<-function(setfolder=getwd(),pattern1='_filtered.csv',pattern2="_p.0.[0-9]{2}_dist.[0-9]{1,}bp_filtered.csv",pValue=pValue,distance=distance){
# library('pvclust')
hitsfiles=list.files(setfolder,pattern1)
print (hitsfiles)
hits_Lst =list()
	for (i in 1:length(hitsfiles)){
	hits_Lst[[i]]=read.csv(file=hitsfiles[i])
	}
	names.temp <- lapply(hitsfiles,function(x)gsub(pattern2,"",x))
	names(hits_Lst) = unlist(names.temp)
	
	complexity=1:length(hits_Lst)
	labels <- names(hits_Lst)
	# allcomb <- lapply(complexity, function(x) combn(labels, m=x, simplify=FALSE))
	allcomb2 <- combn(labels, m=2, simplify=FALSE) #choose only pairwise ones
	# complevels <- sapply(allcombl, length)
	intersect_genes<-list()
	p.intersect<-list()
	pre_dist_matrix=matrix(ncol=length(hits_Lst),nrow=length(hits_Lst))
	rownames(pre_dist_matrix)=colnames(pre_dist_matrix)=names(hits_Lst)
		for (i in 1: length(allcomb2)){
			combi <- which(names(hits_Lst) %in% allcomb2[[i]]) 
			num1<-dim(hits_Lst[[combi[1]]])[1]
			num2<-dim(hits_Lst[[combi[2]]])[1]
			intersect_genes[[i]]<-intersect(hits_Lst[[combi[1]]][,1],hits_Lst[[combi[2]]][,1])
				if (length(intersect_genes[[i]])!=0){
				p.intersect[[i]]=phyper(length(intersect_genes[[i]]),num1,2400-num1,num2,lower.tail=F)
				} else p.intersect[[i]]=1
			pre_dist_matrix[combi[1],combi[2]]=pre_dist_matrix[combi[2],combi[1]]=p.intersect[[i]]
		}
		names(intersect_genes)=sapply(allcomb2, paste, collapse='_')

dist_m<-as.dist(pre_dist_matrix)

#hclust
hc<-hclust(dist_m)
temp<-gsub(" |:","-",substr(Sys.time(),1,19))

pdf(file=paste('hclust based on gene.hits_p.',pValue,'_dist.',distance,'bp_',temp,'.pdf',sep=''))
plclust(hc,hang=0.1,unit=T)
title(paste('hclust based on gene hits_p.',pValue,'_dist.',distance,'bp'))
dev.off()

# pvclust:boostrap hclust
# pvc<-pvclust(pre_dist_matrix, method.hclust = "complete",method.dist="euclidean",nboot=1000)
# pdf(file=paste('hclust w bootstrap(1000) based on gene.hits_p.',pValue,'_dist.',distance,'bp_',temp,'.pdf',sep=''))
# plot(pvc)
# pvrect(pvc)
# title(sub=paste('hclust w/bootstrap(1000) based on gene hits_p.',pValue,'_dist.',distance,'bp'))
# dev.off()


dir.create(paste('intersect.genes_','dis.',distance,'_p.',pValue,'_',temp,sep=''))
setwd(paste(getwd(),'/',paste('intersect.genes_','dis.',distance,'_p.',pValue,'_',temp,sep=''),sep=''))

	for (i in 1:length(intersect_genes)){
	write.csv(intersect_genes[[i]],file=paste(as.name(names(intersect_genes)[i]),'.csv',sep=''),quote=F,row.names=F)
	}
	
	temp.m <- as.matrix (dist_m)
	
	for ( i in 1:dim(temp.m)[1]){
	temp.m [i,i] = NA}
	write.csv (temp.m,file = paste ('Distance matrix_p.',pValue,'_dis.',distance,'.csv',sep=''),quote=F)
return(list(hc=hc,dist_m=dist_m,intersect_genes=intersect_genes))

}

#################################################
# get operons list
#################################################

`get.operons.halo`<-
function() {
  if ( ! exists( "operons" ) ) try( load('U://R/work/operons.RData')  )
  if ( exists( "operons" ) ) return( operons )
  # try( load( "data/halo.coords.RData" ) ) ## loads "halo.coords"
  operons
}


##################################################
## get downstream gene list in the same operon ###
##################################################

`get.operons.list`<-function(genelist){
# genelist: a character vector of genes canonical name
# operons: operon list
    # genelist <- as.character(sapply (genelist,function(x)substr(x,1,8)))
	operons <- get.operons.halo()
	coords <- get.gene.coords.halo()
	gene.matrix <- coords[coords$canonical_Name %in% genelist,]
	is.for <- gene.matrix$Orientation == "For"
	
	operon.list=list()
	for (i in 1:dim(gene.matrix)[1]){
		if (!(gene.matrix[i,1] %in% unlist(operons))) operon.list[[i]]<-as.character(gene.matrix[i,1])
		else{	for (j in 1:length(operons)){
				if (gene.matrix[i,1] %in% operons[[j]]){
				indx = which (operons[[j]] %in% gene.matrix[i,1])
				if (is.for[i]) operon.list[[i]] = operons[[j]][indx:length(operons[[j]])] 
				else operon.list[[i]] = operons[[j]][1:indx]
				} 
			}
			}
	}
	return(operon.list)
}
##################################################
## get all genes  in the same operon ###
##################################################

`get.all.operons.list`<-function(genelist){
# genelist: a character vector of genes canonical name
# operons: operon list
    # genelist <- as.character(sapply (genelist,function(x)substr(x,1,8)))
	operons <- get.operons.halo()
	coords <- get.gene.coords.halo()
	gene.matrix <- coords[coords$canonical_Name %in% genelist,]
	is.for <- gene.matrix$Orientation == "For"
	
	operon.list=list()
	for (i in 1:dim(gene.matrix)[1]){
		if (!(gene.matrix[i,1] %in% unlist(operons))) operon.list[[i]]<-as.character(gene.matrix[i,1])
		else{	for (j in 1:length(operons)){
				if (gene.matrix[i,1] %in% operons[[j]]){
					indx = which (operons[[j]] %in% gene.matrix[i,1])
					operon.list[[i]] = operons[[j]]
				} 
			}}
	}
	return(operon.list)
}
####################################################################
## get all genes  in the same operon  for calculating topGO      ###
####################################################################

`get.all.operons.list_Lrp.tg` <- function(genelist,Lrp,Lrp.ChIP.targets,corr.coeff.matrix,cor.cutoff = 0.2,DirectionOfCorrelation="pos")
{
# Include genes that are (1) on the same operon & downstream of the input gene(s); 
#                        (2) on the same operon & upstream of the input gene(s) & ChIP target;
#                        (3) on the same operon & upstream of the input gene(s) & correlated with lrp ( same direction and pass a cutoff); 
# genelist: a character vector of genes canonical name
# operons: operon list
# DirectionOfCorrelation = "pos" or "neg"
    
	operons <- get.operons.halo()
	coords <- get.gene.coords.halo()
	gene.matrix <- coords[coords$canonical_Name %in% genelist,]
	is.for <- gene.matrix$Orientation == "For"
	targets <- Lrp.ChIP.targets[[Lrp]]
	
	
	operon.list=list()
	for (i in 1:dim(gene.matrix)[1])
	{
		if (!(gene.matrix[i,1] %in% unlist(operons))) operon.list[[i]]<-as.character(gene.matrix[i,1])
		else{	
				for (j in 1:length(operons))
				{
					if (gene.matrix[i,1] %in% operons[[j]])
					{
						indx = which (operons[[j]] %in% gene.matrix[i,1])
						# operon.list[[i]] = operons[[j]]
						if (is.for[i])
						{
							downsreamhits <- operons[[j]][indx:length(operons[[j]])]
							upstreamhits <- character()
							
							if (indx > 1)
							{
								count <- 1
								for (k in 1: (indx-1))
								{
									# if it's a target:included. correlation coeff.
									if (operons[[j]][k] %in% rownames(corr.coeff.matrix))
									{
										if ((operons[[j]][k] %in% targets)|
										((corr.coeff.matrix[operons[[j]][k],Lrp] >= cor.cutoff)&(DirectionOfCorrelation =="pos"))|
										((corr.coeff.matrix[operons[[j]][k],Lrp] <= -cor.cutoff)&(DirectionOfCorrelation =="neg")))
										{
											upstreamhits[count] <- operons[[j]][k];
											count <- count+1;
										}
									} else if (!operons[[j]][k] %in% rownames(corr.coeff.matrix)) 
									{
										if (operons[[j]][k] %in% targets)
										{
											upstreamhits[count] <- operons[[j]][k];
											count <- count+1;
										}
									} 
								}
							}
							operon.list[[i]] <- c(upstreamhits,downsreamhits)
						  
						} else 
						{
							downsreamhits <- operons[[j]][1:indx]
							upstreamhits <- character()
							if (indx <length(operons[[j]]) )
							{
								count <- 1
								for (k in (indx+1): length(operons[[j]]))
								{
									if (operons[[j]][k] %in% rownames(corr.coeff.matrix))
									{
										if ((operons[[j]][k] %in% targets)|
										((corr.coeff.matrix[operons[[j]][k],Lrp] >= cor.cutoff)&(DirectionOfCorrelation =="pos"))|
										((corr.coeff.matrix[operons[[j]][k],Lrp] <= -cor.cutoff)&(DirectionOfCorrelation =="neg")))
										{
											upstreamhits[count] <- operons[[j]][k];
											count <- count+1;
										}
									} else if (!operons[[j]][k] %in% rownames(corr.coeff.matrix)) 
									{
										if (operons[[j]][k] %in% targets)
										{
											upstreamhits[count] <- operons[[j]][k];
											count <- count+1;
										}
									}
								}
							}
							operon.list[[i]] <- c(upstreamhits,downsreamhits)
						} 
					} 
				}
			}
	}
	return(operon.list)
}


`get.operons.list`<-function(genelist){
# genelist: a character vector of genes canonical name
# operons: operon list
    # genelist <- as.character(sapply (genelist,function(x)substr(x,1,8)))
	operons <- get.operons.halo()
	coords <- get.gene.coords.halo()
	gene.matrix <- coords[coords$canonical_Name %in% genelist,]
	is.for <- gene.matrix$Orientation == "For"
	
	operon.list=list()
	for (i in 1:dim(gene.matrix)[1])
	{
		if (!(gene.matrix[i,1] %in% unlist(operons))) operon.list[[i]]<-as.character(gene.matrix[i,1])
		else
		{	for (j in 1:length(operons))
			{
				if (gene.matrix[i,1] %in% operons[[j]])
				{
					indx = which (operons[[j]] %in% gene.matrix[i,1])
					if (is.for[i]) operon.list[[i]] = operons[[j]][indx:length(operons[[j]])] 
					else operon.list[[i]] = operons[[j]][1:indx]
				} 
			}
		}
	}
	return(operon.list)
}
##############################################################
## get fractions of peaks in proximity to gene start & end ###
##############################################################

`get.fraction.close.startoend`<-
function( coeffs, max.npeaks=NA, slop=0, chr=NA,dist2start = 200 ,dist2end = 200 ) {
  ##slop <- 20 ## Allow +/- 20 bp "slop" in gene start/stop sites -- see if this does anything.
  coords <- get.gene.coords.halo()
  if ( is.na( max.npeaks ) || max.npeaks > nrow( coeffs ) ) max.npeaks <- nrow( coeffs )
  if ( is.na( chr ) ) chr <- unique( as.character( coords$where ) )
  
  out <- NULL
  # out.chr <- character()
  hits.in.start <- hits.in.end <- hits.in <- hits.out.start <- hits.out.end <- hits.out <- bg.in.start <- bg.in.end <- bg.in <- bg.out.start <- bg.out.end <- bg.out <- total.len <- coe.len <- 0
  
  for ( w in chr ) {
    if ( "where" %in% colnames( coeffs ) ) coe <- coeffs[ coeffs$where == w, ,drop=F ]
    else coe <- coeffs[ rownames( coeffs ) == w, ,drop=F ]
	
    if ( nrow( coe ) <= 0 ) next
    if ( nrow( coe ) > max.npeaks ) coe <- coe[ order( coe[ ,2 ], decreasing=TRUE ), ][ 1:max.npeaks, ]
    max.np <- nrow( coe )
    cc <- coords[ as.character( coords$where ) == w, ]
	
	coo.start <- rep( FALSE, max( c( cc$Start, cc$Stop ) ) )
	coo.end <- rep( FALSE, max( c( cc$Start, cc$Stop ) ) )
	coo <- rep( FALSE, max( c( cc$Start, cc$Stop ) ) )
	
	for ( i in 1:nrow( cc ) )
		{
			coo.start[ (cc$Start[ i ]-dist2start):(cc$Start[ i ]+dist2start)] <- TRUE # close to start
			coo.end[ (cc$Stop[ i ]-dist2end):(cc$Stop[ i ]+dist2end)] <- TRUE # close to end
			coo[ cc$Start[ i ]:cc$Stop[ i ] ] <- TRUE # in gene
		}
	
	
    ##if ( verbose ) cat( w, sum( coo ) / length( coo ), sum( ! coo ) / length( coo ), "\t" )
    coe[ ,1 ] <- round( coe[ ,1 ] )
	
    hits.in <- hits.in + sum( coo[ coe[ ,1 ] ], na.rm=T )
	hits.in.start <- hits.in.start + sum( coo.start[ coe[ ,1 ] ], na.rm=T )
	hits.in.end <- hits.in.end + sum( coo.end[ coe[ ,1 ] ], na.rm=T )
		
    hits.out <- hits.out + sum( ! coo[ coe[ ,1 ] ], na.rm=T )
	hits.out.start <- hits.out.start + sum( ! coo.start[ coe[ ,1 ] ], na.rm=T )
	hits.out.end <- hits.out.end + sum( ! coo.end[ coe[ ,1 ] ], na.rm=T )
	
    bg.in <- bg.in + sum( coo, na.rm=T )
	bg.in.start <- bg.in.start + sum( coo.start, na.rm=T )
	bg.in.end <- bg.in.end + sum( coo.end, na.rm=T )
    bg.out <- bg.out + sum( ! coo, na.rm=T )
	bg.out.start <- bg.out.start + sum( ! coo.start, na.rm=T )
	bg.out.end <- bg.out.end + sum( ! coo.end, na.rm=T )
	
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

  out <- c( coe.len, slop,  bg.out.start / total.len, bg.out.end / total.len, bg.out / total.len,
			bg.in.start / total.len, bg.in.end / total.len, bg.in / total.len,
			hits.out.start / coe.len, hits.out.end / coe.len,hits.out / coe.len,
			hits.in.start / coe.len, hits.in.end / coe.len,hits.in / coe.len,
			pbinom( hits.in.start, coe.len, bg.in.start / total.len, lower=F )
			# pbinom( hits.out, coe.len, bg.out / total.len, lower=F  ),
			# pbinom( hits.in, coe.len, bg.in / total.len ) 
           )
  
  names( out ) <- c( "N. peaks", "Slop", "Expected Out of start", "Expected Out of end", "Expected Out " ,
					"Expected In Start", "Expected In in End", "Expected In", 
					"Observed Out of start", "Observed Out of end ", "Observed Out", 
					"Observed In of start", "Observed In of end", "Observed In",
					"p-value close to start"
                    # "P-value Out", "P-value In" 
					)
					
	 out1 <- c( coe.len,  bg.out / total.len,
			bg.in.start / total.len,  bg.in / total.len,
			hits.out / coe.len,
			hits.in.start / coe.len, hits.in.end / coe.len,hits.in / coe.len,
			pbinom( hits.in.start, coe.len, bg.in.start / total.len, lower=F )
			# pbinom( hits.out, coe.len, bg.out / total.len, lower=F  ),
			# pbinom( hits.in, coe.len, bg.in / total.len ) 
           )
	names( out1 ) <- c( "N. peaks",  "Expected Out " ,
					"Expected In Start",  "Expected In", 
					"Observed Out", 
					"Observed In of start", "Observed In of end", "Observed In",
					"p-value close to start"
                    # "P-value Out", "P-value In" 
					)
			   
	
  t( t( out1 ) )
}
