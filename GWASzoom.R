#!/usr/bin/env Rscript

## CPCantalapiedra 2017-2018
##
## GWASzoom
##
## An R script to load GWAS data in the surrounding region of a genetic marker
##
## Use case: a GWAS analysis was performed with 2 different sets of markers.
## For example, one set consists of low density SNPs and the other one of high density SNPs.
## The aim is to show both sets of markers, in a given region, to compare the association values
## and, hopefully, identify candidate genes.

# Load R libraries
suppressMessages(library(dplyr))
suppressMessages(library(lazyeval))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(gtable))
suppressMessages(library("bigmemory"))
suppressMessages(library("biganalytics"))

## Read parameters
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=10) {
  stop("10 arguments must be supplied (input file).n", call.=FALSE)
}

markers <- args[1] # a list of markers to analyze
markers_hapmap <- args[2] # the map and genotypic data of no-excap markers
markers_values <- args[3] # association values of no-excap markers
excap_hapmap <- args[4] # the map and genotypic data of exome capture markers
excap_values <- args[5] # association values of excap markers
genes_map <- args[6] # map data of genes
genes_annot <- args[7] # description of genes
interval <- args[8] # region to survey around each marker
outdir <- args[9] # output directory
pheno_file <- args[10] # phenotypic data

cat("Arguments:\n", file=stderr())
cat(paste("\t", markers, "\n"), file=stderr())
cat(paste("\t", markers_hapmap, "\n"), file=stderr())
cat(paste("\t", markers_values, "\n"), file=stderr())
cat(paste("\t", excap_hapmap, "\n"), file=stderr())
cat(paste("\t", excap_values, "\n"), file=stderr())
cat(paste("\t", genes_map, "\n"), file=stderr())
cat(paste("\t", genes_annot, "\n"), file=stderr())
cat(paste("\t", interval, "\n"), file=stderr())
cat(paste("\t", outdir, "\n"), file=stderr())
cat(paste("\t", pheno_file, "\n"), file=stderr())
cat("\n", file=stderr())

########################################### Functions
###########################################

################# Functions used to create a table of markers-genes
#################

# A function to retrieve the genes hit by each marker
genes_marker <- function(x, markers_pos, genes_pos){

    # data of the current marker
    marker <- markers_pos[x,]
    
    # genes hit by the marker
    # (marker$c is the chromosome)
    marker_genes <- genes_pos[genes_pos$chrom == marker$c &
                                genes_pos$start <= marker$pos &
                                genes_pos$end >= marker$pos,]
    
    # If the marker hits NO genes, create a row with empty gene information
    if (nrow(marker_genes)==0){
        marker_genes = data.frame(matrix(nrow = 1, ncol = 8))
        colnames(marker_genes) <- c("markerID", "markerChrom", "markerPos", "markerLogP", "gene", "start", "end", "desc")
        marker_genes[1,] = c(marker$rs, marker$c, marker$pos, marker$var, "-", "-", "-", "-")
    
    # If the marker hits some gene, join the gene and marker information
    } else {
        marker_genes <- marker_genes %>%
                                mutate(markerID = marker$rs,
                                        markerChrom = marker$c,
                                        markerPos = marker$pos,
                                        markerLogP = marker$var) %>%
                                select(markerID, markerChrom, markerPos, markerLogP,
                                        gene, start, end, desc)
    }
    
    return(marker_genes)
}

# A function to create the data for those genes without marker hits
f_genes_no_marker <- function(x, genes_no_markers) {
    
    # data of the current gene
    gene <- genes_no_markers[x,]
    
    # create a row with empty marker information
    gene_no_markers = data.frame(matrix(nrow = 1, ncol = 8))
    
    colnames(gene_no_markers) <- c("markerID", "markerChrom", "markerPos", "markerLogP",
                                    "gene", "start", "end", "desc")
        
    gene_no_markers[1,] = c("-", gene$chrom, gene$start, "-",
                            as.character(gene$gene), gene$start, gene$end, as.character(gene$desc))
    
    return(gene_no_markers)
}

# The function which creates the whole table of markers-genes
genes_table <- function(marker, marker_hapmap_data, other_markers_data, excap_hapmap_data,
                                marker_genes_map_data, genes_annot_data, plottable){
    
    print("Retrieving genes hit by each marker...", file=stderr())
    
    # retrieve marker name ("rs"), chromosome ("c"), and position ("pos")
    # for all the markers
    marker_pos <- marker_hapmap_data %>% select(rs, c, pos) # main marker data
    other_pos <- other_markers_data %>% select(rs, c, pos) # other no-excap markers
    excap_pos <- excap_hapmap_data %>% select(rs, c, pos) # excap markers
    
    # join all the markers and sort by chromosome and position
    markers_pos <- marker_pos %>% rbind(other_pos) %>% rbind(excap_pos) %>%
                    arrange(c, pos) %>%
                    left_join(plottable, by = c("rs" = "marker_ID" )) %>%
                    select(-mbpos, -type)
    
    # retrieve genes data
    genes_pos <- marker_genes_map_data %>% select(gene, chrom, start, end) %>%
                left_join(genes_annot_data, by = "gene")
    
    # join markers and genes
    markers_genes <- do.call("rbind",
                            lapply(seq(nrow(markers_pos)),
                            genes_marker, markers_pos = markers_pos, genes_pos = genes_pos))
    # note: I could use dplyr::bind_rows but it is yielding an error
    
    # genes without markers
    genes_no_markers <- genes_pos[!(genes_pos$gene %in% markers_genes$gene),]
    genes_no_markers <- do.call("rbind",
            lapply(seq(nrow(genes_no_markers)), f_genes_no_marker, genes_no_markers = genes_no_markers))
    
    # join all the markers and genes in a isngle table
    # and sort by chromosome and position
    markers_genes = rbind(markers_genes, genes_no_markers) %>%
                    arrange(markerChrom, markerPos)
    
    outfile=paste(outdir, "/", marker, ".genes.tab", sep="")
    write.table(markers_genes, file = outfile, quote = FALSE, sep = "\t", row.names = TRUE)
    
    print("Table of markers and genes created.", file=stderr())
    
    return();
}

### Auxiliar function to convert ATGC genotypes to numeric ones
numerical_genotype <- function(genotype){
  #print(genotype, file=stderr())
  nts <- unique(genotype)
  nts <- nts[nts %in% c("A", "C", "G", "T")]
  #print(nts, file=stderr())
  genotype[!genotype %in% nts] = -1 # for all which are not in ACGT
  # For those in ACGT:
  for (nt in nts){
    genotype[genotype == nt] = which(nts == nt) - 1
  }
  
  #print(genotype, file=stderr())
  return(genotype)
}

numerical_genotypes <- function(genotypes){
  
  genotypes <- apply(genotypes, 1, numerical_genotype)
  
  return(t(genotypes))
}

################# Graphical genotypes
graphical_genotypes <- function(marker, marker_hapmap_data, other_markers_data, excap_hapmap_data,
                                pheno_data){
  
  print("Generating graphical genotypes...")
  
  ## Keep only those columns to be used in the graphical genotypes
  marker_cols <- colnames(marker_hapmap_data)
  excap_cols <- colnames(excap_hapmap_data)
  common_cols <- intersect(marker_cols, excap_cols)
  #print(common_cols)
  
  genotypes <- rbind(select(marker_hapmap_data, common_cols), 
                     select(other_markers_data, common_cols),
                     select(excap_hapmap_data, common_cols))
  
  genotypes <- select(genotypes, -alleles, -c, -strand, -assembly, -center, -protLSID, -assayLSID, -panelLSID) %>%
    arrange(pos)
  
  rownames(genotypes) <- unlist(select(genotypes, rs))
  genotypes <- select(genotypes, -rs, -pos)
  #print(head(genotypes))
  
  # convert to numerical genotypes
  genotypes <- numerical_genotypes(genotypes)
  
  all <- genotypes
  # Join both phenotypic and genotypic data in a single table
  # all <- rbind(genotypes, phenotypes)
  
  ## Obtain phenotypic data
  entries <- colnames(genotypes)
  phenotypes <- filter(pheno_data, Genotype %in% entries)
  #print(t(phenotypes), file=stderr())
  
  all <- rbind(all, t(phenotypes))
  ## Write the table to a file
  outfile=paste(outdir, "/", marker, ".genotypes.tab", sep="")
  write.table(all, outfile, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
  
  cat("Table of graphical genotypes was written.\n", file=stderr());
  
  gmat <- apply(genotypes, 2, as.numeric)
  rownames(gmat) <- rownames(genotypes)
  #print(head(gmat), file=stderr())
  #gmat = data.matrix(gmat)
  
  library(gplots)
  library(leaflet)
  
  ########## Plot of graphical genotypes
  ##########
  
  envrange <- select(phenotypes, -Genotype)
  minrange = min(envrange, na.rm = T)# - 1
  maxrange = max(envrange, na.rm = T)# + 1
  
  colsidecolors <- colorNumeric(c("red", "gray", "blue"), minrange:maxrange)
  colsidecolors <- colsidecolors(envrange)
  
  mypalette <- colorRampPalette(c("red", "gray", "blue"))(n = 1000)
  
  outfile=paste(outdir, "/", marker, ".genotypes.pdf", sep="")
  pdf(outfile, paper="a4")
  # outfile=paste(outdir, "/", markeralt, ".", env, ".", type, ".genotypes.tiff", sep="")
  # tiff(outfile);#, height = 21, width = 21, units = "in", res = 300)
  # 1 is the column key over the plot
  # 2 is the plot
  lmat = rbind(c(4, 1),c(3,2))
  # lmat = rbind(4:3,2:1)
  lwid = c(0.5, 4)
  lhei = c(0.2, 12)
  
  heatmap.2(x = gmat, Rowv = FALSE, 
            key = FALSE, dendrogram = "none", #Colv = FALSE, 
            trace="none", cexCol=0.5, cexRow=0.3, col=mypalette,
            ColSideColors=colsidecolors, lmat=lmat, lwid=lwid, lhei=lhei,
            margins = c(3,4),
            sepcolor="black", sepwidth=c(0.1,0.1), colsep=c(0,ncol(gmat)), rowsep=c(0,nrow(gmat)))
  
  dev.off()
  
  cat("Plot of graphical genotypes was created.\n", file=stderr());
  
}

#############################################
## Plot
generate_plots <- function(marker, plottable){
  
  print("Generating plots...", file=stderr())
    
  ## Main plot function: data
  bfplot <- ggplot(plottable, aes(x=as.numeric(mbpos), y=as.numeric(var), 
                                  color=type, shape=type, linetype=type, size=type))+
    geom_point()+
    ylim(0, max(as.numeric(plottable$var)))
    theme_bw()
  
  ## Scale, legend and labels
  bfplot = bfplot +
    scale_color_manual("", values = c(alpha("red", 0.8), # GWAS
                                      alpha("black", 0.5), # gene
                                      alpha("blue", 0.8), # other markers
                                      alpha("blue", 1.0)), # marker
                       #breaks=c("B_GWAS", "gene", paste("Z",marker, sep="")),
                       labels=c("Markers", "Genes", "Others", marker))+
    scale_shape_manual("", values=c(16, 17, 5, 18), guide=FALSE)+
    scale_linetype_manual("", values=c(1,1,1,1), guide=FALSE) +
    scale_size_manual("", values=c(1.5,2.5,1.5,3.5), guide=FALSE)+
    guides(color=guide_legend(title=NULL))+
    labs(x = "Mbp", y = "-log10P", title = "")+
    theme(plot.title = element_text(hjust = 0.5), legend.background = element_rect(fill="gray90"),
           legend.position = "bottom")
  
  gt <- bfplot
  
  # As TIFF
  outfile=paste(outdir, "/", marker, ".tiff", sep="")
  tiff(outfile, height = 7, width = 7, units = "in", res = 300)
  grid.arrange(gt, nrow=1, ncol=1, newpage=FALSE)
  
  dev.off()
  cat(paste("Output to ", outfile, "\n", sep=""), file=stderr())
  
  # The same plot but with smooth loess adjusted function
  
  gt_smooth <- bfplot+
    geom_smooth(data=subset(plottable,type=="B_GWAS"),
                aes(x=as.numeric(mbpos), y=as.numeric(var), 
                    color=type, shape=type, linetype=type, size=type),
                method=loess, legend=FALSE, size=.5, se=TRUE, level=0.99)

  outfile=paste(outdir, "/", marker, ".smooth.tiff", sep="")
  tiff(outfile, height = 7, width = 7, units = "in", res = 300)
  grid.arrange(gt_smooth, nrow=1, ncol=1, newpage=FALSE)

  dev.off()
  cat(paste("Output to ", outfile, "\n", sep=""), file=stderr())
}

f_log <- function(x){
    #return(-log10(x)+1)
    return(-log10(x))
}

################################################
### Main function for each single marker
process_marker <- function(x){
  
  marker <- as.character(unlist(x))
  cat("\n******************** Processing marker:", marker, "\n\n", file=stderr())
  
  ### Obtain pvalue of marker
  #print(head(markers_values_data), file=stderr())
  marker_pvalue <- filter(markers_values_data, SNP == marker)[,2]
  print(paste("Previous P.value:", marker_pvalue), file=stderr())
  
  ### Obtain position of the marker
  marker_hapmap_data <- filter(markers_hapmap_data, rs == marker)
  #print(marker_hapmap_data, file=stderr())
  marker_chrom = marker_hapmap_data[1,3]
  marker_pos = marker_hapmap_data[1,4]
  print(paste("Chr:", marker_chrom, ", Pos:", marker_pos), file=stderr())
  
  ### Interval
  pos_start = as.numeric(marker_pos) - as.numeric(interval)
  if (pos_start < 0) pos_start = 0;
  pos_end = as.numeric(marker_pos) + as.numeric(interval)
  print(paste("Interval:", pos_start, "-", pos_end), file=stderr())
  
  ### Obtain other markers from the original set in the interval
  # defined by the position of the original marker plus the interval parameter
  other_markers_hapmap_data <- filter(markers_hapmap_data, c == marker_chrom &
                                       pos >= pos_start & pos <= pos_end & 
                                        rs != marker)
  
  other_markers_data <- select(other_markers_hapmap_data, rs, pos) %>%
                              mutate(mbpos = pos / 1000000) %>%
                              select(-pos) %>%
                              mutate(type = "other_markers")
  
  print(paste("Num other original markers found:", nrow(other_markers_data)), file=stderr())
  
  # obtain their pvalues
  #print(head(markers_values_data))
  other_markers_data = left_join(other_markers_data, markers_values_data, by = c("rs" = "SNP"))
  #print(nrow(other_markers_hapmap_data))
  # -log10(P.value)
  other_markers_data <- mutate(other_markers_data, var = ifelse(is.na(P.value), 0, f_log(P.value))) %>%
                        select(-P.value, marker_ID = rs)
  #print(other_markers_data)
  
  ### Obtain excap markers in the interval
  # defined by the position of the original marker plus the interval parameter
  marker_excap_hapmap_data <- filter(excap_hapmap_data, c == marker_chrom &
                                                        pos >= pos_start &
                                                        pos <= pos_end)
  print(paste("Num excap markers found:", nrow(marker_excap_hapmap_data)), file=stderr())
  #print(t(head(t(head(marker_excap_hapmap_data)))), file=stderr())
  excap_markers_list = select(marker_excap_hapmap_data, rs)[,1]
  #print(head(excap_markers_list), file=stderr())
  
  # Obtain GWAS values of those markers
  marker_excap_values_data <- filter(excap_values_data, SNP %in% excap_markers_list)
  print(paste("Num GWAS results found:", nrow(marker_excap_values_data)), file=stderr())
  #print(head(marker_excap_values_data), file=stderr())
  
  marker_excap_values_no_data <- excap_markers_list[!excap_markers_list %in% select(marker_excap_values_data, SNP)[,1]]
  #print(head(marker_excap_values_no_data), file=stderr())
  print(paste("Num GWAS results NOT found:", length(marker_excap_values_no_data)), file=stderr())
  
  # Assign 0 to P.value of not found markers
  marker_excap_values_no_data = data.frame(marker_excap_values_no_data)
  colnames(marker_excap_values_no_data) = c("SNP")
  marker_excap_values_no_data <- mutate(marker_excap_values_no_data, P.value = NA)
  #print(head(data.frame(marker_excap_values_no_data)))
  
  #### Join all the data of excap markers
  marker_excap_values_data <- rbind(marker_excap_values_data, marker_excap_values_no_data)
  
  # Calculate log10 of association values
  # -log10(P.value)
  marker_excap_values_data <- mutate(marker_excap_values_data, log10P = ifelse(is.na(P.value), 0, f_log(P.value)))
  
  # Join association values with map positions
  ### map positions --> marker_excap_hapmap_data
  ### association values --> marker_excap_values_data
  # print(head(marker_excap_hapmap_data))
  # print(head(marker_excap_values_data))
  marker_excap_data <- inner_join(marker_excap_hapmap_data, marker_excap_values_data, by = c("rs" = "SNP")) %>%
    mutate(mbpos = as.numeric(as.character(pos))/1000000) 
    #select("rs", "c", "pos", "P.value", "log10P") %>%
    
  # the last line is to reduce the magnitude of map positions (so that the value can be plotted easily)
  
  #print(nrow(marker_excap_data), file=stderr())
  #print(head(marker_excap_data), file=stderr())
  
  # Prepare for plot
  #print(marker_excap_data, file=stderr())
  marker_excap_data <- select(marker_excap_data, marker_ID = rs, mbpos = mbpos, var = log10P) %>% 
    mutate(type = "B_GWAS")
  # The "B" was used in the original script to make sure that the series with excap markers
  # was the second after the dummy series, so that the plot colors, etc
  # were applied correctly
  
  #print(head(marker_excap_data), file=stderr())
  ### Therefore, marker_excap_data is the table to plot markers and pvalues
  
  ### Genes data
  
  # Obtain genes in the interval
  print(paste("Num total genes", nrow(genes_map_data)), file=stderr())
  
  marker_genes_map_data <- filter(genes_map_data, chrom == marker_chrom &
                                       ((start >= pos_start & start <= pos_end) |
                                       (end >= pos_start & end <= pos_end) | 
                                         (start <= pos_start & end >= pos_end)))
  print(paste("Num genes in interval:", nrow(marker_genes_map_data)), file=stderr())
  
  # Obtain mid position of gene
  marker_genes_map_data <- mutate(marker_genes_map_data, pos = round(start + (end-start)/2))
  # The position to plot (mbmid) will be the mid position, unless the latter is outside the interval
  marker_genes_map_data <- mutate(marker_genes_map_data, mbmid = (ifelse(pos < pos_start, end, 
                                                                    ifelse(pos > pos_end, start, 
                                                                                            pos)))
                                                                /1000000) # show in Mbp
  #print(head(marker_genes_map_data, file=stderr()))
  
  # Prepare for the plot
  genestable <- select(marker_genes_map_data, mbpos = mbmid, marker_ID = gene) %>% 
    mutate(var = 0.0) %>% 
    mutate(type = "gene")
  
  #print(head(genestable, file=stderr()))
  
  ### Add the original marker to the table of markers
  # -log10(marker_pvalue)
  marker_excap_data <- rbind(marker_excap_data, c(marker_ID = marker, mbpos = marker_pos/1000000, 
                                                  var = f_log(marker_pvalue), type = paste("Z", marker, sep="")))
  # The "Z" is to make sure that the marker is the last series
  # (see above "The B...")
  #print(head(marker_excap_data), file=stderr())
  #print(head(other_markers_data), file=stderr())
  
  ### Join all the tables
  
  plottable <- marker_excap_data %>% rbind(other_markers_data) %>% rbind(genestable) %>% arrange(type)
  
  #print(plottable, file=stderr())
  
  generate_plots(marker, plottable)
  graphical_genotypes(marker, marker_hapmap_data, other_markers_hapmap_data, marker_excap_hapmap_data,
                      pheno_data)
                      
    genes_table(marker, marker_hapmap_data, other_markers_hapmap_data, marker_excap_hapmap_data,
                        marker_genes_map_data, genes_annot_data, plottable)
}

########################################### BEGIN
###########################################

gmat <- NA
# This is done here because add.expr of heatmap.2 only searches variables
# in the global environment, and not in the local function from where heatmap.2 is called
# This was known as a bug already in 2015

# Reading markers list
cat("Reading primary markers list...\n", file=stderr())
markers_list <- read.table(markers, header = FALSE)

# Reading markers data
cat("Reading primary markers hapmap data...\n", file=stderr())
markers_hapmap_data <- read.table(markers_hapmap, header = TRUE, sep="\t")
#print(t(head(t(head(markers_hapmap_data)))), file=stderr())
print(t(head(t(head(markers_hapmap_data)))), file=stderr())

# Reading markers association values
cat("Reading primary markers association values...\n", file=stderr())
markers_values_data <- read.table(markers_values, header = TRUE)
#print(head(markers_values_data), file=stderr())
print(head(markers_values_data), file=stderr())

# Reading excap data
cat("Reading secondary markers hapmap data...\n", file=stderr())
excap_hapmap_data <- read.table(excap_hapmap, header = TRUE, sep="\t")
#print(t(head(t(head(excap_hapmap_data)))), file=stderr())
print(t(head(t(head(excap_hapmap_data)))), file=stderr())

# Reading excap association values
cat("Reading secondary markers association values...\n", file=stderr())
excap_values_data <- read.table(excap_values, header = TRUE)
#print(head(excap_values_data), file=stderr())
print(head(excap_values_data), file=stderr())

# Reading genes map
cat("Reading genes map...\n", file=stderr())
genes_map_data <- read.table(genes_map, header = FALSE)
colnames(genes_map_data) <- c("chrom", "start", "end", "gene")
#print(head(genes_map_data), file=stderr())
print(head(genes_map_data), file=stderr())

# Reading genes annotation
cat("Reading genes annotation...\n", file=stderr())
genes_annot_data <- read.table(genes_annot, header = FALSE, sep = "\t")
colnames(genes_annot_data) <- c("gene", "desc")
#print(head(genes_annot_data), file=stderr())
print(head(genes_annot_data), file=stderr())

# Reading phenotypic data
cat("Reading phenotypic data...\n", file=stderr())
pheno_data <- read.table(pheno_file, header = TRUE, sep = "\t")
colnames(pheno_data) <- c("Genotype", "Phenotype")
#print(head(pheno_data), file=stderr())
print(head(pheno_data), file=stderr())

# Process each marker
devnull <- lapply(markers_list[,1], process_marker)
cat("\nFinished.\n************************\n\n", file=stderr())
q()

## END
