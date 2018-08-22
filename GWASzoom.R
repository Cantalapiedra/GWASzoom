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

write("Running GWASzoom...", file=stderr())

# Load R libraries
suppressMessages(library(dplyr))
suppressMessages(library(lazyeval))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(gtable))
suppressMessages(library("bigmemory"))
suppressMessages(library("biganalytics"))
suppressMessages(library("optparse"))
suppressMessages(library(gplots))
suppressMessages(library(leaflet))

## Parameters specification

option_list = list(
    make_option("--marker", type="character", default=NULL, 
        help="A file with the name of the marker to be processed.
        This will be the main marker; the zoom will be performed to its position.",
        metavar="char"),
    
    make_option("--markers_hapmap", type="character", default=NULL, 
        help="A hapmap file with genotyping of markers. This file should include the markers in --markers file.",
        metavar="FILE"),
    
    make_option("--markers_values", type="character", default=NULL, 
        help="A file with association values for the markers in the --markers_hapmap file.",
        metavar="FILE"),
        
    make_option("--excap_hapmap", type="character", default=NULL, 
        help="A hapmap file with genotyping of markers, different than those in --markers_hapmap file.",
        metavar="FILE"),
        
    make_option("--excap_values", type="character", default=NULL, 
        help="A file with association values for the markers in the --excap_hapmap file.",
        metavar="FILE"),
        
    make_option("--pheno_file", type="character", default=NULL, 
        help="A file with the phenotypic values used for the association.",
        metavar="FILE"),
        
    make_option("--outdir", type="character", default="./", 
        help="Base output directory. Subdirectories will be created for each marker in --markers.",
        metavar="char"),
        
    make_option("--interval", type="numeric", default=-1, 
        help="If this value is present, the zoom will be performed to a region defined
        from upstream --interval distance to the main marker (see --markers),
        to downstream --interval distance to the main marker.
        If --interval is not present, the output will include all the range of positions
        of the input hapmap files (i.e. zoom will be not performed).",
        metavar="int"),
        
    make_option("--gwastype", type="character", default="pvalue", 
        help='Type of GWAS value ("pvalue", from GAPIT for example; or "factor", from Bayenv2 for example).',
        metavar="char"),
        
    make_option("--alleles_format", type="character", default="biallelic", 
        help='Type of alleles in the hapmap files or how they are coded: either "biallelic" or "monoallelic").',
        metavar="char"),
        
    make_option("--ldr2file", type="character", default=NULL, 
        help="A file with LD r2 values to plot along with the association values.",
        metavar="FILE"),
    
    make_option("--ldthres", type="numeric", default=0.5, 
        help="This LD threshold will be used just to reduce the plot of graphical genotypes based on LD decay.
        The other plots and tables will be not affected by this parameter.",
        metavar="float"),
    
    # Removed because the table of genes was not working properly
    #make_option("--genes_annot", type="character", default=NULL, 
    #    help="Annotation features.",
    #    metavar="FILE"),
        
    make_option("--genes_map", type="character", default=NULL, 
        help="Map position of annotation features.",
        metavar="FILE"),
    
    make_option(c("-v", "--verbose"), action="store_true", default=FALSE)
); 

## Read and check parameters

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$marker)){
    message("Missing the --marker parameter, which is mandatory.")
    print_help(opt_parser)
    stop()
} else {
    markers = opt$marker
}

if (is.null(opt$markers_hapmap)){
    message("Missing the --markers_hapmap parameter, which is mandatory.")
    print_help(opt_parser)
    stop()
} else {
    markers_hapmap = opt$markers_hapmap
}

if (is.null(opt$markers_values)){
    message("Missing the --markers_values parameter, which is mandatory.")
    print_help(opt_parser)
    stop()
} else {
    markers_values = opt$markers_values
}

if (is.null(opt$excap_hapmap)){
    message("Missing the --excap_hapmap parameter, which is mandatory.")
    print_help(opt_parser)
    stop()
} else {
    excap_hapmap = opt$excap_hapmap
}

if (is.null(opt$excap_values)){
    message("Missing the --excap_values parameter, which is mandatory.")
    print_help(opt_parser)
    stop()
} else {
    excap_values = opt$excap_values
}

if (is.null(opt$pheno_file)){
    message("Missing the --pheno_file parameter, which is mandatory.")
    print_help(opt_parser)
    stop()
} else {
    pheno_file = opt$pheno_file
}

outdir = opt$outdir
interval = opt$interval
gwastype = opt$gwastype
alleles_format = opt$alleles_format
ldr2file = opt$ldr2file
#genes_annot = opt$genes_annot
genes_map = opt$genes_map
LDTHRES = opt$ldthres

verbose = opt$verbose

if (verbose == TRUE) {
cat("Arguments:\n", file=stderr())
cat(paste("\t", markers, "\n"), file=stderr())
cat(paste("\t", markers_hapmap, "\n"), file=stderr())
cat(paste("\t", markers_values, "\n"), file=stderr())
cat(paste("\t", excap_hapmap, "\n"), file=stderr())
cat(paste("\t", excap_values, "\n"), file=stderr())
cat(paste("\t", pheno_file, "\n"), file=stderr())
cat(paste("\t", outdir, "\n"), file=stderr())
cat(paste("\t", interval, "\n"), file=stderr())
cat(paste("\t", gwastype, "\n"), file=stderr())
cat(paste("\t", alleles_format, "\n"), file=stderr())
cat(paste("\t", ldr2file, "\n"), file=stderr())
#cat(paste("\t", genes_annot, "\n"), file=stderr())
cat(paste("\t", genes_map, "\n"), file=stderr())
cat("\n", file=stderr())
}

########################################### Functions
###########################################

################# Functions used to create a table of markers-genes
#################

# A function to retrieve the genes hit by each marker
genes_marker <- function(x, markers_pos, genes_pos){

    # data of the current marker
    marker <- markers_pos[x,]
    
    marker_chrom = unfactor(marker$c)
    
    # genes hit by the marker
    # (marker$c is the chromosome)
    marker_genes <- genes_pos[genes_pos$chrom == marker_chrom &
                                genes_pos$start <= marker$pos &
                                genes_pos$end >= marker$pos,]
    
    # If the marker hits NO genes, create a row with empty gene information
    if (nrow(marker_genes)==0){
        marker_genes = data.frame(matrix(nrow = 1, ncol = 8))
        colnames(marker_genes) <- c("markerID", "markerChrom", "markerPos", "markerLogP", "gene", "start", "end", "desc")
        marker_genes[1,] = c(marker$rs, marker_chrom, marker$pos, marker$var, "-", "-", "-", "-")
    
    # If the marker hits some gene, join the gene and marker information
    } else {
        marker_genes <- marker_genes %>%
                                mutate(markerID = marker$rs,
                                        markerChrom = marker_chrom,
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
    
    cat("Retrieving genes hit by each marker...\n", file=stderr())
    
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
    
    genes_pos <- marker_genes_map_data %>% select(gene, chrom, start, end)
    
    print(head(genes_pos), file=stderr())
    print(head(genes_annot_data), file=stderr())
    
    genes_pos <- genes_pos %>%
                left_join(genes_annot_data, by = "gene")
    
    print(head(genes_pos, 20), file=stderr())
    
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
    
    cat("Table of markers and genes created.\n", file=stderr())
    
    return();
}

### Auxiliar function to convert ATGC genotypes to numeric ones
numerical_genotype <- function(genotype){
  #print(genotype, file=stderr())
  nts <- unique(genotype)
  
  #print(nts, file=stderr())
  
  if (alleles_format == "biallelic") {
        genotype[genotype %in% "NN"] = -1
        for (nt in nts){
            #print(nt, file=stderr())
            
            unnt = unique(unlist(strsplit(nt, "")[[1]]))
            #print(unnt, file=stderr())
            
            if (nchar(unnt)>1) {
              genotype[genotype == nt] = -1
            } else {
              genotype[genotype == nt] = which(nts == nt) - 1
            }
        }
  } else {
        nts <- nts[nts %in% c("A", "C", "G", "T")]
        
        genotype[!genotype %in% nts] = -1 # for all which are not in ACGT
        # For those in ACGT:
        for (nt in nts){
          genotype[genotype == nt] = which(nts == nt) - 1
        }
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
                                pheno_data, ldr2_data){
  
  cat("Generating graphical genotypes...\n", file=stderr())
  
  #print(head(ldr2_data), file=stderr())
  ldr2_data_filtered <- ldr2_data %>% filter(r2 >= LDTHRES)
  
  #print(head(ldr2_data_filtered), file=stderr())
  
  minpos = min(ldr2_data_filtered$pos)
  maxpos = max(ldr2_data_filtered$pos)
  
  ## Keep only those columns to be used in the graphical genotypes
  marker_cols <- colnames(marker_hapmap_data)
  excap_cols <- colnames(excap_hapmap_data)
  common_cols <- intersect(marker_cols, excap_cols)
  #print(common_cols)
  
  genotypes <- rbind(select(marker_hapmap_data, common_cols), 
                     select(other_markers_data, common_cols),
                     select(excap_hapmap_data, common_cols))
  
  genotypes <- select(genotypes, -alleles, -c, -strand, -assembly,
                                 -center, -protLSID, -assayLSID, -panelLSID,
                                 -QCCode) %>%
    arrange(pos)

    genotypes <- genotypes %>% filter(pos >= minpos) %>% filter(pos <= maxpos)
  
  rownames(genotypes) <- unlist(select(genotypes, rs))
  genotypes <- select(genotypes, -rs, -pos)
  #print(head(genotypes))
  
  # convert to numerical genotypes
  genotypes <- numerical_genotypes(genotypes)
  #print(head(genotypes))
  
  all <- genotypes
  # Join both phenotypic and genotypic data in a single table
  # all <- rbind(genotypes, phenotypes)
  
  ## Obtain phenotypic data
  entries <- colnames(genotypes)
  phenotypes <- filter(pheno_data, Genotype %in% entries)
  
  all <- rbind(all, t(phenotypes))
  ## Write the table to a file
  outfile=paste(outdir, "/", marker, ".genotypes.tab", sep="")
  write.table(all, outfile, quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
  
  cat("Table of graphical genotypes was written.\n", file=stderr());
  
  gmat <- apply(genotypes, 2, as.numeric)
  rownames(gmat) <- rownames(genotypes)
  #print(head(gmat), file=stderr())
  #gmat = data.matrix(gmat)
  
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
  
  cat(paste("Plot of graphical genotypes saved as ", outfile, "\n", sep=""), file=stderr());
  
}

#############################################
## Plot
generate_plots <- function(marker, plottable, ldr2_data, value_colname){
  
  cat("Generating plots...\n", file=stderr())
  
  #print(head(plottable), file=stderr())
  
  ## Main plot function: data
  bfplot <- ggplot(plottable, aes(x=as.numeric(mbpos), y=as.numeric(var),  
                                  color=type, shape=type, linetype=type, size=type))+
    geom_point()+
    ylim(0, max(as.numeric(plottable$var)))#+
    #theme_bw()
  
  ## Adding LD r2
  ldr2_data$r2 <- ldr2_data$r2 * max(as.numeric(plottable$var))
  #print(head(ldr2_data), file=stderr())
  
  bfplot = bfplot +
        geom_line(data = ldr2_data,
                aes(x=mbpos, y=r2, color=type, linetype=type, size=type))+
        geom_point(data = ldr2_data,
                aes(x=mbpos, y=r2), color="black", shape=16, size=0.2, alpha=0.8)

  # To plot LD as a smoothed curve instead:
  #    geom_smooth(data = ldr2_data,
  #        aes(x=mbpos, y=r2, color=type, linetype=type, size=type, se=FALSE, level=1.0))+
  
  # secondary axis (in fact, it is creating a new primary axis and leaving the previous as secondary)
  bfplot = bfplot + scale_y_continuous(sec.axis = sec_axis(~./max(as.numeric(plottable$var)),
                                        name = "LD", breaks=c(0,1)))
  
  ## Scale, legend and labels
  bfplot = bfplot +
    scale_color_manual("",
                        values = c(alpha("red", 0.8), # GWAS
                                      alpha("black", 1.0), # gene
                                      alpha("black", 0.8), # LDr2
                                      alpha("blue", 1.0), # other_markers
                                      alpha("blue", 1.0)), # main marker
                       breaks=c("B_GWAS", "gene", "LDr2", "other_markers", paste("Z",marker, sep="")),
                       labels=c("Markers", "Genes", "LD", "Others", marker))+
    scale_shape_manual("", values=c(16, 2, 16, 16, 18), guide=FALSE)+
    scale_linetype_manual("", values=c(1,1,1,1,1), guide=FALSE)+
    scale_size_manual("", values=c(2.0, 2.5, 0.3, 3.0, 5.0), guide=FALSE)+
    guides(color=guide_legend(title=NULL))
    
  bfplot = bfplot +
    theme_bw()+
    labs(x = "Position in Mbp", y = value_colname, title = "")+
    theme(plot.title = element_text(hjust = 0.5),
            legend.background = element_rect(fill="gray90"),
            legend.position = "bottom",
            legend.text = element_text(size=15),
            axis.text.x = element_text(size=15),
            axis.text.y = element_text(size=15),
            axis.title=element_text(size=14,face="bold"))
            
    #theme(text = element_text(size=20),
    #    axis.text.x = element_text(angle=90, hjust=1)) 

  gt <- bfplot
  
  # As TIFF
  outfile=paste(outdir, "/", marker, ".tiff", sep="")
  tiff(outfile, height = 7, width = 7, units = "in", res = 300)
  grid.arrange(gt, nrow=1, ncol=1, newpage=FALSE)
  
  dev.off()
  cat(paste("Plot saved as ", outfile, "\n", sep=""), file=stderr())
  
  # As table
  outfile=paste(outdir, "/", marker, ".tsv", sep="")
  write.table(plottable, file=outfile, sep="\t")
  cat(paste("Table saved as ", outfile, "\n", sep=""), file=stderr())
  
  # The same plot but with smooth loess adjusted function
  
  gt_smooth <- bfplot+
    geom_smooth(data=subset(plottable,type=="B_GWAS"),
                aes(x=as.numeric(mbpos), y=as.numeric(var), 
                    color=type, linetype=type, size=type),
                method=loess, size=.5, se=TRUE, level=0.99)

  outfile=paste(outdir, "/", marker, ".smooth.tiff", sep="")
  tiff(outfile, height = 7, width = 7, units = "in", res = 300)
  grid.arrange(gt_smooth, nrow=1, ncol=1, newpage=FALSE)
  
  dev.off()
  cat(paste("Plot saved as ", outfile, "\n", sep=""), file=stderr())
}

f_log <- function(x){
    #return(-log10(x)+1)
    return(-log10(x))
}

################################################
### Main function for each single marker
process_marker <- function(x){
  
  marker <- as.character(unlist(x))
  cat("\n**** Processing marker:", marker, "\n\n", file=stderr())
  
  ### Obtain pvalue of marker
  #print(head(markers_values_data), file=stderr())
  marker_pvalue <- filter(markers_values_data, SNP == marker)[,2]
  #print(paste("Previous P.value:", marker_pvalue), file=stderr())
  
  # this is the association value column name
  value_colname = colnames(markers_values_data)[[2]]
  
  ### Obtain position of the marker
  marker_hapmap_data <- filter(markers_hapmap_data, rs == marker)
  #print(marker_hapmap_data, file=stderr())
  marker_chrom = marker_hapmap_data[1,3]
  library(varhandle)
  marker_chrom <- unfactor(marker_chrom)
  marker_pos = marker_hapmap_data[1,4]
  cat(paste("Chr:", marker_chrom, ", Pos:", marker_pos, "\n"), file=stderr())
  
  ### Interval
  if (as.numeric(interval) > 0) {
    pos_start = as.numeric(marker_pos) - as.numeric(interval)
    if (pos_start < 0) pos_start = 0;
    pos_end = as.numeric(marker_pos) + as.numeric(interval)
    
  } else {
    other_markers_hapmap_data <- filter(markers_hapmap_data, rs %in% markers_values_data$SNP &
                                                             c == marker_chrom)
    
    marker_excap_hapmap_data <- filter(excap_hapmap_data, rs %in% excap_values_data$SNP &
                                                          c == marker_chrom)
    
    pos_start_other = min(other_markers_hapmap_data$pos)
    pos_start_excap = min(marker_excap_hapmap_data$pos)
    pos_start = min(c(marker_pos, pos_start_other, pos_start_excap))
    
    pos_end_other = max(other_markers_hapmap_data$pos)
    pos_end_excap = max(marker_excap_hapmap_data$pos)
    pos_end = max(c(marker_pos, pos_end_other, pos_end_excap))
  }
  
  cat(paste("Interval:", pos_start, "-", pos_end, "\n"), file=stderr())
  
  ### Obtain other markers from the original set in the interval
  # defined by the position of the original marker plus the interval parameter
  other_markers_hapmap_data <- filter(markers_hapmap_data, c == marker_chrom &
                                                            pos >= pos_start &
                                                            pos <= pos_end &
                                                            rs != marker)
  
  other_markers_data <- select(other_markers_hapmap_data, rs, pos) %>%
                              mutate(mbpos = pos / 1000000) %>%
                              select(-pos) %>%
                              mutate(type = "other_markers")
    
  if (verbose == TRUE) { print(paste("Num other original markers found:", nrow(other_markers_data)), file=stderr()) }
  
  # obtain their pvalues
  other_markers_data = left_join(other_markers_data, markers_values_data, by = c("rs" = "SNP"))
  
  if (gwastype == "pvalue"){ # -log10 of the 'value_colname' column
  
    other_markers_data <- mutate(other_markers_data,
                                var = ifelse(is.na(other_markers_data[,value_colname]), 0,
                                            -f_log(other_markers_data[,value_colname])))
    
    other_markers_data <- select(other_markers_data, -!!value_colname, marker_ID = rs)
    
    
    #other_markers_data <- mutate(other_markers_data, var = ifelse(is.na(P.value), 0, f_log(P.value))) %>%
    #                    select(-P.value, marker_ID = rs)
        
  } else { # the value of the 'value_colname' column
    other_markers_data <- mutate(other_markers_data,
                                var = ifelse(is.na(other_markers_data[,value_colname]), 0,
                                            other_markers_data[,value_colname]))
    
    other_markers_data <- select(other_markers_data, -!!value_colname, marker_ID = rs)
  }
  
  ### Obtain excap markers in the interval
  # defined by the position of the original marker plus the interval parameter
  
  marker_excap_hapmap_data <- filter(excap_hapmap_data, c == marker_chrom &
                                                        pos >= pos_start &
                                                        pos <= pos_end)
  
  if (verbose == TRUE) { print(paste("Num excap markers found:", nrow(marker_excap_hapmap_data)), file=stderr()) }
  #print(t(head(t(head(marker_excap_hapmap_data)))), file=stderr())
  excap_markers_list = select(marker_excap_hapmap_data, rs)[,1]
  #print(head(excap_markers_list), file=stderr())
  
  # Obtain GWAS values of those markers
  marker_excap_values_data <- filter(excap_values_data, SNP %in% excap_markers_list)
  if (verbose == TRUE) { print(paste("Num GWAS results found:", nrow(marker_excap_values_data)), file=stderr()) }
  #print(head(marker_excap_values_data), file=stderr())
  
  marker_excap_values_no_data <- excap_markers_list[!excap_markers_list %in% select(marker_excap_values_data, SNP)[,1]]
  #print(head(marker_excap_values_no_data), file=stderr())
  if (verbose == TRUE) { print(paste("Num GWAS results NOT found:", length(marker_excap_values_no_data)), file=stderr()) }
  
  # Assign 0 to value_colname column of not found markers
  marker_excap_values_no_data = data.frame(marker_excap_values_no_data)
  colnames(marker_excap_values_no_data) = c("SNP")
  marker_excap_values_no_data <- mutate(marker_excap_values_no_data, !!value_colname := NA)
  #print(head(data.frame(marker_excap_values_no_data)), file=stderr())
  
  #### Join all the data of excap markers
  marker_excap_values_data <- rbind(marker_excap_values_data, marker_excap_values_no_data)
  
  if (gwastype == "pvalue"){  # -log10 of the 'value_colname' column
    #marker_excap_values_data <- mutate(marker_excap_values_data, log10P = ifelse(is.na(P.value), 0, f_log(P.value)))
    
    marker_excap_values_data <- mutate(marker_excap_values_data,
                                var = ifelse(is.na(marker_excap_values_data[,value_colname]), 0,
                                            -f_log(marker_excap_values_data[,value_colname])))
    
    #marker_excap_values_data <- select(marker_excap_values_data, value_colname, marker_ID = rs)
    
  } else { # the value of the 'value_colname' column
    #marker_excap_values_data <- mutate(marker_excap_values_data, log10P = ifelse(is.na(P.value), 0, P.value))
    
    marker_excap_values_data <- mutate(marker_excap_values_data,
                                var = ifelse(is.na(marker_excap_values_data[,value_colname]), 0,
                                            marker_excap_values_data[,value_colname]))
    
    #marker_excap_values_data <- select(marker_excap_values_data, value_colname, marker_ID = SNP)
    
  }
  
  # Join association values with map positions
  ### map positions --> marker_excap_hapmap_data
  ### association values --> marker_excap_values_data
    
  marker_excap_data <- inner_join(marker_excap_hapmap_data, marker_excap_values_data, by = c("rs" = "SNP")) %>%
    mutate(mbpos = as.numeric(as.character(pos))/1000000) 
    
  # the last line is to reduce the magnitude of map positions (so that the value can be plotted easily)
  
  # Prepare for plot
  
  marker_excap_data <- select(marker_excap_data, marker_ID = rs, mbpos = mbpos, var = var) %>% 
    mutate(type = "B_GWAS")
  # The "B" was used in the original script to make sure that the series with excap markers
  # was the second after the dummy series, so that the plot colors, etc
  # were applied correctly
  
  ### Therefore, marker_excap_data is the table to plot markers and pvalues
  
  ### Genes data
  
  # Obtain genes in the interval
  if (verbose == TRUE) { print(paste("Num total genes", nrow(genes_map_data)), file=stderr()) }
  
  marker_genes_map_data <- filter(genes_map_data, chrom == marker_chrom &
                                       ((start >= pos_start & start <= pos_end) |
                                       (end >= pos_start & end <= pos_end) | 
                                         (start <= pos_start & end >= pos_end)))
  
  if (verbose == TRUE) { print(paste("Num genes in interval:", nrow(marker_genes_map_data)), file=stderr()) }
  
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
  
  if (gwastype == "pvalue"){
    marker_excap_data <- rbind(marker_excap_data, c(marker_ID = marker, mbpos = marker_pos/1000000, 
                                                  var = f_log(marker_pvalue),
                                                  type = paste("Z", marker, sep="")))
  } else {
    marker_excap_data <- rbind(marker_excap_data, c(marker_ID = marker, mbpos = marker_pos/1000000, 
                                                  var = marker_pvalue,
                                                  type = paste("Z", marker, sep="")))
  }
  
  # The "Z" is to make sure that the marker is the last series
  # (see above "The B...")
  #print(head(marker_excap_data), file=stderr())
  #print(head(other_markers_data), file=stderr())
  
  ## Prepare LD r2 data
  
  ldr2_data <- mutate(ldr2_data, mbpos = pos / 1000000) %>% mutate(type = "LDr2")
  
  ### Join all the tables
  
  plottable <- marker_excap_data %>% rbind(other_markers_data) %>% rbind(genestable) %>% arrange(type)
  
  #print(plottable, file=stderr())
  
  generate_plots(marker, plottable, ldr2_data, value_colname)
  
  graphical_genotypes(marker, marker_hapmap_data, other_markers_hapmap_data, marker_excap_hapmap_data,
                      pheno_data, ldr2_data)

    # removed because it was not working properly
    #genes_table(marker, marker_hapmap_data, other_markers_hapmap_data, marker_excap_hapmap_data,
    #            marker_genes_map_data, genes_annot_data, plottable)
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
if (verbose == TRUE) { print(t(head(t(head(markers_hapmap_data, 3)))), file=stderr()) }

# Reading markers association values
cat("Reading primary markers association values...\n", file=stderr())
markers_values_data <- read.table(markers_values, header = TRUE)
if (verbose == TRUE) { print(head(markers_values_data, 3), file=stderr()) }

# Reading excap data
cat("Reading secondary markers hapmap data...\n", file=stderr())
excap_hapmap_data <- read.table(excap_hapmap, header = TRUE, sep="\t")
if (verbose == TRUE) { print(t(head(t(head(excap_hapmap_data, 3)))), file=stderr()) }

# Reading excap association values
cat("Reading secondary markers association values...\n", file=stderr())
excap_values_data <- read.table(excap_values, header = TRUE)
if (verbose == TRUE) { print(head(excap_values_data, 3), file=stderr()) }

# Reading genes map
cat("Reading genes map...\n", file=stderr())
genes_map_data <- read.table(genes_map, header = FALSE)
colnames(genes_map_data) <- c("chrom", "start", "end", "gene")
if (verbose == TRUE) { print(head(genes_map_data, 3), file=stderr()) }

## Removed because this was not working properly
# Reading genes annotation
#cat("Reading genes annotation...\n", file=stderr())
#genes_annot_data <- read.table(genes_annot, header = FALSE, sep = "\t")
#colnames(genes_annot_data) <- c("gene", "desc")
#if (verbose == TRUE) { print(head(genes_annot_data, 3), file=stderr()) }

# Reading phenotypic data
cat("Reading phenotypic data...\n", file=stderr())
pheno_data <- read.table(pheno_file, header = TRUE, sep = "\t")
colnames(pheno_data) <- c("Genotype", "Phenotype")
if (verbose == TRUE) { print(head(pheno_data, 3), file=stderr()) }

# Reading phenotypic data
cat("Reading LD r2 values...\n", file=stderr())
ldr2_data <- read.table(ldr2file, header = TRUE, sep = "\t")
if (verbose == TRUE) { print(head(ldr2_data, 3), file=stderr()) }

# Process each marker
devnull <- lapply(markers_list[,1], process_marker)
cat("\nFinished.\n************************\n\n", file=stderr())
q()

## END
