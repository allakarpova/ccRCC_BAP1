#Modified by Alla from Ilya Strunilin's script
### RCC_ATAC_PBRM1_pseudobulk_DEG_calculate1
### 
# Clear plots
if(!is.null(dev.list())) dev.off()
# Clear console
cat("\014") 
# Clean workspace
rm(list=ls())
# Clean memory
gc(verbose = FALSE)

# Set R options

options(error = NULL)
options(digits = 3)
options(scipen = 3)

## Load packages

my0packages <- c(
   "tidyverse",
   "data.table",
   "Seurat",
    'Signac',
   "RColorBrewer",
   "patchwork",
   "circlize",
   "ggrepel",
       "doParallel",
   "Matrix",
   "optparse"
)

for (pkg_name_tmp in my0packages) {
   if (!(pkg_name_tmp %in% installed.packages()[,1])) {
      message(pkg_name_tmp," package not installed")
   }
   library(package = pkg_name_tmp, character.only = T, warn.conflicts = FALSE)
}

### Functions live here ----

dfread <- function(l.path) {
   fread(l.path, data.table = FALSE)
}

rename <- dplyr::rename
count <- dplyr::count
select <- dplyr::select
slice <- dplyr::slice

usort <- function(l.vec) {
  l.vec <- sort(unique(l.vec),na.last = TRUE)
  return(l.vec)
}

nafwrite <- function(l.table,l.path, ...) {
  fwrite(x = l.table,file = l.path,na = NA, quote = FALSE, sep = "\t",...)
}

meltdf <- function(l.df, ...) {
  if (is.data.frame(l.df)) {
    l.df <- as.data.frame(data.table::melt(data = as.data.table(l.df), variable.factor = FALSE, ...))
    return(l.df)
  } else {stop("Not a data frame")}
}

dcastdf <- function(l.df, ...) {
  if (is.data.frame(l.df)) {
    l.df <- as.data.frame(data.table::dcast(data = as.data.table(l.df),...))
    return(l.df)
  } else {stop("Not a data frame")}
}

tabla <- function(...) {
  call <- list(...)
  custom_args <- list(useNA = "ifany",deparse.level = 1) # could extend this list for more customization
  non_overlap_args <- setdiff(names(custom_args),names(call))
  call <- c(call, custom_args[non_overlap_args])
  do.call(table, call) # exectue table() with the custom settings
}

pct.positive.helper <- function(l.obj,l.genes) {
  l.genes <- unique(l.genes)
  l.dat = l.obj[['X500peaksMACS2']]@data
  l.dat <- l.dat[l.genes, ,drop = FALSE]
  l.ncells = ncol(l.dat)
  l.dat@x <- rep(1,length(l.dat@x))
  l.gene_pcr <- data.frame(pct = Matrix::rowSums(l.dat)/l.ncells*100, Gene = rownames(l.dat)) 
  return(l.gene_pcr)
}

pct.positive <- function(l.obj,l.genes,l.groupby) {
  #
  l.obj_list = SplitObject(l.obj, split.by = l.groupby)
  l.results <- lapply(l.obj_list,pct.positive.helper, l.genes = l.genes)
  for (i in names(l.results)) {l.results[[i]]$Group <- i}
  l.results <- bind_rows(l.results)
  l.results_wide <- l.results %>% dcastdf(Gene ~ Group, value.var = 'pct')
  colnames(l.results_wide)[2:3] = c('pct.1', 'pct.2')
  return(l.results_wide)
}


###options###
######################
option_list = list(
  make_option(c("-i", "--input.object"),
              type="character",
              default="/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/data_freeze/V1/snATAC/RCC_snATAC_with_chromvar.rds", 
              help="path to merged object",
              metavar="character"),
  make_option(c("-o", "--output"),
              type="character",
              default="./", 
              help="output folder path",
              metavar="character"),
  make_option(c("-e", "--extra"),
              type="character",
              default="./", 
              help="add unique string identifier for your data",
              metavar="character"),
  make_option(c("-m","--metadata.file"),
              type="character",
              default="/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/data_freeze/V1/snATAC/snATAC_celltypes.tsv",
              help = "path to metadats file with cell types, make cell barcodes in the 1st column",
              metavar="character"),
  make_option(c("-c","--clinical"),
              type="character",
              default='/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/data_freeze/V1/clinical/RCC_2022_clinical_v13_2024.10.15.tsv',
              help = "path to clinical table",
              metavar="character"),
  make_option(c("--mutation"),
              type="character",
              default='BAP1',
              help = "BAP1 or PBRM1",
              metavar="character"),
  make_option(c("--min.pct"),
              type="double",
              default=0.1,
              help = "% of cells expressing a feature",
              metavar="numeric")
  
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

# read in initial arguments
input.path <- opt$input.object
out_path <- opt$output
add_filename <- opt$extra
meta.path <- opt$metadata.file
clinic <- opt$clinical
mut <- opt$mutation
min.pct.cutoff <- opt$min.pct

dir.create(out_path, showWarnings = F)
setwd(out_path)
dir.create('out')


### -----



# load clinical data
clinical <- dfread(clinic)
clinical <- clinical %>% filter(hist_type %in% c("ccRCC","NAT ccRCC"))

# load hgnc
hgnc <- dfread("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/input/ccRCC/hgnc_complete_set.txt")
protein_coding <- hgnc %>% filter(locus_group == "protein_coding_gene") %>% pull(symbol)


annotated_peaks = dfread(paste0('/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/data_freeze/V1/snATAC/PanImmune_peaks_df_v8.0_EnsDb.Hsapiens.v100.annot.tsv'))


### Load files ----
object = read_rds(input.path)

celltypes <- dfread(meta.path)


# load TME colors
tme_colors <- read_rds("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/input/ccRCC/RCC_TME_colors_v1.rds")



## adjust colors
#tme_colors$Group <- tme_colors$PBRM1_mut

celltypes %>% head

## add meta
addmeta <- celltypes %>% 
  left_join(clinical[c('Aliquot.snRNA','Case_comparison')], by = c('orig.ident' = 'Aliquot.snRNA'))
rownames(addmeta) <- addmeta$merged_barcode
object <- AddMetaData(object,addmeta)

#subset only cancer cells
object <- subset(object, celltype_final_short	 == 'ccRCC cc')
#subset out HT293 and HT282
object$meta.data$HTAN.case <- object$meta.data$Case %in% c('HT293H1', 'HT282H1')
object <- subset(object, HTAN.case, invert = TRUE)


#save.signif.deg <- function(l.deg_df,l.log2FC_df,l.coverage_df, l.sex) {
    
    ## convert NA pvals to 1
    l.deg_df$pval[is.na(l.deg_df$pval)] = 1
    
    ## re-calculate p.adjust with fdr
    l.deg_df$fdr = p.adjust(l.deg_df$pval, method = 'fdr')
    
    ## format
    l.deg_df = l.deg_df %>% rename(Peak = Gene)
    #print(head(l.deg_df))
    l.log2FC_df = l.log2FC_df %>% rename(Peak = Gene)
    #print(head(l.deg_df))
    l.coverage_df = l.coverage_df %>% rename(Peak = Gene)


    ## add peak annotation
    l.deg_df = l.deg_df %>% 
        left_join(annotated_peaks[c('new_peak','SYMBOL','transcriptBiotype')], by = c('Peak' = 'new_peak'))
    
  ## add correct log2FC
  l.deg_df <- l.deg_df %>% select(!log2FC) %>% inner_join(l.log2FC_df[c('Peak','log2FC')], by = 'Peak') %>% 
    relocate(log2FC, .after = Peak) %>% select(!c(pct.1,pct.2)) %>% inner_join(l.coverage_df, by = 'Peak')
  
  #cap pval
  l.deg_df$pval <- ifelse(l.deg_df$pval < 10^-300, 10^-300,l.deg_df$pval)
  l.deg_df$fdr <- ifelse(l.deg_df$fdr < 10^-300, 10^-300,l.deg_df$fdr)
  
  # save almost all
  l.deg_df_format <- l.deg_df %>% arrange(pval)
  l.deg_df_format[colnames(l.deg_df_format) %in% c("pval","fdr")] <- format(l.deg_df_format[colnames(l.deg_df_format) %in% c("pval","fdr")], scientific = TRUE, digits = 2)
  l.deg_df_format[colnames(l.deg_df_format) %in% c("log2FC")] <- round(l.deg_df_format[colnames(l.deg_df_format) %in% c("log2FC")], digits = 3)
  #save
  nafwrite(l.deg_df_format, paste0(l.sex,"_ALL_DEG_.tsv"))
  
  
  for (l.log2FC_cutoff in c(0,1, 1.5, 2)) {
    
    for (l.pct_cutoff in c( 0, 5, 10)) {
      # browser()
      # l.abun_Peaks <- l.coverage_df %>% filter() %>% pull(Peak)
      
      l.signif_df <- l.deg_df %>%  
        filter(abs(log2FC) >= l.log2FC_cutoff & (pct.1 >= l.pct_cutoff | pct.2 >= l.pct_cutoff) ) %>%
        filter(fdr < 0.05) 
      
      if (nrow(l.signif_df) > 0) {
        # print("here")
        l.deg_df_format <- l.signif_df
        l.deg_df_format[colnames(l.deg_df_format) %in% c("pval","fdr")] <- format(l.deg_df_format[colnames(l.deg_df_format) %in% c("pval","fdr")], scientific = TRUE, digits = 2)
        l.deg_df_format[colnames(l.deg_df_format) %in% c("log2FC")] <- round(l.deg_df_format[colnames(l.deg_df_format) %in% c("log2FC")], digits = 3)
        l.deg_df_format[colnames(l.deg_df_format) %in% c("pct.1",'pct.2')] <- round(l.deg_df_format[colnames(l.deg_df_format) %in% c("pct.1",'pct.2')], digits = 1)
        #save
        nafwrite(l.deg_df_format, paste0(l.sex,"_log2FC_",l.log2FC_cutoff,"_pct_",l.pct_cutoff,"_.tsv"))
        
      }
      
    }
    
  }
  
}


# for (N in c(1000)){
#     open_peaks = AccessiblePeaks(
#       sct_subset,
#       assay = 'X500peaksMACS2',
#       min.cells = N
#     )
#     message(N, 'cells')
#     message(length(open_peaks), 'peaks')
# }

cnv.table <- readRDS('/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/Alla_help/match_peaks_with_CNV_by_genes/Case2Peak.CNV.20241018.v1.RDS')
colnames(cnv.table)=gsub('-','_',colnames(cnv.table))

meta.cc <- object@meta.data %>% 
    filter(celltype_final_short	 == 'ccRCC cc') %>%
    select(Case, sample)


###Now check how FindMarkers modules will work:

## set ident.use.1
ident.use.1=ifelse(mut=='BAP1', "BAP1_MUT", "PBRM1_MUT")
## set ident.use.2
ident.use.2='WTWT'
## make sure set Idents for the object
Idents(object)=object$Genotype_group

##assay
assay='X500peaksMACS2'

## input the CNV values per peak per cell
### column is named after the peak by rownames(object)
### row is named after the cells by colnames(object)
cnv_per_feature_df=barcode

# set pre-filtering parameters --------------------------------------------
## set min.pct
min.pct=min.pct.cutoff
## set min.diff.pct (lower to 0, since it filtered a lot of peaks)
min.diff.pct=0
## set logfc.threshold
logfc.threshold=0

# feature selection -------------------------------------------------------
fc.results <- FoldChange(
  object = object,
  slot = "data",
  ident.1 = ident.use.1,
  ident.2 = ident.use.2,
  base = 2
)

# feature selection (based on percentages)
alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
names(x = alpha.min) <- rownames(x = fc.results)
features <- names(x = which(x = alpha.min >= min.pct))
if (length(x = features) == 0) {
  warning("No features pass min.pct threshold; returning empty data.frame")
  return(fc.results[features, ])
}
alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
features <- names(
  x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
)
if (length(x = features) == 0) {
  warning("No features pass min.diff.pct threshold; returning empty data.frame")
  return(fc.results[features, ])
}
# feature selection (based on logFC)
#slot
slot='data'
only.pos=FALSE

if (slot != "scale.data") {
  total.diff <- fc.results[, 1] #first column is logFC
  names(total.diff) <- rownames(fc.results)
  features.diff <- if (only.pos) {
    names(x = which(x = total.diff >= logfc.threshold))
  } else {
    names(x = which(x = abs(x = total.diff) >= logfc.threshold))
  }
  features <- intersect(x = features, y = features.diff)
  if (length(x = features) == 0) {
    warning("No features pass logfc.threshold threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
}

filtered.peaks=fc.results[rownames(fc.results) %in% features,]
filtered.peaks$peak=rownames(filtered.peaks)
write.table(filtered.peaks,glue::glue('out/{mut}_comparison_Filtered_peaks_byMinPct{min.pct.cutoff}_MinPctDiff.tsv'),sep='\t',
            row.names=F,quote=F)

## set "features" -- we do it later
features=colnames(barcode)


# set inputs for the below process --------------------------------------------------------------
cells.1 <- WhichCells(object = object, idents = ident.use.1)
cells.2 <- WhichCells(object = object, idents = ident.use.2)
data.use=object[[assay]]
data.use <- data.use[features, c(cells.1, cells.2)]

# process inputs further --------------------------------------------------
## prepare latent.vars data frame
latent.vars <- FetchData(
  object = object,
  vars = "peak_RF_500MACS2",
  cells = c(cells.1, cells.2)
)
## prepare group.info data frame
group.info <- data.frame(row.names = c(cells.1, cells.2))
group.info[cells.1, "group"] <- "Group1"
group.info[cells.2, "group"] <- "Group2"
group.info[, "group"] <- factor(x = group.info[, "group"])

## prepare data.use object
data.use <- data.use[, rownames(group.info), drop = FALSE]
rownames(data.use)=gsub('-','_',rownames(data.use))
latent.vars <- latent.vars[rownames(group.info), , drop = FALSE]
#colnames(latent.vars.full)=gsub('-','_',colnames(latent.vars.full))



  
# run test ----------------------------------------------------------------
my.sapply <- ifelse(
  nbrOfWorkers() == 1,
  #  test = verbose && nbrOfWorkers() == 1,
  yes = pbsapply,
  no = future_sapply
)
p_val <- my.sapply(
  X = rownames(x = data.use),
  FUN = function(x) {
    cnv_per_feature_df <- meta.cc %>% 
      rownames_to_column('R') %>% 
      left_join(cnv.table[,c('Case', x)], by = 'Case') %>%
      column_to_rownames('R') %>%
      select(all_of(x))
    
    latent.vars.full <- cbind(latent.vars, cnv_per_feature_df[rownames(group.info),])
    
    model.data <- cbind(GENE = data.use[x,], group.info, latent.vars.full)
    fmla <- as.formula(object = paste(
      "group ~ GENE +",
      paste(c(x, 'peak_RF_500MACS2'), collapse = "+")
    ))
    fmla2 <- as.formula(object = paste(
      "group ~",
      paste(c(x, 'peak_RF_500MACS2'), collapse = "+")
    ))
    model1 <- glm(formula = fmla, data = model.data, family = "binomial")
    model2 <- glm(formula = fmla2, data = model.data, family = "binomial")
    lrtest <- lrtest(model1, model2)
    return(lrtest$Pr[2])
  }
)
to.return <- data.frame(p_val, row.names = rownames(data.use))
to.return$p_val=as.numeric(as.character(unlist(to.return$p_val)))
to.return$FDR=p.adjust(to.return$p_val,method='fdr')
to.return$p_adjust_bonf=p.adjust(to.return$p_val,method='bonferroni')

###Now merge with fc.results
to.return$peak=row.names(to.return)
fc.results$peak=row.names(fc.results)
to.return$peak=gsub('_','-',to.return$peak)
to.return=merge(to.return,fc.results,all.x=TRUE)


to.return$chr_peak=gsub('(.*)-.*-.*','\\1',to.return$peak)

version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

write.table(to.return,glue:glue("out/DA_peaks_{mut}mutants_vs_nonMutants_correctedbyCNV.{run_id}.tsv"),
            sep='\t',quote=FALSE,row.names=F)


#
#save.signif.deg(pval_df,log2FC_df,coverage_df, "All")
  











