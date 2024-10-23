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
  "optparse",
  "future",
  "pbapply",
  "lmtest"
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
              metavar="numeric"), 
  make_option(c("--min.diff.pct"),
              type="double",
              default=0,
              help = "minimum % difference between 2 groups",
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
min.diff.pct.cutoff <- opt$min.diff.pct

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
  left_join(clinical[c('Aliquot.snRNA','Case_comparison', 'Case', 'Genotype_group')], by = c('orig.ident' = 'Aliquot.snRNA'))
rownames(addmeta) <- addmeta$merged_barcode
object <- AddMetaData(object,addmeta)

head(object@meta.data)

#subset only cancer cells
print(dim(object))
object <- subset(object, celltype_final_short=='ccRCC cc')
print(dim(object))
#subset out HT293 and HT282
object@meta.data$HTAN.case <- object@meta.data$Case %in% c('HT293N1', 'HT282N1')
object <- subset(object, HTAN.case, invert = TRUE)
print(dim(object))

#cnv.table <- readRDS('/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/Alla_help/match_peaks_with_CNV_by_genes/Case2Peak.CNV.20241018.v1.RDS')
#remove peaks with NA CNVs
#peaks.withNA <- apply(cnv.table,2, FUN = function(x) all(is.na(x)))
#cnv.table <- cnv.table[,colnames(cnv.table)[!peaks.withNA]]

meta.cc <- object@meta.data %>% 
  filter(celltype_final_short	 == 'ccRCC cc') %>%
  select(Case, Sample)


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
#cnv_per_feature_df=barcode

# set pre-filtering parameters --------------------------------------------
## set min.pct
min.pct=min.pct.cutoff
## set min.diff.pct (lower to 0, since it filtered a lot of peaks)
min.diff.pct=min.diff.pct.cutoff
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
write.table(filtered.peaks,glue::glue('out/{mut}_comparison_Filtered_peaks_byMinPct{min.pct.cutoff}_MinPctDiff{min.diff.pct.cutoff}.tsv'),sep='\t',
            row.names=F,quote=F)

## set "features" -- we do it later
#features=colnames(barcode)
#keep only peaks that have CNV information for them
#filtered.peaks <- filtered.peaks %>%
#  filter(peak %in% colnames(cnv.table))
features=filtered.peaks$peak

#colnames(cnv.table)=gsub('-','_',colnames(cnv.table))


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
#rownames(data.use)=gsub('-','_',rownames(data.use))
latent.vars <- latent.vars[rownames(group.info), , drop = FALSE]
#colnames(latent.vars.full)=gsub('-','_',colnames(latent.vars.full))




# run test ----------------------------------------------------------------
my.lapply <- ifelse(
  nbrOfWorkers() == 1,
  #  test = verbose && nbrOfWorkers() == 1,
  yes = pblapply,
  no = future_lapply
)
fc_p_val <- my.lapply(
  X = rownames(x = data.use),
  FUN = function(x) {
    
    model.data <- cbind(PEAK = data.use[x,], group.info, latent.vars)
    model.data <- cbind(model.data, PEAK.log2 = log2(exp(model.data[,'PEAK'])))
    model.data <- cbind(model.data, group.binary = ifelse(model.data[,'group']=='Group1',1,0))
    
    fmla <- as.formula(object = "group.binary ~ PEAK.log2 + peak_RF_500MACS2")
    fmla2 <- as.formula(object = "group.binary ~ peak_RF_500MACS2")
    
    model1 <- glm(formula = fmla, data = model.data, family = "binomial")
    model2 <- glm(formula = fmla2, data = model.data, family = "binomial")
    fc.glm <- coef(summary(model1))['PEAK.log2','Estimate']
    lrtest <- lrtest(model1, model2)
    pv <- lrtest$Pr[2]
    return(c(fc.glm, pv))
  }
)
mat <- do.call('rbind', fc_p_val)
colnames(mat) <- c('log2FC_glm', 'p_val')
to.return <- data.frame(mat, row.names = rownames(data.use))
to.return$p_val=as.numeric(as.character(unlist(to.return$p_val)))
to.return$FDR=p.adjust(to.return$p_val,method='fdr')
to.return$p_adjust_bonf=p.adjust(to.return$p_val,method='bonferroni')

###Now merge with fc.results
to.return$peak=row.names(to.return)
fc.results$peak=row.names(fc.results)
to.return$peak=gsub('_','-',to.return$peak)
to.return=merge(to.return,fc.results,all.x=TRUE)


to.return$chr_peak=gsub('(.*)-.*-.*','\\1',to.return$peak)
to.return <- to.return %>% arrange(p_val)

to.return <- to.return %>% left_join((annotated_peaks %>% select(new_peak, peak.position, SYMBOL)), by=c('peak'='new_peak'))


version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

write.table(to.return,glue::glue("out/DA_peaks_{mut}mutants_vs_nonMutants_NOTcorrectedbyCNV.{run_id}.tsv"),
            sep='\t',quote=FALSE,row.names=F)


#
#save.signif.deg(pval_df,log2FC_df,coverage_df, "All")












