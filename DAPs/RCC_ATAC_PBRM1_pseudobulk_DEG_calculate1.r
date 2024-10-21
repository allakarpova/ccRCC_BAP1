### RCC_ATAC_PBRM1_pseudobulk_DEG_calculate1
###
### v1
### - sn_PBRM1_WT_ids_clinical_v12_.tsv and sn_PBRM1_MUT_ids_clinical_v12_.tsv
### - combine CD8 T and EX CD8 T
### - combine macrophages
### - v12 clinical
### - uses RCC_final_celltype_transfer4 celltypes
### desc: For each celltype find DEGs between PBRM1 WT and MUT samples using downsampling. On a case-level.
### regress sex and grade with lm
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
   "Matrix"
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


### -----

obj_merge_id <- "ATACaug2024"
celltype_id <- "ATACaug2024"
res <- 0.5

# create output directory
output_dir <- paste0("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/output/RCC_ATAC_PBRM1_pseudobulk_DEG_calculate1/",
                     obj_merge_id,"_",celltype_id,"/")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# set working directory
setwd(output_dir)
normalizePath('../')

# load clinical data
clinical <- dfread("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/input/ccRCC/RCC_2022_clinical_v12_2024.08.30.tsv")
clinical <- clinical %>% filter(hist_type %in% c("ccRCC","NAT ccRCC"))
# load hgnc
hgnc <- dfread("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/input/ccRCC/hgnc_complete_set.txt")
protein_coding <- hgnc %>% filter(locus_group == "protein_coding_gene") %>% pull(symbol)
## load sample groups
group1_sample_ids <- dfread("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/input/ccRCC/sn_PBRM1_WT_ids_clinical_v12_.tsv") %>% 
pull()
group2_sample_ids <- dfread("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/input/ccRCC/sn_PBRM1_MUT_ids_clinical_v12_.tsv") %>% 
pull()
### -----

annotated_peaks = dfread(paste0('/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/output/RCC_ATAC_peak_annotation1/',obj_merge_id,'/PanImmune_peaks_df_v8.0_EnsDb.Hsapiens.v100.annot.tsv'))


### Load files ----
sct_raw = read_rds(list.files(paste0('/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/output/ATAC_merging_1/',obj_merge_id), 
                              pattern = 'final_.rds', full.names = TRUE))

celltypes <- dfread(paste0("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/output/RCC_ATAC_celltype_assignment_v1/",
                       celltype_id,"/",res,"_resolution/",celltype_id,"_addmeta.tsv"))

# load celltype matching table
celltype_matching <- dfread(paste0("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/output/RCC_final_celltype_transfer4/",
                       'feb2024_harmony1',"/","cell_type_name_matching.tsv"))
# load TME colors
tme_colors <- read_rds("/diskmnt/Projects/ccRCC_scratch/RCC_snRNA_2022/Results/input/ccRCC/RCC_TME_colors_v1.rds")

## adjust clinical
clinical <- clinical %>% 
    filter(Aliquot.snRNA %in% sct_raw$orig.ident) %>% 
    mutate(Case_comparison = case_when(Aliquot.snRNA %in% group1_sample_ids ~ paste0(Case,"_WT"),
                                       Aliquot.snRNA %in% group2_sample_ids ~ paste0(Case,"_MUT")))
# head(clinical2[c(1:6,ncol(clinical2))])
#
group1_sample_cases <- clinical %>% filter(Aliquot.snRNA %in% group1_sample_ids) %>% pull(Case_comparison) %>% 
    unique()
length(group1_sample_cases)
group2_sample_cases <- clinical %>% filter(Aliquot.snRNA %in% group2_sample_ids) %>% pull(Case_comparison) %>% 
    unique()
length(group2_sample_cases)

## Adjust clinical again
clinical <- clinical %>% 
  mutate(Group = case_when(Aliquot.snRNA %in% group1_sample_ids ~ "WT",
                           Aliquot.snRNA %in% group2_sample_ids ~ "MUT"
  )) %>% 
  filter(!is.na(Group))
#
clinical$Group <- ordered(clinical$Group,levels =(c("WT","MUT")))

#
case_clinical <- clinical %>% select(Case_comparison,Group) %>% unique()


## adjust colors
tme_colors$Group <- tme_colors$PBRM1_mut

celltypes %>% head

## add meta
addmeta <- celltypes %>% left_join(clinical[c('Aliquot.snRNA','Case_comparison')], by = c('orig.ident' = 'Aliquot.snRNA'))
rownames(addmeta) <- addmeta$merged_barcode
sct <- AddMetaData(sct_raw,addmeta)

## Keep only cancer cells
sct <- subset(sct, celltype_final == "ccRCC cancer cell" & hist_type == 'ccRCC')
sct
## subset
sct_subset <- subset(sct, orig.ident %in% c(group1_sample_ids,group2_sample_ids))


# ## combine cell types
# sct_subset$celltype_final <- ifelse(sct_subset$celltype_final %in% c("CD8 T","EX CD8 T"),"CD8 T",sct_subset$celltype_final)







save.signif.deg <- function(l.deg_df,l.log2FC_df,l.coverage_df, l.sex) {
    
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


  ## filter DEG by consistency
  pseudo_ifnb <- AggregateExpression(sct_subset, assays = "X500peaksMACS2", 
                                     return.seurat = T, group.by = c("Case_comparison"), slot = "counts")
#   pseudo_ifnb <- NormalizeData(pseudo_ifnb,normalization.method = "LogNormalize")
  pseudo_addmeta <- pseudo_ifnb@meta.data
  pseudo_addmeta$Case_comparison <- rownames(pseudo_addmeta)
  pseudo_addmeta <- pseudo_addmeta %>% left_join(case_clinical, by = "Case_comparison")
  rownames(pseudo_addmeta) <- pseudo_addmeta$Case_comparison
  pseudo_ifnb <- AddMetaData(pseudo_ifnb, pseudo_addmeta)
  
  

  coverage_df <- pct.positive(sct_subset,rownames(sct_subset),"PBRM1_mut")

Idents(pseudo_ifnb) <- "Group"

pval_df <- FindMarkers(object = pseudo_ifnb,
#                        features = open_peaks,
                               # slot = "counts",
                               ident.1 = 'MUT', 
                               ident.2 = 'WT',
                               logfc.threshold = 0,
                               min.cells.group = 1,
                               pseudocount.use = 10,
                               min.pct = 0,
                               test.use = "DESeq2")  %>% rownames_to_column("Gene") %>%
    rename(pval = p_val, log2FC = avg_log2FC, fdr = p_val_adj) %>% arrange(desc(log2FC))
  #
  log2FC_df <- FoldChange(object = pseudo_ifnb,
#                           features = open_peaks,
                                 slot = "data",
                                 ident.1 = 'MUT', 
                                 ident.2 = 'WT',
                                 logfc.threshold = 0,
                                 min.cells.group = 1,
                                 pseudocount.use = 1,
                                 min.pct = 0,
                                 test.use = "wilcox")  %>% rownames_to_column("Gene") %>%
    rename(log2FC = avg_log2FC) %>% arrange(desc(log2FC))

  #
  save.signif.deg(pval_df,log2FC_df,coverage_df, "All")
  



deg_df = dfread('All_ALL_DEG_.tsv')
nrow(deg_df)

deg_df %>% filter( fdr < 0.05) %>%  nrow

deg_df %>% filter(log2FC > 1, fdr < 0.05)  %>% arrange(pval) %>% head(n=10)

deg_df %>% filter( fdr < 0.05) %>% 
    filter(SYMBOL %in% c('LIX1')) %>% arrange(pval)

deg_df %>% filter( Peak == 'chr5-97140161-97140661') %>% 
    filter(SYMBOL %in% c('LIX1')) %>% arrange(pval)



print(paste(Sys.time(),"Done!"))








