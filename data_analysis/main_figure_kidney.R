library(ggplot2)
library(RColorBrewer)
library(plyr)
library(stats)
library(ggpubr)
library(viridis)
library(scales)
library(gprofiler2)

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "palegreen", "gold1",
  "skyblue2", "orchid1", # lt pink
  "palegreen4",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "#FB9A99", "deeppink1", "steelblue4","blue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)

get_gprofiler_enrich <- function(markers, model_animal_name){
  gostres <- gost(query = markers,
                  ordered_query = TRUE, exclude_iea =TRUE, 
                  sources=c('GO:BP' ,'REAC', 'GO:MF', 'GO:CC', 'KEGG', 'CORUM', 'HP', 'WP'), #'TF',
                  organism = model_animal_name)
  return(gostres)
}

model_animal_name ='hsapiens'
head(PT_male_df_sort,30)
num_genes = 200

########################################################
############## data exploration and visualization
########################################################
umap_coord = read.table('~/scLMM/LMM-scRNAseq/Data/kidney_atlas/UMAP.coords.tsv.gz')
cell_meta = read.csv('~/scLMM/LMM-scRNAseq/Data/kidney_atlas/meta.tsv', sep = '\t')

sum(cell_meta$Cell != umap_coord$V1)
cell_meta = cbind(cell_meta, umap_coord)

ggplot(cell_meta, aes(x=V2, y=V3, color=Cell_Types_Broad))+geom_point(alpha=1, size=2)+
  theme_classic()+scale_color_manual(values = c25)+xlab('UMAP 1')+ylab('UMAP 2')+
  theme(text = element_text(size=16),legend.title = element_blank())#

ggplot(cell_meta, aes(x=V2, y=V3, color=sex))+geom_point(alpha=0.4, size=2)+
  theme_classic()+scale_color_manual(values = c("orchid1", 'lightskyblue'))+xlab('UMAP 1')+ylab('UMAP 2')+
  theme(text = element_text(size=16),legend.title = element_blank())#
a_cell_type = 'CNT'
cell_meta$a_cell_type = ifelse(cell_meta$Cell_Types_Broad==a_cell_type, a_cell_type, '')

ggplot(cell_meta, aes(x=V2, y=V3, color=a_cell_type))+geom_point(alpha=0.4, size=1.5)+
  theme_classic()+scale_color_manual(values = c('grey85', "green4"))+xlab('UMAP 1')+ylab('UMAP 2')+
  theme(text = element_text(size=16),legend.title = element_blank())#


unique_cell_types <- unique(cell_meta$Cell_Types_Broad)  # Get unique cell types
new_palette <- setNames(
  muted(rainbow(length(unique_cell_types)), l = 94, c = 30),  # More muted colors
  unique_cell_types
)
new_palette["CNT"] <- "#6A3D9A"  # Assign bold purple to cTAL

# Plot using the new palette
ggplot(cell_meta, aes(x=V2, y=V3, color=Cell_Types_Broad)) +
  geom_point(alpha=1, size=2) +
  theme_classic() +
  scale_color_manual(values = new_palette) +
  xlab('UMAP 1') +
  ylab('UMAP 2') +
  theme(text = element_text(size=17), legend.title = element_blank())


df = data.frame(table(cell_meta$Cell_Types_Broad))
df[order(df$Freq, decreasing = T),]
#### based on sample-type ##### 
cell_meta_counts <- ddply(cell_meta, .(cell_meta$sex, cell_meta$Cell_Types_Broad), nrow)
names(cell_meta_counts) <- c("Sex", "CellType", "Freq")
#names(c25) = cell_meta_counts$CellType
#cell_meta_counts$CellType= factor(cell_meta_counts$CellType, levels = as.character(c25) ) 

ggplot(data=cell_meta_counts, aes(x=CellType, y=Freq, fill=Sex)) +
  geom_bar(stat="identity",color='black')+theme_classic()+#+scale_fill_brewer(palette = "Blues")+
  ylab('Counts')+xlab('Clusters')+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'),
        legend.title = element_blank())+
  xlab('')

cell_meta_counts_split <- split( cell_meta_counts , f = cell_meta_counts$CellType )
cell_meta_counts_split_norm <- lapply(cell_meta_counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,cell_meta_counts_split_norm )
head(counts_norm)

ggplot(data=counts_norm, aes(x=CellType, y=Freq, fill=Sex)) +
  geom_bar(stat="identity",color='black',alpha=0.9)+theme_classic()+
  ylab('Fraction per cell type (%)')+
  scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=15,angle=90,color='black'),
        legend.title = element_blank()) +  
  xlab('')


######### scMM lmmfit functions
source("~/FLASH-MM/R/lmmfit.nt.R")
source("~/FLASH-MM/R/lmmfitSS.R")
source("~/FLASH-MM/R/lmmtest.R")
source("~/FLASH-MM/R/qqpvalue.R")

######### importing kidney data  
#load('LMM-scRNAseq-jan2024/Kidney_reanalysis_CC/data/kidney-counts-lmmfit.RData')
load("~/FLASH-MM/Results_data//kidney-counts-lmm.beta.RData")

##running time
rtlmm
##t-values
tvlmm <- t(fit$t)
##p-values
pvlmm <- t(fit$p)
dim(pvlmm)
sum(apply(is.na(pvlmm), 1, any))
felmm <- t(fit$coef)
slmm <- fit$theta

##LMM tests
#test <- lmmtest(fit) ##t-values
#tvlmm <- test[, grep("_t", colnames(test)), drop = F]


##p-values
#plmm <- test[, grep("_pvalue", colnames(test)), drop = F]
hist(pvlmm)
pvalues_df=data.frame(pvalues=as.vector(pvlmm))
# basic histogram
ggplot(pvalues_df, aes(x=pvalues))+
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9)+theme_classic()+
  xlab('p-values')+ylab('Counts')+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=15,color='black'),
        legend.title = element_blank()) 


pvlmm_adj = sapply(1:ncol(pvlmm), function(i)p.adjust(pvlmm[,i], method = "fdr"))
colnames(pvlmm_adj) = colnames(pvlmm)
colnames(felmm) == colnames(pvlmm)
sum(apply(is.na(pvlmm), 1, any))
dim(pvlmm)
dim(felmm)

#### counts extracted from data object shared via paper 
dirData = '~/scLMM/LMM-scRNAseq/Data/kidney_atlas/'
datafile <- paste0(dirData, "/Human_Kidney_data.rds")
data = readRDS(file = datafile)
coldata = data@meta.data
Idents(data) = data$Cell_Types_Broad

                                      
########################################################
############## cell type marker identification ##########
########################################################
COEF_THR = 0.2
cov_marker_list = list()
for(i in 20:ncol(pvlmm_adj)){
  df = data.frame(genes=rownames(pvlmm_adj),
                  tvalue=tvlmm[,i],
                  pvalue=pvlmm[,i],
                  coef = felmm[,i],
                  pvalue_adj=pvlmm_adj[,i],
                  score = -log(pvlmm_adj[,i]+1e-20)*felmm[,i])
  df_ord = df[order(df$score, decreasing = TRUE),]
  df_ord$is_sig = df_ord$pvalue_adj<P_VAL_THR & (df_ord$coef > COEF_THR | df_ord$coef < -COEF_THR)
  cov_marker_list[[colnames(pvlmm)[i]]] = df_ord
}
lapply(cov_marker_list, head)
res= lapply(cov_marker_list, function(x) sum(x$is_sig))
res= data.frame(numDEG=t(data.frame(res)))
res$cellType.cov= rownames(res)
res=res[order(res$numDEG, decreasing = T),]
res

#col_name = 'Proximal.Tubule:Male'
col_name ='CNT:Male' #'cTAL:Male'
cell_type_name = 'CNT'#'cTAL'

PT_male_df = data.frame(pvalue=pvlmm[,col_name],
                        pvalue_adj = pvlmm_adj[,col_name],
                        tvalue = tvlmm[,col_name],
                        coef = felmm[,col_name],
                        gene=rownames(pvlmm))

PT_male_df$pvalue_adj_log = -log(PT_male_df$pvalue_adj+1e-800)
hist(PT_male_df$tvalue+1e-10)
hist(PT_male_df$coef)
hist(PT_male_df$pvalue)
hist(PT_male_df$pvalue_adj)

PT_male_df[is.na(PT_male_df$pvalue_adj),]
# colname="cTAL:Male" ; 1848
sum(PT_male_df$pvalue_adj<P_VAL_THR & (PT_male_df$coef > COEF_THR | PT_male_df$coef < -COEF_THR))
PT_male_df[PT_male_df$pvalue_adj<P_VAL_THR & (PT_male_df$coef > COEF_THR | PT_male_df$coef < -COEF_THR),]


data_sub = data[,data$Cell_Types_Broad %in% cell_type_name]

PT_male_df$score = -log(PT_male_df$pvalue_adj+1e-600)*PT_male_df$coef
PT_male_df_ord = PT_male_df[order(PT_male_df$score, decreasing = TRUE),]

for(i in 1:10){
  gene_name = PT_male_df_ord$gene[i]
  gene_name
  #hist(counts[gene_name,], main = gene_name)
  df = data.frame(gene=counts[gene_name,], status=coldata$sex)
  head(df)
  df = df[coldata$Cell_Types_Broad==names(table(coldata$Cell_Types_Broad))[15],]
  #ggplot(df, aes(x=status, y=gene))+geom_boxplot()+ggtitle(paste0(gene_name,' ' ,col_name))
  p=ggplot(df, aes(x=gene, color=status))+geom_density()+ggtitle(paste0(gene_name,' ' ,col_name))
  print(p)
}

############### Visualizing top male specific genes
PT_male_df_sort = PT_male_df[order(PT_male_df$score, decreasing = T),]
PT_male_df_sort_m = head(PT_male_df_sort,35)
PT_male_df_sort_m$gene = factor(PT_male_df_sort_m$gene, levels = rev(PT_male_df_sort_m$gene))
ggplot(PT_male_df_sort_m, aes(y=gene, x=score,color=pvalue_adj_log))+
  scale_color_continuous('')+
  geom_point(size=3)+theme_classic()+ 
  theme(axis.text.x = element_text(angle = -90, size = 12.5), 
        axis.text.y = element_text(size = 12.5), 
        axis.title.x = element_text(size = 17))+xlab('Coefficient')+ylab('')


PT_female_df_sort = PT_male_df[order(PT_male_df$score, decreasing = F),]
PT_female_df_sort_m = head(PT_female_df_sort,35)
PT_female_df_sort_m = PT_female_df_sort_m[nrow(PT_female_df_sort_m):1,]
PT_female_df_sort_m$gene = factor(PT_female_df_sort_m$gene, levels = (PT_female_df_sort_m$gene))
ggplot(PT_female_df_sort_m, aes(y=gene, x=score,color=pvalue_adj_log))+
  scale_color_continuous('')+
  geom_point(size=3)+theme_classic()+ 
  scale_color_gradient(name='', low = "pink1",
                        high = "maroon", space = "Lab" )+
  theme(axis.text.x = element_text(angle = -90, size = 12.5), 
        axis.text.y = element_text(size = 12.5), 
        axis.title.x = element_text(size = 17))+
  xlab('')+ylab('')



################################
#### Pathway analysis 
############################
enrich_res_male = get_gprofiler_enrich(markers=PT_male_df_sort$gene[1:num_genes], model_animal_name)
enrich_res_female = get_gprofiler_enrich(markers=PT_female_df_sort$gene[1:num_genes], model_animal_name)

View(enrich_res_female$result)
View(enrich_res_male$result)

enrich_res = enrich_res_female# enrich_res_female

head(enrich_res$result,30)
enrich_res_pos = data.frame(enrich_res$result)
enrich_res_pos = enrich_res_pos[1:20,]
enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'p_value')]
enrich_res_pos$log_p = -log(enrich_res_pos$p_value)
title = paste0(a_cell_type)
ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(paste0(title))


################################
enrich_res_pos = data.frame(enrich_res$result)
enrich_res_pos = enrich_res_pos[,!colnames(enrich_res_pos)%in%'evidence_codes']
enrich_res_pos$log_p = -log(as.numeric(enrich_res_pos$p_value))
enrich_res_pos = enrich_res_pos[order(enrich_res_pos$log_p, decreasing = T),]

enrich_res_pos_vis = enrich_res_pos
enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'log_p')]
enrich_res_pos$term_name = gsub('metabolic process', 'metabolism',enrich_res_pos$term_name)
enrich_res_pos_vis  = enrich_res_pos[!1:nrow(enrich_res_pos) %in% c(2,10,16:20),]
rownames(enrich_res_pos_vis) = NULL
enrich_res_pos_vis$term_name[12] = 'Transport of ions and amino acids'

enrich_res_pos_vis$term_name <- factor(enrich_res_pos_vis$term_name, 
                                   levels =  enrich_res_pos_vis$term_name[length(enrich_res_pos_vis$term_name):1])

title = ''#'stim'#'Male'
ggplot(enrich_res_pos_vis, aes(y=term_name,x=log_p))+
  geom_bar(stat = 'identity',fill='skyblue',color='grey10')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(title)+
  scale_fill_manual(values = c('skyblue'))+
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = .5, face = "plain"))


ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+
  geom_bar(stat = 'identity',fill='pink',color='grey10')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(title)+
  scale_fill_manual(values = c('pink'))+
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = .5, face = "plain"))
################################
    
    
    


