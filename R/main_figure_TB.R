library(ggplot2)
library(dplyr)
library(ggplot2)
library(dplyr)

counts = readRDS('~/sciFA/Data/Nathan_NatImm_2021.rds')
coldata = counts@meta.data

coldata = readRDS('~/scLMM/LMM-scRNAseq/Data/TB_immune_df.rds')
cell_count_df = data.frame(table(coldata$cluster_name))
cell_count_df = cell_count_df[order(cell_count_df$Freq, decreasing = T),]
summary(coldata)
colnames(coldata)
length(table(coldata$batch))
length(table(coldata$donor))
table(coldata$sex)
table(coldata$age)
table(coldata$season)


c25.2 = c(
    "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0", "#F0027F", "#BF5B17", "#666666", 
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", 
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", 
    "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#DECBE4", 
    "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC", "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", 
    "#F4CAE4", "#E6F5C9", "#FFF2AE", "#F1E2CC", "#CCCCCC", "#E41A1C", "#377EB8", "#4DAF4A", 
    "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", 
    "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#8DD3C7", "#FFFFB3", 
    "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", 
    "#CCEBC5", "#FFED6F"
  )

library(RColorBrewer)
head(coldata)
ggplot(coldata, aes(x=UMAP_1, y=UMAP_2, color=cluster_name))+geom_point(alpha=1, size=2)+
  theme_classic()+scale_color_manual(values = c25.2)+xlab('UMAP 1')+ylab('UMAP 2')+
  theme(text = element_text(size=16),legend.title = element_blank())#


a_cell_type = as.character(cell_count_df$Var1)[1]#'CNT' "CD4+ CD27+CD161+"
coldata$a_cell_type = ifelse(coldata$cluster_name==a_cell_type, a_cell_type, '')
cell_types = c("CD4+ activated" , "CD8+ activated")
coldata$selected_cell_type = ifelse(coldata$cluster_name %in% cell_types, coldata$cluster_name, '')
table(coldata$a_cell_type)
table(coldata$selected_cell_type)

ggplot(coldata, aes(x=UMAP_1, y=UMAP_2, color=a_cell_type))+geom_point(alpha=0.4, size=1.5)+
  theme_classic()+scale_color_manual(values = c('grey85', "green4"))+xlab('UMAP 1')+ylab('UMAP 2')+
  theme(text = element_text(size=16),legend.title = element_blank())

ggplot(coldata, aes(x=UMAP_1, y=UMAP_2, color=selected_cell_type))+geom_point(alpha=1, size=2.2)+
  theme_classic()+scale_color_manual(values = c('grey85', "maroon", "green4"))+xlab('UMAP 1')+ylab('UMAP 2')+
  theme(text = element_text(size=16),legend.title = element_blank())




# Define bold colors for selected cell types
selected_colors <- c("CD4+ activated" = "deepskyblue3",
                     "CD8+ activated" = "coral3")

# Define a muted palette for other cell types
muted_palette <- c(
  "#B3CDE3", "#CCEBC5", "#DECBE4", "#FED9A6", "#FFFFCC", "#E5D8BD", "#FDDAEC",
  "#F2F2F2", "#B3E2CD", "#FDCDAC", "#CBD5E8", "#F4CAE4", "#E6F5C9", "#FFF2AE"
)

# Assign muted colors to non-selected cell types
other_cell_types <- setdiff(unique(coldata$cluster_name), names(selected_colors))
other_colors <- setNames(rep(muted_palette, length.out = length(other_cell_types)), other_cell_types)

# Combine all colors
all_colors <- c(selected_colors, other_colors)

# Create the UMAP plot
ggplot(coldata, aes(x = UMAP_1, y = UMAP_2)) +
  # Plot non-selected cells with lighter colors and slightly smaller size
  geom_point(data = coldata %>% filter(!cluster_name %in% names(selected_colors)),
             aes(color = cluster_name), size = 2, alpha = 0.9) + #size = 1.6, alpha = 0.4
  # Plot selected cells with bold colors and slightly larger size
  geom_point(data = coldata %>% filter(cluster_name %in% names(selected_colors)),
             aes(color = cluster_name), size = 2, alpha = 1) + #size = 1, alpha = 0.6
  # Apply custom color scale
  scale_color_manual(values = all_colors) +
  theme_classic() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(
    text = element_text(size = 16),
    legend.title = element_blank()
  )

#### based on sample-type ##### 
coldata$cluster_name[is.na(coldata$cluster_name)] = 'unknown'
table(coldata$cluster_name)
coldata_counts <- ddply(coldata, .(coldata$TB_status, coldata$cluster_name), nrow)
names(coldata_counts) <- c("TB_status", "CellType", "Freq") #Sex, season, TB_status

ggplot(data=coldata_counts, aes(x=CellType, y=Freq, fill=TB_status)) +
  geom_bar(stat="identity",color='black')+theme_classic()+
  #+scale_fill_brewer(palette = "Blues")+
  ylab('Counts')+xlab('Clusters')+
  scale_fill_manual(values = c("gray70", '#FFD92F'))+
  #scale_fill_manual(values = c("pink1", "skyblue1"))+
  theme(text = element_text(size=15),
        axis.text.x = element_text(size=10,angle=90,color='black'),
        legend.title = element_blank())+
  xlab('')




coldata_counts_split <- split( coldata_counts , f = coldata_counts$CellType )
coldata_counts_split_norm <- lapply(coldata_counts_split, function(x) {x$Freq=x$Freq/sum(x$Freq);x})
counts_norm <- do.call(rbind,coldata_counts_split_norm )
head(counts_norm)

ggplot(data=counts_norm, aes(x=CellType, y=Freq, fill=TB_status)) +
  geom_bar(stat="identity",color='black',alpha=0.9)+theme_classic()+
  ylab('Proportion (%)')+
  #scale_fill_manual(values = c("#E69F00", "#009E73"))+
  #scale_fill_brewer(palette = 'Set2')+
  #scale_fill_manual(values = c("pink1", "skyblue1"))+
  scale_fill_manual(values = c("gray70", '#FFD92F'))+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=15,angle=-90,color='black'),
        legend.title = element_blank()) +  coord_flip()+
  xlab('')

ggplot(coldata, aes(x=Cell_Types_Broad, fill=sex))+
  theme_classic()+geom_bar(stat="identity")+
  theme(text = element_text(size=16),legend.title = element_blank())



age_df = data.frame(donor=coldata$donor,age=coldata$age)
age_df_unique = unique.data.frame(age_df)
age_df_unique = age_df_unique[order(age_df_unique$age, decreasing = F),]
dim(age_df_unique)
hist(age_df_unique$age)

ggplot(age_df_unique, aes(x=age)) + theme_classic()+
  geom_histogram(color="black", fill="palegreen")+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16,angle=-90,color='black'),
        axis.text.y = element_text(size=16,color='black'),
        legend.title = element_blank()) +coord_flip()+xlab('Age')+ylab('Donor count')


ggplot(age_df_unique, aes(x=age)) + theme_classic()+
  geom_histogram(color="black", fill="grey89")+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=16,color='black'),
        axis.text.y = element_text(size=16,color='black'),
        legend.title = element_blank())+xlab('Age')+ylab('Donor count')


######### scMM lmmfit functions
source("~/FLASH-MM/R/lmmfit.nt.R")
source("~/FLASH-MM/R/lmmfitSS.R")
source("~/FLASH-MM/R/lmmtest.R")
source("~/FLASH-MM/R/qqpvalue.R")
######### importing TB data  
#load("~/FLASH-MM/Results_data/TB_lmer_default_intercept.RData") 
load("~/FLASH-MM/Results_data/TB_lmm_intercept.RData") 
#fit <- readRDS('~/scLMM/LMM-scRNAseq-jan2024/lmmfitSS_Nathan_NatImm_2021_X_sexAgeSeason_Z_donor.rds') # Time difference of 1.200262 hours
#fit<- readRDS('~/scLMM/LMM-scRNAseq-jan2024/lmmfitSS_Nathan_NatImm_2021.rds') ## Time difference of 6.184403 hours
tvlmm <- t(fit$t)
pvlmm <- t(fit$p)
sum(apply(is.na(pvlmm), 1, any))
felmm <- t(fit$coef)
slmm <- t(fit$theta)

pvlmm_adj = sapply(1:ncol(pvlmm), function(i) p.adjust(pvlmm[,i], method = "fdr"))
colnames(pvlmm_adj) = colnames(pvlmm)

sum(apply(is.na(pvlmm), 1, any))
dim(pvlmm)
colnames(pvlmm)

P_VAL_THR = 0.05
COEF_THR = 0.2
cov_marker_list = list()
for(i in 31:ncol(pvlmm_adj)){
  df = data.frame(genes=rownames(pvlmm_adj),
                  tvalue=tvlmm[,i],
                  pvalue=pvlmm[,i],
                  coef = felmm[,i],
                  pvalue_adj=pvlmm_adj[,i],
                  score = -log(pvlmm_adj[,i]+1e-20)*felmm[,i])
  df_ord = df[order(df$score, decreasing = TRUE),]
  #df_ord$is_sig = df_ord$pvalue_adj<P_VAL_THR & (df_ord$coef > COEF_THR | df_ord$coef < -COEF_THR)
  #df_ord$is_sig = df_ord$pvalue_adj<P_VAL_THR & (df_ord$coef > COEF_THR | df_ord$coef < -COEF_THR)
  df_ord$is_sig = df_ord$pvalue_adj<P_VAL_THR & df_ord$coef > 0
  cov_marker_list[[colnames(pvlmm)[i]]] = df_ord
  
}
lapply(cov_marker_list, head)
res= lapply(cov_marker_list, function(x) sum(x$is_sig))
res= data.frame(numDEG=t(data.frame(res)))
res$cellType.cov= rownames(res)
res=res[order(res$numDEG, decreasing = T),]
res

num_de_df = res
num_de_df$cell_type = gsub('.trt', '',num_de_df$cellType.cov)
num_de_df <- num_de_df %>%
  mutate(cell_type = reorder(cell_type, numDEG))
# Create the bar plot with flipped axes, sorted cell types, text labels, and increased x-axis tick size
ggplot(num_de_df, aes(x = cell_type, y = numDEG)) +
  geom_bar(stat = "identity", fill = "skyblue", color = "black") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Number of TB-specific DE genes per cell type",
       x = "",
       y = "Number of DE genes") +
  geom_text(aes(label = numDEG), hjust = -0.2, color = "black", size = 3.5) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size = 11))  # Increase the size of x-axis tick labels (cell types)




colnames(pvlmm)
col_name = 'CD4p.activated:trt'#'CD8p.activated:trt' 
cell_trt_df = data.frame(pvalue=pvlmm[,col_name],
                        pvalue_adj = pvlmm_adj[,col_name],
                        tvalue = tvlmm[,col_name],
                        coef = felmm[,col_name],
                        gene=rownames(pvlmm))
cell_trt_df$pvalue_adj_log = -log(cell_trt_df$pvalue_adj+1e-800)
hist(cell_trt_df$tvalue+1e-10)
hist(cell_trt_df$coef)
hist(cell_trt_df$pvalue)
hist(cell_trt_df$pvalue_adj)

cell_trt_df[is.na(cell_trt_df$pvalue_adj),]
sum(cell_trt_df$pvalue_adj<P_VAL_THR & (cell_trt_df$coef > COEF_THR | cell_trt_df$coef < -COEF_THR))
cell_trt_df[cell_trt_df$pvalue_adj<P_VAL_THR & (cell_trt_df$coef > COEF_THR | cell_trt_df$coef < -COEF_THR),]
dim(cell_trt_df)
cell_trt_df = cell_trt_df[cell_trt_df$coef>0,]
dim(cell_trt_df)

cell_trt_df$score = -log(cell_trt_df$pvalue_adj+1e-600)*cell_trt_df$coef
cell_trt_df_ord = cell_trt_df[order(cell_trt_df$score, decreasing = TRUE),]
head(cell_trt_df_ord, 20)


############### Visualizing top TB specific genes
cell_trt_df_ord_m = head(cell_trt_df_ord, 25)
cell_trt_df_ord_m$gene = factor(cell_trt_df_ord_m$gene, levels = rev(cell_trt_df_ord_m$gene))

ggplot(cell_trt_df_ord_m, aes(x=gene, y=score,color=pvalue_adj_log))+
  scale_color_continuous('')+
  geom_point(size=3)+theme_classic()+ 
  scale_color_gradient(name='', low = "pink",
                       high = "brown4", space = "Lab" )+
  theme(axis.text.x = element_text(angle = -90, size = 12.5), 
        axis.text.y = element_text(size = 12.5), 
        axis.title.y = element_text(size = 17))+ylab('Score')+xlab('')


ggplot(cell_trt_df_ord_m, aes(x=gene, y=score,color=pvalue_adj_log))+
  scale_color_continuous('')+
  geom_point(size=3)+theme_classic()+ 
  scale_color_gradient(name='', low = "skyblue", #"#B2DF8A"
                       high = "royalblue4", space = "Lab" )+ #"darkgreen"
  theme(axis.text.x = element_text(angle = -90, size = 12.5), 
        axis.text.y = element_text(size = 12.5), 
        axis.title.y = element_text(size = 17))+ylab('Score')+xlab('')





#source('~/RatLiver/Codes/Functions.R')
#Initialize()
library(gprofiler2)
library(ggplot2)

get_gprofiler_enrich <- function(markers, model_animal_name){
  gostres <- gost(query = markers,
                  ordered_query = TRUE, exclude_iea =TRUE, 
                  sources=c('GO:BP' ,'REAC'),
                  organism = model_animal_name)
  return(gostres)
}

model_animal_name ='hsapiens'
num_genes = 300

################################
enrich_res_cd8 = get_gprofiler_enrich(markers=cell_tb_df_sort_m$gene[1:num_genes], model_animal_name)
enrich_res_cd4 = get_gprofiler_enrich(markers=cell_tb_df_sort_m$gene[1:num_genes], model_animal_name)
enrich_res = enrich_res_cd8#enrich_res_male

head(enrich_res$result,30)
enrich_res_pos = data.frame(enrich_res$result)

enrich_res_pos = enrich_res_pos[,!colnames(enrich_res_pos)%in%'evidence_codes']
enrich_res_pos$log_p = -log(as.numeric(enrich_res_pos$p_value))
enrich_res_pos = enrich_res_pos[order(enrich_res_pos$log_p, decreasing = T),]
View(data.frame(1:nrow(enrich_res_pos),enrich_res_pos$term_name, enrich_res_pos$intersection))
selected_terms=c(2,5,7:9,12,13,14,16) #cd4
selected_terms=c(1,3,4:8) #cd8
enrich_res_pos = enrich_res_pos[selected_terms,] 
enrich_res_pos$term_name[3] = 'T cell activation via TCR contact with MHC molecule on APC'
enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'log_p')]
enrich_res_pos$term_name <- factor(enrich_res_pos$term_name, levels =  enrich_res_pos$term_name[length(enrich_res_pos$term_name):1])

enrich_res_pos
title = ''
ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity',fill='pink3',color='grey10')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(title)+
  scale_fill_manual(values = c('pink3'))+
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = .5, face = "plain"))

ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity',fill="#B2DF8A",color='grey10')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(title)+
  scale_fill_manual(values = c("#B2DF8A"))+
  theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 0, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 15, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 17, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 17, angle = 90, hjust = .5, vjust = .5, face = "plain"))


################################
num_genes = 300
enrich_res = get_gprofiler_enrich(markers=a_cell_type_sex_df_female$genes[1:num_genes], model_animal_name)
head(enrich_res$result,30)
enrich_res_pos = data.frame(enrich_res$result)
enrich_res_pos = enrich_res_pos[1:20,]
enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'p_value')]
enrich_res_pos$log_p = -log(enrich_res_pos$p_value)
enrich_res_pos = enrich_res_pos[!is.na(enrich_res_pos$log_p),]
title = gsub(pattern = 'Male', 'Female', a_cell_type_sex)
ggplot(enrich_res_pos, aes(y=term_name,x=log_p))+geom_bar(stat = 'identity')+xlab('-log(p value)')+
  theme_classic()+ylab('')+ggtitle(paste0(title))

################################
enrich_res_pos = data.frame(enrich_res$result)
enrich_res_pos = enrich_res_pos[,!colnames(enrich_res_pos)%in%'evidence_codes']
enrich_res_pos$log_p = -log(as.numeric(enrich_res_pos$p_value))
enrich_res_pos = enrich_res_pos[order(enrich_res_pos$log_p, decreasing = T),]
View(data.frame(1:nrow(enrich_res_pos),enrich_res_pos$term_name, enrich_res_pos$intersection))
selected_terms=c(2:5,11,12,13,16,28) #male
selected_terms = c(1)
enrich_res_pos = enrich_res_pos[selected_terms,] ## positive - F18 - kidney dataset
enrich_res_pos = enrich_res_pos[,colnames(enrich_res_pos) %in% c('term_name', 'log_p')]




###############################################
###############################################
###############################################
hist(pvlmm)
pvalues_df=data.frame(pvalues=as.vector(pvlmm))
# basic histogram
ggplot(pvalues_df, aes(x=pvalues))+
  geom_histogram(fill="#69b3a2", color="#e9ecef", alpha=0.9)+theme_classic()+
  xlab('p-values')+ylab('Counts')+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=15,color='black'),
        legend.title = element_blank()) 


colnames(pvlmm_adj) = colnames(pvlmm)
colnames(pvlmm_adj)
head(pvlmm_adj)[,7:35]
ncol(head(pvlmm_adj))
Idents(counts) = counts$cluster_name

markers_list = list()
for(i in 7:35){
  df = data.frame(genes=rownames(pvlmm_adj),
                  tvalue=tvlmm[,i],
                  pvalue=pvlmm[,i],
                  pvalue_adj=pvlmm_adj[,i],
                  score = -log(pvlmm_adj[,i]+1e-20)*tvlmm[,i])
  df_ord = df[order(df$score, decreasing = TRUE),]
  markers_list[[colnames(pvlmm)[i]]] = df_ord
}

lapply(markers_list, head)
top_marker_genes = lapply(markers_list, function(x) rownames(x)[1:4])
### identifying top genes
all_markers_list = lapply(markers_list, function(x) rownames(x)[1:30])
pool_rep_genes = table(unlist(all_markers_list))>3 ### TRUE: rep
pool_rep_genes = names(pool_rep_genes)[pool_rep_genes] ### selecting genes which are repeated


#all_markers_list = unlist(lapply(markers_list, function(x) rownames(x)[1:20]))
all_markers_vec = unlist(all_markers_list)
all_markers_vec = all_markers_vec[!all_markers_vec %in% pool_rep_genes]
all_markers_vec = all_markers_vec[all_markers_vec %in% unique(unlist(top_marker_genes))]





names <- make.unique(cancer.rna$X)
rownames(cancer.rna) <- names
cancer.rna <- cancer.rna[,-1] # get rid of old names
cancer <- CreateSeuratObject(counts = cancer.rna)

markers_to_vis = unique(all_markers_vec)
markers_to_vis = c(T.cell,MNP.cell,NK.cell,CCD_like.cell,CNT.cell,cTAL.cell,
                   DCT.cell, Endothelial.cell, IC_A.cell, IC_B.cell, LOH_like.cell,
                   Mesangial.cell,PC.cell,PEC.cell,Podocyte.cell,Proximal.Tubule.cell,
                   U1.cell,U2.cell )

idents_order = c("B cell",  "T cell" ,  "MNP", "NK cell" ,"CCD-like",
                 "CNT",  "cTAL" , "DCT", "Endothelial", "IC-A" , "IC-B",
                 "LOH-like","Mesangial",'PC','PEC', "Podocyte",
                 "Proximal Tubule","U1" ,"U2"  )
Idents(data) = factor(data@active.ident, 
                      levels=idents_order)
unique_rownames = make.unique(rownames(counts))

rownames(counts@assays$originalexp) = unique_rownames
rownames(counts[['originalexp']]$counts) = unique_rownames
rownames(counts[['originalexp']]$data) = unique_rownames


counts_sub = counts[rownames(counts) %in% unique(make.unique(markers_to_vis)),]
test = CreateSeuratObject(counts_sub[['originalexp']]$data)
test@meta.data = counts@meta.data
Idents(test) = counts$cluster_name
test = SCTransform(object=test, 
                   assay = "RNA",
                   new.assay.name = "SCT",)
DotPlot(test, assay='RNA',features = unique(markers_to_vis))+
  theme(axis.text.x = element_text(angle = -90, size = 13),
        axis.text.y = element_text(size = 13))+
  ylab('')+xlab('')

