# Setup
# Bioconductor and CRAN libraries used
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(tximport)
library(ggplot2)
library(ggrepel)

# 读取
samples <- dir("./data/raw/GSE163165_RAW/", full.names = TRUE)
sample_list <- lapply(samples, readr::read_delim, col_names = TRUE)

head(sample_list[[1]])
head(sample_list[[2]])
head(sample_list[[4]])
# 检查每张原始表中有多少个重复名称
sapply(sample_list, function(x) {
  sum(duplicated(x$gene_name))
})

x <- sample_list[[1]]

# 找出所有涉及重复名称的行，包括第一次出现的行
dup <- duplicated(x$gene_name) |
  duplicated(x$gene_name, fromLast = TRUE)

x[dup, ] |>
  dplyr::arrange(gene_name)

sapply(sample_list, function(x) {
  sum(duplicated(x$Geneid))
})

# 合并
merged_counts <- purrr::reduce(
  sample_list,
  dplyr::inner_join,
  by = c("Geneid", "gene_name")
) 

# 保存名称对应关系，方便后续注释和绘图
gene_annotation <- merged_counts |>
  dplyr::select(Geneid, gene_name)

# 使用 Geneid 作为行名
counts_table <- merged_counts |>
  dplyr::select(-gene_name) |>
  tibble::column_to_rownames(var = "Geneid")

# 修改原列名
newcolnames <- paste0(
  rep(c("Control_","LPS_"), times = 3),
  rep(1:3, each = 2))
colnames(counts_table) <- newcolnames

# 创建metadata
sampletype <- paste0(rep(c("Control", "LPS"), 3))
donor <- rep(paste0("D", 1:3), each = 2)
                     
meta <- data.frame(sampletype, donor, row.names = colnames(counts_table))
meta$donor <- factor(meta$donor)
# 创建 dds 对象
all(colnames(counts_table) %in% rownames(meta))
dds <- DESeqDataSetFromMatrix(counts_table, colData = meta, design = ~ donor + sampletype)

# 无监督聚类分析
rld <- rlog(dds, 
            blind = TRUE)
# PCA
pca_data <- DESeq2::plotPCA(
  rld,
  intgroup = c("sampletype", "donor"),
  returnData = TRUE
)
percent_var <- round(100 * attr(pca_data, "percentVar"))

ggplot2::ggplot(
  pca_data,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    colour = sampletype,
    shape = donor
  )
) +
  ggplot2::geom_point(size = 4) +
  ggplot2::labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    colour = "Condition",
    shape = "Donor"
  ) +
  ggplot2::coord_fixed() +
  ggplot2::theme_bw()

ggplot2::ggsave(
  "results/figures/PCA.pdf",
  width = 7,
  height = 6,
)
ggplot2::ggsave(
  "results/figures/PCA.png",
  width = 7,
  height = 6,
)
# 层次聚类
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)
pheatmap::pheatmap(
  rld_cor,
  annotation_col = meta,
  filename = "results/figures/sample_correlation_heatmap.png",
  width = 8,
  height = 7
)

# 运行 DESeq2
dds <- DESeq(dds)

# 检查离散度估计值的拟合情况
plotDispEsts(dds)

# wald假设检验
#由于只有一个分组变量和两个水平所以可以省略contrast
res <- results(dds, contrast = c("sampletype", "LPS", "Control"), alpha = 0.05)
png(
  filename = "results/figures/MA_before_shrinkage.png",
  width = 7,
  height = 5,
  units = "in",
  res = 300
)
DESeq2::plotMA(res, main = "Before LFC shrinkage")
dev.off()

# 效应大小收缩
res <- lfcShrink(dds, 
                 coef = "sampletype_LPS_vs_Control", 
                 type = "apeglm")
png(
  filename = "results/figures/MA_after_shrinkage.png",
  width = 7,
  height = 5,
  units = "in",
  res = 300
)
DESeq2::plotMA(res, main = "After LFC shrinkage")
dev.off()
# 差异基因结果与可视化
# set thresholds
padj.cutoff <- 0.05
# 将 results 对象转换为 tibble，以便与 tidyverse 函数一起使用
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

sig_res <- dplyr::filter(res_tbl, 
                         padj < padj.cutoff)

normalized_counts <- counts(dds, normalized = TRUE) %>% 
  data.frame() %>%
  rownames_to_column(var = "gene")

# DGE pheatmap
norm_sig <- normalized_counts %>% 
  dplyr::filter(gene %in% sig_res$gene)

# 设置颜色
heat_colors <- RColorBrewer::brewer.pal(6, "YlGnBu")

pheatmap::pheatmap(
  rld_cor,
  annotation_col = meta,
  
  width = 8,
  height = 7
)


pheatmap::pheatmap(
  norm_sig[2:7],
  color = heat_colors,
  cluster_rows = TRUE,
  show_rownames = FALSE,
  annotation_col = meta,
  filename = "results/figures/sig_heatmap.png",
  border_color = NA,
  fontsize = 10,
  scale = "row",
  width = 8,
  height = 10
)

# DGE volcano plot
# 添加列
res_tbl <- res_tbl %>% 
  dplyr::mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= log2(1.5))

# 对padj进行排序，标记前十基因
res_tbl <- res_tbl |> 
  dplyr::arrange(res_tbl$padj) |> 
  dplyr::mutate(
    genelabels = dplyr::if_else(
    dplyr::row_number() <= 10, 
    as.character(gene),
    ''
    )
  )


ggplot(res_tbl,
       aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(colour = threshold)) +
  geom_text_repel(aes(label = genelabels)) +
  labs(title = "LPS vs Control",
       x = expression(log[2]~Fold~Change),
       y = expression(-log[10]~Adjusted~P~Value)) +
  theme_bw() +
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

ggplot2::ggsave("results/DGE_volcano.png",
                width = 7,
                height = 6)
# 过度富集分析ORA
# Load libraries
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
# 进行功能分析需要 Ensembl 和 Entrez ID，而元数据中只有Gene Symbol
library(AnnotationHub)
library(ensembldb)
ah <- AnnotationHub()
# 将gene name添加到计结果中
res_new_tbl <- left_join(res_tbl, gene_annotation, by = c('gene' = 'Geneid'))
# 获取orgdb数据库
human_orgdb <- query(ah, c("org.Hs.eg.db.sqlite"))
human_orgdb <- human_orgdb[["AH121953"]]

# 映射Ensembl 和 Entrez ID
annotations_orgdb <- AnnotationDbi::select(
  human_orgdb, res_new_tbl$gene_name, c("SYMBOL", "ENSEMBL", "GENENAME", "ENTREZID"), "SYMBOL")
# 去重
entrez_map <- annotations_orgdb |>
  dplyr::select(SYMBOL, ENTREZID) |>
  dplyr::filter(!is.na(ENTREZID)) |>
  dplyr::distinct() |>
  dplyr::group_by(SYMBOL) |>
  # 保守处理：仅保留对应一个 Entrez ID 的名称
  dplyr::filter(dplyr::n() == 1) |>
  dplyr::ungroup()
# 映射添加到结果表
res_ora<- res_new_tbl |> 
  left_join(entrez_map, by = c('gene_name' = 'SYMBOL')) |>
  dplyr::filter(!is.na(padj), !is.na(ENTREZID)) |>
  dplyr::group_by(ENTREZID) |>
  # gene 保存原始 Geneid：
  # 仅保留一个 Entrez ID 对应一个有效原始记录的情况
  dplyr::filter(dplyr::n_distinct(gene) == 1) |>
  dplyr::ungroup()

# 将上调基因和下调基因分组
up_ids <- res_ora |>
  dplyr::filter(
    padj < 0.05,
    log2FoldChange >= log2(1.5)
  ) |>
  dplyr::pull(ENTREZID) |>
  unique()

down_ids <- res_ora |>
  dplyr::filter(
    padj < 0.05,
    log2FoldChange <= -log2(1.5)
  ) |>
  dplyr::pull(ENTREZID) |>
  unique()
# GO
background_ids <- unique(res_ora$ENTREZID)
de_ids <- res_ora |>
  dplyr::filter(
    padj < 0.05,
    abs(log2FoldChange) >= log2(1.5)
  ) |>
  dplyr::pull(ENTREZID) |>
  unique()

ego <- enrichGO(gene = de_ids, 
                universe = background_ids,
                keyType = "ENTREZID",
                OrgDb = human_orgdb, 
                ont = "ALL", 
                pAdjustMethod = "BH", 
                pool = FALSE, 
                readable = TRUE)

cluster_summary <- data.frame(ego)
write.csv(cluster_summary, "results/clusterProfilerGO.csv")
# 气泡图
ego_bp <- ego |>
  dplyr::filter(ONTOLOGY == "BP")

dotplot(ego, showCategory = 30)

ggplot2::ggsave(
  "results/figures/GO_BPdotplot.pdf",
  width = 8,
  height = 20
)
ggplot2::ggsave(
  "results/figures/GO_BPdotplot.png",
  width = 8,
  height = 20
)
# 富集图
ego_bp <- enrichplot::pairwise_termsim(ego_bp)
enrichplot::emapplot(
  ego_bp,
  showCategory = 30
)
ggplot2::ggsave('results/figures/enrichment_GO.pdf', width = 10, height = 10)
# 类别网络图
foldchanges <- res_ora$log2FoldChange
names(foldchanges) <- res_ora$gene

cnetplot(ego_bp, 
         showCategory = 5, 
         foldChange = foldchanges)

ggplot2::ggsave('results/figures/category_netplot.pdf', width = 10, height = 10)

# GO的上调分析
ego_up <- enrichGO(gene = up_ids, 
                universe = background_ids,
                keyType = "ENTREZID",
                OrgDb = human_orgdb, 
                ont = "ALL", 
                pAdjustMethod = "BH", 
                pool = FALSE, 
                readable = TRUE)
cluster_summary_up <- data.frame(ego_up)
write.csv(cluster_summary_up, "results/clusterProfilerGOup.csv")
# 上调的bp气泡图
ego_up_bp <- ego_up |>
  dplyr::filter(ONTOLOGY == "BP")

dotplot(ego_up_bp, showCategory = 30) +
  ggplot2::ggtitle("GO BP enrichment: upregulated genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave(
  "results/figures/GO_UP_BPdotplot.pdf",
  width = 8,
  height = 20
)
ggplot2::ggsave(
  "results/figures/GO_UP_BPdotplot.png",
  width = 8,
  height = 20
)
# 上调的bp网络富集图
ego_up_bp <- enrichplot::pairwise_termsim(ego_up_bp)
enrichplot::emapplot(
  ego_up_bp,
  showCategory = 30
) + 
  ggplot2::ggtitle("GO BP enrichment: upregulated genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave('results/figures/enrichment_GOup.pdf', width = 10, height = 10)
ggplot2::ggsave('results/figures/enrichment_GOup.png', width = 10, height = 10)
# 上调的类别网络图
foldchanges <- res_ora$log2FoldChange
names(foldchanges) <- res_ora$gene

cnetplot(ego_up_bp, 
         showCategory = 5, 
         foldChange = foldchanges)+ 
  ggplot2::ggtitle("GO BP enrichment: upregulated genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave('results/figures/category_netplot_GOup.png', width = 10, height = 10)

# GO的下调分析
ego_down <- enrichGO(gene = down_ids, 
                   universe = background_ids,
                   keyType = "ENTREZID",
                   OrgDb = human_orgdb, 
                   ont = "ALL", 
                   pAdjustMethod = "BH", 
                   pool = FALSE, 
                   readable = TRUE)
cluster_summary_down <- data.frame(ego_down)
write.csv(cluster_summary_down, "results/clusterProfilerGOdown.csv")
# 下调的bp气泡图
ego_down_bp <- ego_down |>
  dplyr::filter(ONTOLOGY == "BP")

dotplot(ego_down_bp, showCategory = 30) +
  ggplot2::ggtitle("GO BP enrichment: downregulated genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave(
  "results/figures/GO_DOWN_BPdotplot.png",
  width = 8,
  height = 20
)
# 下调的bp网络富集图
ego_down_bp <- enrichplot::pairwise_termsim(ego_down_bp)
enrichplot::emapplot(
  ego_down_bp,
  showCategory = 30
) + 
  ggplot2::ggtitle("GO BP enrichment: downregulated genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave('results/figures/enrichment_GOdown.png', width = 10, height = 10)
# 下调的类别网络图
cnetplot(ego_down_bp, 
         showCategory = 5, 
         foldChange = foldchanges)+ 
  ggplot2::ggtitle("GO BP enrichment: downregulated genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave('results/figures/category_netplot_GOdown.png', width = 10, height = 10)

# KEgg的上调分析
ekegg_up <- clusterProfiler::enrichKEGG(
  gene = up_ids,
  universe = background_ids,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 1
)

# 将结果中用于展示的 Entrez ID 转换为基因符号
ekegg_up <- clusterProfiler::setReadable(
  ekegg_up,
  OrgDb = human_orgdb,
  keyType = "ENTREZID"
)
kegg_summary_up <- data.frame(ekegg_up)
write.csv(kegg_summary_up, "results/clusterProfilerKEGGup.csv")
# 上调的气泡图
enrichplot::dotplot(
  ekegg_up,
  showCategory = 30
) +
  ggplot2::ggtitle("KEGG enrichment: upregulated genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave(
  "results/figures/KEGG_UP_BPdotplot.png",
  width = 8,
  height = 20
)
# 上调的网络图
cnetplot(ekegg_up, 
         showCategory = 5, 
         foldChange = foldchanges)+ 
  ggplot2::ggtitle("KEGG enrichment: upregulated genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave('results/figures/category_netplot_KEGGup.png', width = 10, height = 10)
# 上调通路图，TNF，hsa04668
BiocManager::install("pathview")

fc_tbl <- res_ora |>
  dplyr::filter(
    !is.na(ENTREZID),
    is.finite(log2FoldChange)
  )

# 确保每个 Entrez ID 只对应一条记录
stopifnot(!anyDuplicated(fc_tbl$ENTREZID))

gene_fc <- setNames(
  fc_tbl$log2FoldChange,
  as.character(fc_tbl$ENTREZID)
)

head(gene_fc)

as.data.frame(ekegg_up) |>
  dplyr::select(ID, Description, p.adjust, Count) |>
  head(20)

dir.create(
  "results/figures/KEGG_pathview",
  recursive = TRUE,
  showWarnings = FALSE
)

pv <- pathview::pathview(
  gene.data = gene_fc,
  pathway.id = "04668",  # 替换为你选择的通路编号
  species = "hsa",
  gene.idtype = "entrez",
  kegg.native = TRUE, 
  limit = list(gene = 2, cpd = 1), # 设置 log2FC 色阶范围：-2 到 2
  kegg.dir = "results/figures/KEGG_pathview",
  out.suffix = "LPS_vs_Control_TNF",
  pdf.size = c(12, 10)  # PDF 页面宽、高，单位为英寸
)

# KEgg的下调分析
ekegg_down <- clusterProfiler::enrichKEGG(
  gene = down_ids,
  universe = background_ids,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 1
)

# 将结果中用于展示的 Entrez ID 转换为基因符号
ekegg_down <- clusterProfiler::setReadable(
  ekegg_down,
  OrgDb = human_orgdb,
  keyType = "ENTREZID"
)
kegg_summary_down <- data.frame(ekegg_down)
write.csv(kegg_summary_down, "results/clusterProfilerKEGGdown.csv")

# 下调的气泡图
enrichplot::dotplot(
  ekegg_down,
  showCategory = 30
) +
  ggplot2::ggtitle("KEGG enrichment: downregulated genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave(
  "results/figures/KEGG_down_dotplot.png",
  width = 8,
  height = 20
)
# 下调的网络图
cnetplot(ekegg_down, 
         showCategory = 5, 
         foldChange = foldchanges)+ 
  ggplot2::ggtitle("KEGG enrichment: downregulated genes") +
  ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

ggplot2::ggsave('results/figures/category_netplot_KEGGdown.png', width = 10, height = 10)

# 下调通路图，氧化磷酸化（hsa00190）
pv_oxphos <- pathview::pathview(
  gene.data = gene_fc,
  pathway.id = "00190",  # 替换为你选择的通路编号
  species = "hsa",
  gene.idtype = "entrez",
  kegg.native = TRUE, 
  limit = list(gene = 2, cpd = 1), # 设置 log2FC 色阶范围：-2 到 2
  kegg.dir = "results/figures/KEGG_pathview",
  out.suffix = "LPS_vs_Control_00190",
  pdf.size = c(12, 10)  # PDF 页面宽、高，单位为英寸
)

# 功能类别评分FCS
# GSEA
remove.packages("curl")
install.packages("curl")
install.packages("msigdbr")
# 获取人类 Hallmark 基因集
hallmark <- msigdbr::msigdbr(
  species = "Homo sapiens",
  collection = "H"
)
hallmark_t2g <- hallmark |>
  dplyr::distinct(gs_name, ncbi_gene) |>
  dplyr::filter(!is.na(ncbi_gene))

# 由于Hallmark分析要传入res结果中的stat列，但是由于之前res表在lfcShrink()收缩后就不在保留stat的列
# ，所以这里要重新获取res结果
res_raw <- DESeq2::results(
  dds,
  contrast = c("sampletype", "LPS", "Control"),
  alpha = 0.05
)
library(tibble)

colnames(res_raw)
res_tbl_1 <- res_raw |> 
  data.frame() |> 
  rownames_to_column(var="gene") |>  
  as_tibble()

#将stat映射到去重的res_ora中
res_ora_stat <- left_join(res_ora, res_tbl_1[, c('gene', 'stat')], by = 'gene')
gene_list <- res_ora_stat$stat
names(gene_list) <- res_ora_stat$ENTREZID
sum(duplicated(res_ora_stat$ENTREZID))

gene_list <- sort(gene_list, decreasing = TRUE)

set.seed(123)

gsea_hallmark <- clusterProfiler::GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_t2g,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  minGSSize = 15,
  maxGSSize = 500,
  verbose = FALSE
)

# 提取调整后 P 值 < 0.05 的结果
hallmark_sig <- as.data.frame(gsea_hallmark) |>
  dplyr::filter(p.adjust < 0.05) |>
  dplyr::arrange(p.adjust)

write.csv(hallmark_sig, "results/gsea_hallmark.csv", row.names = FALSE, quote = TRUE)
view(hallmark_sig)
# NES排序图
plot_tbl <- hallmark_sig |>
  dplyr::mutate(
    direction = ifelse(NES > 0, "Positive NES", "Negative NES")
  ) |>
  dplyr::group_by(direction) |>
  dplyr::slice_min(p.adjust, n = 10, with_ties = FALSE) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    label = sub("^HALLMARK_", "", Description),
    label = gsub("_", " ", label)
  )

p_nes <- ggplot2::ggplot(
  plot_tbl,
  ggplot2::aes(
    x = NES,
    y = reorder(label, NES), # 按 NES 大小排列通路名称
    fill = direction # 根据方向给图形内部填色
  )
) +
  ggplot2::geom_col(width = 0.7) +
  ggplot2::geom_vline(xintercept = 0, colour = "grey50") +
  ggplot2::scale_fill_manual(
    values = c("Positive NES" = "#D73027",
               "Negative NES" = "#4575B4")
  ) +
  ggplot2::labs(
    title = "Hallmark GSEA",
    x = "Normalized enrichment score (NES)",
    y = NULL
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "none",
                 plot.title = element_text(hjust = 0.5, face = "bold"))

p_nes

ggplot2::ggsave(
  filename = "results/figures/Hallmark_GSEA_NES.png",
  plot = p_nes,
  width = 10,
  height = 8,
  units = "in"
)

# 绘制GSEA富集曲线
selected_sets <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)

# 查看所选基因集的统计结果
as.data.frame(gsea_hallmark) |>
  dplyr::filter(ID %in% selected_sets) |>
  dplyr::select(ID, NES, p.adjust, setSize)

# 同一张图绘制三条富集曲线
p_gsea <- enrichplot::gseaplot2(
  gsea_hallmark,
  geneSetID = selected_sets,
  color = c("#D73027", "#E69F00", "#4575B4")
)

p_gsea
ggplot2::ggsave(
  filename = "results/figures/Hallmark_GSEA_NES_enrich.png",
  plot = p_gsea,
  width = 10,
  height = 8,
  units = "in"
)

# 想展示的基因
genes_to_plot <- c("IL1B", "IL6", "TNF", "CCL2")

# 创建保存目录
out_dir <- "results/figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (gene_name in genes_to_plot) {
  
  # 检查基因是否存在于 dds 行名中
  if (!gene_name %in% rownames(dds)) {
    warning(paste("跳过：dds 中找不到基因", gene_name))
    next
  }
  
  # 提取当前基因的归一化计数
  plot_data <- DESeq2::plotCounts(
    dds,
    gene = gene_name,
    intgroup = c("sampletype", "donor"),
    returnData = TRUE
  )
  
  # 绘图
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = sampletype,
      y = count,
      group = donor,
      colour = donor
    )
  ) +
    ggplot2::geom_line(alpha = 0.6) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = gene_name,
      x = NULL,
      y = "Normalized count (log10 scale)"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold"
      )
    )
  
  # 循环中显示图片，需要显式 print()
  print(p)
  
  # 每个基因保存一张 PNG
  ggplot2::ggsave(
    filename = file.path(out_dir, paste0(gene_name, ".png")),
    plot = p,
    width = 5,
    height = 4,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}

infla <- sig_res |> dplyr::filter(gene %in% c('IL1B','IL6','TNF', 'CCL2'))
infla

# Print sessionInfo()
sessionInfo()
