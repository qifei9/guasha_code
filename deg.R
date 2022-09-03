library("tidyverse")
library("DESeq2")
library("limma")
library("edgeR")
library("biomaRt")
library("openxlsx")
library("tximport")
source("./function.R")



tx2gene <- read.csv("./tx2gene.csv", header = F)

dir <- "./salmon_output/"
files <- list.files(path = dir, pattern = ".*sf", recursive = T, full.names = T)
txi <- tximport(files = files, type = "salmon", tx2gene = tx2gene)


sample_id <- str_c(c("before", "after"), "_", sort(rep(1:9, 2)))
group <- as.factor(rep(c("before", "after"), 9))
participant <- sort(rep(1:9, 2))

colnames(txi$abundance) <- sample_id
colnames(txi$counts) <- sample_id
colnames(txi$length) <- sample_id

meta <- data.frame(group = group, participant = participant, row.names = sample_id)

save(txi, file = "txi.RData")
save(meta, file = "meta.RData")


gene_id <- txi$counts %>% rownames()
mart <- useDataset(
    dataset = "hsapiens_gene_ensembl",
    useMart(biomart = "ensembl"))

ids <- getBM(
    attributes = c("ensembl_gene_id", "external_gene_name", "description"),
    filters = "ensembl_gene_id", values = gene_id, mart = mart)


dir.create("./diff_exp_analysis")
#-------------------------------------------------------------------------------
# DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ participant + group)
dds <- dds[rowSums(counts(dds)) > 18, ]
dds2 <- DESeq(dds)
res <- results(object = dds2, contrast = c("group", "after", "before"))
res_order <- res[order(res$padj), ] %>% as.data.frame()
res_order$ensembl_gene_id <- rownames(res_order)
deseq_result <- res_order %>%
    left_join(ids) %>%
     dplyr::select(ensembl_gene_id, everything())

write.xlsx(deseq_result, file = "./diff_exp_analysis/DESeq2_result.xlsx")


# edgeR
dat <- txi$counts
group <- meta$group
group <- relevel(group, ref = "before")  
participant <- meta$participant

y <- DGEList(counts = dat,group = group)
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ participant + group)
colnames(design) <- gsub("group", "", colnames(design))
rownames(design) <- colnames(y)

y <- estimateDisp(y, design = design, robust = T)
fit <- glmFit(y, design, robust = TRUE)
lrt <- glmLRT(fit)

topTags_res <- topTags(lrt, n = nrow(lrt), adjust.method = "BH") %>%
    as.data.frame()

edger_result <- topTags_res %>%
    mutate(ensembl_gene_id = rownames(.)) %>%
    left_join(ids) %>%
    dplyr::select(ensembl_gene_id, everything())

write.xlsx(edger_result, file = "./diff_exp_analysis/edgeR.xlsx")



# limma
dat <- txi$counts
group <- meta$group
participant <- meta$participant
sample <- rownames(meta)

dge <- DGEList(counts = dat)

dge$samples$group <- meta$group
dge$samples$participant <- meta$participant
keep.exprs <- filterByExpr(dge, group = group)
dge <- dge[keep.exprs, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge, method = "TMM")

design <- model.matrix(~ 0 + group + participant)
colnames(design) <- gsub("group", "", colnames(design))

constrasts <- paste(levels(meta$group), collapse = "-")
contrast_matrix <- makeContrasts(contrasts = constrasts, levels = design)

v <- voom(dge, design, plot = TRUE)
vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contrast_matrix)
efit <- eBayes(vfit)
deg <- topTable(efit, coef = 1, n = Inf)

limma_result <- deg %>%
    mutate(ensembl_gene_id = rownames(.)) %>%
    left_join(ids) %>%
    dplyr::select(ensembl_gene_id, everything())

write.xlsx(limma_result, file = "./diff_exp_analysis/limma.xlsx")
#-------------------------------------------------------------------------------


# GO reactome
deg_gene <- c(
  deseq_result %>% filter(padj < 0.1) %>% pull(ensembl_gene_id),
  edger_result %>% filter(FDR < 0.1) %>% pull(ensembl_gene_id),
  limma_result %>% filter(adj.P.Val < 0.1) %>% pull(ensembl_gene_id)) %>%
    table() %>%
    as_tibble() %>%
    filter(n == 3) %>%
    pull(1)

deg_gene <- list(DEG = deg_gene)
run_go(genelist = deg_gene, xlsxfilename = "DEG", outputpath = "./GO_DEG/")
run_reactome(genelist = deg_gene, xlsxfilename = "DEG", outputpath = "./reactome_DEG/")
