library("tidyverse")
library("DESeq2")
library("limma")
library("psych")
library("openxlsx")
source("./function.R")


load("txi.RData")
load("meta.RData")

dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ participant + group)
dds <- dds[rowSums(counts(dds)) > 18, ]
vsd <- varianceStabilizingTransformation(dds, blind = T)
design <- model.matrix(~ 0 + meta$group)
colnames(design) <- levels(meta$group)
rownames(design) <- rownames(meta)


# removeBatchEffect
assay(vsd) <- limma::removeBatchEffect(assay(vsd), batch = vsd$participant, design = design)
dat <- assay(vsd) %>% as.data.frame()

h2_id <- "ENSG00000187837"
h3_id <- "ENSG00000124575"
h4_id <- "ENSG00000168298"
id <- c(h2_id, h3_id, h4_id)

h_dat <-  as.data.frame(t(dat[id, ]))
dat <- dat[!(rownames(dat) %in% id), ] %>% t() %>% as.data.frame()


# Calculating the correlation
cor_result <- corr.test(h_dat, dat, method = "spearman", adjust = "none")
cor_result_r <- as.data.frame(t(cor_result$r)) %>% set_names(paste0("cor_", id))
cor_result_p <- as.data.frame(t(cor_result$p)) %>% set_names(paste0("pvalue_", id))

adjp1 <- p.adjust(cor_result_p$pvalue_ENSG00000187837, method = "BH", n = nrow(cor_result_p))
adjp2 <- p.adjust(cor_result_p$pvalue_ENSG00000124575, method = "BH", n = nrow(cor_result_p))
adjp3 <- p.adjust(cor_result_p$pvalue_ENSG00000168298, method = "BH", n = nrow(cor_result_p))

cor_result_adjp <- cbind(adjp1, adjp2, adjp3)
rownames(cor_result_adjp) <- rownames(cor_result_p)
colnames(cor_result_adjp) <- paste0("padj_", id)


cor_result_r_padj <- cbind(cor_result_r, cor_result_adjp)
write.csv(cor_result_r_padj,file = "cor.csv", row.names = T)

cor_result_r_padj_filter_postive <- cor_result_r_padj %>% filter((cor_ENSG00000187837 > 0.6 & padj_ENSG00000187837 < 0.05) |
                                                                 (cor_ENSG00000124575 > 0.6 & padj_ENSG00000124575 < 0.05) |
                                                                 (cor_ENSG00000168298 > 0.6 & padj_ENSG00000168298 < 0.05))


cor_result_r_padj_filter_negative <- cor_result_r_padj %>% filter((cor_ENSG00000187837 < -0.6 & padj_ENSG00000187837 < 0.05) |
                                                                  (cor_ENSG00000124575 < -0.6 & padj_ENSG00000124575 < 0.05) |
                                                                  (cor_ENSG00000168298 < -0.6 & padj_ENSG00000168298 < 0.05))



# GO reactome
postive_cor_gene <- list(postively = rownames(cor_result_r_padj_filter_postive))
run_go(genelist = postive_cor_gene, xlsxfilename = "postively", outputpath = "./GO_co_expression_postively/")
run_reactome(genelist = postive_cor_gene, xlsxfilename = "postively", outputpath = "./reactome_co_expression_postively/")


negative_cor_gene <- list(negatively = rownames(cor_result_r_padj_filter_negative))
run_go(genelist = negative_cor_gene , xlsxfilename = "negatively", outputpath = "./GO_co_expression_negatively/")
run_reactome(genelist = negative_cor_gene, xlsxfilename = "negatively", outputpath = "./reactome_co_expression_negatively/")
