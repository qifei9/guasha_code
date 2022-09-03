run_go <- function(genelist, xlsxfilename = "", outputpath, keyType = "ENSEMBL"){  
  library(openxlsx)
  library(clusterProfiler)
  library(stringr)
  
  wb <- createWorkbook()
  wb2 <- createWorkbook()
  
  if (!dir.exists(outputpath)) {
    dir.create(outputpath, recursive = T)
  }
  
  dir.create(paste0(outputpath, "/fig"))

  
  for (i in 1 : length(genelist)) {
    print(i)
    name <- names(genelist[i])
    gene <- as.data.frame(genelist[i])
    gene <- gene[, 1]
    gene <- grep("^ENSG", gene, value = T)
    
    enrich_go <- enrichGO(gene         = gene,
                         OrgDb         = "org.Hs.eg.db",
                         keyType       = keyType,
                         ont           = "all",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 1,
                         qvalueCutoff  = 1
    )
    
    if (is.null(enrich_go)) {
      addWorksheet(wb, sheetName = name)
      addWorksheet(wb2, sheetName = name)
      next
    }
    
    go_result <- enrich_go@result
    
    addWorksheet(wb, sheetName = name)
    writeData(wb, i, go_result, colNames = T, rowNames = T)
    
    filter_go_all <- enrich_go[enrich_go@result$p.adjust < 0.05, asis=T]
    addWorksheet(wb2, sheetName = name)
    writeData(wb2, i, filter_go_all, colNames = T, rowNames = T)

    filter_go_BP <- enrich_go[enrich_go@result$p.adjust < 0.05 & enrich_go@result$ONTOLOGY == "BP", asis=T]
    filter_go_MF <- enrich_go[enrich_go@result$p.adjust < 0.05 & enrich_go@result$ONTOLOGY == "MF", asis=T]
    filter_go_CC <- enrich_go[enrich_go@result$p.adjust < 0.05 & enrich_go@result$ONTOLOGY == "CC", asis=T]
    
    if (length(filter_go_BP@result$ID) >= 1) {
      pdf(file = paste0(outputpath, "/fig/GO_enrichment_BP_dotplot_", name, ".pdf"), width = 9, height = 5 + (length(filter_go_BP@result$ID) * 0.2))
      p <- dotplot(filter_go_BP, showCategory = length(filter_go_BP@result$ID), x = "Count") + scale_size(range = c(4, 8)) + scale_y_discrete(labels = function(x) str_wrap(x, width = 45))
      print(p)
      dev.off()
    }
    
    
    if (length(filter_go_MF@result$ID) >= 1) {
      pdf(file = paste0(outputpath, "/fig/GO_enrichment_MF_dotplot_", name, ".pdf"), width = 9, height = 5 + (length(filter_go_BP@result$ID) * 0.2))
      p <- dotplot(filter_go_MF, showCategory = length(filter_go_BP@result$ID), x = "Count") + scale_size(range = c(4, 8)) + scale_y_discrete(labels = function(x) str_wrap(x, width = 45))
      print(p)
      dev.off()
    }
    
    
    if (length(filter_go_CC@result$ID) >= 1) {
      
      p <- dotplot(filter_go_CC, showCategory = length(filter_go_BP@result$ID), x = "Count") + scale_size(range = c(4, 8)) + scale_y_discrete(labels = function(x) str_wrap(x, width = 45))
      pdf(file = paste0(outputpath, "/fig/GO_enrichment_CC_dotplot_", name, ".pdf"), width = 9, height = 5 + (length(filter_go_BP@result$ID) * 0.2))
      print(p)
      dev.off()
    }
    
    
  }
  
  saveWorkbook(wb, file = paste0(outputpath, "/GO_enrichment_result_all_", xlsxfilename, ".xlsx"), overwrite = T)
  saveWorkbook(wb2, file = paste0(outputpath, "/GO_enrichment_result_padj005_", xlsxfilename, ".xlsx"), overwrite = T)
  
}

run_reactome <- function(genelist, xlsxfilename = "", outputpath){ 
  library(ReactomePA)
  library(openxlsx)
  library(clusterProfiler)
  library(stringr)
  
  if (!file.exists(outputpath)) {
    dir.create(outputpath, recursive = T)
  }

  wb <- createWorkbook()
  wb2 <- createWorkbook()

  for (i in 1:length(genelist)) {
    print(i)
    name <- names(genelist[i])
    gene <- as.data.frame(genelist[i])
    gene <- gene[, 1]
    gene <- grep("^ENSG", gene, value = T)
    
    possibleError <- tryCatch(expr = {entrez <- bitr(gene, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")},
             error = function(e){
               "ENTREZ_id is NULL"
             })
    
    if (possibleError == "ENTREZ_id is NULL") {
      addWorksheet(wb,sheetName = name)
      addWorksheet(wb2,sheetName = name)
      print(str_c(name, ":ENTREZ_id is NULL"))
      next
    }

    entrez_gene <- entrez$ENTREZID
    enrich_reactome <-enrichPathway(gene = entrez_gene, organism = "human", pvalueCutoff = 1, qvalueCutoff = 1, readable = T)

    
    
    if (is.null(enrich_reactome)) {
      addWorksheet(wb, sheetName = name)
      addWorksheet(wb2, sheetName = name)
      next
    }
    
    enrich_reactome_result <- enrich_reactome@result
  
    addWorksheet(wb, sheetName = name)
    writeData(wb, i, enrich_reactome_result, colNames = T, rowNames = T)
    
    
    filter_enrich_reactome <- enrich_reactome[enrich_reactome@result$p.adjust < 0.05, asis=T] %>% as.data.frame()
    addWorksheet(wb2, sheetName = name)
    writeData(wb2, i, filter_enrich_reactome, colNames = T, rowNames = T)
     
    filter_enrich_reactome <- enrich_reactome[enrich_reactome@result$p.adjust < 0.05, asis=T]
    if (length(filter_enrich_reactome@result$ID) >= 1) {
      pdf(file = paste0(outputpath, "Reactome_enrichment_dotplot_", name, ".pdf"), width = 9, height = 5+length(filter_enrich_reactome@result$ID) * 0.2)
      p <- dotplot(filter_enrich_reactome, showCategory = length(filter_enrich_reactome@result$ID), x = "Count") + scale_size(range = c(4, 8)) + scale_y_discrete(labels = function(x) str_wrap(x, width = 45))
      print(p)
      dev.off()
    }
    
    
  }
  
  saveWorkbook(wb, file = paste0(outputpath, "reactome_enrichment_result_all_", xlsxfilename, ".xlsx"), overwrite = T)
  saveWorkbook(wb2, file = paste0(outputpath, "reactome_enrichment_result_padj005_", xlsxfilename, ".xlsx"), overwrite = T)
}
