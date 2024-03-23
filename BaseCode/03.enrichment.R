


my.breaks <- c(seq(-2.5, -0.01, by = 0.001), seq(0.01, 2.5, by = 0.001) ) 
my.colors <- c(colorRampPalette(colors = c("#00599F","#287dba","#62aee5","#a8dcff","white"))(length(my.breaks)/2), 
               colorRampPalette(colors = c("white","#ffca79","#ef6e51","#db1b18","#b50600"))(length(my.breaks)/2))



plot.data <- data.frame(Gene = union(gene.enrich$SHAM_UCCAO_up, gene.enrich$Rep_UCCAO_up),
                        GeneFC = NA,
                        Type = "SHAM_UCCAO_up")
plot.data$GeneFC <- plot.data.2$logFC[match(plot.data$Gene, plot.data.2$Gene)]
plot.data$Type[plot.data$Gene %in% gene.enrich$Rep_UCCAO_up] <- "Rep_UCCAO_up"
plot.data$Type[plot.data$Gene %in% gene.enrich$SHAM_UCCAO_up & plot.data$Gene %in% gene.enrich$Rep_UCCAO_up] <- 'Common_up'
plot.data$Type <- factor(as.character(plot.data$Type), levels = c("Rep_UCCAO_up","Common_up","SHAM_UCCAO_up"))
plot.data <- plot.data[order(plot.data$Type, plot.data$GeneFC, decreasing = T), ]
rownames(plot.data) <- plot.data$Gene

plot.mat <- exp.fpkm.log[plot.data$Gene, ]

p <- pheatmap(as.matrix(plot.mat), scale = "row",
              color = my.colors, breaks = my.breaks,
              cluster_row = F, cluster_col = T, border_color = NA,
              annotation_col = anno.col,
              annotation_row = plot.data[, c(2,3)],
              annotation_colors = anno.color,
              clustering_method = "ward.D",
              clustering_distance_rows = "manhattan",
              clustering_distance_cols = "manhattan",
              fontsize_col = 10,
              fontsize_row = 1)
p




gmt.file <- rbind(gsea.hallmark, gsea.kegg, gsea.go )
deg.data.sub <- deg.data[match(unique(deg.data$GeneHuman), deg.data$GeneHuman), ]
foldChange <- deg.data.sub$logFC 
names(foldChange) <- deg.data.sub$GeneHuman
foldChange <- foldChange[order(foldChange, decreasing = T)]
set.seed(1)
edo.gsea <- GSEA(foldChange, TERM2GENE = gmt.file, minGSSize = 10, maxGSSize = 1000,
                 pvalueCutoff = 1, nPerm = 1000, exponent = 0.01)







