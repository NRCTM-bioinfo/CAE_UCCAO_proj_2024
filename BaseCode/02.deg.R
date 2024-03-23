

####### Read GSEA file
gsea.path <- c("msigdb_v2023.1.Hs_GMTs")
gsea.hallmark <- read.gmt(paste0(gsea.path, "/h.all.v2023.1.Hs.symbols.gmt"))
gsea.kegg <- read.gmt(paste0(gsea.path, "/c2.cp.kegg.v2023.1.Hs.symbols.gmt"))
gsea.go <- read.gmt(paste0(gsea.path, "/c5.go.bp.v2023.1.Hs.symbols.gmt"))
gsea.all <- read.gmt(paste0(gsea.path, "/msigdb.v2023.1.Hs.symbols.gmt"))


gene.id.mm <- read.table("gene.id.vM34.txt")
gene.id.mm.coding <- gene.id.mm[which(gene.id.mm$V6 == "protein_coding"), ]

gene.id <- read.table("gene.id.v43.txt")
gene.id.coding <- gene.id[which(gene.id$V6 == "protein_coding"), ]






############# PCA
res.pca.comp <- prcomp(exp.fpkm.log[rownames(exp.fpkm.log) %in% gene.id.mm.coding$V7, ], scale = TRUE)

plot.data <- as.data.frame(res.pca.comp$rotation)
plot.data$Name = rownames(plot.data)
plot.data$Group <- meta.data$Group


sub.1 <- signif(unlist(summary(aov(plot.data$PC1 ~ as.factor(plot.data$Group) )))[9], 5)
sub.2 <- signif(unlist(summary(aov(plot.data$PC2 ~ as.factor(plot.data$Group) )))[9], 5)
sub.3 <- signif(unlist(summary(aov(plot.data$PC3 ~ as.factor(plot.data$Group) )))[9], 5)

p <- ggscatter(plot.data, x = "PC1", y = "PC2",
               fill = "Group", color = "Group", size = 2,
               palette = color.group, 
               ellipse = F, ellipse.level = 0.9,
               label = plot.data$Name, font.label = 6,
               xlab = paste0("PC1 (P = ", sub.1, " )"), 
               ylab = paste0("PC2 (P = ", sub.2, " )")) + theme_base() 
p <- p + scale_x_continuous(limits = c(-0.3345,-0.332)) + scale_y_continuous(limits = c(-0.8,0.8))
ggsave(paste0("Sham.pca12.label.pdf"), p, width = 7, height = 5)

p <- ggscatter(plot.data, x = "PC2", y = "PC3",
               fill = "Group", color = "Group", size = 2,
               palette = color.group, 
               ellipse = F, ellipse.level = 0.9,
               label = plot.data$Name, font.label = 6,
               xlab = paste0("PC2 (P = ", sub.2, " )"), 
               ylab = paste0("PC3 (P = ", sub.3, " )")) + theme_base() 
ggsave(paste0("Sham.pca23.label.pdf"), p, width = 7, height = 5)

p <- ggscatter(plot.data, x = "PC1", y = "PC2",
               fill = "Group", color = "Group", size = 2,
               palette = color.group, 
               ellipse = F, ellipse.level = 0.9,
               #label = plot.data$Name, font.label = 6,
               xlab = paste0("PC1 (P = ", sub.1, " )"), 
               ylab = paste0("PC2 (P = ", sub.2, " )")) + theme_base() 
p <- p + scale_x_continuous(limits = c(-0.3345,-0.332)) + scale_y_continuous(limits = c(-0.8,0.8))
ggsave(paste0("Sham.pca12.pdf"), p, width = 7, height = 5)

p <- ggscatter(plot.data, x = "PC2", y = "PC3",
               fill = "Group", color = "Group", size = 2,
               palette = color.group, 
               ellipse = F, ellipse.level = 0.9,
               #label = plot.data$Name, font.label = 6,
               xlab = paste0("PC2 (P = ", sub.2, " )"), 
               ylab = paste0("PC3 (P = ", sub.3, " )")) + theme_base() 
ggsave(paste0("Sham.pca23.pdf"), p, width = 7, height = 5)

write.xlsx(plot.data, "PCA.xlsx")





############# DEG
pdata <- meta.data
pdata$contrast <- pdata$Group
design <- model.matrix(~ 0 + as.factor(contrast), data = pdata)
colnames(design) <- str_replace_all(colnames(design), fixed("as.factor(contrast)"), "")
fit <- lmFit(exp.fpkm.log[, pdata$Sample], design)
contrast <- makeContrasts(SHAM_UCCAO = SHAM - UCCAO, 
                          Rep_UCCAO = Reperfusion - UCCAO, 
                          Rep_SHAM = Reperfusion - SHAM, 
                          levels = design)
fits <- contrasts.fit(fit, contrast)
ebFit <- eBayes(fits)

for (kkk in 1:3) {
  deg_sig_list <- topTable(ebFit, coef = kkk, adjust.method = 'fdr', number = Inf)
  deg.data <- deg_sig_list[which(!is.na(deg_sig_list$adj.P.Val)), ]
  deg.data$logP <- -log10(deg.data$P.Value)
  deg.data$group = "not_sig"
  deg.data$group[which( (deg.data$P.Value < 0.05) & (deg.data$logFC > 0.58) )] = "up"
  deg.data$group[which( (deg.data$P.Value < 0.05) & (deg.data$logFC < -0.58) )] = "down"
  deg.data$tag <- colnames(contrast)[kkk]
  deg.data$Gene <- rownames(deg.data)
  deg.data$GeneType <- gene.id.mm$V6[match(deg.data$Gene, gene.id.mm$V7)]
  deg.data$GeneHuman <- out$HGNC.symbol[match(deg.data$Gene, out$MGI.symbol)]
  deg.data$GeneHuman[which(is.na(deg.data$GeneHuman))] <- toupper(deg.data$Gene[which(is.na(deg.data$GeneHuman))])
  deg.data$GeneHuman[! deg.data$GeneHuman %in% gene.id$V7] <- NA

  write.xlsx(deg.data, paste0( colnames(contrast)[kkk], ".DEGs.xlsx" ), rowNames = T, overwrite = T)
 }



