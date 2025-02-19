setwd("C:/Users/matub/OneDrive/Documentos/TMF/Análisis")
list.of.packages <- c("Seurat","tidyverse", "clusterProfiler", "org.Hs.eg.db",
                      "biomaRt", "speckle", "limma", "ggvenn", "hdWGCNA")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages())]
if(length(new.packages)> 0) install.packages(new.packages)

invisible(lapply(list.of.packages, FUN=library, character.only=TRUE))
set.seed(123)

############################################### PROCESAMIENTO DE CONTROL DE CALIDAD ###############################################
# Importar objeto Seurat
data_TP53 <- readRDS("C:/Users/matub/OneDrive/Documentos/TMF/Datos/so_TP53_ensembl.rds")
data_TP53 <- readRDS("C:/Users/matub/OneDrive/Documentos/TMF/Datos/so_TP53_symbol.rds") # repetir con los datos en formato symbol

# Visualizar métricas de control de calidad
summary(data_TP53@meta.data$percent.mt)
summary(data_TP53@meta.data$nFeature_RNA)
summary(data_TP53@meta.data$nCount_RNA)

vln1 <- VlnPlot(data_TP53, features = "nFeature_RNA", group.by = "sample", split.by = "TP53_mutation")
vln2 <- VlnPlot(data_TP53, features = "nCount_RNA", group.by = "sample", split.by = "TP53_mutation")
vln1 + vln2

RidgePlot(data_TP53, features="nFeature_RNA", group.by = "TP53_mutation") +
  scale_fill_viridis_d(option = "mako") +
  theme(axis.title.y = element_blank()) +
  NoLegend()

RidgePlot(data_TP53, features="nFeature_RNA", group.by = "sample") +
  scale_fill_viridis_d(option = "mako") +
  theme(axis.title.y = element_blank()) +
  NoLegend()

FeatureScatter(data_TP53, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
               group.by = "TP53_mutation", split.by = "sample", ncol = 4, shuffle = T)

# Filtrar células
data_TP53 <- subset(data_TP53, subset = nFeature_RNA >= 500 & nFeature_RNA <= 7500 & percent.mt <= 5)

data_TP53 <- subset(data_TP53, subset = nCount_RNA >= 1000 & nCount_RNA <= 25000)

# Normalizar los datos
data_TP53 <- NormalizeData(data_TP53, normalization.method = "LogNormalize", scale.factor = 10000) # Parámetros por defecto

# Identificación de características altamente variables (selección de características)
data_TP53 <- FindVariableFeatures(data_TP53, selection.method = "vst", nfeatures = 2000)

## Identificar los 10 genes más variables
top10 <- head(VariableFeatures(data_TP53), 10)

## Representar características variables con y sin etiquetas
var.feat.plot <- VariableFeaturePlot(data_TP53)
var.feat.plot <- LabelPoints(plot = var.feat.plot, points = top10, repel = TRUE)
var.feat.plot

# Escalar los datos
all.genes <- rownames(data_TP53)
data_TP53 <- ScaleData(data_TP53, features = all.genes)

# Realizar reducción dimensional lineal
data_TP53 <- RunPCA(data_TP53, features = VariableFeatures(object = data_TP53))
VizDimLoadings(data_TP53, dims = 1:2, reduction = "pca")
DimPlot(data_TP53, reduction = "pca", group.by = "ann_coarse", split.by = "TP53_mutation", pt.size = 1) +
  theme(plot.title = element_blank()) +
  scale_color_brewer(palette = "Paired")

DimHeatmap(data_TP53, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(data_TP53, dims = 1:30, cells = 500, balanced = TRUE)

# Determinar la ‘dimensionalidad’ del conjunto de datos
ElbowPlot(data_TP53, ndims = 50)

# Agrupar las células en clústeres
data_TP53 <- FindNeighbors(data_TP53, dims = 1:30)
data_TP53 <- FindClusters(data_TP53, resolution = 0.5)

head(Idents(data_TP53), 5) # Ver los IDs de clúster de las primeras 5 células

# Ejecutar reducción dimensional no lineal (UMAP/tSNE)
data_TP53 <- RunUMAP(data_TP53, dims = 1:30)
DimPlot(data_TP53, reduction = "umap", label = F, label.size = 5, group.by = "ann_coarse", split.by = "TP53_mutation") +
  theme(plot.title = element_blank()) +
  scale_color_brewer(palette = "Paired")

# Guardar el objeto
saveRDS(data_TP53, file = "C:/Users/matub/OneDrive/Documentos/TMF/Datos/so_TP53_ensembl_processed.rds")
saveRDS(data_TP53, file = "C:/Users/matub/OneDrive/Documentos/TMF/Datos/so_TP53_symbol_processed.rds")

########################################## PROPORCIONES DE TIPOS CELULARES ###########################################

# Representar proporciones de tipos celulares
meta.data_df <- data_TP53@meta.data %>%
  group_by(sample, TP53_mutation) %>%
  mutate(total_per_mutation = n()) %>%
  group_by(ann_coarse, sample, TP53_mutation) %>%
  reframe(Proportion = n() * 100 / unique(total_per_mutation))

means <- plyr::ddply(meta.data_df, ~ ann_coarse, function(x) c(mean=mean(x$Proportion)))

meta.data_df$ann_coarse <- factor(meta.data_df$ann_coarse, levels = c(means[order(-means$mean), "ann_coarse"]))

ggplot(meta.data_df, aes(x=sample, y=Proportion, fill=ann_coarse)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_brewer(palette="Paired") +
  theme(text=element_text(size=15),
        panel.background = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_wrap(~ TP53_mutation, scales = "free_x") +
  labs(x = "Muestras", fill = "Tipo celular")

# Ejecutar la prueba de propeller para evaluar diferencias en la proporción de tipos celulares entre los dos grupos
celltype_stats_coarse <- propeller(clusters = data_TP53@meta.data$ann_coarse, sample = data_TP53@meta.data$sample, 
                                   group = data_TP53@meta.data$TP53_mutation)

celltype_stats_fine <- propeller(clusters = data_TP53@meta.data$ann_fine, sample = data_TP53@meta.data$sample, 
                                 group = data_TP53@meta.data$TP53_mutation)

write.csv(celltype_stats_coarse, "celltype_stats_coarse.csv", row.names = F)
write.csv(celltype_stats_fine, "celltype_stats_fine.csv", row.names = F)

############################################### ANÁLISIS DE EXPRESIÓN DIFERENCIAL ###############################################

# Importar objeto Seurat procesado
data_TP53 <- readRDS("C:/Users/matub/OneDrive/Documentos/TMF/Datos/so_TP53_ensembl_processed.rds") # Genes en formato Ensembl
data_TP53 <- readRDS("C:/Users/matub/OneDrive/Documentos/TMF/Datos/so_TP53_symbol_processed.rds") # Genes en formato símbolo

## Subconjuntar objeto Seurat
data_TP53_wt <- subset(data_TP53, subset = TP53_mutation == "not mutated")
data_TP53_mt <- subset(data_TP53, subset = TP53_mutation == "mutated")

# Identificación de marcadores para los grupos de mutación de TP53
TP53.markers <- FindAllMarkers(data_TP53, only.pos = TRUE, group.by = "TP53_mutation") %>%
  filter(avg_log2FC > 1 & p_val_adj <= 0.05)

############################################## COEXPRESIÓN DEL GEN TP53 ##############################################

# Prepara los objetos de Seurat para el análisis de WGCNA (Weighted Gene Co-expression Network Analysis)
data_TP53_mt <- SetupForWGCNA(
  data_TP53_mt,
  gene_select = "fraction", # el método de selección de genes
  fraction = 0.05, # fracción mínima de células que deben expresar un gen para ser incluido
  wgcna_name = "TFM" # el nombre del experimento de hdWGCNA
)
data_TP53_wt <- SetupForWGCNA(
  data_TP53_wt,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "TFM"
)

# Construir metacélulas en cada grupo
data_TP53_mt <- MetacellsByGroups(
  seurat_obj = data_TP53_mt,
  group.by = c("ann_coarse", "TP53_mutation"), # columnas en data_TP53@meta.data para agrupar
  reduction = 'pca', # reducción dimensional para realizar KNN
  k = 25, # parámetro de KNN
  min_cells = 100, # número mínimo de células en un grupo
  max_shared = 15, # número máximo de células compartidas entre dos metacélulas
  ident.group = 'ann_coarse' # columnas que definen los Idents de las metacélulas
)

data_TP53_wt <- MetacellsByGroups(
  seurat_obj = data_TP53_wt,
  group.by = c("ann_coarse", "TP53_mutation"),
  reduction = 'pca',
  k = 25,
  min_cells = 100,
  max_shared = 15,
  ident.group = 'ann_coarse'
)

# Normalizar la matriz de expresión de las metacélulas
data_TP53_mt <- NormalizeMetacells(data_TP53_mt)
data_TP53_wt <- NormalizeMetacells(data_TP53_wt)

# Preparar la matriz de expresión para WGCNA
data_TP53_mt <- SetDatExpr(
  data_TP53_mt,
  group_name = "mutated", # grupo de interés
  group.by='TP53_mutation', # columna de metadatos con la información de mutación
  assay = 'RNA', # usar el ensayo de ARN
  layer = 'data' # usar los datos normalizados
)

data_TP53_wt <- SetDatExpr(
  data_TP53_wt,
  group_name = "not mutated",
  group.by='TP53_mutation',
  assay = 'RNA',
  layer = 'data'
)

# Seleccionar el umbral de soft-power
data_TP53_mt <- TestSoftPowers(
  data_TP53_mt,
  networkType = 'signed' # tipo de red
)
data_TP53_wt <- TestSoftPowers(
  data_TP53_wt,
  networkType = 'signed'
)

# Representar los resultados de las pruebas de soft-power:
plot_list <- PlotSoftPowers(data_TP53_mt)
patchwork::wrap_plots(plot_list, ncol=2)

# Construir la red de coexpresión
data_TP53_mt <- ConstructNetwork(
  data_TP53_mt,
  tom_name = 'TP53mt', # nombre de la matriz de solapamiento topológico
  overwrite_tom = TRUE # sobrescribir la matriz existente
)
data_TP53_wt <- ConstructNetwork(
  data_TP53_wt,
  tom_name = 'TP53wt',
  overwrite_tom = TRUE
)

# Representar el dendrograma para TP53mt
PlotDendrogram(data_TP53_mt, main='TP53mt hdWGCNA Dendrogram')

# Obtener los módulos de genes
gene_modules_mt <- GetModules(data_TP53_mt)
gene_modules_wt <- GetModules(data_TP53_wt)

# Extraer la expresión de TP53
tp53_expression_mt <- FetchData(data_TP53_mt, vars = "ENSG00000141510")
tp53_expression_wt <- FetchData(data_TP53_wt, vars = "ENSG00000141510")

# Encontrar los módulos de genes coexpresados con TP53
tp53_module_mt <- gene_modules_mt[gene_modules_mt$gene_name == "ENSG00000141510", "module"]
coexpressed_genes_mt <- gene_modules_mt[gene_modules_mt$module == tp53_module_mt, "gene_name"]

tp53_module_wt <- gene_modules_wt[gene_modules_wt$gene_name == "ENSG00000141510", "module"]
coexpressed_genes_wt <- gene_modules_wt[gene_modules_wt$module == tp53_module_wt, "gene_name"]

# Calcular correlación de Pearson entre TP53 y los genes coexpresados
correlation_results_mt <- lapply(coexpressed_genes_mt, function(gene) {
  gene_expression <- FetchData(data_TP53_mt, vars = gene)  # Obtener expresión del gen
  cor_test <- cor.test(tp53_expression_mt[,1], gene_expression[,1], method = "pearson")  # Correlación de Pearson
  return(data.frame(gene = gene, cor = cor_test$estimate, p.value = cor_test$p.value))
})

correlation_results_wt <- lapply(coexpressed_genes_wt, function(gene) {
  gene_expression <- FetchData(data_TP53_wt, vars = gene)  # Obtener expresión
  cor_test <- cor.test(tp53_expression_wt[,1], gene_expression[,1], method = "pearson")
  return(data.frame(gene = gene, cor = cor_test$estimate, p.value = cor_test$p.value))
})

# Convertir la lista de resultados en un dataframe
correlation_results_mt <- do.call(rbind, correlation_results_mt)
correlation_results_wt <- do.call(rbind, correlation_results_wt)

# Ajustar los valores p utilizando el método de Benjamini-Hochberg (FDR)
correlation_results_mt$adj.p.value <- p.adjust(correlation_results_mt$p.value, method = "BH")
correlation_results_wt$adj.p.value <- p.adjust(correlation_results_wt$p.value, method = "BH")

# Filtrar genes con p-value ajustado < 0.05
significant_genes_mt <- correlation_results_mt[correlation_results_mt$adj.p.value <= 0.01 & correlation_results_mt$cor >= 0, ]
significant_genes_wt <- correlation_results_wt[correlation_results_wt$adj.p.value <= 0.01 & correlation_results_wt$cor >= 0.04558091, ] # valor del ultimo gen en significant_genes_mt

# Ordenar por la correlación (valores absolutos mayores primero)
significant_genes_mt <- significant_genes_mt[order(significant_genes_mt$cor, decreasing = TRUE), ]
significant_genes_wt <- significant_genes_wt[order(significant_genes_wt$cor, decreasing = TRUE), ]

# Visualizar la ubicación de la expresión de los genes significativos
gene_symbols_mt <- mapIds(org.Hs.eg.db,
                          keys = significant_genes_mt$gene[1:12],
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")
gene_symbols_wt <- mapIds(org.Hs.eg.db,
                          keys = significant_genes_wt$gene[1:12],
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

# Representar la expresión de los genes significativos en las metacélulas
data_TP53_sym <- readRDS("C:/Users/matub/OneDrive/Documentos/TMF/Datos/so_TP53_symbol_processed.rds")
plots_mt <- FeaturePlot(data_TP53_sym, features = gene_symbols_mt[1:6])
plots_wt <- FeaturePlot(data_TP53_sym, features = gene_symbols_wt[1:6])

# Añadir etiquetas de tipo celular a los gráficos
labeled_plots_mt <- lapply(plots_mt, function(p) {
  p[["data"]] <- merge(p[["data"]], data_TP53@meta.data["ann_coarse"], by= 'row.names') %>% column_to_rownames("Row.names")
  LabelClusters(p, id = "ann_coarse", repel = TRUE)  # Etiquetar tipos celulares
})

labeled_plots_wt <- lapply(plots_wt, function(p) {
  p[["data"]] <- merge(p[["data"]], data_TP53@meta.data["ann_coarse"], by= 'row.names') %>% column_to_rownames("Row.names")
  LabelClusters(p, id = "ann_coarse", repel = TRUE)  # Etiquetar tipos celulares
})

# Combinar y visualizar los gráficos etiquetados
combined_plot_mt <- patchwork::wrap_plots(labeled_plots_mt, ncol = 3)  # Ajusta el número de columnas
print(combined_plot_mt)
combined_plot_wt <- patchwork::wrap_plots(labeled_plots_wt, ncol = 3)
print(combined_plot_wt)

venn <- list(`TP53mt` = significant_genes_mt$gene,
             `TP53wt` = significant_genes_wt$gene)
ggvenn(venn,
       fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
       set_name_size = 10, text_size = 10)

############################################### ENRIQUECIMIENTO DE GENES ###############################################

# Convertir los identificadores de Ensembl a ENTREZ
# Se utiliza la función `bitr` para convertir los identificadores de genes en el formato de Ensembl a ENTREZID
mt_entrez_cor <- bitr(significant_genes_mt$gene[-1], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
wt_entrez_cor <- bitr(significant_genes_wt$gene[2:223], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
mt_entrez <- bitr(TP53.markers$gene[TP53.markers$cluster == "mutated"], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
wt_entrez <- bitr(TP53.markers$gene[TP53.markers$cluster == "not mutated"], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Análisis de enriquecimiento de genes (GO) utilizando la base de datos de Gene Ontology (org.Hs.eg.db)
ego_mt_cor <- enrichGO(gene = mt_entrez_cor$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
ego_wt_cor <- enrichGO(gene = wt_entrez_cor$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
ego_mt <- enrichGO(gene = mt_entrez$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)
ego_wt <- enrichGO(gene = wt_entrez$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05)

# Visualización de los resultados de enriquecimiento de genes con un gráfico de puntos (dotplot)
# Para cada conjunto de genes enriquecidos (TP53mt y TP53wt, con y sin correlación), se realiza un gráfico para mostrar los 10 términos más significativos
plot_mt <- dotplot(ego_mt_cor, showCategory = 10, title = "GO Enrichment - TP53mt correlated genes")
plot_wt <- dotplot(ego_wt_cor, showCategory = 10, title = "GO Enrichment - TP53wt correlated genes")

# Mostrar los resultados de enriquecimiento de los genes diferencialmente expresados (DE)
dotplot(ego_mt, showCategory = 10, title = "GO Enrichment - TP53mt DE genes")
dotplot(ego_wt, showCategory = 10, title = "GO Enrichment - TP53wt DE genes")
