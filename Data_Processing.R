# Conjunto de datos descargado de https://cellxgene.cziscience.com/collections/edb893ee-4066-4128-9aec-5eb2b03f8287
setwd("C:/Users/matub/OneDrive/Documentos/TMF/Datos")

list.of.packages <- c("Seurat","dplyr", "AnnotationDbi", "org.Hs.eg.db")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages())]
if(length(new.packages)> 0) install.packages(new.packages)

invisible(lapply(list.of.packages, FUN=library, character.only=TRUE))

# Importar datos del atlas
data <- readRDS("so_Atlas.rds")

# Exportar estadísticas de los datos
result <- lapply(names(dt), function(col) {
  # Calcular frecuencias y ordenar en orden descendente
  res <- dt[, .N, by = .(value = get(col))][order(-N)]
  
  # Identificar filas fuera del top 10
  if (nrow(res) > 10) {
    others <- sum(res$N[11:nrow(res)]) # Sumar frecuencias fuera del top 10
    res <- res[1:10]                   # Mantener solo las 10 primeras filas
    res <- rbind(res, data.table(value = "Otros", N = others)) # Agregar "Otros"
  }
  
  return(res)
})

names(result) <- names(dt) # Asignar nombres a la lista
write_xlsx(result, "atlas_stats.xlsx", format_headers = TRUE)

# Filtrar datos
data_TP53 <- subset(data, subset = dataset == "Laughney_Massague_2020")

cat("Número de células con TP53 mutado:", ncol(data_TP53), "\n")

# Extraer la matriz de conteo y los metadatos
counts_df <- data_TP53@assays[["RNA"]]@counts %>%
  as.data.frame()

metadata_cells <- data_TP53@meta.data %>% as.data.frame()

# Crear un nuevo objeto Seurat
so_TP53 <- CreateSeuratObject(counts = counts_df, meta.data = metadata_cells)

# Simplificar etiquetas
so_TP53@meta.data$sample <- sub("Laughney_Massague_2020_", "", so_TP53@meta.data$sample)

# Verificar los porcentajes mitocondriales y agregar una nueva columna
summary(so_TP53@meta.data$pct_counts_mito)
so_TP53[["percent.mt"]] <- PercentageFeatureSet(so_TP53, pattern = "^MT-")
so_TP53@meta.data$percent.mt <- as.numeric(so_TP53@meta.data$percent.mt)

# Guardar el nuevo objeto Seurat sin procesar
saveRDS(so_TP53, file = "so_TP53_ensembl.rds")

# Traducir genes ensembl a genes con símbolos
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(counts_df),
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

## Manejar duplicados agregando sufijos únicos
gene_symbols[is.na(gene_symbols)] <- rownames(counts_df)[is.na(gene_symbols)]
gene_symbols <- make.unique(gene_symbols)

## Asignar nuevos nombres de fila sin errores
rownames(counts_df) <- gene_symbols

# Crear un nuevo objeto Seurat
so_TP53 <- CreateSeuratObject(counts = counts_df, meta.data = metadata_cells)

# Guardar el nuevo objeto Seurat sin procesar
saveRDS(so_TP53, file = "so_TP53_symbol.rds")

