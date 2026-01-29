library(ggplot2)
library(ggrepel)

pca_plot(final,
         pheno,
         group = "group",
         point_shape = "batch",
         pal = "auto",
         ellipse_conf = 0.75,
         point_size = 5,
         base_size = 18,
         title = "PCA — Group vs Batch",show_anova = TRUE)




pca_plot <- function(data,
                     pheno,
                     group,
                     names = FALSE,
                     x = 1,
                     y = 2,
                     alpha = 0.7,
                     point_size = 4,
                     base_size = 16,
                     title = NULL,
                     pal = "auto",
                     scale = TRUE,
                     ellipse_conf = 0.95,
                     show_ellipses = TRUE,
                     legend_position = "right",
                     point_shape = NULL,
                     show_anova = FALSE,
                     return_pca = FALSE){

  library(ggplot2)
  library(dplyr)
  library(ggrepel)

  # --- Palette auto ---
  .auto_palette <- function(groups) {
    n <- length(groups)
    if (n == 1) return("#1f77b4")
    if (n == 2) return(c("#1f77b4", "#ff7f0e"))
    if (n <= 8) return(RColorBrewer::brewer.pal(n, "Dark2"))
    viridisLite::viridis(n)
  }

  # --- Vérifications ---
  if (!is.data.frame(data)) stop("'data' doit être un data.frame numérique.")
  if (!all(sapply(data, is.numeric))) stop("Toutes les colonnes de 'data' doivent être numériques.")

  mat <- as.matrix(data)

  # --- PCA ---
  pca <- prcomp(t(mat), scale. = scale)

  # --- Variance expliquée ---
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  xlab <- paste0("PC", x, " (", round(var_expl[x] * 100, 1), "%)")
  ylab <- paste0("PC", y, " (", round(var_expl[y] * 100, 1), "%)")

  # --- Data PCA ---
  df <- pheno
  df$PC1 <- pca$x[, x]
  df$PC2 <- pca$x[, y]

  # --- ANOVA ---
  subtitle_text <- NULL
  if (isTRUE(show_anova)) {
    df$group_var <- df[[group]]
    anova_res <- summary(aov(PC1 ~ group_var, data = df))
    pval <- anova_res[[1]][["Pr(>F)"]][1]
    subtitle_text <- paste0("ANOVA PC1: p = ", signif(pval, 3))
  }

  # --- Aesthetics ---
  aes_points <- aes(
    x = PC1,
    y = PC2,
    color = !!as.name(group)
  )

  # --- Shapes optionnels ---
  if (!is.null(point_shape)) {
    df[[point_shape]] <- as.factor(df[[point_shape]])
    aes_points$shape <- as.name(point_shape)
  }

  # --- Plot ---
  p <- ggplot(df, aes_points) +
    geom_point(size = point_size, alpha = alpha) +
    theme_classic(base_size = base_size) +
    labs(x = xlab, y = ylab, title = title, subtitle = subtitle_text) +
    theme(legend.position = legend_position)

  # --- Palette couleurs ---
  if (is.character(pal) && pal == "auto") {
    pal_vals <- .auto_palette(unique(df[[group]]))
    p <- p + scale_color_manual(values = pal_vals)
  } else if (is.character(pal) && length(pal) == 1) {
    p <- p + scale_color_brewer(palette = pal)
  } else {
    p <- p + scale_color_manual(values = pal)
  }

  # --- Palette shapes ---
  if (!is.null(point_shape)) {
    p <- p + scale_shape_manual(values = c(0:25))
  }

  # --- Ellipses (corrigées : basées uniquement sur group) ---
  if (isTRUE(show_ellipses)) {
    p <- p +
      stat_ellipse(
        data = df,
        mapping = aes(
          x = PC1,
          y = PC2,
          group = !!as.name(group),
          color = !!as.name(group)
        ),
        level = ellipse_conf,
        linewidth = 1,
        linetype = "dashed",
        inherit.aes = FALSE
      )
  }

  # --- Labels ---
  if (isTRUE(names)) {
    p <- p +
      geom_text_repel(
        label = rownames(df),
        size = base_size / 4
      )
  }

  # --- Retour optionnel ---
  if (return_pca) {
    return(list(
      plot = p,
      pca = pca,
      anova_p = if (show_anova) pval else NULL
    ))
  }

  return(p)
}

