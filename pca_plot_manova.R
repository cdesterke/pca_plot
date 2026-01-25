pca_plot(process,
         pheno,
         group = "group",
         point_shape = "batch",
         pal = "auto",
         ellipse_conf = 0.75,
         point_size = 5,
         base_size = 18,
         title = "PCA Group vs Batch",
         show_manova = TRUE)





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
                     show_manova = FALSE,
                     return_pca = FALSE){

  # --- Palette auto robuste ---
  .auto_palette <- function(groups) {
    n <- length(groups)
    if (n == 1) return("#1f77b4")
    if (n == 2) return(c("#1f77b4", "#ff7f0e"))
    if (n <= 8) return(RColorBrewer::brewer.pal(n, "Dark2"))
    return(viridisLite::viridis(n))
  }

  # --- VÃ©rification + conversion ---
  if (!is.data.frame(data)) stop("'data' doit Ãªtre un data.frame numÃ©rique.")
  if (!all(sapply(data, is.numeric))) stop("Toutes les colonnes de 'data' doivent Ãªtre numÃ©riques.")
  mat <- as.matrix(data)

  # --- PCA ---
  trans <- t(mat)
  pca <- prcomp(trans, scale. = scale)

  # --- Variance expliquÃ©e ---
  var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
  xlab <- paste0("PC", x, " (", round(var_expl[x] * 100, 1), "%)")
  ylab <- paste0("PC", y, " (", round(var_expl[y] * 100, 1), "%)")

  # --- Data PCA ---
  df <- pheno
  df$PC1 <- pca$x[, x]
  df$PC2 <- pca$x[, y]
  df$group_var <- df[[group]]

  # --- MANOVA sur PC1 + PC2 ---
  subtitle_text <- NULL
  if (isTRUE(show_manova)) {
    manova_res <- summary(manova(cbind(PC1, PC2) ~ group_var, data = df))
    pval <- manova_res$stats[1, "Pr(>F)"]
    subtitle_text <- paste0("MANOVA (PC1+PC2): p = ", signif(pval, 3))
  }

  # --- Aesthetics pour les points ---
  aes_points <- aes(
    x = PC1,
    y = PC2,
    color = !!as.name(group)
  )

  if (!is.null(point_shape)) {
    aes_points$shape <- as.name(point_shape)
  }

  # --- Plot de base ---
  p <- ggplot(df, aes_points) +
       geom_point(size = point_size, alpha = alpha) +
       theme_classic(base_size = base_size) +
       labs(x = xlab, y = ylab, title = title, subtitle = subtitle_text) +
       theme(legend.position = legend_position)

  # --- Palette ---
  if (is.character(pal) && length(pal) == 1 && pal == "auto") {
    pal_vals <- .auto_palette(unique(df[[group]]))
    p <- p + scale_color_manual(values = pal_vals)
  } else if (is.character(pal) && length(pal) == 1) {
    p <- p + scale_color_brewer(palette = pal)
  } else {
    p <- p + scale_color_manual(values = pal)
  }

  # --- Ellipses uniquement selon group ---
  if (isTRUE(show_ellipses)) {
    p <- p +
      stat_ellipse(
        data = df,
        aes(x = PC1, y = PC2, color = !!as.name(group)),
        level = ellipse_conf,
        linewidth = 1,
        linetype = "dashed",
        inherit.aes = FALSE
      )
  }

  # --- Labels ---
  if (isTRUE(names)) {
    p <- p +
      ggrepel::geom_text_repel(
        label = rownames(df),
        size = base_size / 4
      )
  }

  # --- Retour optionnel ---
  if (return_pca) {
    return(list(plot = p, pca = pca, manova_p = if (show_manova) pval else NULL))
  }

  return(p)
}
