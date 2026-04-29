library(tibble)

outdir <- "Dynamic"
q_levels <- c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6")
long_list <- list()
for(q in q_levels){
  infile <- file.path(outdir, paste0(q, ".inv.residual.tsv"))
  mat_dt <- fread(infile)
  mat <- as.data.frame(mat_dt)
  apatag <- mat[[1]]
  mat[[1]] <- NULL

  dt_long <- as.data.table(mat)
  dt_long[, APATag := apatag]
  dt_long <- melt(dt_long, id.vars = "APATag", variable.name = "inv", value.name = "value")
  dt_long[, Q := q]

  long_list[[q]] <- dt_long
}
dt_dyn <- rbindlist(long_list)

dt_dyn[, Stage_num := as.numeric(factor(Q, levels = paste0("Q",1:6)))]
test_trend_lm <- function(df){
  model_linear <- lm(value ~ Stage_num, data = df)
  model_quad <- lm(value ~ Stage_num + I(Stage_num^2), data = df)

  p_quad <- anova(model_linear, model_quad)[["Pr(>F)"]][2]
  p_linear <- summary(model_linear)$coefficients["Stage_num", "Pr(>|t|)"]

  data.table(p_linear = p_linear, p_quad = p_quad)
}

res_trend <- dt_dyn[, test_trend_lm(.SD), by = APATag]
res_trend <- res_trend %>%
    mutate(dynamic_type = case_when(
      !is.na(p_linear) & p_linear < 0.05 & !is.na(p_quad) & p_quad < 0.05 ~ "Linear & Quadratic",
      !is.na(p_quad) & p_quad < 0.05 ~ "Quadratic",
      !is.na(p_linear) & p_linear < 0.05 ~ "Linear",
      TRUE ~ "Static"
    ))
res_trend <- res_trend %>%
    mutate(
      fdr_linear = p.adjust(p_linear, method = "BH"),
      fdr_quad   = p.adjust(p_quad, method = "BH")
    )
res_trend <- res_trend %>%
    mutate(dynamic_type = case_when(
      !is.na(fdr_linear) & fdr_linear < 0.05 & !is.na(fdr_quad) & fdr_quad < 0.05 ~ "Linear & Quadratic",
      !is.na(fdr_quad) & fdr_quad < 0.05 ~ "Quadratic",
      !is.na(fdr_linear) & fdr_linear < 0.05 ~ "Linear",
      TRUE ~ "Static"
    ))
table(res_trend$dynamic_type)

fwrite(res_trend, "dynamic_PAS_lm.res", sep="\t")


# plot
pas_mean <- fread("mean.data.txt")
res_trend <- fread("dynamic_PAS_lm.res")
res_trend <- subset(res_trend, dynamic_type != "Static")
table(res_trend$dynamic_type)
sigModel <- setNames(res_trend[["dynamic_type"]], res_trend[["APATag"]])

dataPlot = copy(pas_mean) %>%
  column_to_rownames("APATag") %>%
  as.matrix() %>%
  .[names(sigModel), ] %>%
  t() %>%
  scale() %>%
  t()

method_vec <- sigModel[rownames(dataPlot)]
method_vec <- factor(
  method_vec,
  levels = c("Linear", "Quadratic", "Linear & Quadratic")
)
method_col <- c(
  "Linear" = "#1f78b4",
  "Quadratic" = "#e31a1c",
  "Linear & Quadratic" = "#ff7f00"
)

row_ha <- rowAnnotation(
  Method = method_vec,
  col = list(Method = method_col),
  show_annotation_name = FALSE,
  width = unit(3, "mm")
)

colorHeatmap <- list(low = "#5E3C99", mid = "white", high = "#BA2E1F")
p <- Heatmap(
  dataPlot,
  #  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  col=circlize::colorRamp2(c(min(dataPlot, na.rm=T), 0, max(dataPlot, na.rm=T)),
                           c(colorHeatmap[["low"]], colorHeatmap[["mid"]], colorHeatmap[["high"]])),
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_row_dend = FALSE,
  show_column_names = TRUE,
  row_names_rot  = 90,
  column_title = "Dynamic PAS usage (3'UTR)",
  row_title = NULL,
  heatmap_legend_param = list(title = "PAS usage"),
  left_annotation = row_ha,
  row_split = method_vec,
  gap = unit(2, "mm")
)
