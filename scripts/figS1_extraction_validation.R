# =============================================================================
# Supplementary Figure S1 for Royal Society journal  -  single panel
#   Agreement between AI-assisted figure extraction and author-provided raw data
#   (Table S1, "Author data comparison" tab)
#
# House style matches Figures 2-4: uniform 8 pt text, theme_bw, cairo output.
# Ratios (treatment / control) are log2-transformed to match Figure 2.
#
# install.packages(c("readxl","dplyr","ggplot2","scales","broom"))
# =============================================================================

library(readxl)
library(dplyr)
library(ggplot2)
library(scales)
library(broom)

# ---------------------------
# user settings
# ---------------------------
data_file  <- file.path("data", "Table_S1_RNAi.xlsx")
data_sheet <- "Author data comparison"   # header on row 4 -> skip = 3

png_file <- file.path("figures", "figS1_extraction_validation.png")
pdf_file <- file.path("figures", "figS1_extraction_validation.pdf")

# Royal Society single-column width
fig_width_in  <- 3.3     # 8.4 cm
fig_height_in <- 3.45
base_size     <- 8       # uniform printed point size

# shared theme fragment (same as Figure 2)
theme_fig <- function() {
  theme_bw(base_size = base_size) +
    theme(
      text             = element_text(size = base_size),
      plot.title       = element_text(size = base_size, face = "bold", hjust = 0),
      panel.grid.minor = element_blank(),
      axis.text        = element_text(size = base_size),
      axis.title       = element_text(size = base_size),
      legend.text      = element_text(size = base_size),
      legend.title     = element_text(size = base_size),
      legend.key.size  = unit(0.35, "lines"),
      legend.margin    = margin(t = 0, b = 0, l = 0, r = 0),
      plot.margin      = margin(3, 4, 3, 6)
    )
}

# =============================================================================
# DATA
# =============================================================================
cmp <- read_excel(data_file, sheet = data_sheet, skip = 3)

cmp <- cmp %>%
  transmute(
    Paper           = as.character(Paper),
    Assay           = as.character(Assay),
    Target          = as.character(Target),
    author_ratio    = as.numeric(`Author ratio`),
    extracted_ratio = as.numeric(`Extracted ratio`)
  ) %>%
  filter(!is.na(Paper), !is.na(author_ratio), !is.na(extracted_ratio),
         author_ratio > 0, extracted_ratio > 0) %>%
  mutate(
    log2_author    = log2(author_ratio),
    log2_extracted = log2(extracted_ratio),
    Assay = case_when(
      tolower(Assay) == "mrna expression" ~ "mRNA expression",
      TRUE ~ Assay
    ),
    Assay = factor(Assay, levels = c("mRNA expression", "Spore load"))
  )

# order papers by year, then author
paper_year <- suppressWarnings(as.integer(sub(".*\\s(\\d{4})$", "\\1", cmp$Paper)))
paper_levels <- unique(cmp$Paper[order(paper_year, cmp$Paper)])
cmp <- cmp %>% mutate(Paper = factor(Paper, levels = paper_levels))

# =============================================================================
# STATISTICS  (linear regression of extracted on author log2 ratios)
# =============================================================================
fit  <- lm(log2_extracted ~ log2_author, data = cmp)
g    <- glance(fit)
sl   <- tidy(fit) %>% filter(term == "log2_author")
pear <- cor.test(cmp$log2_author, cmp$log2_extracted, method = "pearson")
mean_abs_diff <- mean(abs(cmp$log2_extracted - cmp$log2_author))

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001",
                            paste0("p = ", format.pval(p, digits = 2)))
stats_label <- paste0(
  "R² = ", sprintf("%.3f", g$r.squared), ", ", fmt_p(sl$p.value), "\n",
  "slope = ", sprintf("%.2f", sl$estimate), "\n",
  "n = ", nrow(cmp)
)

cat("\n--- Figure S1 statistics ---\n")
cat("n pairs              :", nrow(cmp), "\n")
cat("n papers             :", n_distinct(cmp$Paper), "\n")
cat("R^2                  :", signif(g$r.squared, 4), "\n")
cat("slope (95% CI)       :", signif(sl$estimate, 4),
    paste0("(", paste(signif(confint(fit)["log2_author", ], 4), collapse = " to "), ")"), "\n")
cat("intercept            :", signif(coef(fit)[1], 4), "\n")
cat("slope p-value        :", signif(sl$p.value, 4), "\n")
cat("Pearson r            :", signif(pear$estimate, 4), "\n")
cat("mean |log2 diff|     :", signif(mean_abs_diff, 4), "\n\n")

# =============================================================================
# PLOT
# =============================================================================
paper_cols <- setNames(
  c("#4e79a7", "#e15759", "#59a14f", "#b07aa1", "#f28e2b", "#76b7b2")[seq_along(paper_levels)],
  paper_levels
)
assay_shapes <- c("mRNA expression" = 16, "Spore load" = 17)

# common square axis range so the 1:1 line runs corner to corner
rng <- range(c(cmp$log2_author, cmp$log2_extracted))
rng <- rng + c(-0.3, 0.3)

fig_S1 <- ggplot(cmp, aes(x = log2_author, y = log2_extracted)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.3,
              colour = "grey40") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.5, colour = "black") +
  geom_point(aes(colour = Paper, shape = Assay), size = 1.6, alpha = 0.9) +
  annotate("text", x = rng[1] + 0.1, y = rng[2] - 0.1, label = stats_label,
           hjust = 0, vjust = 1, size = (base_size - 1) / .pt,
           lineheight = 0.95) +
  scale_colour_manual(values = paper_cols, name = "Study") +
  scale_shape_manual(values = assay_shapes, name = "Assay") +
  scale_x_continuous(limits = rng, breaks = pretty_breaks(5)) +
  scale_y_continuous(limits = rng, breaks = pretty_breaks(5)) +
  coord_equal() +
  labs(x = "Author raw data\nFold change (log2)",
       y = "AI-assisted extraction\nFold change (log2)") +
  theme_fig() +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.spacing = unit(2, "pt"),
        legend.spacing.x = unit(4, "pt"),
        legend.justification = "center") +
  guides(colour = guide_legend(ncol = 2, order = 1, title.position = "top",
                               override.aes = list(size = 1.6)),
         shape  = guide_legend(ncol = 1, order = 2, title.position = "top",
                               override.aes = list(size = 1.6)))

ggsave(png_file, fig_S1, width = fig_width_in, height = fig_height_in,
       dpi = 600, type = "cairo")
ggsave(pdf_file, fig_S1, width = fig_width_in, height = fig_height_in,
       device = cairo_pdf)

print(fig_S1)
