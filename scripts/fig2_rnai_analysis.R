# =============================================================================
# Figure 2 for Royal Society journal  -  SIX panels, full two-column width
#   Panel A : Timeline of microsporidia RNAi papers          (publications.txt)
#   Panel B : RNAi targets by protein category               (rnai target script)
#   Panel C : Phenotype effect, best (most-reduced) timepoint (combined PLOT 4)
#   Panel D : mRNA vs phenotype, best timepoint               (combined PLOT 3)
#   Panel E : Off-target mRNA correlation                     (combined PLOT 5)
#   Panel F : Growth vs non-growth phenotypes examined        (combined PLOT 6)
#
# House style matches Figures 1-2: uniform 8 pt text, abbreviated species,
# compact legends, cairo output devices, patchwork composition.
#
# To keep the figure short, the original per-plot "Summary statistics" text
# strips are removed; R^2 / p are shown as small in-panel annotations instead.
#
# install.packages(c("readxl","dplyr","tidyr","ggplot2","scales","stringr",
#                    "forcats","broom","patchwork"))
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(stringr)
library(forcats)
library(broom)
library(patchwork)

# ---------------------------
# user settings
# ---------------------------
bib_file     <- file.path("data", "Table_S1_RNAi.xlsx")
target_sheet <- "Categories of RNAi targets"
strict_only <- TRUE

png_file <- file.path("figures", "fig2_rnai_analysis.png")
pdf_file <- file.path("figures", "fig2_rnai_analysis.pdf")

# Royal Society full two-column width
fig_width_in  <- 6.85    # 17.4 cm
fig_height_in <- 7.1     # three rows; tightened to minimise white space
base_size     <- 8       # uniform printed point size

# shared theme fragment so all panels match
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

# species abbreviation helper (same as Figs 1-2)
abbreviate_species <- function(x) {
  abbr_one <- function(s) {
    parts <- strsplit(s, "\\s+")[[1]]
    if (length(parts) >= 2 && nchar(parts[1]) > 1) {
      paste0(substr(parts[1], 1, 1), ". ", paste(parts[-1], collapse = " "))
    } else s
  }
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    if (grepl(" and ", s)) {
      pieces <- strsplit(s, " and ")[[1]]
      return(paste(vapply(trimws(pieces), abbr_one, character(1)), collapse = " & "))
    }
    if (grepl(" / |,", s)) return(s)
    abbr_one(s)
  }, character(1))
}

# consistent phenotype colour scale across C and D
phenotype_levels <- c("Spore number", "18S rRNA", "EhSWP1", "Beta-tubulin")
phenotype_cols   <- c(
  "Spore number"        = "#1b9e77",
  "18S rRNA"               = "#d95f02",
  "EhSWP1"              = "#7570b3",
  "Beta-tubulin" = "#e7298a"
)

# =============================================================================
# PANEL A : RNAi PUBLICATION TIMELINE
# =============================================================================
bib <- read_excel(bib_file, sheet = "Master_Bibliography", skip = 1)
if (strict_only && "Strict_Inclusion" %in% names(bib)) {
  bib <- bib %>% filter(Strict_Inclusion == "Yes")
}
bib <- bib %>% distinct(Title, Year, .keep_all = TRUE)

year_species <- bib %>%
  mutate(
    Year = suppressWarnings(as.integer(Year)),
    Parasite_Species = ifelse(is.na(Parasite_Species) | Parasite_Species == "",
                              "Unspecified", Parasite_Species)
  ) %>%
  filter(!is.na(Year)) %>%
  count(Year, Parasite_Species, name = "n")

all_years <- seq(min(year_species$Year, na.rm = TRUE),
                 max(year_species$Year, na.rm = TRUE))
all_species <- sort(unique(year_species$Parasite_Species))

year_species_full <- expand_grid(Year = all_years, Parasite_Species = all_species) %>%
  left_join(year_species, by = c("Year", "Parasite_Species")) %>%
  mutate(n = replace_na(n, 0))

year_totals <- year_species_full %>%
  group_by(Year) %>%
  summarise(papers_per_year = sum(n), .groups = "drop") %>%
  arrange(Year) %>%
  mutate(cumulative_papers = cumsum(papers_per_year))

scale_factor <- max(year_totals$papers_per_year) / max(year_totals$cumulative_papers)
year_totals  <- year_totals %>% mutate(cumulative_scaled = cumulative_papers * scale_factor)

year_species_full <- year_species_full %>%
  mutate(Parasite_Species = abbreviate_species(Parasite_Species))
species_order <- year_species_full %>%
  group_by(Parasite_Species) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>% pull(Parasite_Species)
year_species_full <- year_species_full %>%
  mutate(Parasite_Species = factor(Parasite_Species, levels = species_order))

panel_A <- ggplot() +
  geom_col(data = year_species_full,
           aes(x = Year, y = n, fill = Parasite_Species), width = 0.8) +
  geom_line(data = year_totals, aes(x = Year, y = cumulative_scaled, group = 1),
            linewidth = 0.5) +
  geom_point(data = year_totals, aes(x = Year, y = cumulative_scaled), size = 0.9) +
  scale_x_continuous(breaks = all_years[all_years %% 4 == 0]) +
  scale_y_continuous(
    name = "Papers per year",
    breaks = scales::pretty_breaks(4),
    sec.axis = sec_axis(trans = ~ . / scale_factor, name = "Cumulative"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(title = "RNAi papers over time", x = "Year", fill = "Species") +
  theme_fig() +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = base_size - 1),
    legend.position = "right",
    legend.text     = element_text(face = "italic", size = base_size - 1)
  ) +
  guides(fill = guide_legend(ncol = 1, title.position = "top"))

# =============================================================================
# PANEL B : RNAi TARGETS BY PROTEIN CATEGORY
# =============================================================================
tgt <- read_excel(bib_file, sheet = target_sheet)
count_df <- tgt %>%
  rename(Target = 1, Category = 2) %>%
  mutate(Target = trimws(as.character(Target)),
         Category = trimws(as.character(Category))) %>%
  filter(!is.na(Category), Category != "") %>%
  mutate(Category = ifelse(tolower(Category) == "metabolic enzyme",
                           "Metabolic enzyme", Category)) %>%
  count(Category, name = "Count") %>%
  arrange(desc(Count), Category)
count_df$Category <- factor(count_df$Category, levels = rev(count_df$Category))

panel_B <- ggplot(count_df, aes(x = Count, y = Category)) +
  geom_col(width = 0.72, fill = "grey55") +
  geom_text(aes(label = Count), hjust = -0.25, size = (base_size - 2) / .pt) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.20)),
                     breaks = scales::pretty_breaks(4)) +
  labs(title = "RNAi targets by protein category",
       x = "Number of targets", y = NULL) +
  theme_fig() +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = base_size - 1),
        axis.title.x = element_text(margin = margin(t = -6)))

# =============================================================================
# SHARED DATA PREP for Panels C, D, E (from combined_scatter_and_bar_6.R)
# =============================================================================
dat <- read_excel(bib_file, sheet = "Quant_Timepoint_Data")
dat_on_target <- dat %>% filter(Is_On_Target == TRUE)

clean_assay_label <- function(x) {
  x <- as.character(x); x <- str_trim(x)
  str_replace(x, regex("\\s*qPCR$", ignore_case = TRUE), "")
}
assign_phenotype_label <- function(measurement_type, assay) {
  assay_clean <- clean_assay_label(assay)
  case_when(
    measurement_type == "spore_number" ~ "Spore number",
    measurement_type == "beta_tubulin_copy_number" ~ "Beta-tubulin",
    measurement_type %in% c("ssRNA_copy_number", "ssRNA_relative_expression") &
      str_detect(assay_clean, regex("^EhSWP1$", ignore_case = TRUE)) ~ "EhSWP1",
    measurement_type %in% c("ssRNA_copy_number", "ssRNA_relative_expression") ~ "18S rRNA",
    TRUE ~ measurement_type %>% str_replace_all("_", " ") %>% str_to_sentence()
  )
}
lm_r2_p <- function(yvar, xvar, df) {
  if (nrow(df) >= 2 && n_distinct(df[[xvar]]) >= 2) {
    fit <- lm(reformulate(xvar, yvar), data = df)
    g <- glance(fit); t <- tidy(fit)
    sl <- t %>% filter(term == xvar)
    list(r2 = g$r.squared[1], p = if (nrow(sl) == 1) sl$p.value[1] else NA_real_)
  } else list(r2 = NA_real_, p = NA_real_)
}
fmt_r2p <- function(r2, p) {
  r2s <- ifelse(is.na(r2), "NA", sprintf("%.2f", r2))
  ps  <- ifelse(is.na(p), "NA",
                ifelse(p < 0.001, "p < 0.001",
                       paste0("p = ", format.pval(p, digits = 2))))
  paste0("R\u00B2 = ", r2s, ", ", ps)
}

# Per-group regression stats, returning ONE row per group with a formatted
# label and whether the slope is significant (used to keep only significant
# correlations and to colour each label by its group).
SIG_ALPHA <- 0.05
group_sig_stats <- function(df, group_var, yvar, xvar) {
  df %>%
    group_by(.data[[group_var]]) %>%
    group_modify(~ {
      s <- lm_r2_p(yvar, xvar, .x)
      tibble(r2 = s$r2, p = s$p, n = nrow(.x))
    }) %>%
    ungroup() %>%
    rename(grp = !!group_var) %>%
    filter(!is.na(p), p < SIG_ALPHA) %>%
    arrange(p) %>%
    mutate(label = paste0(grp, ":\n", fmt_r2p(r2, p)))
}

mrna_df <- dat_on_target %>%
  filter(Measurement_Type == "mRNA_expression") %>%
  select(Paper_ID, First_Author, Year, Short_Title, Parasite_Species,
         Gene, Time_h, mrna_ratio = Target_over_Control_Ratio)

phenotype_df <- dat_on_target %>%
  filter(Measurement_Type != "mRNA_expression") %>%
  mutate(Phenotype = assign_phenotype_label(Measurement_Type, Assay)) %>%
  select(Paper_ID, First_Author, Year, Short_Title, Parasite_Species,
         Gene, Time_h, Measurement_Type, Assay, Phenotype,
         phenotype_ratio = Target_over_Control_Ratio) %>%
  filter(!is.na(phenotype_ratio), phenotype_ratio > 0)

known_order  <- c("Spore number", "18S rRNA", "EhSWP1", "Beta-tubulin")
extra_levels <- setdiff(unique(phenotype_df$Phenotype), known_order)
all_levels   <- c(known_order, sort(extra_levels))
phenotype_df <- phenotype_df %>% mutate(Phenotype = factor(Phenotype, levels = all_levels))

join_cols <- c("Paper_ID", "First_Author", "Year", "Short_Title",
               "Parasite_Species", "Gene", "Time_h")

# best (most-reduced) timepoint per paper/gene/phenotype
peak_phenotype_df <- phenotype_df %>%
  mutate(log2_phenotype = log2(phenotype_ratio)) %>%
  arrange(Paper_ID, Gene, Phenotype, phenotype_ratio, Time_h) %>%
  group_by(Paper_ID, First_Author, Year, Short_Title, Parasite_Species, Gene, Phenotype) %>%
  slice(1) %>% ungroup() %>%
  mutate(phenotype_percent = phenotype_ratio * 100,
         Phenotype = factor(Phenotype, levels = all_levels))

# =============================================================================
# PANEL C : phenotype effect at best timepoint (bar + jitter)   [PLOT 4]
# =============================================================================
bar_peak_df <- peak_phenotype_df %>% mutate(Phenotype = factor(Phenotype, levels = all_levels))
bar_peak_summary <- bar_peak_df %>%
  group_by(Phenotype) %>%
  summarise(mean_log2 = mean(log2_phenotype, na.rm = TRUE),
            sd_log2 = sd(log2_phenotype, na.rm = TRUE),
            n_points = n(), n_genes = n_distinct(Gene),
            n_papers = n_distinct(Paper_ID), .groups = "drop") %>%
  mutate(sd_log2 = ifelse(is.na(sd_log2), 0, sd_log2),
         ymin = mean_log2 - sd_log2, ymax = mean_log2 + sd_log2) %>%
  filter(!is.na(Phenotype))

bar_peak_df  <- bar_peak_df  %>% filter(Phenotype %in% bar_peak_summary$Phenotype) %>% droplevels()
bar_peak_summary <- bar_peak_summary %>% droplevels()

panel_C <- ggplot() +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_col(data = bar_peak_summary,
           aes(x = Phenotype, y = mean_log2, fill = Phenotype),
           width = 0.62, alpha = 0.7, show.legend = FALSE) +
  geom_errorbar(data = bar_peak_summary,
                aes(x = Phenotype, ymin = ymin, ymax = ymax, color = Phenotype),
                width = 0.15, linewidth = 0.5, show.legend = FALSE) +
  geom_jitter(data = bar_peak_df,
              aes(x = Phenotype, y = log2_phenotype, color = Phenotype),
              width = 0.12, height = 0, size = 1.1, alpha = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = phenotype_cols, drop = FALSE) +
  scale_color_manual(values = phenotype_cols, drop = FALSE) +
  labs(title = "Growth phenotype effect",
       x = NULL, y = "Growth phenotype\nFold change (log2)") +
  theme_fig() +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = base_size - 1))

# =============================================================================
# PANEL D : mRNA vs phenotype at best timepoint (scatter)        [PLOT 3]
# =============================================================================
scatter_peak_df <- mrna_df %>%
  inner_join(peak_phenotype_df %>%
               select(Paper_ID, First_Author, Year, Short_Title, Parasite_Species,
                      Gene, Time_h, Phenotype, phenotype_ratio, log2_phenotype),
             by = join_cols) %>%
  filter(mrna_ratio > 0, phenotype_ratio > 0) %>%
  mutate(log2_mrna = log2(mrna_ratio),
         Phenotype = factor(Phenotype, levels = all_levels))

d_sig <- group_sig_stats(scatter_peak_df, "Phenotype", "log2_phenotype", "log2_mrna")
# stacked label positions (top-left); stats shown only for significant groups
d_sig <- d_sig %>% mutate(
  Phenotype = factor(grp, levels = all_levels),
  vjust_i   = 1.6 + (row_number() - 1) * 2.6
)

# Dummy rows (fully transparent) guarantee every phenotype appears in the
# legend even when a level has no real points (e.g. EhSWP1 here). This is
# version-proof, unlike relying on scale limits to surface unused levels.
d_dummy <- tibble(
  log2_mrna = mean(scatter_peak_df$log2_mrna, na.rm = TRUE),
  log2_phenotype = mean(scatter_peak_df$log2_phenotype, na.rm = TRUE),
  Phenotype = factor(phenotype_levels, levels = all_levels)
)

panel_D <- ggplot(scatter_peak_df,
                  aes(x = log2_mrna, y = log2_phenotype, color = Phenotype)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  # invisible points to force all phenotype keys into the legend
  geom_point(data = d_dummy, alpha = 0, size = 1.2, show.legend = TRUE) +
  geom_point(size = 1.2, alpha = 0.85) +
  # regression line for ALL phenotypes (groups with too few points are skipped)
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  # coloured stats text ONLY for significant phenotypes
  geom_text(data = d_sig, inherit.aes = FALSE,
            aes(x = -Inf, y = Inf, label = label, color = Phenotype, vjust = vjust_i),
            hjust = -0.03, size = (base_size - 3) / .pt, show.legend = FALSE) +
  scale_color_manual(values = phenotype_cols, drop = FALSE, name = "Growth phenotype",
                     limits = phenotype_levels, breaks = phenotype_levels) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  labs(title = "Target mRNA vs growth phenotype",
       x = "Target mRNA Fold change (log2)",
       y = "Growth phenotype\nFold change (log2)") +
  theme_fig() +
  theme(legend.position = "bottom",
        legend.box.spacing = unit(2, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.key.size = unit(0.3, "lines"),
        legend.title = element_text(size = base_size - 1),
        legend.text  = element_text(size = base_size - 1),
        axis.title.x = element_text(margin = margin(t = -10))) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                              override.aes = list(size = 1.6, alpha = 1),
                              title.position = "left"))

# =============================================================================
# PANEL E : off-target mRNA correlation (scatter)                [PLOT 5]
# =============================================================================
off_join <- c("Paper_ID", "First_Author", "Year", "Short_Title",
              "Parasite_Species", "Time_h")
knockdown_target_mrna_df <- dat %>%
  filter(Is_On_Target == TRUE, Measurement_Type == "mRNA_expression") %>%
  select(Paper_ID, First_Author, Year, Short_Title, Parasite_Species,
         Time_h, Knockdown_Gene = Gene, knockdown_mrna_ratio = Target_over_Control_Ratio)
off_target_mrna_df <- dat %>%
  filter(Is_On_Target == FALSE, Measurement_Type == "mRNA_expression") %>%
  select(Paper_ID, First_Author, Year, Short_Title, Parasite_Species,
         Time_h, Off_Target_Gene = Gene, off_target_mrna_ratio = Target_over_Control_Ratio)

off_target_scatter_df <- knockdown_target_mrna_df %>%
  inner_join(off_target_mrna_df, by = off_join) %>%
  filter(knockdown_mrna_ratio > 0, off_target_mrna_ratio > 0,
         Knockdown_Gene != Off_Target_Gene) %>%
  mutate(log2_knockdown_mrna = log2(knockdown_mrna_ratio),
         log2_off_target_mrna = log2(off_target_mrna_ratio),
         Knockdown_Gene = factor(Knockdown_Gene, levels = sort(unique(Knockdown_Gene))))

# fixed colour map for knockdown genes so the stats labels can match the points
kd_genes <- levels(off_target_scatter_df$Knockdown_Gene)
kd_cols  <- setNames(
  c("#e15759", "#b07aa1", "#59a14f", "#4e79a7", "#e377c2",
    "#9c755f", "#edc948", "#76b7b2")[seq_along(kd_genes)],
  kd_genes
)

e_sig <- group_sig_stats(off_target_scatter_df, "Knockdown_Gene",
                         "log2_off_target_mrna", "log2_knockdown_mrna")
e_sig <- e_sig %>% mutate(
  Knockdown_Gene = factor(grp, levels = kd_genes),
  vjust_i = 1.3 + (row_number() - 1) * 1.9
)

panel_E <- ggplot(off_target_scatter_df,
                  aes(x = log2_knockdown_mrna, y = log2_off_target_mrna,
                      color = Knockdown_Gene)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_point(size = 1.2, alpha = 0.85) +
  # regression line for ALL knockdown genes (groups with too few points skipped)
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
  # coloured stats text ONLY for significant knockdown genes
  geom_text(data = e_sig, inherit.aes = FALSE,
            aes(x = -Inf, y = Inf, label = label, color = Knockdown_Gene, vjust = vjust_i),
            hjust = -0.03, size = (base_size - 3) / .pt, show.legend = FALSE) +
  scale_color_manual(values = kd_cols, name = "Knockdown target") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(title = "Off-target expression",
       x = "Target mRNA\nFold change (log2)",
       y = "Off-target mRNA\nFold change (log2)") +
  theme_fig() +
  theme(legend.position = "right",
        legend.box.spacing = unit(2, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        legend.text = element_text(size = base_size - 1)) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 1.6),
                              title.position = "top"))

# =============================================================================
# PANEL F : growth vs non-growth phenotypes examined             [PLOT 6]
# =============================================================================
master_bib <- read_excel(bib_file, sheet = "Master_Bibliography", skip = 1) %>%
  filter(!is.na(Ref_No), !is.na(Target))

phenotype_measure_df <- master_bib %>%
  select(Ref_No, First_Author, Year, Target, Phenotype_Assay) %>%
  filter(!is.na(Phenotype_Assay)) %>%
  separate_rows(Phenotype_Assay, sep = ";") %>%
  mutate(
    Phenotype_Assay = str_squish(Phenotype_Assay),
    # Target is one cell per paper in Table S1, so each ID corresponds to one paper
    Study_Target_ID = paste(Ref_No, Target, sep = "__"),
    Phenotype_Assay_Std = case_when(
      str_to_lower(Phenotype_Assay) %in% c("spores", "spore number") ~ "spore number",
      str_to_lower(Phenotype_Assay) %in% c("beta-tubulin", "beta-tubulin copies") ~ "beta-tubulin",
      str_to_lower(Phenotype_Assay) == "ssrna" ~ "18S rRNA",
      Phenotype_Assay == "EhSWP1" ~ "EhSWP1",
      TRUE ~ Phenotype_Assay
    ),
    Phenotype_Class = case_when(
      Phenotype_Assay_Std %in% c("spore number", "beta-tubulin", "18S rRNA", "EhSWP1") ~ "Growth",
      str_to_lower(Phenotype_Assay_Std) %in% c("none", "na", "") ~ NA_character_,
      TRUE ~ "Non-growth"
    )
  ) %>%
  filter(!is.na(Phenotype_Class))

phenotype_measure_summary <- phenotype_measure_df %>%
  group_by(Phenotype_Class, Phenotype_Assay = Phenotype_Assay_Std) %>%
  summarise(n_targets = n_distinct(Study_Target_ID), .groups = "drop") %>%
  # merge related phenotypes into broader categories
  mutate(Phenotype_Assay = case_when(
    str_to_lower(Phenotype_Assay) %in% c(
      "host prophenoloxidase activation", "host immunity",
      "host immune expression") ~ "Host immunity",
    str_detect(str_to_lower(Phenotype_Assay), "actin") &
      str_detect(str_to_lower(Phenotype_Assay), "localization|localisation|loc\\.") ~
      "Cytoskeleton protein localization",
    str_detect(str_to_lower(Phenotype_Assay), "^microsporidia actin") ~
      "Cytoskeleton protein localization",
    str_detect(str_to_lower(Phenotype_Assay), "actin") &
      str_detect(str_to_lower(Phenotype_Assay), "tubulin") ~
      "Cytoskeleton protein localization",
    TRUE ~ Phenotype_Assay
  )) %>%
  # re-aggregate after merging so combined categories sum their counts
  group_by(Phenotype_Class, Phenotype_Assay) %>%
  summarise(n_targets = sum(n_targets), .groups = "drop") %>%
  # capitalize the first letter of every category name, but keep established
  # gene/RNA casing (18S rRNA, miRNA, EhSWP1) intact
  mutate(Phenotype_Assay = ifelse(
    str_detect(Phenotype_Assay, "^(18S rRNA|ssRNA|miRNA|EhSWP1)"),
    Phenotype_Assay,
    str_replace(Phenotype_Assay, "^(\\w)", function(m) toupper(m))
  )) %>%
  arrange(Phenotype_Class, n_targets, Phenotype_Assay) %>%
  mutate(Phenotype_Assay = factor(Phenotype_Assay, levels = unique(Phenotype_Assay)))

class_cols <- c("Growth" = "#e15759", "Non-growth" = "#4db6ac")

panel_F <- ggplot(phenotype_measure_summary,
                  aes(x = Phenotype_Assay, y = n_targets, fill = Phenotype_Class)) +
  geom_col(width = 0.7, alpha = 0.85, show.legend = FALSE) +
  facet_grid(Phenotype_Class ~ ., scales = "free_y", space = "free_y") +
  coord_flip(clip = "off") +
  scale_fill_manual(values = class_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     breaks = scales::pretty_breaks(4)) +
  labs(title = "Phenotypes examined", x = NULL, y = "Number of papers") +
  theme_fig() +
  theme(
    panel.grid.major.y = element_blank(),
    strip.background   = element_rect(fill = "grey90", color = NA),
    strip.text.y       = element_text(size = base_size - 1, angle = 0),
    axis.text.y        = element_text(size = base_size - 2)
  )

# =============================================================================
# COMPOSE: 3 rows x 2 cols  (A B / C D / E F)
# =============================================================================
combined <- (panel_A | panel_B) /
            (panel_C | panel_D) /
            (panel_E | panel_F) +
  plot_layout(heights = c(1, 1, 1.25)) +   # bottom row taller for Panel F
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = base_size + 1))

ggsave(png_file, combined, width = fig_width_in, height = fig_height_in,
       dpi = 600, type = "cairo")
ggsave(pdf_file, combined, width = fig_width_in, height = fig_height_in,
       device = cairo_pdf)

print(combined)