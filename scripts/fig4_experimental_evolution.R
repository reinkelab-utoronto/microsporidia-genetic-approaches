# =============================================================================
# Figure 4 for Royal Society journal (full two-column width)
#   Panel A : Timeline of microsporidia experimental evolution papers
#   Panel B : Selection phenotype properties chart
#
# Royal Society single-column width = 8.4 cm (~3.3 in).
# Panels are stacked vertically (A over B) and labelled A / B.
#
# install.packages(c("readxl","dplyr","tidyr","ggplot2","scales",
#                    "stringr","forcats","janitor","cellranger","patchwork"))
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(stringr)
library(forcats)
# janitor is only used for clean_names(); fall back to a base-R version if absent
if (requireNamespace("janitor", quietly = TRUE)) {
  library(janitor)
} else {
  clean_names <- function(df) {
    nm <- names(df)
    nm <- tolower(trimws(nm))
    nm <- gsub("[^a-z0-9]+", "_", nm)
    nm <- gsub("^_+|_+$", "", nm)
    names(df) <- make.unique(nm, sep = "_")
    df
  }
}
library(cellranger)
library(patchwork)
library(ggtext)

# ---------------------------
# user settings
# ---------------------------
file_path  <- file.path("data", "Table_S3_experimental_evolution.xlsx")
sheet_name <- "Experimental_Evolution_Biblio"

# Table S3: title on row 1, column headers on row 2
# Panel A reads with skip = 1; Panel B reads with header_row = 2
header_row  <- 2
strict_only <- FALSE   # set TRUE to keep only Strict_Inclusion == "Yes"

png_file <- file.path("figures", "fig4_experimental_evolution.png")
pdf_file <- file.path("figures", "fig4_experimental_evolution.pdf")

# Royal Society FULL two-column dimensions (~17.4 cm = 6.85 in wide).
# NOTE: the RStudio plot preview rescales to fit your screen, so it is NOT a
# reliable proxy for the exported file. Always judge by the saved PNG/PDF.
#
# Panels are placed side by side (A left, B right). Legends sit BELOW each
# panel so they do not steal horizontal width and squash the plots.
# All text is forced to a uniform 8 pt.
fig_width_in  <- 6.85   # 17.4 cm full two-column width
fig_height_in <- 3.4    # single row; legends sit to the right of each panel
base_size     <- 8      # uniform printed point size for ALL text

# =============================================================================
# PANEL A : TIMELINE
# =============================================================================

bib <- read_excel(
  file_path,
  sheet = sheet_name,
  skip  = 1
)

# optional strict filter
if (strict_only && "Strict_Inclusion" %in% names(bib)) {
  bib <- bib %>% filter(Strict_Inclusion == "Yes")
}

# deduplicate by Title + Year
bib <- bib %>% distinct(Title, Year, .keep_all = TRUE)

# counts by year and species
year_species <- bib %>%
  mutate(
    Year = suppressWarnings(as.integer(Year)),
    Parasite_Species = ifelse(
      is.na(Parasite_Species) | Parasite_Species == "",
      "Unspecified",
      Parasite_Species
    )
  ) %>%
  filter(!is.na(Year)) %>%
  count(Year, Parasite_Species, name = "n")

# fill missing year/species combos with zero
all_years <- seq(
  min(year_species$Year, na.rm = TRUE),
  max(year_species$Year, na.rm = TRUE)
)
all_species <- sort(unique(year_species$Parasite_Species))

year_species_full <- expand_grid(
  Year = all_years,
  Parasite_Species = all_species
) %>%
  left_join(year_species, by = c("Year", "Parasite_Species")) %>%
  mutate(n = replace_na(n, 0))

# yearly totals and cumulative
year_totals <- year_species_full %>%
  group_by(Year) %>%
  summarise(papers_per_year = sum(n), .groups = "drop") %>%
  arrange(Year) %>%
  mutate(cumulative_papers = cumsum(papers_per_year))

# dual-axis scale factor
scale_factor <- max(year_totals$papers_per_year) / max(year_totals$cumulative_papers)
year_totals <- year_totals %>%
  mutate(cumulative_scaled = cumulative_papers * scale_factor)

# Abbreviate species to "G. species" form for a compact legend.
# Multi-species / compound entries are left intact.
abbreviate_species <- function(x) {
  abbr_one <- function(s) {
    parts <- strsplit(s, "\\s+")[[1]]
    if (length(parts) >= 2 && nchar(parts[1]) > 1) {
      paste0(substr(parts[1], 1, 1), ". ", paste(parts[-1], collapse = " "))
    } else s
  }
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    # compound entries like "Genus species and Genus2 species2": abbreviate each
    if (grepl(" and ", s)) {
      pieces <- strsplit(s, " and ")[[1]]
      return(paste(vapply(trimws(pieces), abbr_one, character(1)), collapse = " & "))
    }
    if (grepl(" / |,", s)) return(s)
    abbr_one(s)
  }, character(1))
}

year_species_full <- year_species_full %>%
  mutate(Parasite_Species = abbreviate_species(Parasite_Species))

# order species by overall abundance
species_order <- year_species_full %>%
  group_by(Parasite_Species) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(Parasite_Species)

year_species_full <- year_species_full %>%
  mutate(Parasite_Species = factor(Parasite_Species, levels = species_order))

panel_A <- ggplot() +
  geom_col(
    data = year_species_full,
    aes(x = Year, y = n, fill = Parasite_Species),
    width = 0.8
  ) +
  geom_line(
    data = year_totals,
    aes(x = Year, y = cumulative_scaled, group = 1),
    linewidth = 0.6
  ) +
  geom_point(
    data = year_totals,
    aes(x = Year, y = cumulative_scaled),
    size = 1.1
  ) +
  scale_x_continuous(
    # decade labels avoid overlap now that the panel shares width with a legend
    breaks = all_years[all_years %% 10 == 0]
  ) +
  scale_y_continuous(
    name = "Papers per year",
    breaks = seq(0, ceiling(max(year_totals$papers_per_year)), by = 1),
    sec.axis = sec_axis(
      trans = ~ . / scale_factor,
      name = "Cumulative papers"
    ),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Experimental evolution papers over time",
    x = "Year",
    fill = "Microsporidia species"
  ) +
  theme_bw(base_size = base_size) +
  theme(
    text             = element_text(size = base_size),
    plot.title       = element_text(size = base_size, face = "bold", hjust = 0),
    panel.grid.minor = element_blank(),
    axis.text        = element_text(size = base_size),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = base_size),
    axis.title       = element_text(size = base_size),
    legend.position  = "right",
    legend.box       = "vertical",
    legend.key.size  = unit(0.4, "lines"),
    legend.text      = element_text(size = base_size, face = "italic"),
    legend.title     = element_text(size = base_size),
    legend.margin    = margin(t = 0, b = 0, l = 0, r = 0),
    plot.margin      = margin(2, 4, 2, 2)
  ) +
  guides(fill = guide_legend(ncol = 1, byrow = TRUE, title.position = "top"))

# =============================================================================
# PANEL B : SELECTION PHENOTYPE PROPERTIES CHART
# =============================================================================

clean_blank <- function(x) {
  x <- as.character(x)
  x <- str_squish(x)
  x[x %in% c("", "NA", "Na", "N/A", "na")] <- NA_character_
  x
}

standardize_assay <- function(x) {
  x <- clean_blank(x)
  x_low <- str_to_lower(coalesce(x, ""))
  case_when(
    str_detect(x_low, "phenotype") & str_detect(x_low, "geno|genom") ~ "Phenotype + genotype",
    str_detect(x_low, "transcript") ~ "Transcriptome",
    str_detect(x_low, "geno|genom") ~ "Genotype",
    str_detect(x_low, "phenotype") ~ "Phenotype",
    str_detect(x_low, "^none$|not assayed|no assay") ~ "None",
    is.na(x) ~ "Unspecified",
    TRUE ~ str_to_sentence(x)
  )
}

standardize_change <- function(x) {
  x <- clean_blank(x)
  x_low <- str_to_lower(coalesce(x, ""))
  case_when(
    str_detect(x_low, "^yes$|change|changed|observed") ~ "Yes",
    str_detect(x_low, "^no$|none|not observed|unchanged") ~ "No",
    is.na(x) ~ "Unspecified",
    TRUE ~ str_to_sentence(x)
  )
}

# read workbook robustly
raw <- read_excel(
  path = file_path,
  sheet = sheet_name,
  col_names = FALSE,
  range = cell_cols(1:30)
)

header_values <- raw %>%
  slice(header_row) %>%
  unlist(use.names = FALSE) %>%
  as.character()

last_col <- max(which(!is.na(header_values) & str_squish(header_values) != ""))
headers  <- header_values[1:last_col]

dat <- raw %>%
  slice((header_row + 1):n()) %>%
  select(1:last_col)

names(dat) <- headers
dat <- dat %>% clean_names()

# Table S3 names this column "Selection condition"; older versions used
# "Selection phenotype". Standardize to one name.
if ("selection_condition" %in% names(dat) && !"selection_phenotype" %in% names(dat)) {
  dat <- dat %>% rename(selection_phenotype = selection_condition)
}

# drop fully empty rows
if (nrow(dat) > 0) {
  dat <- dat %>%
    filter(!if_all(everything(), ~ is.na(.x) | str_squish(as.character(.x)) == ""))
}

# keep real entries
if ("ref_no" %in% names(dat)) {
  dat <- dat %>% filter(!is.na(ref_no))
} else if ("title" %in% names(dat)) {
  dat <- dat %>% filter(!is.na(title) & str_squish(title) != "")
}

# remove duplicates
if ("duplicate_or_possible_match" %in% names(dat)) {
  dat <- dat %>%
    filter(is.na(duplicate_or_possible_match) | str_squish(as.character(duplicate_or_possible_match)) == "")
}

# strict inclusion only (optional)
if (strict_only && "strict_inclusion" %in% names(dat)) {
  dat <- dat %>% filter(strict_inclusion == "Yes")
}

# required columns
required_cols <- c("selection_phenotype", "assay_after_passage", "change_observed")
missing_cols  <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

# clean / standardize
plot_df <- dat %>%
  mutate(
    selection_phenotype = clean_blank(selection_phenotype),
    assay_after_passage = standardize_assay(assay_after_passage),
    change_observed     = standardize_change(change_observed)
  ) %>%
  filter(!is.na(selection_phenotype)) %>%
  # wrap long phenotype labels so they need less horizontal room in a 2-col layout
  mutate(selection_phenotype = str_wrap(selection_phenotype, width = 20))

# bar data: total per phenotype
bar_df <- plot_df %>%
  count(selection_phenotype, name = "total") %>%
  mutate(selection_phenotype = fct_reorder(selection_phenotype, total))

# point data: one point per study
point_df <- plot_df %>%
  mutate(
    selection_phenotype = factor(selection_phenotype, levels = levels(bar_df$selection_phenotype)),
    change_observed     = factor(change_observed,     levels = c("Yes", "No", "Unspecified")),
    assay_after_passage = factor(assay_after_passage)
  ) %>%
  left_join(bar_df, by = "selection_phenotype") %>%
  group_by(selection_phenotype) %>%
  mutate(jitter_x = total + 0.4 + (row_number() - 1) * 0.55) %>%
  ungroup()

# shapes / colours
assay_levels  <- levels(point_df$assay_after_passage)
shape_palette <- c(21, 22, 23, 24, 25, 8, 3, 4, 7, 11)[seq_along(assay_levels)]
names(shape_palette) <- assay_levels

change_colours <- c("Yes" = "#1b9e77", "No" = "#d95f02", "Unspecified" = "grey60")
# colours used to tint the legend label text for Change observed (Yes/No only)
change_label_cols <- c("#1b9e77", "#d95f02")

x_max <- max(point_df$jitter_x, na.rm = TRUE) + 1

panel_B <- ggplot() +
  geom_col(
    data  = bar_df,
    aes(x = total, y = selection_phenotype),
    fill  = "grey85",
    color = "black",
    linewidth = 0.3,
    width = 0.6
  ) +
  geom_text(
    data = bar_df,
    aes(x = total / 2, y = selection_phenotype, label = total),
    size = 2.2,
    color = "grey30"
  ) +
  geom_point(
    data  = point_df,
    aes(
      x     = jitter_x,
      y     = selection_phenotype,
      shape = assay_after_passage,
      fill  = change_observed,
      color = change_observed
    ),
    size   = 1.8,
    stroke = 0.4,
    alpha  = 0.9
  ) +
  scale_shape_manual(values = shape_palette, name = "Assay after passage") +
  scale_fill_manual(
    values = change_colours, name = "Change observed",
    breaks = c("Yes", "No"), limits = c("Yes", "No"),
    labels = c(
      Yes = paste0("<span style='color:", change_colours[["Yes"]], "'>**Yes**</span>"),
      No  = paste0("<span style='color:", change_colours[["No"]],  "'>**No**</span>")
    )
  ) +
  scale_color_manual(values = change_colours, name = "Change observed",
                     breaks = c("Yes", "No"), limits = c("Yes", "No")) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(0, x_max),
    breaks = seq(0, floor(x_max), by = 2)
  ) +
  labs(
    title = "Selection conditions and outcomes",
    x = "Number of studies",
    y = "Selection condition"
  ) +
  theme_bw(base_size = base_size) +
  theme(
    text               = element_text(size = base_size),
    plot.title         = element_text(size = base_size, face = "bold", hjust = 0),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    legend.box         = "vertical",
    legend.key.size    = unit(0.4, "lines"),
    legend.text        = element_markdown(size = base_size),
    legend.title       = element_text(size = base_size),
    legend.margin      = margin(t = 0, b = 0),
    axis.title         = element_text(size = base_size),
    axis.text          = element_text(size = base_size),
    plot.margin        = margin(2, 2, 2, 4)
  ) +
  guides(
    shape = guide_legend(override.aes = list(size = 2), ncol = 1,
                         title.position = "top"),
    fill  = guide_legend(
      ncol = 1, title.position = "top",
      # blank key glyph: the colour lives in the label text, not a symbol,
      # because the plotted points use several shapes (a circle key misleads)
      override.aes = list(shape = NA)
    ),
    color = "none"
  )

# =============================================================================
# COMBINE PANELS (A over B) AND SAVE
# =============================================================================

combined <- (panel_A | panel_B) +
  plot_layout(widths = c(1.1, 1)) &
  theme(plot.tag = element_text(face = "bold", size = base_size))
combined <- combined + plot_annotation(tag_levels = "A")

# Use a high-DPI cairo device for the raster preview (crisp anti-aliased text)
# and keep the PDF as the true vector deliverable for submission.
ggsave(png_file, combined, width = fig_width_in, height = fig_height_in,
       dpi = 600, type = "cairo")
ggsave(pdf_file, combined, width = fig_width_in, height = fig_height_in,
       device = cairo_pdf)

print(combined)