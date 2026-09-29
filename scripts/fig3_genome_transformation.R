# =============================================================================
# Figure 3 for Royal Society journal
#   Panel A : Timeline of microsporidia genome-transformation papers
#   Panel B : Summary table of transformation studies (from the "Table" tab)
#
# Same house style as Figure 1: uniform 8 pt text, patchwork composition,
# cairo output devices, abbreviated species names.
#
# This script builds BOTH a stacked single-column version and a side-by-side
# two-column version, so they can be compared. Set LAYOUT below.
#
# install.packages(c("readxl","dplyr","tidyr","ggplot2","scales",
#                    "stringr","forcats","patchwork"))
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(stringr)
library(forcats)
library(patchwork)
library(gridExtra)
library(grid)

# ---------------------------
# user settings
# ---------------------------
file_path  <- file.path("data", "Table_S2_transformation.xlsx")
bib_sheet  <- "Master_Bibliography"   # title on row 1, header on row 2 -> skip = 1
# Panel B summary of the three transformation studies (first author, year,
# DNA delivery approach, fluorescent marker, selection, result)
tab_file   <- file.path("data", "fig3B_transformation_summary.csv")
strict_only <- TRUE                   # keep only Strict_Inclusion == "Yes"

# LAYOUT: "stacked" (single column) or "side" (two column)
LAYOUT <- "side"

# Royal Society widths
col1_width_in <- 3.3    # single column (8.4 cm)
col2_width_in <- 6.85   # full two columns (17.4 cm)
base_size     <- 8      # uniform printed point size

# =============================================================================
# PANEL A : TIMELINE  (adapted from the user's script, restyled to match Fig 1)
# =============================================================================

bib <- read_excel(file_path, sheet = bib_sheet, skip = 1)

if (strict_only && "Strict_Inclusion" %in% names(bib)) {
  bib <- bib %>% filter(Strict_Inclusion == "Yes")
}

bib <- bib %>% distinct(Title, Year, .keep_all = TRUE)

year_species <- bib %>%
  mutate(
    Year = suppressWarnings(as.integer(Year)),
    Parasite_Species = ifelse(
      is.na(Parasite_Species) | Parasite_Species == "",
      "Unspecified", Parasite_Species
    )
  ) %>%
  filter(!is.na(Year)) %>%
  count(Year, Parasite_Species, name = "n")

all_years <- seq(min(year_species$Year, na.rm = TRUE),
                 max(year_species$Year, na.rm = TRUE))
all_species <- sort(unique(year_species$Parasite_Species))

year_species_full <- expand_grid(
  Year = all_years,
  Parasite_Species = all_species
) %>%
  left_join(year_species, by = c("Year", "Parasite_Species")) %>%
  mutate(n = replace_na(n, 0))

year_totals <- year_species_full %>%
  group_by(Year) %>%
  summarise(papers_per_year = sum(n), .groups = "drop") %>%
  arrange(Year) %>%
  mutate(cumulative_papers = cumsum(papers_per_year))

scale_factor <- max(year_totals$papers_per_year) / max(year_totals$cumulative_papers)
year_totals <- year_totals %>%
  mutate(cumulative_scaled = cumulative_papers * scale_factor)

# abbreviate species to "G. species" (same helper as Fig 1)
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

year_species_full <- year_species_full %>%
  mutate(Parasite_Species = abbreviate_species(Parasite_Species))

species_order <- year_species_full %>%
  group_by(Parasite_Species) %>%
  summarise(total = sum(n), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(Parasite_Species)

year_species_full <- year_species_full %>%
  mutate(Parasite_Species = factor(Parasite_Species, levels = species_order))

# legend placement depends on layout (right for side-by-side, bottom for stacked)
legend_pos_A <- if (LAYOUT == "side") "right" else "bottom"

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
  scale_x_continuous(breaks = all_years[all_years %% 2 == 0]) +
  scale_y_continuous(
    name = "Papers per year",
    breaks = seq(0, ceiling(max(year_totals$papers_per_year)), by = 1),
    sec.axis = sec_axis(trans = ~ . / scale_factor, name = "Cumulative papers"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(title = "Genome-transformation papers over time",
       x = "Year", fill = "Microsporidia species") +
  theme_bw(base_size = base_size) +
  theme(
    text             = element_text(size = base_size),
    plot.title       = element_text(size = base_size, face = "bold", hjust = 0),
    panel.grid.minor = element_blank(),
    axis.text        = element_text(size = base_size),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = base_size),
    axis.title       = element_text(size = base_size),
    legend.position  = legend_pos_A,
    legend.box       = "vertical",
    legend.key.size  = unit(0.4, "lines"),
    legend.text      = element_text(size = base_size, face = "italic"),
    legend.title     = element_text(size = base_size),
    legend.margin    = margin(t = 0, b = 0, l = 0, r = 0),
    plot.margin      = margin(2, 4, 2, 2)
  ) +
  guides(fill = guide_legend(
    ncol = if (LAYOUT == "side") 1 else 2,
    byrow = TRUE, title.position = "top"
  ))

# =============================================================================
# PANEL B : SUMMARY TABLE  (rendered as a ggplot grob to match the figure)
# =============================================================================

tab <- read.csv(tab_file, check.names = FALSE, stringsAsFactors = FALSE,
                na.strings = c("", "NA"))

# rename columns for display
display_names <- c(
  "First_Author"                    = "First\nauthor",
  "Year"                            = "Year",
  "Approaches used to delivery DNA" = "DNA delivery\napproach",
  "Fluorescence"                    = "Fluorescence\nmarker",
  "Selection"                       = "Selection",
  "Result"                          = "Transformation\nresult"
)
names(tab) <- ifelse(names(tab) %in% names(display_names),
                     display_names[names(tab)], names(tab))

# fix a small typo in the source data for the figure
tab[] <- lapply(tab, function(col) {
  if (is.character(col)) col <- gsub("flourcent", "fluorescent", col, ignore.case = TRUE)
  col <- as.character(col)
  col[is.na(col) | col == ""] <- "\u2013"   # show empty cells as an en dash, not "NA"
  col
})

# per-column wrap widths (characters) so long cells break instead of overflowing
wrap_w <- c("First\nauthor" = 8, "Year" = 4, "DNA delivery\napproach" = 16,
            "Fluorescence\nmarker" = 8, "Selection" = 11,
            "Transformation\nresult" = 14)

tab_wrapped <- tab
for (cn in names(tab_wrapped)) {
  w <- if (cn %in% names(wrap_w)) wrap_w[[cn]] else 18
  tab_wrapped[[cn]] <- str_wrap(as.character(tab_wrapped[[cn]]), width = w)
}

# table theme matching the figure house style (8 pt, zebra rows, bold header)
tt <- ttheme_minimal(
  base_size = base_size - 2,
  core = list(
    fg_params = list(hjust = 0, x = 0.03, lineheight = 0.9),
    bg_params = list(fill = c("grey95", "white"), col = NA)
  ),
  colhead = list(
    fg_params = list(fontface = "bold", hjust = 0, x = 0.03, lineheight = 0.9),
    bg_params = list(fill = "grey80", col = NA)
  )
)

tbl_grob <- tableGrob(tab_wrapped, rows = NULL, theme = tt)

# set relative column widths so the long-text columns get more room
rel_w <- c(0.13, 0.07, 0.20, 0.19, 0.17, 0.24)
rel_w <- rel_w / sum(rel_w)
tbl_grob$widths <- unit(rel_w, "npc")

# add thin horizontal rules: top, below header, bottom
tbl_grob <- gtable::gtable_add_grob(
  tbl_grob,
  grobs = grid::segmentsGrob(x0 = 0, x1 = 1, y0 = 1, y1 = 1,
                             gp = grid::gpar(lwd = 1.2)),
  t = 1, b = 1, l = 1, r = ncol(tbl_grob))
tbl_grob <- gtable::gtable_add_grob(
  tbl_grob,
  grobs = grid::segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0,
                             gp = grid::gpar(lwd = 0.8)),
  t = 1, b = 1, l = 1, r = ncol(tbl_grob))
tbl_grob <- gtable::gtable_add_grob(
  tbl_grob,
  grobs = grid::segmentsGrob(x0 = 0, x1 = 1, y0 = 0, y1 = 0,
                             gp = grid::gpar(lwd = 1.2)),
  t = nrow(tbl_grob), b = nrow(tbl_grob), l = 1, r = ncol(tbl_grob))

# bold panel title above the table, matching the other panels' style.
# Indent it from the left so it does not collide with the bold "B" panel tag
# that patchwork places in the top-left corner of this panel.
tbl_title <- grid::textGrob(
  "Transformation studies",
  x = grid::unit(2.2, "lines"), hjust = 0,
  gp = grid::gpar(fontface = "bold", fontsize = base_size)
)
tbl_with_title <- gridExtra::arrangeGrob(
  grobs = list(tbl_title, tbl_grob),
  ncol = 1,
  heights = grid::unit.c(grid::unit(1.4, "lines"), grid::unit(1, "null"))
)

# wrap the grob so patchwork can compose it as a panel
panel_B <- wrap_elements(full = tbl_with_title)

# =============================================================================
# COMBINE AND SAVE
# =============================================================================

if (LAYOUT == "side") {
  combined  <- (panel_A | panel_B) + plot_layout(widths = c(0.7, 1.55))
  fig_w <- col2_width_in
  fig_h <- 3.5
} else {
  combined  <- (panel_A / panel_B) + plot_layout(heights = c(1.6, 1))
  fig_w <- col1_width_in
  fig_h <- 5.4
}

combined <- combined +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = base_size))

png_file <- file.path("figures", "fig3_genome_transformation.png")
pdf_file <- file.path("figures", "fig3_genome_transformation.pdf")

ggsave(png_file, combined, width = fig_w, height = fig_h, dpi = 600, type = "cairo")
ggsave(pdf_file, combined, width = fig_w, height = fig_h, device = cairo_pdf)

print(combined)