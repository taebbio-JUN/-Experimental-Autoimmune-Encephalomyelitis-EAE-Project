EAE / CFA Microbiome Bioinformatics Pipeline (Family & Genus)
==============================================================

Overview
--------
This repository contains a reproducible bioinformatics workflow for microbiome count tables comparing
three groups: Control, CFA/PTX, and EAE. The pipeline is organized into three main analyses:

1) Beta-diversity (Aitchison distance) PCoA
   - Count -> relative abundance -> CLR transform (with pseudocount) -> Euclidean distance = Aitchison distance
   - PCoA visualization
   - PERMANOVA (adonis2), betadisper (homogeneity of dispersion), and pairwise PERMANOVA
   - Statistical results are exported to Excel ONLY (not printed on the plot)

2) Family-level summary: log2FC heatmap from DESeq2 normalized counts (Excel Sheet2)
   - Input: DESeq2 normCounts stored in an Excel sheet (usually “Family_DESeq2_normCounts”)
   - log2FC = log2((mean_group + pseudo)/(mean_control + pseudo))
   - Significance stars inside heatmap cells are computed via t-test on log2(x + pseudo)
   - Both Down (CFA & EAE both decreased) taxa are grouped at the TOP and exported as a separate CSV

3) Genus-level focused analysis (Family-selected)
   - Input: the Family-selected list produced in step (2)
   - Subset the Genus sheet to genera belonging to selected families
   - DESeq2(sfType="poscounts") normalization + contrasts:
       (1) CFA vs Control
       (2) EAE vs Control
   - Plot values: per-sample log2FC vs Control mean (for visualization only)
   - Effect filter: drop if |median log2FC| < 1 (default)
   - Keep genera if p <= 0.10 in EITHER contrast; stars only for p < 0.05
   - Split half-violin + half-boxplot: CFA on the LEFT, EAE on the RIGHT


Input Data Requirements
-----------------------
Primary input:
- Taxonomy_abundance_count.xlsx

Required sheets:
- “Family” sheet: family-level counts
- “Genus” sheet: genus-level counts

Sample column naming (strongly recommended):
- Control1, Control2, Control3, Control4  (or “Control 1” etc. — spaces are tolerated)
- CFA1, CFA2, CFA3, CFA4
- EAE1, EAE2, EAE3, EAE4

Notes:
- Group labels are inferred from column names. Keep them consistent.
- Some minor typos (e.g., “Contol”) may be corrected in code, but do not rely on this.


Recommended Repository Structure
-------------------------------
.
├─ data/
│  └─ Taxonomy_abundance_count.xlsx
├─ scripts/
│  ├─ 01_pcoa_aitchison_permanova.R
│  ├─ 02_family_log2fc_heatmap_from_normcounts.R
│  ├─ 03_genus_familySelected_deseq2_splitHalfViolin.R
│  ├─ install_dependencies.R
│  └─ session_versions.R
├─ results/
│  ├─ 01_pcoa/
│  ├─ 02_family_heatmap/
│  └─ 03_genus_selected/
└─ README.txt


Environment and Dependencies
----------------------------
R version:
- R 4.x is recommended (Bioconductor/DESeq2 compatibility)

Required packages:

Base R:
- grid  (included with R; no installation needed)

CRAN:
- readxl
- dplyr
- tidyr
- stringr
- tibble
- ggplot2
- vegan
- openxlsx
- pheatmap

Bioconductor:
- DESeq2

Other:
- gghalves (for split half-violin plots)

Important:
- gghalves installation may fail depending on your R version and platform. If CRAN installation fails,
  install it from GitHub using remotes (see below).


One-shot Installation
---------------------
Save the following as scripts/install_dependencies.R and run it once.

----- BEGIN scripts/install_dependencies.R -----
# 1) CRAN packages
cran_pkgs <- c(
  "readxl","dplyr","tidyr","stringr","tibble",
  "ggplot2","vegan","openxlsx","pheatmap"
)
install.packages(cran_pkgs)

# 2) Bioconductor: DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")

# 3) gghalves
# Try CRAN first (may fail on some setups)
try(install.packages("gghalves"), silent = TRUE)

# If not installed, use GitHub
if (!requireNamespace("gghalves", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("erocoar/gghalves")
}
----- END scripts/install_dependencies.R -----


Reproducibility: Capture Package Versions
-----------------------------------------
To make results reproducible across machines, record package versions and session info.
Save the following as scripts/session_versions.R and run it.

----- BEGIN scripts/session_versions.R -----
pkgs <- c(
  "readxl","dplyr","tidyr","stringr","tibble",
  "ggplot2","vegan","openxlsx","pheatmap",
  "DESeq2","gghalves","grid"
)

get_ver <- function(p) {
  if (p %in% rownames(installed.packages())) as.character(packageVersion(p)) else NA_character_
}

ver_df <- data.frame(
  Package = pkgs,
  Version = vapply(pkgs, get_ver, character(1)),
  stringsAsFactors = FALSE
)

write.csv(ver_df, "package_versions.csv", row.names = FALSE)

sink("sessionInfo.txt")
cat("R version:", as.character(getRversion()), "\n")
cat("\n--- sessionInfo() ---\n")
print(sessionInfo())
sink()

message("Wrote: package_versions.csv")
message("Wrote: sessionInfo.txt")
----- END scripts/session_versions.R -----


Workflow (Execution Order)
--------------------------

Step 1) PCoA (Aitchison distance) + PERMANOVA
- Script: scripts/01_pcoa_aitchison_permanova.R
- Purpose:
  - Convert counts to relative abundance, apply CLR transform with pseudocount, compute Aitchison distance,
    and run PCoA.
  - Test group differences using PERMANOVA (adonis2).
  - Evaluate homogeneity of dispersion using betadisper + permutest.
  - Run pairwise PERMANOVA and BH-adjust p-values.
- Outputs:
  - PCoA_<Rank>_Aitchison_PERMANOVA_stats.xlsx
    * overall_PERMANOVA
    * pairwise_PERMANOVA
    * PCoA_points
  - PCoA_<Rank>_Aitchison_SDellipse_WHITEBG.png

Interpretation notes:
- If PERMANOVA is significant but betadisper is also significant, part of the apparent group separation
  may be driven by differences in within-group dispersion rather than shifts in centroid location.
- This pipeline exports betadisper results alongside PERMANOVA to support careful interpretation.

Plot-only group separation (sep_factor):
- The script optionally increases spacing between groups ONLY for plotting (visual readability).
- Statistical tests (PERMANOVA, betadisper) are computed on the original Aitchison distance, not on the
  shifted coordinates. This prevents visual adjustments from affecting inferential results.


Step 2) Family log2FC Heatmap from DESeq2 normCounts (Excel Sheet2)
- Script: scripts/02_family_log2fc_heatmap_from_normcounts.R
- Input:
  - An Excel file containing DESeq2 normalized counts (normCounts) on Sheet2.
    Typically the sheet is named “Family_DESeq2_normCounts”.
- Core calculations:
  - Control_mean_norm, CFA_mean_norm, EAE_mean_norm
  - log2FC_CFA = log2((CFA_mean_norm + pseudo)/(Control_mean_norm + pseudo))
  - log2FC_EAE = log2((EAE_mean_norm + pseudo)/(Control_mean_norm + pseudo))
  - Stars for visualization:
    t-test on log2(x + pseudo) comparing group vs control
      * p < 0.05  -> *
      * p < 0.01  -> **
      * p < 0.001 -> ***
- Features:
  - Stars are placed inside heatmap cells.
  - Heatmap cell borders are black.
  - “Both Down” (CFA & EAE both decreased) is grouped at the TOP and exported separately.

Outputs:
- Family_normalized_log2FC_stats.xlsx
- Family_log2FC_pval.all.csv
- Family_log2FC_pval.filtered.ordered.DOWNTOP.csv
- Family_BothDown_CFA_and_EAE.csv
- Family_log2FC_left_CFA_right_EAE.DOWNTOP.tiff

Interpretation notes:
- Heatmap log2FC is an effect-size summary based on group means from DESeq2-normalized counts.
- The star annotation is a fast visual significance marker (t-test on log2 transformed normalized counts),
  not the same as DESeq2’s model-based padj. For strict differential abundance, refer to the DESeq2
  contrast results in the Genus analysis step.


Step 3) Genus Analysis (Family-selected) + DESeq2 contrasts + split half-violin
- Script: scripts/03_genus_familySelected_deseq2_splitHalfViolin.R
- Inputs:
  - Selected Family CSV from Step 2 (auto-detected):
      * Family_log2FC_pval.filtered.ordered.DOWNTOP.csv
        or
      * Family_BothDown_CFA_and_EAE.csv
  - Genus counts from Taxonomy_abundance_count.xlsx (“Genus” sheet)
- Process:
  1) Load selected families
  2) Filter genera to those belonging to selected families
  3) DESeq2 run with sfType="poscounts" (recommended for sparse microbiome counts)
  4) DESeq2 contrasts:
       - CFA vs Control
       - EAE vs Control
  5) Plot values:
       per-sample log2FC vs Control mean using normalized counts (+pseudo)
  6) Effect filter:
       drop if |median log2FC| < 1 (default behavior)
  7) Statistical keep rule:
       keep genus if p <= 0.10 in EITHER contrast
       stars only for p < 0.05 (per contrast)

Outputs:
- Genus_selected_DESeq2_normCounts.xlsx
- Genus_stats_DESeq2.CFA_vs_Control__EAE_vs_Control.p_le_0.1.csv
- Genus_split_half_violin.log2FCvsControl.CFA_left_EAE_right.DESeq2_vsControl.p_le_0.1.tiff
- Additional helper outputs:
  - Genus_effect_filter_table.csv
  - Genus_kept_by_effect_cut.csv
  - Genus_kept_by_p_le_0.1.either_contrast.csv
  - Genus_plot_input.long.effectFiltered.p_le_0.1.csv

Interpretation notes:
- The plot shows per-sample deviations from the Control mean (visualization), while DESeq2 provides
  model-based contrast statistics (inference). Using both helps interpret variability and effect size.
- Keeping p <= 0.10 (either contrast) is a “broad candidate retention” strategy, while stars at p < 0.05
  highlight stronger evidence in the figure.


How to Run (Examples)
---------------------
From the repository root:

1) Install dependencies:
   Rscript scripts/install_dependencies.R

2) Run PCoA/PERMANOVA:
   Rscript scripts/01_pcoa_aitchison_permanova.R

3) Run Family log2FC heatmap:
   Rscript scripts/02_family_log2fc_heatmap_from_normcounts.R

4) Run Genus family-selected DESeq2 + split half-violin:
   Rscript scripts/03_genus_familySelected_deseq2_splitHalfViolin.R

5) Capture versions/session info:
   Rscript scripts/session_versions.R


Troubleshooting
---------------
1) Package installation errors / compilation issues
   - On Windows, some packages may require Rtools if installing from source.
   - Prefer CRAN binaries when available.

2) gghalves installation fails
   - Install from GitHub:
     remotes::install_github("erocoar/gghalves")

3) Sample columns not detected
   - Ensure columns are named consistently: Control1.., CFA1.., EAE1..
   - Spaces are tolerated, but consistent naming is best.

4) DESeq2 errors with input data
   - DESeq2 expects non-negative integer counts; if your input is not integer counts,
     review how counts were generated and confirm they are raw counts before DESeq2.

5) Interpretation caution
   - Significant PERMANOVA together with significant betadisper indicates dispersion differences may confound
     centroid-based group separation; interpret with care.


License
-------
Set an appropriate license for your project (e.g., MIT, CC BY-NC, internal-use only).


Contact / Notes
---------------
This README is designed to describe the full workflow and the rationale behind key analysis choices:
- Aitchison/CLR for compositional microbiome data
- PERMANOVA + betadisper pairing for robust beta-diversity interpretation
- DESeq2(poscounts) for sparse count stability
- Separate roles of visualization (per-sample log2FC) and inference (DESeq2 contrasts)
