# peptidomicsR

## Tools for peptidomics analysis of digesta from Protein Digestion

peptidomicsR provides functions to process, filter, analyze, and visualize peptidomics data, especially from MaxQuant protein digestion studies. This Readme file will introduce this package with the following three aspects:

✨**Functions**

🚀**Example analysis workflow**

📊**Example plots**

## ✨ Functions

-   **Data Processing (processPeptides())**

1.  Import MaxQuant output "peptides" file, *`intensity column meatadata, and protein mapping`*`.`

2.  Automatically remove contaminants and reverse sequences.

3.  Compute replicate- and group-level mean intensities and peptide counts.

4.  Map parental protein name and protein group for each peptide

5.  Calculate peptide GRAVY scores.

-   **Filtering & Statistics**

1.  PCA analysis on the clustering of sample groups for pre-checking the difference of samples groups (pcaPeptides()).

2.  Subset peptides by sequence, regex pattern, or grouping variables (filterPeptides()).

3.  Compare specific groups with t-tests or limma::treat (ttestPeptides()) to get significantly different peptides.

-   **Visualization**

1.  Stacked bar plots of peptide intensities with their parental proteins (plot_int()).

2.  Stacked bar plots of number of peptide types with their parental proteins (plot_count()).

3.  Length of peptides distribution with intensity plots (plot_length_distribution()).

4.  The intensity of each amino acid occurring as peptide cleavage sites (plot_cleavage_site())

5.  GRAVY score distribution with intensity plots (plot_gravy_vs_intensity()).

6.  PCA clustering plot for pre-checking the difference of samples groups (plot_pcaPeptides())

7.  Volcano plot(s) for peptide differential tests (t.test or limma) (plot_volcano())

## 🚀 Example Analysis Workflow

1.  Import and process data → processPeptides()

2.  Explore distributions:

-   Intensities → plot_int()

-   Counts → plot_count()

-   Length distribution → plot_length_distribution()

-   Cleavage sites → plot_cleavage_site()

-   Hydrophobicity distribution → plot_gravy_vs_intensity()

3.  Filter subsets of peptides → filterPeptides()

4.  Perform statistical comparisons → ttestPeptides()

5.  Volcano plot(s) → (plot_volcano())

## 📦 Installation

```{r}
# Install from GitHub
if (!require("devtools")) install.packages("devtools")
devtools::install_github("xchuam/peptidomicsR")
```

## 📊 Example Plots

This is only a small sample of the functions included in the package! For a complete list, please check the package vignette.

### **1. Import and process data**

```{r}

library(peptidomicsR)

result <- processPeptides(
  peptides_file          = "Data/peptides.txt",
  intensity_columns_file = "Data/Intensity_columns.csv",
  protein_mapping_file   = "Data/protein_mapping.csv"
)


#Place peptides.txt, Intensity_columns.csv, and protein_mapping.csv under Data/ as shown above.
```

### **2.** Plots of peptide mean intensities with parental proteins

```{r}
plot_int(
  result,
  x_var         = "Lipid",
  type          = "mean",
  facet_cols    = "Casein.ratio",
  filter_params = list(Digest.stage = "I"),
  color_by      = "Protein.name"
)

```

### 3. Length of peptides distribution with intensity plots

```{r}
plot_length_distribution(
  result,
  facet_cols    = "Casein.ratio",
  facet_rows    = "Lipid",
  filter_params = list(Digest.stage = "I")
)
```

### 4. The intensity of peptide cleavage sites

```{r}

plot_cleavage_site(
result,
terminal = c("both"),
measure = c("intensity"),
replicate_mode = c("mean"),
filter_params = list(Digest.stage = "I"),
scientific_10_y   = TRUE,
drop_constant_groups = TRUE)

```

### 5. GRAVY score distribution with intensity plots

```{r}

plot_gravy_vs_intensity(
  result,
  facet_cols    = "Casein.ratio",
  facet_rows    = "Lipid",
  filter_params = list(Digest.stage = "G"),
)
```

### 6. T-test on the two specific groups

```{r}
ttest_result <- ttestPeptides(
    result,
    comparisons,
    pseudocount       = 1,
    test_method       = c("treat", "plain"),
    lfc_thresh        = 1,
    alpha             = 0.05,
    equal_var         = FALSE,
    adjust            = "BH",
    min_reps_per_side = 2
)
```

### 7. Volcano plot(s) for peptide differential tests

```{r}
plot_volcano(
    ttest_result,
    comparisons,
    test_method        = c("treat","plain"),
    show_threshold     = TRUE,
    lfc_thresh         = 1,
    alpha              = 0.05,
    fill_values        = c(no = "grey75", yes = "#FFC010"),
    point_size         = 2,
    point_alpha        = 0.85,
    label_seqs         = NULL,
    label_size         = 3,
    label_col          = "black",
    highlight_seqs     = NULL,
    highlight_size     = NULL,
    highlight_stroke   = 1.2,
    y_log_scale        = FALSE)
```

## 🛠️ Dependencies

-   **data.table**, **ggplot2**, **ggpubr**, **scales** (plots & helpers)
-   **limma** (optional; for `treat`)

## 📖 Citation

If you use peptidomicsR in your work, please cite this repository.

## 📜 License

This project is licensed under the MIT License -- see the LICENSE file for details.




