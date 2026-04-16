# Bacterial Stress Responses Lower mRNA-Protein Level Correlations

## Requirements

- **R ≥ 4.4** — download from <https://cran.r-project.org>
- **Rtools** (Windows only) — must match your R version exactly.
  Download from <https://cran.r-project.org/bin/windows/Rtools/>

## 1 — Clone the repository

```bash
git clone https://github.com/avicanlab/RNAProtCorr
cd RNAProtCorr
```

## 2 — Restore the R environment with renv

This project uses [renv](https://rstudio.github.io/renv/) to pin all package
versions. Open R (or RStudio) in the project root and run:

```r
install.packages("renv")
renv::restore()
```

**Common issues**

| Problem                                          | Solution                                                                                                                      |
|--------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------|
| `BiocParallel` fails to compile                  | Update Rtools to match your R version, then re-run `renv::restore()`                                                          |
| A single package cannot be installed             | Exclude it: `renv::restore(exclude = "pkgName")`, then install manually with `install.packages()` or `BiocManager::install()` |
| `renv::restore()` asks to activate — say **yes** | This links the project library so packages are isolated                                                                       |

## 3 — Download data from Zenodo

All raw data files are archived on Zenodo at
**[https://zenodo.org/records/19488690](https://zenodo.org/records/19488690)**.

### Option A — R (recommended, reproducible)

```{r}
# install.packages("zen4R")   # run once if not already installed
library(zen4R)

zenodo <- ZenodoManager$new()
record <- zenodo$getRecordByDOI("10.5281/zenodo.19488689")
# List available files in the record (to find the exact filename)
files <- record$listFiles()
print(files)


dir.create(file.path("DATA"), showWarnings = FALSE)

# Download the DATA archive
record$downloadFiles(
  files = "DATA_notebook.zip",  # exact filename
  path = "DATA"
)

unzip(
  zipfile = file.path("DATA", "DATA_notebook.zip"),
  exdir   = "DATA"
)

# Optional: remove the zip after extraction to save space
file.remove(file.path("DATA", "DATA_notebook.zip"))
```

### Option B — Command line

```bash
# requires curl ≥ 7.x
data_notebook = "https://zenodo.org/records/19488690/files/DATA_notebook.zip?download=1"
curl -L "$data_notebook" -o $data_notebook.zip
unzip $data_notebook.zip -d DATA/
```

### Option C — Browser

1. Open [https://zenodo.org/records/19488690](https://zenodo.org/records/19488690)
2. Click **Download** for the DATA_notebook.zip archive.
3. Unzip into the `DATA/` folder at the project root

### Expected DATA folder structure

```
DATA/
├── Read_counts/          # RNA-seq TPM Excel files  (*.xlsx)
├── Proteome/             # TMT result TSV files, annotation, sample mapping
│   ├── sample_mapping.xlsx
│   ├── <Species>/
│   │   ├── annotation.xlsx
│   │   └── <Species>_R<n>.tsv
├── Fasta_files/          # Proteome FASTA files     (*.fa)
├── Essential genes/      # Essential gene lists     (*.xlsx)
├── Enrichment/           # GO term files            (*.xlsx / *.csv)
├── Stimulons/            # Stimulon gene lists      (*.xlsx)
└── ORA/                  # MOBILE framework output  (*.csv / *.txt)
```

## 4 — Render the notebook

```r
rmarkdown::render("main.Rmd")
```

Or use the **Knit** button in RStudio.

---
