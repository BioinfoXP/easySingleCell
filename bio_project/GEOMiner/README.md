<div align="center">

<img src="https://wandering.oss-cn-hangzhou.aliyuncs.com/OB_Zotero/20260119103543.png" width="256" alt="GEOMiner Logo">

# GEOMiner: AI-Driven Semantic Curation of NCBI GEO Datasets

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/Language-R-276DC3.svg)](https://www.r-project.org/)
[![GEO Mining](https://img.shields.io/badge/GEO-Mining-green.svg)](https://www.ncbi.nlm.nih.gov/geo/)

</div>

> **GEOMiner** is a sophisticated R package designed to automate the retrieval, curation, and analysis of transcriptomic datasets from the NCBI Gene Expression Omnibus (GEO). Leveraging **Large Language Models (LLMs)** and **Retrieval-Augmented Generation (RAG)**, it moves beyond traditional keyword searches to provide context-aware semantic curation.

---

## 🚀 Key Features

* **🤖 LLM-Powered Query Translation**
  Converts natural language queries into precise NCBI Boolean search strings, automatically handling synonyms, MeSH terms, and complex logic.
* **🎯 Context-Aware Scoring**
  Scores datasets on a **0-100 scale** based on user-defined inclusion/exclusion criteria, ensuring relevance to specific contexts (e.g., tissue types, disease states).
* **⚡ Parallel Processing**
  Optimizes large-scale analysis with efficient, multi-threaded metadata parsing and screening.
* **📊 Automated Extraction**
  Automatically extracts critical metadata such as Sample Size ($N$), Organism, and Study Design.
* **📑 Smart Export**
  Generates color-coded Excel reports with relevance scores and direct hyperlinks for rapid manual verification.

---

## 📦 Installation

Install the development version of **GEOMiner** directly from GitHub:

```r
# Install devtools if not already available
if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}

# Install GEOMiner
devtools::install_github("BioinfoXP/GEOMiner")
```

## 🔑 Configuration

**GEOMiner** requires an LLM API key (OpenAI, DeepSeek, or compatible). For security, set this in your R environment:

R

```
Sys.setenv(OPENAI_API_KEY = "sk-...")
```

*Alternatively, you can pass the API key directly via function arguments.*

------

## 📖 Usage

The core function `geo_mine()` handles the end-to-end workflow.

### 1. Basic Search

*Standard search for a specific disease or topic.*

R

```
library(GEOMiner)

# Example: Hepatocellular Carcinoma (HCC) immunotherapy studies
results <- geo_mine(
  keyword = "HCC immunotherapy",
  limit = 10,
  output_file = "hcc_results.xlsx"
)
```

### 2. Context-Specific Mining (RAG-Powered)

*Refine results by defining specific inclusion/exclusion criteria.*

R

```
# Example: Lung cancer, focusing on EGFR mutations in human tissue
# Explicitly excluding cell lines and PDX models
results <- geo_mine(
  keyword = "Lung Cancer EGFR",
  context = "Only human tissue samples with EGFR mutation status. Exclude cell lines and PDX models.",
  limit = 20,
  output_file = "lung_egfr_clinical.xlsx"
)
```

### 3. Integration with DeepSeek / Local LLMs

*Switch providers by changing the `base_url`.*

R

```
results <- geo_mine(
  keyword = "Single cell RNA-seq breast cancer",
  limit = 10,
  base_url = "[https://api.deepseek.com/v1](https://api.deepseek.com/v1)",  # DeepSeek Endpoint
  model = "deepseek-chat",                   # Model Name
  api_key = "sk-..."                         # Your API Key
)
```

------

## 🏆 Scoring System

GEOMiner assigns a **Match Score** to categorize dataset relevance:

| **Score** | **Level**          | **Description**                                              |
| --------- | ------------------ | ------------------------------------------------------------ |
| 🟢         | **90-100 (High)**  | **Strong Match.** Dataset contains detailed clinical/phenotype data, correct organism, and study design aligning with criteria. |
| 🟡         | **70-89 (Medium)** | **Moderate Match.** Topic is relevant, but key details (e.g., survival data) may be missing or sample size ($N$) is small. |
| ⚪         | **0-69 (Low)**     | **Low Relevance.** Mismatched organisms, inappropriate sample types (e.g., cell lines instead of tissue), or sparse metadata. |

------

## 🛠 Methodology

The pipeline follows a five-step automated process:

1. **Query Generation**: Processes natural language keywords into valid NCBI Entrez query strings.
2. **Data Retrieval**: Fetches metadata for potential GEO Series (GSE) IDs via `rentrez`.
3. **Parallel Parsing**: Scrapes and analyzes abstracts using parallel processing for speed.
4. **Semantic Evaluation**: Scores dataset relevance against the research context using LLM logic.
5. **Structured Output**: Standardizes data into an Excel report for review.

------

## ⚠️ Disclaimer

While **GEOMiner** streamlines metadata curation, **results should be verified** against original GEO records. Automated scoring and entity extraction assist effectively, but rigorous research requires manual validation of the selected datasets.

------

## 👤 Author

Peng Xia, PhD Candidate

School of Basic Medical Sciences, Lanzhou University

- **Focus**: Cancer Genomics, Single-Cell Analysis, AI-Driven Data Mining
- 📧 **Email**: [xp294053@163.com](mailto:xp294053@163.com)
- 🐈 **GitHub**: [BioinfoXP](https://github.com/BioinfoXP)

------

## 📄 License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).
