# Practical ChIP-seq Analysis Tutorial

A simplified, "Zero-to-Hero" guide to ChIP-seq analysis using **mdbook**. This repository contains the source code for the interactive book.

**[Read the Tutorial Book](./book/index.html)** *(If viewing locally after build)*

## About This Project
This tutorial walks through a complete ChIP-seq analysis pipeline using real **ENCODE** data (CEBPA, H3K27me3, H3K9ac). It is structured using the "Tiered Learning" method:
1.  **Concept:** Understanding the biology.
2.  **Execution:** Running the code.
3.  **Interpretation:** Analyzing the results.

## Requirements
To build this book locally, you need:
*   **mdbook** (The static site generator)
*   **Git**

## How to Run Locally

### 1. Clone the Repository
```bash
git clone https://github.com/your-username/Chipseq_analysis_tutorial.git
cd Chipseq_analysis_tutorial
```

### 2. Install mdbook
If you don't have it installed:
```bash
# macOS
brew install mdbook

# Or via Cargo
cargo install mdbook
```

### 3. Build & Serve
Run the server to view the book in your browser (at `http://localhost:3000`):
```bash
mdbook serve --open
```

## Project Structure
*   `src/`: Contains all the Markdown content (Chapters).
*   `book/`: The generated HTML site (excluded from git).
*   `book.toml`: Configuration file.

## Data Availability
The analysis uses public ENCODE data. See `src/00_introduction.md` for specific accession numbers and DOIs.
