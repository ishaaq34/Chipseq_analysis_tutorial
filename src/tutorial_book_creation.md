# Appendix: How We Built This Book

This tutorial documents the exact steps taken to turn our loose collection of tutorial files into this structured **mdbook**.

## 1. Getting Started
The very first step is to get the code and set up the book tool.

### Step 1: Clone the Repository
We started by downloading the tutorial material from GitHub:
```bash
git clone https://github.com/your-username/Chipseq_analysis_tutorial.git
cd Chipseq_analysis_tutorial
```

### Step 2: Install mdbook
On macOS (using Homebrew):
```bash
brew install mdbook
```

### Step 3: Initialize the Book
We initialized a new book structure inside our folder:
```bash
mdbook init my-first-book
cd my-first-book
```
*   This created the `book.toml` configuration file and the `src/` directory.

---

## 2. Organizing the Content
We moved our existing tutorial files into this new structure.

### File Migration
*   **Moved:** All markdown files (`00_setup_environment.md` through `pipeline_tutorial_simple.md`) were moved into the `src/` directory.
*   **Excluded:** `02_bash_automation.ipynb` and the RMarkdown file remained in the root (for now) as they needed conversion.

### Configuration (`SUMMARY.md`)
We updated `src/SUMMARY.md` to define the order of chapters in the sidebar:
```markdown
# Summary
- [Setup Environment](./00_setup_environment.md)
...
```

---

## 3. The "Notebook" Problem (.ipynb)
**Issue:** `mdbook` does not natively support Jupyter Notebooks.
**Solution:** Convert them to standard Markdown.

**Command Used:**
```bash
# Convert the notebook to markdown
jupyter nbconvert --to markdown 02_bash_automation.ipynb

# Move the result to src
mv 02_bash_automation.md src/
```

---

## 4. The "ChIPseeker" Problem (HTML Reports)
**Issue:** We had an RMarkdown report (`ChIP-seq_analysis_with_Chipseeker.html`) with interactive plots. Markdown cannot display dynamic HTML directly, but we wanted to keep the interactivity.

**Solution:** The Hybrid Wrapper (Iframe).

1.  **Move the HTML:** The `.html` file **must** be inside `src/` so it gets copied to the final website.
    ```bash
    mv ChIP-seq_analysis_with_Chipseeker.html src/
    ```

2.  **Create a Wrapper File:** We created `src/14_chipseeker_annotation.md`.

3.  **Embed it:** We used an HTML `<iframe>` tag inside the markdown file:
    ```html
    # Annotation Tutorial
    
    (Normal text explaining the science goes here...)
    
    <!-- This window shows the other file -->
    <iframe src="ChIP-seq_analysis_with_Chipseeker.html" width="100%" height="800px" style="border:none;"></iframe>
    ```

---

## 5. Building & Serving
To see the book, we run:
```bash
mdbook serve --open
```
This serves the website at `http://localhost:3000` and auto-updates whenever we save a file.
