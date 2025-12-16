# 01 : Setting Up Your Digital Lab Bench (Conda Environment)

##  Basic Concept 

Imagine you are about to start a complex experiment in a wet lab. You wouldn't just dump all your chemicals and tools onto a cluttered desk, right? You would set up a specific **Lab Bench** with exactly the pipettes, reagents, and machines you need for *that specific experiment*.

In bioinformatics, **Conda** allows us to do the exact same thing on a computer.
*   **The Problem:** Different software tools often conflict with each other (like needing different versions of Python). Installing them all on your main system is like mixing all your chemicals in one bucket—messy and dangerous.
*   **The Solution:** We create a "Virtual Environment" (our digital lab bench). Inside this environment, we install only the tools we need for ChIP-seq. When we are done, we can "close" the environment and go back to a clean computer.

In this tutorial, we will build a simplified environment named `chip` that contains all the tools for our analysis.

---

##  Execution 

### Step 1: Get the Environment Manager (Anaconda)
First, you need the software that builds these environments. We recommend **Anaconda** (or the lighter version, **Miniconda**).

*   **Check if you already have it:**
    Open your terminal and type:
    ```bash
    conda --version
    ```
    *(If you see a version number like `conda 24.7.1`, skip to Step 2.)*

*   **If not, download and install it:**
    *   **Anaconda (Recommended):** [https://www.anaconda.com/download](https://www.anaconda.com/download)
    *   **Miniconda (Lightweight):** [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)

### Step 2: Creates the "Recipe" File
To build our lab bench, we give Conda a "recipe" list called a YAML file. This tells Conda exactly which tools to fetch.

1.  Create a new file named `chip_env.yml`.
2.  Copy and paste the exact code below into that file:

```yaml
name: chip   # Name of our environment
channels:    # The "App Stores" where we find tools
  - bioconda
  - conda-forge
  - defaults
dependencies: # The tools we want to install
  - multiqc=1.31
  - fastqc=0.12.1
  
prefix: /opt/anaconda3/envs/chip 
```

> **Note:** If you downloaded the file [chip.yml](https://github.com/ishaaq34/Chipseq_analysis_tutorial/blob/main/Chip.yml) . You can use that to install all the tools required for the Chip-seq analysis!

### Step 3: Build the Environment
Now, let's look at the instructions to build the bench.

1.  **Create the environment:**
    Run this command in the same folder where you saved your file:
    ```bash
    conda env create -f chip_env.yml
    ```
    *(This takes a few minutes as it downloads all the tools.)*

2.  **Enter the environment:**
    To start working, you must "step into" your new lab bench:
    ```bash
    conda activate chip
    ```
    *(You should see `(chip)` appear next to your cursor in the terminal.)*

### Step 4: Verify Your Tools
Let's make sure our tools are actually there. Run these commands:

```bash
# Check if key tools are reachable
which fastqc
which bowtie2
which macs3

# Check versions to ensure successful installation
fastqc --version
bowtie2 --version
```
If these commands print paths (like `/Users/.../envs/chip/bin/fastqc`) and version numbers, **congratulations!** Your digital lab bench is ready.

---


### Understanding the YAML "Recipe"
Let's break down the `chip_env.yml` file we just used:

*   **`channels`**: These are repositories (like App Stores).
    *   **`bioconda`**: The community hub for bioinformatics software.
    *   **`conda-forge`**: A massive collection of general-purpose tools.
    *   *The order matters!* Conda looks in the first channel on the list, then the second.
*   **`dependencies`**: This is your shopping list.
    *   **Version Pinned (e.g., `bowtie2=2.5.4`)**: We specify the *exact* version. Why? So that if you run this analysis 2 years from now, you get the exact same results. This is key for **Reproducibility**.

### Managing Your Environment
Here are some useful "Housekeeping" commands:

| Command | Explanation |
| :--- | :--- |
| `conda deactivate` | "Step out" of the environment. Returns you to your base system. |
| `conda env list` | Lists all environments you have created on your computer. |
| `conda env remove -n chip` | Deletes the `chip` environment mostly used if you made a mistake and want to start over. |
| `conda env update -f chip_env.yml --prune` | Updates the environment if you change the YAML file. The `--prune` flag removes old tools you don't need anymore. |

---

## **Summary**
1.  **Analogy:** We built a dedicated "Lab Bench" (Environment) to keep our work clean.
2.  **Action:** We used `conda env create` with a YAML recipe to install tools like Bowtie2 and MACS2.
3.  **Result:** We now have a `(chip)` environment ready for the actual ChIP-seq analysis.

