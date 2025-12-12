# Tutorial 02: Bash Automation (Your Digital Robot)

## Level 1: Basic Concept (The Robot Assistant)

### Manual vs. Automated
Imagine you have 100 samples to analyze. You have two options:
1.  **Manual:** Type the command for Sample 1. Wait. Type for Sample 2. Wait... (This takes days and is prone to typos).
2.  **Automated:** Write a small instruction sheet (Script) and give it to a robot. The robot does all 100 samples while you drink coffee.

In this tutorial, we will learn how to build that "Robot" using **Loops** and **Variables**.

### Reproducibility (Setting the Stage)
Before we program our robot, we need data. We will create a `raw` folder and some dummy files. This ensures that everyone following this tutorial has the exact same setup.


```bash
%%bash

# --- 1. Setup: Creating Dummy Data ---
# We start by creating a directory and some dummy files to work with.
ls 

```

    [34manaconda_projects[m[m
    bash_tutorial.ipynb



```bash
%%bash

mkdir -p raw
ls raw
```


```bash
%%bash
touch raw/Input1.fastq.gz
touch raw/Input2.fastq.gz
touch raw/SampleA.fastq.gz
 ls raw
```

    Input1.fastq.gz
    Input2.fastq.gz
    SampleA.fastq.gz


**Explanation**: We now have a folder named `raw` containing three empty files. This mimics a real sequencing run.

## Level 2: Execution (The Loop)

### 2.1 The Basic Loop
The `for` loop is the heart of automation. It says: **"For every file you see, do this action."**

**The Syntax:**
```bash
for variable in pattern
do
    Action
done
```

Here, we use the wildcard `*` (which means "anything") to find all files ending in `.fastq.gz`.


```bash
%%bash
echo "--- Listing files ---"
for fq in raw/*.fastq.gz
do
    echo "Found file: $fq"
done

```

    --- Listing files ---
    Found file: raw/Input1.fastq.gz
    Found file: raw/Input2.fastq.gz
    Found file: raw/SampleA.fastq.gz


**Start Simple**: Notice we just `echo` (print) the filename first. **Always print your list before running real commands!** It prevents the robot from accidentally deleting your files.

### 2.2 Cleaning Names (`basename`)
The variable `$fq` contains the entire path: `raw/Input1.fastq.gz`.
But we often just want the name: `Input1`.

We use a tool called `basename` to strip away the folder (`raw/`) and the extension (`.fastq.gz`).

This allows us to create nice output filenames like `Input1.bam` instead of `raw_Input1.fastq.gz.bam`.


```bash
%%bash
echo -e "\n--- Extracting Sample Names ---"
for fq in raw/*.fastq.gz
do
    # syntax: basename path suffix
    sample=$(basename "$fq" .fastq.gz)
    echo "Sample name is: $sample"
done

```

    
    --- Extracting Sample Names ---
    Sample name is: Input1
    Sample name is: Input2
    Sample name is: SampleA


**Logic Check**: See how `Input1.fastq.gz` became just `Input1`? This is perfect for naming our upcoming alignment files.

### 2.3 The Full Command Generator
Now we combine everything. Instead of running the alignment, we will `echo` the exact command the robot WOULD run.

**Goal**: `bowtie2 -U raw/Input1.fastq.gz | samtools sort -o mapped/Input1.sorted.bam`

We construct this string using our variables:
*   Input: `$fq` (Full path)
*   Output: `${sample}.sorted.bam` (Clean name)


```bash
 %%bash

for fq in raw/*.fastq.gz
do
    # syntax: basename path suffix
    sample=$(basename "$fq" .fastq.gz)
    
    echo "  Input file: $fq"

    echo "Sample name is: $sample"

     echo "Processing sample: $sample"

    echo "  Command: bowtie2 -U $fq | samtools sort -o mapped/${sample}.sorted.bam"
done


```

      Input file: raw/Input1.fastq.gz
    Sample name is: Input1
    Processing sample: Input1
      Command: bowtie2 -U raw/Input1.fastq.gz | samtools sort -o mapped/Input1.sorted.bam
      Input file: raw/Input2.fastq.gz
    Sample name is: Input2
    Processing sample: Input2
      Command: bowtie2 -U raw/Input2.fastq.gz | samtools sort -o mapped/Input2.sorted.bam
      Input file: raw/SampleA.fastq.gz
    Sample name is: SampleA
    Processing sample: SampleA
      Command: bowtie2 -U raw/SampleA.fastq.gz | samtools sort -o mapped/SampleA.sorted.bam


**Success**: The loop generated three unique, correct alignment commands. If this were a real run, we would simply remove the word `echo` to execute them.

## Level 3: Advanced Automation (The List)

### 3.1 Arrays (The Shopping List)
Using `*` (Globbing) grabs *everything*. But what if you only want to process specific samples (e.g., just the Controls)?

We use an **Array** (a specific list of items).

Instead of saying "Get everything in the fridge", we say "Get only Milk and Eggs".


```bash
%%bash
echo -e "\n--- Looping over Array ---"
SAMPLES=("Input1" "Input2" "TumorA")

for s in "${SAMPLES[@]}"
do
    echo "Processing specific sample: $s"
done


echo " total sample files : ${SAMPLES[@]}"
```

    
    --- Looping over Array ---
    Processing specific sample: Input1
    Processing specific sample: Input2
    Processing specific sample: TumorA
     total sample files : Input1 Input2 TumorA


**Why use this?** It gives you total control. You don't accidentally process temporary files or old data.

### 3.2 Collecting Outputs (`multiBamSummary`)
Some tools (like `multiBamSummary` for correlation plots) need **ALL** the BAM files at once.
We can't just process them one by one. We need to collect them into a "basket" (Array) as we go.


```bash
%%bash


samp_files=() #“Create an empty list that we will fill with BAM file paths.”
for fq in raw/*.fastq.gz    #Append BAM paths to the array
do
    samp_files+=("$fq")       #Append BAM paths to the array
   echo "samples files Added: $fq" # optional print to see what is going 
   
  sample_name=$(basename "$fq" .fastq.gz)
  echo "sample_name of the file used : $sample_name"

done


echo
echo "--- Final FASTQ array --- "
printf "%s\n" "${samp_files[@]}"
```

    samples files Added: raw/Input1.fastq.gz
    sample_name of the file used : Input1
    samples files Added: raw/Input2.fastq.gz
    sample_name of the file used : Input2
    samples files Added: raw/SampleA.fastq.gz
    sample_name of the file used : SampleA
    
    --- Final FASTQ array --- : raw/Input1.fastq.gz raw/Input2.fastq.gz raw/SampleA.fastq.gz
    raw/Input1.fastq.gz
    raw/Input2.fastq.gz
    raw/SampleA.fastq.gz


**Logic**: `samp_files+=("$fq")` adds the file to our basket. At the end, we have the full basket to use.


```bash
%%bash


samp_files=() #“Create an empty list that we will fill with fastq file paths.”
for fq in raw/*.fastq.gz    
do
    samp_files+=("$fq")       #Append FG paths to the array
   echo "samples files Added: $fq" # optional print to see what is going 
   
  sample_name=$(basename "$fq" .fastq.gz)
  echo "sample_name of the file used : $sample_name"

   echo "  Command: bowtie2 -U $fq | samtools sort -o mapped/${sample}.sorted.bam"


done


 echo "total fastq files : ${samp_files[@]}" 
```

    samples files Added: raw/Input1.fastq.gz
    sample_name of the file used : Input1
      Command: bowtie2 -U raw/Input1.fastq.gz | samtools sort -o mapped/.sorted.bam
    samples files Added: raw/Input2.fastq.gz
    sample_name of the file used : Input2
      Command: bowtie2 -U raw/Input2.fastq.gz | samtools sort -o mapped/.sorted.bam
    samples files Added: raw/SampleA.fastq.gz
    sample_name of the file used : SampleA
      Command: bowtie2 -U raw/SampleA.fastq.gz | samtools sort -o mapped/.sorted.bam
    
    --- Final FASTQ array --- 
    total fastq files : raw/Input1.fastq.gz raw/Input2.fastq.gz raw/SampleA.fastq.gz



```bash
%%bash


samp_files=() #“Create an empty list that we will fill with fastq file paths.”
bam_files=()    # BAM array

for fq in raw/*.fastq.gz    
do
    samp_files+=("$fq")       #Append FG paths to the array
   echo "samples files Added: $fq" # optional print to see what is going 
   
    sample_name=$(basename "$fq" .fastq.gz)
   echo "sample_name of the file used : $sample_name"

     echo "  Command: bowtie2 -U $fq | samtools sort -o mapped/${sample_name}.sorted.bam"


    echo "  Command: bamCovergae -bam mapped/${sample_name}.sorted.bam -o mapped/${sample_name}.bw"

  

  


done


 echo "total fastq files : ${samp_files[@]}" 
  echo "total bamfiles : ${bam_files[@]}" 

echo "Command : multibamSummary -bam ${bam_files[@]} -o out.npz"

```

    samples files Added: raw/Input1.fastq.gz
    sample_name of the file used : Input1
      Command: bowtie2 -U raw/Input1.fastq.gz | samtools sort -o mapped/Input1.sorted.bam
      Command: bamCovergae -bam mapped/Input1.sorted.bam -o mapped/Input1.bw
    samples files Added: raw/Input2.fastq.gz
    sample_name of the file used : Input2
      Command: bowtie2 -U raw/Input2.fastq.gz | samtools sort -o mapped/Input2.sorted.bam
      Command: bamCovergae -bam mapped/Input2.sorted.bam -o mapped/Input2.bw
    samples files Added: raw/SampleA.fastq.gz
    sample_name of the file used : SampleA
      Command: bowtie2 -U raw/SampleA.fastq.gz | samtools sort -o mapped/SampleA.sorted.bam
      Command: bamCovergae -bam mapped/SampleA.sorted.bam -o mapped/SampleA.bw
    total fastq files : raw/Input1.fastq.gz raw/Input2.fastq.gz raw/SampleA.fastq.gz
    total bamfiles : mapped/Input1.sorted.bam mapped/Input2.sorted.bam mapped/SampleA.sorted.bam
    Command : multibamSummary -bam mapped/Input1.sorted.bam mapped/Input2.sorted.bam mapped/SampleA.sorted.bam -o out.npz


### 3.3 Adding Labels
We can also collect a separate list of labels (e.g., "Treatment" vs "Control"). The concept is the same: one basket for files, another basket for labels.


```bash
%%bash


samp_files=() #“Create an empty list that we will fill with fastq file paths.”
bam_files=()    # BAM array

for fq in raw/*.fastq.gz    
do
    samp_files+=("$fq")       #Append FG paths to the array
   echo "samples files Added: $fq" # optional print to see what is going 
   
    sample_name=$(basename "$fq" .fastq.gz)
   echo "sample_name of the file used : $sample_name"

     echo "  Command: bowtie2 -U $fq | samtools sort -o mapped/${sample_name}.sorted.bam"


    echo "  Command: bamCovergae -bam mapped/${sample_name}.sorted.bam -o mapped/${sample_name}.bw"

  
    bam="mapped/${sample_name}.sorted.bam"
    bam_files+=("$bam")

  


done


 echo "total fastq files : ${samp_files[@]}" 
  echo "total bamfiles : ${bam_files[@]}" 

echo "Command : multibamSummary -bam ${bam_files[@]} -o out.npz"

```

    samples files Added: raw/Input1.fastq.gz
    sample_name of the file used : Input1
      Command: bowtie2 -U raw/Input1.fastq.gz | samtools sort -o mapped/Input1.sorted.bam
      Command: bamCovergae -bam mapped/Input1.sorted.bam -o mapped/Input1.bw
    samples files Added: raw/Input2.fastq.gz
    sample_name of the file used : Input2
      Command: bowtie2 -U raw/Input2.fastq.gz | samtools sort -o mapped/Input2.sorted.bam
      Command: bamCovergae -bam mapped/Input2.sorted.bam -o mapped/Input2.bw
    samples files Added: raw/SampleA.fastq.gz
    sample_name of the file used : SampleA
      Command: bowtie2 -U raw/SampleA.fastq.gz | samtools sort -o mapped/SampleA.sorted.bam
      Command: bamCovergae -bam mapped/SampleA.sorted.bam -o mapped/SampleA.bw
    total fastq files : raw/Input1.fastq.gz raw/Input2.fastq.gz raw/SampleA.fastq.gz
    total bamfiles : mapped/Input1.sorted.bam mapped/Input2.sorted.bam mapped/SampleA.sorted.bam
    Command : multibamSummary -bam mapped/Input1.sorted.bam mapped/Input2.sorted.bam mapped/SampleA.sorted.bam -o out.npz


**Result**: We have a single command at the end that includes ALL our files. This is essential for summary plots.


```bash
%%bash


samp_files=() #“Create an empty list that we will fill with fastq file paths.”
bam_files=()    # BAM array

for fq in raw/*.fastq.gz    
do
   
    sample_name=$(basename "$fq" .fastq.gz)
    echo "sample_name of the file used : $sample_name"

     echo "  Command: bowtie2 -U $fq | samtools sort -o mapped/${sample_name}.sorted.bam"


    echo "  Command: bamCovergae -bam mapped/${sample_name}.sorted.bam -o mapped/${sample_name}.bw"

  
    bam="mapped/${sample_name}.sorted.bam"
    bam_files+=("$bam")

    labels+=("$sample_name") # add label to array


done

echo "total bamfiles : ${bam_files[@]}" 



echo "Command : multiBamSummary bins \
    --bamfiles "${bam_files[@]}" \
    --labels "${labels[@]}" \
    -o multiBam.npz "

```

    sample_name of the file used : Input1
      Command: bowtie2 -U raw/Input1.fastq.gz | samtools sort -o mapped/Input1.sorted.bam
      Command: bamCovergae -bam mapped/Input1.sorted.bam -o mapped/Input1.bw
    sample_name of the file used : Input2
      Command: bowtie2 -U raw/Input2.fastq.gz | samtools sort -o mapped/Input2.sorted.bam
      Command: bamCovergae -bam mapped/Input2.sorted.bam -o mapped/Input2.bw
    sample_name of the file used : SampleA
      Command: bowtie2 -U raw/SampleA.fastq.gz | samtools sort -o mapped/SampleA.sorted.bam
      Command: bamCovergae -bam mapped/SampleA.sorted.bam -o mapped/SampleA.bw
    total bamfiles : mapped/Input1.sorted.bam mapped/Input2.sorted.bam mapped/SampleA.sorted.bam
    Command : multiBamSummary bins     --bamfiles mapped/Input1.sorted.bam mapped/Input2.sorted.bam mapped/SampleA.sorted.bam     --labels Input1 Input2 SampleA     -o multiBam.npz 


### 3.4 Production Scale (Reading from Files)

If you have 500 samples, writing them in a script is messy. It's better to keep a separate text file (e.g., `samples.txt`) and have your robot read that.


```bash

%%bash

# Create a dummy samples.txt for demonstration
ls raw/ > samples.txt

echo "samples"
cat samples.txt

sed 's/.fastq.gz//' samples.txt> sample_id.txt 

echo "sample_id"
cat sample_id.txt

```

    samples
    Input1.fastq.gz
    Input2.fastq.gz
    SampleA.fastq.gz
    sample_id
    Input1
    Input2
    SampleA


Now we use a `while read` loop. This is safer than a `for` loop because it handles spaces and special characters correctly.


```bash
%%bash

echo -e "\n--- Looping over File List ---"
while read sample
do
    echo "Reading from file: $sample"
done < sample_id.txt

```

    
    --- Looping over File List ---
    Reading from file: Input1
    Reading from file: Input2
    Reading from file: SampleA


**Logic**: The code `< sample_id.txt` feeds the file into the loop, one line at a time.

### 3.5 Example: Running the List
See how the loop reads `Input1`, then `Input2`, then `SampleA`? We can just plug our `bowtie2` command right inside this loop.


```bash
%%bash

echo -e "\n--- Looping over File List ---"

while read -r sample
do
    echo "sample_id: $sample"

    fq="${sample}.fastq.gz"

     echo "sample: $fq"

    echo " Command: bowtie2 -U $fq | samtools sort -o mapped/${sample}.sorted.bam"
     
done < sample_id.txt



```

    
    --- Looping over File List ---
    sample_id: Input1
    sample: Input1.fastq.gz
     Command: bowtie2 -U Input1.fastq.gz | samtools sort -o mapped/Input1.sorted.bam
    sample_id: Input2
    sample: Input2.fastq.gz
     Command: bowtie2 -U Input2.fastq.gz | samtools sort -o mapped/Input2.sorted.bam
    sample_id: SampleA
    sample: SampleA.fastq.gz
     Command: bowtie2 -U SampleA.fastq.gz | samtools sort -o mapped/SampleA.sorted.bam


### 3.6 Paired-End Data (R1 and R2)
For paired-end sequencing, you don't just have one file. You have two (`_R1` and `_R2`).

Since we know the **Sample ID** (e.g., `Input1`), we can just tell the robot: "Look for `Input1` plus `_R1` and `Input1` plus `_R2`."


```bash
%%bash

while read -r sample
do
    echo "sample_id: $sample"

        fq1="${sample}_R1.fastq.gz"
        fq2="${sample}_R2.fastq.gz"


 echo " paired end -  $fq1 : $fq2"



    echo " Command: bowtie2 -1 $fq1 -2 $fq2 | samtools sort -o mapped/${sample}.sorted.bam"

done < sample_id.txt

```

    sample_id: Input1
     paired end -  Input1_R1.fastq.gz : Input1_R2.fastq.gz
     Command: bowtie2 -1 Input1_R1.fastq.gz -2 Input1_R2.fastq.gz | samtools sort -o mapped/Input1.sorted.bam
    sample_id: Input2
     paired end -  Input2_R1.fastq.gz : Input2_R2.fastq.gz
     Command: bowtie2 -1 Input2_R1.fastq.gz -2 Input2_R2.fastq.gz | samtools sort -o mapped/Input2.sorted.bam
    sample_id: SampleA
     paired end -  SampleA_R1.fastq.gz : SampleA_R2.fastq.gz
     Command: bowtie2 -1 SampleA_R1.fastq.gz -2 SampleA_R2.fastq.gz | samtools sort -o mapped/SampleA.sorted.bam


**Scientific Note**: This logic works for almost all modern Illumina data. Just check the file suffix!


```bash
%%bash

for sample in $(cat sample_id.txt); do

  echo "sample_id: $sample"

  fq1="${sample}_R1_val_1.fastq.gz"
  fq2="${sample}_R2_val_2.fastq.gz"

  echo "paired end:  $fq1  :  $fq2"

  echo "Command: bowtie2 -x index \
    -1 $fq1 \
    -2 $fq2 \
    -p 6 --no-unal \
    2> bowalign/${sample}.log | samtools sort -@ 6 -o bowalign/${sample}.sorted.bam"

  echo ""   # blank line for readability

done 
```

    sample_id: Input1
    paired end:  Input1_R1_val_1.fastq.gz  :  Input1_R2_val_2.fastq.gz
    Command: bowtie2 -x index     -1 Input1_R1_val_1.fastq.gz     -2 Input1_R2_val_2.fastq.gz     -p 6 --no-unal     2> bowalign/Input1.log | samtools sort -@ 6 -o bowalign/Input1.sorted.bam
    
    sample_id: Input2
    paired end:  Input2_R1_val_1.fastq.gz  :  Input2_R2_val_2.fastq.gz
    Command: bowtie2 -x index     -1 Input2_R1_val_1.fastq.gz     -2 Input2_R2_val_2.fastq.gz     -p 6 --no-unal     2> bowalign/Input2.log | samtools sort -@ 6 -o bowalign/Input2.sorted.bam
    
    sample_id: SampleA
    paired end:  SampleA_R1_val_1.fastq.gz  :  SampleA_R2_val_2.fastq.gz
    Command: bowtie2 -x index     -1 SampleA_R1_val_1.fastq.gz     -2 SampleA_R2_val_2.fastq.gz     -p 6 --no-unal     2> bowalign/SampleA.log | samtools sort -@ 6 -o bowalign/SampleA.sorted.bam
    


### 3.7 Pro Tip: `mapfile`
If you want to be extra fancy (and efficient), you can use `mapfile` to load the entire list into memory instantly. It's faster for huge lists.


```bash
%%bash

# create array from file
mapfile -t samples < sample_id.txt

for sample in "${samples[@]}"; do

  echo "sample_id: $sample"

  fq1="${sample}_R1_val_1.fastq.gz"
  fq2="${sample}_R2_val_2.fastq.gz"

  echo "paired end:  $fq1  :  $fq2"

  echo "Command: bowtie2 -x index \
    -1 $fq1 \
    -2 $fq2 \
    -p 6 --no-unal \
    2> bowalign/${sample}.log | samtools sort -@ 6 -o bowalign/${sample}.sorted.bam"

  echo ""

done

```

    sample_id: Input1
    paired end:  Input1_R1_val_1.fastq.gz  :  Input1_R2_val_2.fastq.gz
    Command: bowtie2 -x index     -1 Input1_R1_val_1.fastq.gz     -2 Input1_R2_val_2.fastq.gz     -p 6 --no-unal     2> bowalign/Input1.log | samtools sort -@ 6 -o bowalign/Input1.sorted.bam
    
    sample_id: Input2
    paired end:  Input2_R1_val_1.fastq.gz  :  Input2_R2_val_2.fastq.gz
    Command: bowtie2 -x index     -1 Input2_R1_val_1.fastq.gz     -2 Input2_R2_val_2.fastq.gz     -p 6 --no-unal     2> bowalign/Input2.log | samtools sort -@ 6 -o bowalign/Input2.sorted.bam
    
    sample_id: SampleA
    paired end:  SampleA_R1_val_1.fastq.gz  :  SampleA_R2_val_2.fastq.gz
    Command: bowtie2 -x index     -1 SampleA_R1_val_1.fastq.gz     -2 SampleA_R2_val_2.fastq.gz     -p 6 --no-unal     2> bowalign/SampleA.log | samtools sort -@ 6 -o bowalign/SampleA.sorted.bam
    


**Efficiency**: This is how the pros do it for cleaner code.


```bash
%%bash

# create array from file
mapfile -t samples < sample_id.txt

bams=()


for sample in "${samples[@]}"; do

  echo "sample_id: $sample"

  fq1="${sample}_R1_val_1.fastq.gz"
  fq2="${sample}_R2_val_2.fastq.gz"

  echo "paired end:  $fq1  :  $fq2"

  echo "Command: bowtie2 -x index \
    -1 $fq1 \
    -2 $fq2 \
    -p 6 --no-unal \
    2> bowalign/${sample}.log | samtools sort -@ 6 -o bowalign/${sample}.sorted.bam"


  bam="bowalign/${sample}.sorted.bam"
  bams+=( "$bam" )




done


echo "BAM array:"
printf '%s\n' "${bams[@]}"



echo "Command : multibamSummary -bam "${bams[@]}" -o out.npz"
```

    sample_id: Input1
    paired end:  Input1_R1_val_1.fastq.gz  :  Input1_R2_val_2.fastq.gz
    Command: bowtie2 -x index     -1 Input1_R1_val_1.fastq.gz     -2 Input1_R2_val_2.fastq.gz     -p 6 --no-unal     2> bowalign/Input1.log | samtools sort -@ 6 -o bowalign/Input1.sorted.bam
    sample_id: Input2
    paired end:  Input2_R1_val_1.fastq.gz  :  Input2_R2_val_2.fastq.gz
    Command: bowtie2 -x index     -1 Input2_R1_val_1.fastq.gz     -2 Input2_R2_val_2.fastq.gz     -p 6 --no-unal     2> bowalign/Input2.log | samtools sort -@ 6 -o bowalign/Input2.sorted.bam
    sample_id: SampleA
    paired end:  SampleA_R1_val_1.fastq.gz  :  SampleA_R2_val_2.fastq.gz
    Command: bowtie2 -x index     -1 SampleA_R1_val_1.fastq.gz     -2 SampleA_R2_val_2.fastq.gz     -p 6 --no-unal     2> bowalign/SampleA.log | samtools sort -@ 6 -o bowalign/SampleA.sorted.bam
    BAM array:
    bowalign/Input1.sorted.bam
    bowalign/Input2.sorted.bam
    bowalign/SampleA.sorted.bam
    Command : multibamSummary -bam bowalign/Input1.sorted.bam bowalign/Input2.sorted.bam bowalign/SampleA.sorted.bam -o out.npz


## Conclusion

You have now mastered the art of the **Digital Robot**. 
1.  **Loops** (`for`, `while`) repeat actions.
2.  **Variables** (`$file`, `$sample`) hold data.
3.  **Arrays** store lists of files for group analysis.

Use these powers to automate your pipelines!
