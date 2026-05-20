# 🚀 Step-by-Step Guide: Analyzing Non-poly(A) RNAs from CLIP-seq and Ribosome Footprinting Data

Because standard RNA-seq libraries often rely on oligo-dT enrichment, RNAs lacking a poly(A) tail (such as replication-dependent histones, rRNAs, tRNAs, and snoRNAs) are severely underrepresented. Consequently, standard translation efficiency calculations (RPF / RNA-seq) fail or produce artifacts for these transcripts.

This guide outlines a bioinformatic pipeline to analyze non-poly(A) RNAs directly using your Ribosome Footprinting (RPF) and CLIP-seq BAM files, bypassing the flawed RNA-seq baseline.

---

## **Step 1: Data Preparation & Quality Control**
Before diving into transcript-specific analysis, ensure your BAM files are properly sorted and indexed.

* **Objective:** Prepare BAM files for rapid coordinate-based querying and ensure high mapping quality.
* **Bioinformatic Tools:** * `samtools` (for sorting and indexing)
    * `MultiQC` (to aggregate QC metrics)
* **Commands:**
    ```bash
    samtools sort -o CLIP-35L33G.sorted.bam CLIP-35L33G.bam
    samtools index CLIP-35L33G.sorted.bam
    # Repeat for RPF-siLuc.bam and RPF-siLin28a.bam
    ```
* **Visualization:** Use `MultiQC` HTML reports to visualize mapping rates and ensure the files are comparable.

---

## **Step 2: Custom Reference Annotation**
Standard GTF/GFF files often prioritize protein-coding genes. You need to ensure non-coding and non-poly(A) architectural RNAs are properly annotated.

* **Objective:** Create or filter a reference annotation file (.gtf/.bed) that explicitly contains tRNAs, rRNAs, snoRNAs, and histone genes.
* **Bioinformatic Tools:**
    * `BEDTools` (to manipulate genomic intervals)
    * UCSC Table Browser or Ensembl (to download specific biotypes)
* **Action:** Extract coordinates for histone gene clusters (e.g., `HIST1` cluster) and structural RNAs into a custom `non_polyA_targets.bed` file.

---

## **Step 3: Direct Quantification of Histone Translation**
Instead of dividing RPF reads by RNA-seq reads (which are depleted of histones), we quantify absolute ribosome occupancy and normalize against global translation metrics.

* **Objective:** Determine if LIN28A depletion changes the translational output of histone mRNAs.
* **Bioinformatic Tools:**
    * `featureCounts` (from the Subread package)
    * `DESeq2` (in R, for normalization)
* **Commands:**
    ```bash
    featureCounts -T 4 -a non_polyA_targets.gtf -o histone_rpf_counts.txt RPF-siLuc.sorted.bam RPF-siLin28a.sorted.bam
    ```
* **Analysis:** Import `histone_rpf_counts.txt` into R. Instead of standard library size normalization, calculate size factors using a set of stable, invariant housekeeping genes to act as an internal control.
* **Visualization:** * **Boxplots / Violin Plots:** Generated via `ggplot2` in R or `seaborn` in Python, showing the log2 normalized RPF counts for histone genes in siLuc vs. siLin28a.
    * **MA-plots:** To show the fold-change of translation specifically for histone genes compared to the rest of the transcriptome.

---

## **Step 4: Peak Calling on Structural RNAs (tRNAs, rRNAs)**
LIN28A binds non-coding architectural RNAs. We will identify precise binding sites by looking for crosslinking events (peaks) in the CLIP-seq data.

* **Objective:** Identify direct LIN28A binding sites on rRNAs, tRNAs, and snoRNAs without relying on RNA-seq for background normalization.
* **Bioinformatic Tools:**
    * `PureCLIP` or `Piranha` (Peak callers designed for CLIP-seq data)
    * `bedtools getfasta` (to extract sequences under the peaks)
    * `MEME Suite` (for motif discovery)
* **Commands:**
    ```bash
    # Call peaks using PureCLIP
    pureclip -i CLIP-35L33G.sorted.bam -bai CLIP-35L33G.sorted.bam.bai -g mm9.fasta -iv non_polyA_targets.bed -o lin28a_nonpolyA_peaks.bed
    
    # Extract sequences for motif analysis
    bedtools getfasta -fi mm9.fasta -bed lin28a_nonpolyA_peaks.bed -fo lin28a_peaks.fasta
    ```
* **Visualization:**
    * **Motif Logos:** Use `WebLogo` or `ggseqlogo` (R) to visualize if the identified peaks on tRNAs/rRNAs contain the canonical AAGNNG or UGUG motifs.
    * **Genome Browser Tracks:** Use **IGV (Integrative Genomics Viewer)**. Load the BAM files and the peak `.bed` files. Zoom in on tRNA and rRNA loci to visually inspect read pileups and crosslinking-induced mutations (G-to-A or C-to-T transitions).

---

## **Step 5: Meta-Gene Mapping of Histone Transcripts**
Histone mRNAs lack a poly(A) tail and instead terminate in a conserved 3' stem-loop structure. We want to see if LIN28A specifically targets this structure.

* **Objective:** Map the spatial distribution of LIN28A CLIP tags along the body of histone transcripts.
* **Bioinformatic Tools:**
    * `deepTools` (`computeMatrix` and `plotProfile`)
* **Commands:**
    ```bash
    # Convert BAM to BigWig for coverage visualization
    bamCoverage -b CLIP-35L33G.sorted.bam -o CLIP-35L33G.bw --normalizeUsing CPM
    
    # Compute the matrix over histone genes
    computeMatrix scale-regions -S CLIP-35L33G.bw -R histones.bed -b 200 -a 200 -o histone_matrix.gz
    
    # Plot the metagene profile
    plotProfile -m histone_matrix.gz -out histone_profile.pdf --plotTitle "LIN28A Binding on Histone mRNAs"
    ```
* **Visualization:**
    * **Metagene Profile Plots:** A line chart (generated by `plotProfile`) showing the average read density across scaled histone transcripts (from 5' UTR to the 3' stem-loop). A sharp peak near the 3' end would indicate targeted binding to the stem-loop.
