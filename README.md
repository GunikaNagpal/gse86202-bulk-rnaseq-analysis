🧬 Bulk RNA-seq Mini Analysis — GSE86202 (Glioblastoma)

This project performs a clean and minimal bulk RNA-Seq differential expression analysis using the publicly available dataset GSE86202 (Glioblastoma).
It includes the essential steps normally shown in real RNA-seq projects:
	•	Loading differential expression (DEG) results
	•	Cleaning and preprocessing
	•	Volcano Plot
	•	Lollipop Plot of most significant genes
	•	Distribution plots of log2FC & p-values
	•	Ranked gene list for GSEA

The goal is to create a simple, replicable, visually clear mini-project that demonstrates bioinformatics analysis skills in Python.
🧪 Dataset Information
	•	Accession ID: GSE86202
	•	Organism: Homo sapiens
	•	Study Type: Glioblastoma – Differential Gene Expression
	•	Source: NCBI GEO database

We use the differential expression file that contains:
	•	Gene identifiers
	•	log2 fold changes
	•	p-values


🚀 Analysis Steps Covered

1️⃣ Load & Clean Differential Expression Table
	•	Read Excel
	•	Normalize column names
	•	Detect columns: log2FC, p-value, gene

2️⃣ Volcano Plot

A volcano plot helps visualize the global DEG pattern:
	•	X-axis: log2 Fold Change
	•	Y-axis: −log10(p-value)
	•	Highlight biological significance cutoff lines

📌 Output: volcano_gse86202.png

3️⃣ Top Up & Downregulated Genes

Creates two barplots:
	•	Top 10 Upregulated Genes
	•	Top 10 Downregulated Genes

📌 Output:
	•	top10_upregulated.png
	•	top10_downregulated.png

4️⃣ Lollipop Plot — Top 20 Most Significant Genes

Shows the 20 most significant genes ranked by p-value.

📌 Output: lollipop_top20.png

5️⃣ Distribution Plots

Essential QC-style plots showing:
	•	Distribution of log2 fold change
	•	Distribution of −log10(p-value)

📌 Output: distributions.png

6️⃣ Ranked Gene List for GSEA

This is compatible with:
	•	GSEA Desktop
	•	GSEApy
	•	Enrichr preranked analysis

📌 Output: ranked_gene_list.rnk


 How to Run the Script

Prerequisites

Install required packages:
pip install pandas numpy matplotlib seaborn openpyxl

Run the script

Navigate to the project folder and run:
python scripts/analyse_gse86202.py

All results will appear inside the outputs/ folder.

📊 Example Outputs

✔ Volcano Plot

Shows overall up/down regulation trends.

✔ Lollipop Plot

Highlights top 20 most significant genes.

✔ Distribution Plots

Shows fold change & p-value distributions.

These plots together give a complete quick summary of the dataset.
