---

# Trichosanthes pilosa Sex Determination Analysis Pipeline

This repository provides the full analytical workflow used to investigate the sex determination mechanism in *Trichosanthes pilosa*.
It includes deep learning–based candidate gene prioritization, WGCNA module analysis, sex-linked marker identification, and small RNA (especially hc-siRNA) investigation. All steps are fully reproducible.

---

## 📌 Pipeline Overview

---

The full workflow consists of the following four components (in the exact order they were executed):

1. **Deep learning** for prioritizing candidate sex-determining genes
2. **WGCNA analyses** identifying sex-associated modules
3. **Identification of sex-linked molecular markers**
4. **Small RNA mapping and hc-siRNA analysis** revealing potential regulatory roles

---

## 🔬 1. Deep Learning Framework for Candidate Gene Prioritization

---

This module predicts the probability that each SDR-region gene participates in sex determination.

### **Key Features**

* Uses transcriptome expression matrix (M1/M2/F1/F2) 
* Multi-layer perceptron (MLP) architecture
* Outputs a **Sex Determination Probability (SDP)** per gene
* Includes loss curve, ROC metrics, PR AUC

### **Example command**

```
python train.py --input combined_training.tsv
python predict.py --input tpi.txt
```

---

## 🧬 2. WGCNA: Identification of Sex-Associated Co-expression Modules

---

The WGCNA workflow covers all **14 analytical steps** (described here in a concise summary).

### **Major Steps (condensed summary)**

* Filtering, normalization, and batch correction
* Soft-threshold selection
* Construction of weighted co-expression network (TOM)
* Module detection via dynamic tree cutting
* Module–trait correlation analysis (M/F phenotype)
* Extraction of sex-associated modules (p < 0.05)
* Hub gene selection and module visualization

### **Example command**

```
Rscript WGCNA.R 
```

---

## 🧪 3. Identification of Sex-Linked Molecular Markers

---

Sex-linked markers were identified based on SDR-region SNPs and population resequencing data.

### **Key Steps**

* SNP calling within the SDR
* Screening male-specific or female-absent variants
* Primer design with dimer-free criteria
* PCR verification in wild-collected individuals
* Final selection of **male-specific marker combinations**

### **Example command**

```
bash coverage_SDR.sh
```

---

## 🧫 4. Small RNA Mapping and hc-siRNA Characterization

---

This module uncovers potential regulatory small RNAs—especially **high-confidence siRNAs (hc-siRNAs)**—involved in sex differentiation.

### **Analytical Components**

* Adapter trimming and quality control
* Size-based separation of 21-, 22-, and 24-nt sRNAs
* Coverage profiling of sRNAs over SDR genes
* Identification of hc-siRNAs inversely correlated with candidate gene expression
* Construction of a *candidate gene–hc-siRNA* regulatory network

### **Example command**

Please refer to the software manual for details: [Shortstack](https://github.com/MikeAxtell/ShortStack),[sRNAminer](https://github.com/kli28/sRNAminer)
```
ShortStack --genomefile female.fa --known_miRNAs mature.fa --readfile F1_cutadapt.fq F2_cutadapt.fq M1_cutadapt.fq M2_cutadapt.fq --outdir shortstack --threads 64 --dn_mirna
sRNAminer One_step_sRNAminer -p 32 -g female.fa -i F1_cutadapt.fq,F2_cutadapt.fq,M1_cutadapt.fq,M2_cutadapt.fq --nd Rfam_withoutMIR.fa --od Organelle.genomic.fa --miRNA --PHAS21 --PHAS24 --hc_siRNA -o srnaminer_result/
```

---

## 📄 Citation

---

If you use this repository or its methods in your research, please cite our forthcoming manuscript.

---

## 📬 Contact

---

For questions, missing scripts, or data requests, please contact:

**Trichosanthes pilosa Sex Determination Project Team**
Email: *pengzhen22@mails.ucas.ac.cn*

---
