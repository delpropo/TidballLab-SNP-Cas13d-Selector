
# Cas13d gRNA/SNP Pipeline

## Overview
This pipeline identifies heterozygous SNPs in coding and 3' UTR regions and designs optimal Cas13d gRNAs for each SNP, outputting a filterable Excel table.

## Objective

- **Input:**
    *100x WGS (short and linked reads), VCF, human reference (GRCh38), comprehensive GTF annotation.*
- **Task:**
    *Identify all heterozygous SNPs in CDS and 3’ UTR regions genome-wide. For each SNP, output the highest-scoring Cas13d gRNA that overlaps the SNP (with the SNP ideally in the last 10 bp).*
- **Output:**
    *Single Excel with one row per SNP/gRNA, filterable by gene.*

## Dependencies
- Python ≥3.8, pandas, pysam, biopython, openpyxl
- fastp, bwa-mem2, bedtools, Cas13design CLI
- Reference genome (GRCh38), GTF annotation, VCF file

## Inputs
- 100x WGS (short/linked reads)
- VCF file
- Human reference genome (GRCh38)
- GTF annotation

## Outputs
- Excel file: Cas13d_CDS_UTR3_HetSNP_gRNAs.xlsx

## Customization
- Adjust quality/depth thresholds in variant parsing
- Change gRNA scoring parameters as needed

## Limitations
- Assumes VCF is indexed and compatible with pysam
- GTF annotation must match reference genome

## Pipeline Steps

## **Pipeline Steps**

### 1. **Preprocessing**

> *Standard steps as before:*
- Quality Control: `fastp`
- Alignment: `bwa-mem2`
- Variant Calling
- Phasing (with linked reads for best accuracy)
- QC as previously described

---


### 2. **Annotation: Extract CDS and 3' UTR Intervals**

<!-- Example code: adjust as needed for your annotation source and requirements -->
```python
import pandas as pd

# Load GTF annotation
GTF = pd.read_csv("gencode.v38.annotation.gtf", sep='\t', comment='#', header=None,
                  names=['chr','source','feature','start','end','score','strand','frame','info'])

# Extract CDS and 3' UTR features
cds = GTF[GTF['feature'] == 'CDS']
utr3 = GTF[GTF['feature'] == 'three_prime_utr']
cds_utr3 = pd.concat([cds, utr3])

# Export BED for sequence extraction
cds_utr3[['chr','start','end','info']].to_csv('cds_utr3.bed', sep='\t', header=False, index=False)
```

**Extract sequences (BEDTools):**
```bash
bedtools getfasta -fi hg38.fa -bed cds_utr3.bed -fo cds_utr3.fa -name
```

---


### 3. **Identify Heterozygous SNPs in CDS/3’ UTR (from VCF)**

<!-- Example code: adjust thresholds and logic as needed -->
```python
import pysam

vcf = pysam.VariantFile('WTC11.vcf.gz')

def parse_variants(chrom, start, end):
    for rec in vcf.fetch(chrom, start, end):
        # Het, quality, depth
        if rec.filter.keys()[0] == 'PASS' and rec.qual >= 30 and rec.info.get('DP',0) >= 20:
            gt = rec.samples[0]['GT']
            if gt == (0,1) or gt == (1,0):
                yield rec.pos, rec.ref, rec.alts[0], rec.id
```
Loop over all exons (from BED/CDS_UTR3) and collect any overlapping het SNPs.

---


### 4. **Cas13d gRNA Extraction & Scoring**

#### **A. For each SNP (in CDS or 3' UTR):**
- **Extract context sequence:** Get ~50bp window around the SNP (exonic only).
- **Enumerate all 30bp gRNAs** overlapping the SNP, with SNP in positions 21–30 (last 10 nt).
- **Filter gRNAs:** No long homopolymers; GC 35–60%.
- **Batch SCORE:** Use `Cas13design` (CLI) for activity scores.

<!-- Example code: adjust logic and parameters as needed -->
```python
from Bio import SeqIO
import subprocess

def get_candidate_gRNAs(seq, snp_rel_pos):
    window_size = 30
    grnas = []
    for i in range(max(0, snp_rel_pos-29), min(len(seq)-window_size+1, snp_rel_pos+1)):
        window = seq[i:i+window_size]
        snp_in_gRNA_pos = snp_rel_pos - i
        if snp_in_gRNA_pos >= 20:
            if 0.35 <= (window.count('G')+window.count('C'))/window_size <= 0.6 and not any(base*6 in window for base in "ATGC"):
                grnas.append((window, snp_in_gRNA_pos, i))
    return grnas

def cas13design_score_batch(grna_list, outfile="batch.txt"):
    # Write candidate gRNAs to file
    with open("cas13d_candidates.txt", "w") as fout:
        for idx, (gseq,_,_) in enumerate(grna_list):
            fout.write(f">{idx}\n{gseq}\n")
    # Run Cas13design CLI
    subprocess.run("cas13design batch --input cas13d_candidates.txt --output scores.txt", shell=True)
    # Parse scores
    scores = {}
    with open("scores.txt") as fin:
        for line in fin:
            if line.startswith("Score"): continue
            idx, score = line.strip().split()
            scores[int(idx)] = float(score)
    return [scores[i] for i in range(len(grna_list))]
```


#### **B. Output: Highest-scoring gRNA per SNP**
- For each SNP, choose the candidate with the top `Cas13design` score.
- Output all relevant metadata.

---


### 5. **Output Excel Design**

**Columns:**

| Chrom | Position | Ref | Alt | Type (CDS/UTR) | Gene | Transcript | Exon | gRNA Sequence | SNP pos in gRNA | gRNA start (genome) | Score | Context ±10bp | GC% | Homopolymer (Y/N) |

Each row = 1 SNP, 1 best-scored gRNA.

<!-- Example code: adjust output columns and logic as needed -->
```python
import pandas as pd

df = pd.DataFrame(results_list)
df.to_excel("Cas13d_CDS_UTR3_HetSNP_gRNAs.xlsx")
```

---


### 6. **Summary Table: Workflow**

| Step           | Tool/Library      | Purpose                  |
|----------------|-------------------|--------------------------|
| Preprocessing  | fastp, bwa-mem2   | QC, mapping, variant calling |
| Annotation     | GENCODE, pandas   | CDS, 3’ UTR intervals    |
| SNP parsing    | pysam             | Extract all relevant SNPs|
| gRNA finding   | biopython, pandas | Enumerate/filter gRNAs   |
| gRNA scoring   | Cas13design       | Score/edit gRNAs         |
| Output         | pandas, openpyxl  | Final Excel              |

---


## Searchability and Scaling
- Output Excel can be filtered for any gene, transcript, or region.
- For large genome-scale results, output to CSV and partition by chromosome as needed.
- Pipeline supports batch/cloud/HPC/hybrid processing.

---


## In Summary

You will produce a comprehensive genome-wide gRNA/SNP lookup table:
- **All heterozygous SNPs in CDS/3’UTR**
- **Best available Cas13d gRNA overlapping each SNP (with full annotations for downstream querying)**

## References
- [Cas13design](https://github.com/your-link)
- [GENCODE](https://www.gencodegenes.org/)
- [Snakemake](https://snakemake.readthedocs.io/en/stable/)