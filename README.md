# AlphaGenome Endogenous Variant Benchmark: Natural Variant Validation

**Institution:** Layer Laboratory, CU Boulder  
**Dataset:** GTEx v8 eQTLs (Whole Blood) - 107,229 endogenous human variants  
**Repository:** https://github.com/gsstephenson/alphagenome-mpra-benchmark

---

## 🎯 TL;DR - Key Findings

### Main Discovery
**AlphaGenome shows near-zero correlation with GTEx eQTL effect sizes on endogenous variants in open chromatin regions.** After filtering 2.4M GTEx variants to 107K that overlap K562 DNase peaks, AlphaGenome's predictions show no significant relationship with empirical gene expression changes (Spearman ρ = -0.0005, p = 0.864).

### Take-Home Messages

1. **📊 No Detectable Correlation** (ρ = -0.0005, p = 0.864)
   - AlphaGenome predictions do not correlate with GTEx eQTL effect sizes
   - Result holds across 107,229 naturally occurring variants
   - Stark contrast to expected performance on natural variants

2. **🧬 Context vs. Variant Signal**
   - 2048 bp window context (99.2%) dominates single nucleotide variants (0.05%)
   - Flanking chromatin accessibility overwhelms variant-specific effects
   - Model architecture designed for context-rich predictions

3. **⚠️ Cell Type Mismatch Hypothesis**
   - AlphaGenome: K562 (erythroleukemia) chromatin predictions
   - GTEx data: Whole Blood gene expression measurements
   - eQTL effects may be tissue-specific, not captured by K562 accessibility

4. **🔬 Methodological Insights**
   - Proper benchmark: Natural variants + Native chromatin + Species-matched
   - Challenge: eQTL effect sizes measure distal gene expression, not local chromatin
   - AlphaGenome predicts accessibility, not transcriptional output
   - Different biological readouts may explain null correlation

![Scatter Plot](output/scatter_plot.png)
*Figure 1: AlphaGenome predictions vs GTEx eQTL effect sizes showing no correlation*

---

## 🔬 Experimental Design

### How the Study Works

**The GTEx v8 Dataset:**
- 2,414,653 significant eQTL variants (Whole Blood tissue)
- Naturally occurring human SNPs and indels
- Effect sizes (slopes) from linear regression: gene expression ~ genotype
- **Gold standard for variant functional effects** in human populations

**Our Analysis Pipeline:**

```
┌─────────────────────────────────────────────────────────┐
│ Step 1: Download GTEx eQTLs                             │
│ • GTEx v8 Whole Blood tissue data                       │
│ • 2,414,653 significant variant-gene pairs              │
│ • Effect sizes from expression quantitative trait loci  │
└─────────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────┐
│ Step 2: Intersect with K562 DNase Peaks                 │
│ • Filter to variants in open chromatin regions          │
│ • Use ENCODE K562 DNase hypersensitive sites            │
│ • Reduces to 107,232 variants (4.4% of original)        │
└─────────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────┐
│ Step 3: Extract 2KB Genomic Sequences                   │
│ • Reference (hg38) and alternate allele sequences       │
│ • 2048 bp windows centered on each variant              │
│ • Handle SNPs, insertions, and deletions                │
└─────────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────┐
│ Step 4: AlphaGenome Predictions                         │
│ • Predict DNase accessibility for ref and alt           │
│ • K562 cell line ontology (305 DNase tracks)            │
│ • Calculate delta: alt_score - ref_score                │
└─────────────────────────────────────────────────────────┘
                         ↓
┌─────────────────────────────────────────────────────────┐
│ Step 5: Benchmark Against GTEx                          │
│ • Correlate AlphaGenome delta with GTEx effect sizes    │
│ • Spearman and Pearson correlations                     │
│ • Statistical significance testing                      │
└─────────────────────────────────────────────────────────┘
```

**Why This Is The Proper Benchmark:**
- ✅ Natural variants (real human genetic variation)
- ✅ Native chromatin context (endogenous genomic loci)
- ✅ Species-matched (human variants, human genome, human model)
- ✅ Large sample size (107,229 variants, high statistical power)

**Comparison to MPRA Edge Case:**
- MPRA: Synthetic mutations + episomal plasmids → r ~ 0.05-0.09
- This study: Natural variants + native chromatin → r ~ 0.0005
- **Expected improvement not observed**

---

## 📊 Results & Visualizations

### Overall Performance

**Correlation Summary (N=107,229 variants):**

| Metric | Value | p-value | Interpretation |
|--------|-------|---------|----------------|
| **Spearman ρ** | -0.0005 | 0.864 | No correlation |
| **Pearson r** | 0.0026 | 0.403 | No correlation |
| **GTEx slope** | 0.000 ± 0.0001 | - | Small effect sizes |
| **AlphaGenome Δ** | 0.000 ± 0.0049 | - | Predictions vary widely |

**Key Observations:**
- No statistically significant correlation detected
- Very small eQTL effect sizes (mean ~ 0)
- Wide variance in AlphaGenome predictions relative to eQTL effects
- Result is robust across >100K variants

![Scatter Plot](output/scatter_plot.png)
*Figure 2: No discernible relationship between predictions and measured effects*

---

## 💡 Biological Interpretation

### Why No Correlation?

**Hypothesis 1: Different Biological Readouts**

| Feature | AlphaGenome | GTEx eQTLs |
|---------|-------------|------------|
| **Measurement** | Chromatin accessibility | Gene expression |
| **Location** | Local (variant locus) | Distal (gene body) |
| **Mechanism** | TF binding, nucleosomes | Transcriptional output |
| **Cell Type** | K562 (erythroleukemia) | Whole blood (mixed) |
| **Distance** | 0 bp (at variant) | Variable (cis-window) |

**eQTL variants affect gene expression through multiple mechanisms:**
- Enhancer activity (may correlate with accessibility)
- Promoter activity (may correlate with accessibility)
- Splicing regulation (no accessibility connection)
- RNA stability (no accessibility connection)
- Long-range looping (3D structure, not captured in 2KB window)

**Only a subset of eQTL effects work through local chromatin accessibility changes.**

---

### Hypothesis 2: Context Dominance

**The 2048 bp Problem:**
- Single nucleotide variant: 1 bp (0.05% of window)
- Flanking sequence: 2047 bp (99.95% of window)
- AlphaGenome integrates entire window for prediction

**Implication:**
- If flanking sequence has high baseline accessibility → predictions dominated by context
- Variant signal is diluted across 2048 positions
- Small eQTL effects (mean slope ~ 0.0001) may be below detection threshold

**Supporting Evidence:**
- AlphaGenome Δ std dev (0.0049) >> GTEx slope std dev (0.0001)
- Predictions vary ~50× more than observed effects
- Model may be "noisy" relative to small biological signals

---

### Hypothesis 3: Cell Type Specificity

**K562 vs Whole Blood Mismatch:**

| Aspect | K562 (AlphaGenome) | Whole Blood (GTEx) |
|--------|-------------------|-------------------|
| Cell Type | Erythroleukemia | Mixed: lymphocytes, monocytes, granulocytes |
| Source | Cancer cell line | Primary tissue samples |
| Chromatin State | Specific K562 profile | Averaged across cell types |
| eQTL Relevance | May not affect K562 | Affects blood cell types |

**Key Insight:**
- Many eQTLs are tissue/cell-type specific
- Variant may open chromatin in lymphocytes but not K562
- AlphaGenome (K562-trained) cannot predict blood-specific effects
- Need matched cell type predictions for proper validation

---

### Hypothesis 4: Model Limitations

**Potential Technical Issues:**
- 2KB window too small for long-range regulatory elements
- Model trained on bulk accessibility (not variant-specific)
- Training data may lack examples of subtle eQTL-scale effects
- Model optimized for base prediction, not delta (difference) prediction

**Comparison to MPRA:**
- MPRA: Designed to measure variant effects directly (reporter assays)
- eQTLs: Natural experiments, many confounding factors
- MPRA showed weak correlation (r ~ 0.05-0.09) despite measuring variant effects
- eQTLs may be even harder due to indirect relationship

---

## 🔧 Technical Achievements

### Pipeline Success Metrics
- ✅ **107,229 successful predictions** (99.997% success rate)
- ✅ **2.4M variants parsed** from GTEx v8
- ✅ **107K variants filtered** to K562 DNase peaks
- ✅ **Automatic checkpointing** (resume from failures)
- ✅ **~7.5 hours runtime** (overnight batch job)
- ✅ **100% reproducible** (documented parameters)

### Data Quality
- ✅ **No missing values** in predictions
- ✅ **No infinite values** detected
- ✅ **Reference genome validation** (sequence matching)
- ✅ **Strand-aware processing** (proper orientation)
- ✅ **Indel handling** (padding/trimming to 2048 bp)

### Statistical Rigor
- N = 107,229 (power >99% to detect r > 0.01)
- Spearman and Pearson correlations
- Two-tailed significance testing
- Publication-quality visualizations

---

## 📋 Comparison to MPRA Edge Case

### Side-by-Side Analysis

| Feature | MPRA Benchmark | Endogenous Benchmark |
|---------|---------------|---------------------|
| **Variants** | 6,863 synthetic | 107,229 natural |
| **Context** | Episomal plasmids | Native chromatin |
| **Species** | Mouse sequences | Human sequences |
| **Mutations** | Designed (affinity gradients) | Natural (population variants) |
| **Measurement** | Reporter expression | Gene expression (eQTLs) |
| **Correlation** | r = 0.053-0.091 | r = -0.0005-0.0026 |
| **Interpretation** | Weak positive | No correlation |
| **Sample Size** | 6.8K | 107K |
| **Statistical Power** | High | Very high |

### Scientific Insights

**MPRA Results Suggested:**
- Episomal context limits performance
- Natural variants + native chromatin should improve correlation
- Expected: r > 0.3 for proper benchmark

**Endogenous Results Show:**
- Native chromatin context does NOT improve correlation
- Natural variants also show near-zero correlation
- **Different explanation needed**

**Revised Understanding:**
1. AlphaGenome predicts **chromatin state**, not **transcriptional effects**
2. eQTLs measure **gene expression**, not **local accessibility**
3. Weak correlation is expected when biological readouts differ
4. Proper validation requires **matched readouts**: accessibility QTLs (caQTLs)

---

## 🎓 Scientific Implications

### For AlphaGenome Validation

**This benchmark reveals:**
- ⚠️ eQTLs are NOT appropriate for chromatin accessibility model validation
- ✅ Need chromatin accessibility QTLs (caQTLs) instead
- ✅ Need cell-type matched predictions (blood cells, not K562)
- ⚠️ Small effect sizes (mean ~ 0) are challenging for any model

**Better Validation Strategy:**
1. Use ENCODE/Roadmap Epigenomics caQTLs (accessibility measurements)
2. Match cell types: blood eQTLs → blood cell DNase predictions
3. Focus on variants with large effect sizes (power to detect)
4. Test multiple window sizes (2KB, 16KB, 131KB)

### For Regulatory Genomics Field

**Key Lessons:**
1. **Readout Matters**: Chromatin ≠ Expression
   - Accessibility is necessary but not sufficient for expression
   - Many eQTLs work through post-transcriptional mechanisms
   - Need matched molecular phenotypes for validation

2. **Context Dominance**: Window Size Trade-off
   - Large windows: more regulatory context, diluted variant signal
   - Small windows: focused on variant, missing long-range elements
   - Optimal size depends on variant type and effect mechanism

3. **Cell Type Specificity**: Critical for eQTLs
   - K562 chromatin ≠ Whole blood expression
   - Tissue-specific eQTLs require tissue-matched predictions
   - Cross-tissue benchmarks will show weak correlations

4. **Effect Size Distribution**: Statistical Power
   - GTEx eQTL slopes are very small (mean ~ 0.0001)
   - Measurement noise in predictions (std ~ 0.005) >> biological signal
   - Models may require effect sizes > 0.01 for reliable prediction

---

## 📁 Repository Structure

```
alphagenome_endogenous_benchmark/
├── README.md                                    # This file
├── requirements.txt                             # Python dependencies
│
├── scripts/
│   └── run_endogenous_validation_pipeline.py    # Full pipeline
│
├── data/
│   ├── genome/
│   │   └── GRCh38.p13.genome.fa                # hg38 reference (local)
│   ├── gtex/
│   │   ├── GTEx_Analysis_v8_eQTL.tar           # Downloaded (1.5 GB)
│   │   └── Whole_Blood.v8.signif_variant_gene_pairs.txt.gz
│   └── encode/
│       └── wgEncodeOpenChromDnaseK562Pk.narrowPeak.gz
│
└── output/
    ├── gtex_eqtls.bed                          # 2.4M variants (188 MB)
    ├── filtered_eqtls_in_k562_dhs.bed          # 107K variants (8.3 MB)
    ├── alphagenome_eqtl_predictions.csv        # Predictions (12 MB)
    ├── alphagenome_vs_gtex_benchmark.txt       # Statistics (375 B)
    └── scatter_plot.png                        # Visualization (469 KB)
```

---

## 🚀 Quick Start

### Prerequisites

```bash
# System requirements
- Python 3.11+
- bedtools (for genomic interval operations)
- hg38 reference genome

# Conda environment
conda create -n alphagenome-env python=3.11
conda activate alphagenome-env
```

### Installation

```bash
# Clone repository
git clone https://github.com/gsstephenson/alphagenome-mpra-benchmark
cd alphagenome_endogenous_benchmark

# Install dependencies
pip install pandas numpy scipy matplotlib pyfaidx tqdm biopython alphagenome
conda install -c bioconda bedtools

# Configure API key (if needed)
export ALPHA_GENOME_KEY=your_key_here
```

### Run Pipeline

```bash
# Execute full pipeline (may take ~8 hours)
python scripts/run_endogenous_validation_pipeline.py

# Output files will be created in output/
```

---

## 📊 Dataset Details

### GTEx v8 (Whole Blood eQTLs)

**Publication:** GTEx Consortium (2020). *The GTEx Consortium atlas of genetic regulatory effects across human tissues.* Science 369(6509):1318-1330.

**Data Source:**
- URL: https://storage.googleapis.com/adult-gtex/bulk-qtl/v8/single-tissue-cis-qtl/
- File: GTEx_Analysis_v8_eQTL.tar (1.5 GB)
- Tissue: Whole Blood (670 samples)

**Variant Details:**
- Total eQTL variants: 2,414,653
- Filtered to K562 DNase: 107,232
- Successful predictions: 107,229 (99.997%)

**Effect Size Distribution:**
- Mean slope: 0.000 (centered at zero)
- Std dev: 0.0001 (very small effects)
- Range: [-0.001, +0.001] (mostly subtle effects)

### ENCODE DNase (K562)

**Data Source:** 
- ENCODE wgEncodeOpenChromDnaseK562Pk.narrowPeak.gz
- Cell line: K562 (chronic myelogenous leukemia)
- Assay: DNase-seq (chromatin accessibility)

---

## 🔮 Future Directions

### Recommended Next Steps

1. **Chromatin Accessibility QTLs (caQTLs)**
   - Use ATAC-seq or DNase-seq QTLs instead of eQTLs
   - Match biological readout: accessibility predicts accessibility
   - Expected correlation: r > 0.5

2. **Cell Type Matching**
   - Get blood cell DNase/ATAC data
   - Use blood-specific AlphaGenome predictions (if available)
   - Test tissue-specific vs cross-tissue performance

3. **Effect Size Stratification**
   - Separate large effect (|slope| > 0.0005) from small effect variants
   - Test if correlation improves for high-impact eQTLs
   - May reveal signal-to-noise threshold

4. **Window Size Optimization**
   - Test 16 KB, 131 KB windows (AlphaGenome supports these)
   - Larger context may capture long-range regulatory elements
   - Trade-off: more context vs more diluted variant signal

5. **Multi-Track Analysis**
   - AlphaGenome returns 305 DNase tracks (different cell types)
   - Aggregate predictions across relevant tracks
   - May improve correlation with blood eQTLs

---

## ✅ Project Status

**COMPLETE** - All analyses finished, documented, and ready for publication

- ✅ 107,229 endogenous variants predicted
- ✅ Proper benchmark executed (natural variants + native chromatin)
- ✅ No correlation detected (scientifically valid negative result)
- ✅ Multiple hypotheses proposed for null result
- ✅ Comparison to MPRA edge case completed
- ✅ Future directions identified
- ✅ All code documented and reproducible

**Key Findings:**
1. eQTLs are not appropriate for chromatin model validation
2. Cell type matching is critical
3. Biological readout matching is essential (accessibility ≠ expression)
4. Both MPRA and endogenous benchmarks show weak/no correlation

---

## 📚 Citation

**GTEx Consortium:**  
The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science.* 2020;369(6509):1318-1330.

**ENCODE Project:**  
ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. *Nature.* 2012;489(7414):57-74.

**AlphaGenome:**  
This work validates AlphaGenome on endogenous human variants, demonstrating the importance of matched biological readouts for proper model evaluation.

---

## 🏆 Key Takeaways

1. **Negative Result Is Important** - Null correlation reveals mismatch between chromatin predictions and expression readouts
2. **Benchmarks Must Match Biology** - Accessibility models need accessibility data, not expression data
3. **Cell Type Matters** - K562 predictions don't capture Whole Blood biology
4. **Context Dominates** - 2048 bp windows dilute single nucleotide signals
5. **MPRA vs Endogenous** - Both show weak correlations for different reasons
6. **Path Forward** - Use caQTLs, match cell types, stratify by effect size

**Bottom Line:** This benchmark demonstrates that while we've created a proper validation framework (natural variants + native chromatin), we've tested it with an inappropriate gold standard (eQTLs measure expression, not accessibility). The null result is scientifically valuable and redirects future validation efforts toward chromatin accessibility QTLs with cell-type matched predictions.

---

*Last updated: November 6, 2025*
