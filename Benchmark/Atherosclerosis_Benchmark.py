#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from upsetplot import from_contents, plot as upset_plot

# ==========================================
# 1. SETUP: Define The Expanded Benchmark Key
# ==========================================

benchmark_key = {
    'Tier_1': [
        # Lipid Lowering (Statins, PCSK9i, Ezetimibe, Bempedoic Acid, etc.)
        'HMGCR', 'PCSK9', 'APOB', 'NPC1L1', 'ACLY', 'MTTP', 'LDLR', 'ANGPTL3', 'PPARA',
        # Anti-Inflammatory (Canakinumab, Colchicine targets)
        'IL1B', 'TUBA1A', 'TUBB', # Colchicine binds tubulin
        # Anti-Thrombotic (Aspirin, P2Y12 inhibitors - relevant to plaque rupture)
        'PTGS1', 'P2RY12', 'F2R'
    ],
    'Tier_2': [
        # Major GWAS Hits (CARDIoGRAMplusC4D, etc.)
        'SORT1', 'LPA', 'PHACTR1', 'CDKN2A', 'CDKN2B', 'ANRIL', # 9p21
        'SH2B3', 'NOS3', 'GUCY1A3', 'COL4A1', 'COL4A2',
        'LPL', 'LIPA', 'CELSR2', 'PSRC1', 'MIA3', 'SMAD3', 'TCF21', 'ADAMTS7',
        'EDNRA', 'HDAC9', 'KCNK5', 'ZC3HC1'
    ],
    'Tier_3': [
        # Inflammatory/Regulatory/Mechanistic (Literature Validated)
        'NLRP3', 'IL6', 'IL18', 'CCL2', 'VCAM1', 'ICAM1', 'SELE', 'MMP9', 'MMP2',
        'ABCA1', 'ABCG1', # Cholesterol Efflux
        'SREBF1', 'SREBF2', # SREBP Transcription factors
        'MIR155', 'MIR33A', 'MIR33B', 'MIR143', 'MIR145', 'MIR92A1', 'MIR126'
    ]
}

# Define Pathway Mapping for PCM Score
# Maps genes to lists of pathways: ['Lipid'], ['Inflammation'], ['Vascular']
# Adhesion molecules map to BOTH ['Inflammation', 'Vascular']
# Thrombosis targets map to ['Vascular']
gene_pathway_map = {
    # Lipid Metabolism
    'HMGCR': ['Lipid'], 'PCSK9': ['Lipid'], 'APOB': ['Lipid'], 'NPC1L1': ['Lipid'], 'ACLY': ['Lipid'],
    'MTTP': ['Lipid'], 'LDLR': ['Lipid'], 'ANGPTL3': ['Lipid'], 'PPARA': ['Lipid'],
    'SORT1': ['Lipid'], 'LPA': ['Lipid'], 'LPL': ['Lipid'], 'LIPA': ['Lipid'], 'CELSR2': ['Lipid'], 'PSRC1': ['Lipid'],
    'ABCA1': ['Lipid'], 'ABCG1': ['Lipid'], 'SREBF1': ['Lipid'], 'SREBF2': ['Lipid'], 'MIR33A': ['Lipid'], 'MIR33B': ['Lipid'],

    # Inflammation & Immune Response
    'IL1B': ['Inflammation'], 'NLRP3': ['Inflammation'], 'IL6': ['Inflammation'], 'IL18': ['Inflammation'],
    'SH2B3': ['Inflammation'], 'CCL2': ['Inflammation'], 'MMP9': ['Inflammation'], 'MMP2': ['Inflammation'],
    'MIR155': ['Inflammation'],

    # Vascular Wall Biology (EC/VSMC function, Matrix, cell cycle)
    'PHACTR1': ['Vascular'], 'CDKN2A': ['Vascular'], 'CDKN2B': ['Vascular'], 'ANRIL': ['Vascular'],
    'NOS3': ['Vascular'], 'GUCY1A3': ['Vascular'], 'COL4A1': ['Vascular'], 'COL4A2': ['Vascular'],
    'EDNRA': ['Vascular'], 'HDAC9': ['Vascular'], 'TCF21': ['Vascular'], 'ADAMTS7': ['Vascular'], 'MIA3': ['Vascular'],
    'TUBA1A': ['Vascular'], 'TUBB': ['Vascular'], # Cytoskeleton/Colchicine
    'MIR143': ['Vascular'], 'MIR145': ['Vascular'], 'MIR92A1': ['Vascular'], 'MIR126': ['Vascular'],

    # Dual Pathway: Adhesion Molecules (Inflammation + Vascular)
    'VCAM1': ['Inflammation', 'Vascular'],
    'ICAM1': ['Inflammation', 'Vascular'],
    'SELE': ['Inflammation', 'Vascular'],

    # Thrombosis & Coagulation (Merged into Vascular)
    'PTGS1': ['Vascular'], 'P2RY12': ['Vascular'], 'F2R': ['Vascular']
}

# ==========================================
# 2. LOAD METHOD RESULTS
# ==========================================

def load_method_results(filepath, method_name, gene_col='GeneSymbol'):
    try:
        df = pd.read_csv(filepath)
        # Handle variations in column naming
        if gene_col not in df.columns:
             cols = df.columns
             if 'Gene' in cols: gene_col = 'Gene'
             elif 'ID' in cols: gene_col = 'ID'
             elif 'GeneSymbol' in cols: gene_col = 'GeneSymbol'
             
        genes = df[gene_col].dropna().unique().tolist()
        # Clean up: remove empty strings, convert to uppercase for matching
        genes = [str(g).upper() for g in genes if str(g).strip() != '']
        print(f"Loaded {len(genes)} targets for {method_name}")
        return set(genes)
    except Exception as e:
        print(f"Error loading {method_name} from {filepath}: {e}")
        return set()

# Load all 6 methods
methods = {
    'MAGMA (DNA)': load_method_results('magma_significant_HGNC.csv', 'MAGMA', gene_col='GeneSymbol'),
    'TWAS (DNA/RNA)': load_method_results('TWAS_targets_mapped.csv', 'TWAS', gene_col='GeneSymbol'),
    'DESeq2 (RNA)': load_method_results('deseq2_HGNC_genes.csv', 'DESeq2', gene_col='Gene'),
    'WGCNA (RNA)': load_method_results('wgcna_significant_genes.csv', 'WGCNA', gene_col='Gene'),
    'Limma (Protein)': load_method_results('limma_suggestive_HGNC.csv', 'Limma', gene_col='GeneSymbol'),
    'DIABLO (Multi)': load_method_results('diablo_targets_mapped.csv', 'DIABLO', gene_col='GeneSymbol')
}

# ==========================================
# 3. DEFINE BASELINES FOR DTS
# ==========================================

# Define "Signal" sets based on the standard method for that datatype
baseline_dna = methods['MAGMA (DNA)']
baseline_rna = methods['DESeq2 (RNA)']
baseline_prot = methods['Limma (Protein)']

# Combined Baselines
dts_dna_rna = baseline_dna.union(baseline_rna)
dts_rna_prot = baseline_rna.union(baseline_prot) # Baseline for RNA+Protein integration
dts_all = baseline_dna.union(baseline_rna).union(baseline_prot)

# Map methods to their appropriate Datatype Specific DTS
method_dts_map = {
    'MAGMA (DNA)': ('DNA Baseline', baseline_dna),
    'TWAS (DNA/RNA)': ('DNA+RNA Baseline', dts_dna_rna),
    'DESeq2 (RNA)': ('RNA Baseline', baseline_rna),
    'WGCNA (RNA)': ('RNA Baseline', baseline_rna),
    'Limma (Protein)': ('Protein Baseline', baseline_prot),
    'DIABLO (Multi)': ('RNA+Protein Baseline', dts_rna_prot) 
}

# ==========================================
# 4. TEST 2: CALCULATE ABS SCORES
# ==========================================

results_table = []

for name, genes in methods.items():
    # Calculate TVS (Target Validation Score)
    t1_hits = [g for g in genes if g in benchmark_key['Tier_1']]
    t2_hits = [g for g in genes if g in benchmark_key['Tier_2']]
    t3_hits = [g for g in genes if g in benchmark_key['Tier_3']]
    
    tvs = (len(t1_hits) * 10) + (len(t2_hits) * 5) + (len(t3_hits) * 3)
    
    # Calculate PCM (Pathway Coverage Multiplier)
    all_hits = t1_hits + t2_hits + t3_hits
    pathways_hit = set()
    for g in all_hits:
        if g in gene_pathway_map:
            # Add all pathways listed for this gene
            for pathway in gene_pathway_map[g]:
                pathways_hit.add(pathway)
            
    pcm = len(pathways_hit)
    if pcm == 0: pcm = 1 
    
    abs_score = tvs * pcm
    
    # Store Stats
    results_table.append({
        'Method': name,
        'Total Targets': len(genes),
        'Tier 1 Hits': len(t1_hits),
        'Tier 2 Hits': len(t2_hits),
        'Tier 3 Hits': len(t3_hits),
        'Pathways Hit': pcm,
        'TVS': tvs,
        'ABS Score': abs_score
    })

scores_df = pd.DataFrame(results_table).sort_values('ABS Score', ascending=False)

# ==========================================
# 5. TEST 3: DATASET SIGNAL RECOVERY
# ==========================================

dts_results = []

for name, genes in methods.items():
    # 1. Datatype Specific DTS
    dts_name, dts_set = method_dts_map[name]
    
    # Intersection with Specific DTS
    recovered_specific = genes.intersection(dts_set)
    recall_specific = len(recovered_specific) / len(dts_set) if len(dts_set) > 0 else 0
    
    # 2. Overall DTS (Union of all baselines)
    recovered_overall = genes.intersection(dts_all)
    recall_overall = len(recovered_overall) / len(dts_all) if len(dts_all) > 0 else 0

    dts_results.append({
        'Method': name,
        'Specific Baseline Used': dts_name,
        'Specific DTS Size': len(dts_set),
        'Specific Recovery (%)': round(recall_specific * 100, 2),
        'Overall Multi-Omic Recovery (%)': round(recall_overall * 100, 2)
    })

dts_df = pd.DataFrame(dts_results)

# ==========================================
# 6. TEST 4: OVERLAP & CONSENSUS
# ==========================================

# Calculate Consensus Targets (Found by at least 3 methods)
all_genes_list = []
for genes in methods.values():
    all_genes_list.extend(list(genes))

from collections import Counter
gene_counts = Counter(all_genes_list)
consensus_genes_3plus = [g for g, count in gene_counts.items() if count >= 3]
consensus_genes_all = set.intersection(*methods.values())

print(f"\nConsensus Targets (In >= 3 methods): {len(consensus_genes_3plus)}")
print(f"Strict Consensus (In ALL methods): {len(consensus_genes_all)}")

# ==========================================
# 7. VISUALIZATION & REPORTING
# ==========================================

print("\n=== FINAL BENCHMARK REPORT ===")
print(scores_df.to_string(index=False))
print("\n=== DTS SIGNAL RECOVERY ===")
print(dts_df.to_string(index=False))

# Plot 1: ABS Scores
plt.figure(figsize=(12, 6))
sns.barplot(data=scores_df, x='Method', y='ABS Score', palette='viridis')
plt.title('Atherosclerosis Benchmark Score (ABS) Comparison')
plt.xticks(rotation=45)
plt.ylabel('ABS Score')
plt.tight_layout()
plt.savefig('benchmark_abs_score.png')
print("Saved ABS plot to benchmark_abs_score.png")

# Plot 2: Tier Breakdown
tier_df = scores_df[['Method', 'Tier 1 Hits', 'Tier 2 Hits', 'Tier 3 Hits']].set_index('Method')
tier_df.plot(kind='bar', stacked=True, figsize=(12, 6), color=['gold', 'silver', '#cd7f32'])
plt.title('Validated Targets Identified by Tier')
plt.xticks(rotation=45)
plt.ylabel('Count of Targets')
plt.tight_layout()
plt.savefig('benchmark_tier_breakdown.png')
print("Saved Tier plot to benchmark_tier_breakdown.png")

# Plot 3: Upset Plot
# Creates a dictionary formatted for upsetplot
try:
    upset_data = from_contents(methods)
    plt.figure(figsize=(14, 7))
    upset_plot(upset_data, subset_size='count', show_counts=True)
    plt.title('Intersection of Targets Across Methods')
    plt.savefig('benchmark_upset.png')
    print("Saved Upset plot to benchmark_upset.png")
except Exception as e:
    print(f"Could not generate UpSet plot: {e}")

# Save tables to CSV
scores_df.to_csv("final_benchmark_scores.csv", index=False)
dts_df.to_csv("dts_recovery_scores.csv", index=False)