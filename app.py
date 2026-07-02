import sys
print("MY_PYTHON_PATH:", sys.executable)

import sys
import subprocess

# مصفوفة تحتوي على المكتبات الحرج وجودها لتشغيل التطبيق
required_libraries = ['seaborn', 'biopython', 'pandas', 'matplotlib']

for lib in required_libraries:
    try:
        __import__(lib)
    except ImportError:
        # إذا لم يجد المكتبة، سيقوم بتثبيتها ديناميكياً داخل البيئة النشطة حالياً
        subprocess.check_call([sys.executable, "-m", "pip", "install", lib])

import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio import Align
import io

# Set page configuration
st.set_page_config(
    page_title="Microbial Gene Similarity Analyzer",
    page_icon="🧬",
    layout="wide"
)

# App Title & Description
st.title("🧬 Microbial Gene Similarity Analyzer")
st.markdown("""
Compare gene sequences across multiple microorganism strains using global sequence alignment. 
Upload your own GenBank and FASTA files, or explore the built-in sample data below to see it in action.
""")

# Initialize the pair-wise aligner
aligner = Align.PairwiseAligner()
aligner.mode = 'global'

# Sidebar for settings WITH DEFAULT VALUES
st.sidebar.header("Alignment Parameters")
match = st.sidebar.number_input("Match Score", value=1.0, step=0.5)
mismatch = st.sidebar.number_input("Mismatch Penalty", value=-1.0, step=0.5)
gap_open = st.sidebar.number_input("Gap Open Penalty", value=-2.0, step=0.5)
gap_extend = st.sidebar.number_input("Gap Extend Penalty", value=-0.5, step=0.5)

# Update aligner settings based on user input
aligner.match_score = match
aligner.mismatch_score = mismatch
aligner.open_gap_score = gap_open
aligner.extend_gap_score = gap_extend

# Sample FASTA data for demonstration
SAMPLE_FASTA = """>Strain_A_WildType
ATGCGATCGATCGATCGATCGATCGATC
>Strain_B_Mutant1
ATGCGATCGATCGATCGATTGATCGATC
>Strain_C_Mutant2
ATGCGATCGATC---CGATCGATCGATC
>Strain_D_Divergent
ATGCGACCCCCCGATCGATCGATCGATC
"""

st.header("📂 Upload Sequence Files")

col_files1, col_files2 = st.columns(2)

with col_files1:
    st.subheader("1. Reference Strain (GenBank)")
    gb_file = st.file_uploader("Upload GenBank File (.gb, .gbk)", type=["gb", "gbk"])

with col_files2:
    st.subheader("2. Query Strains (FASTA)")
    fasta_files = st.file_uploader("Upload one or multiple FASTA files (.fasta, .fa)", type=["fasta", "fa"], accept_multiple_files=True)

# تفعيل منطق البيانات الافتراضية إذا لم يتم رفع ملفات بعد
if gb_file is not None and fasta_files:
    # ----------------------------------------------------
    # حالة قيام المستخدم برفع ملفاته الخاصة
    # ----------------------------------------------------
    is_sample = False
    try:
        gb_string = gb_file.getvalue().decode("utf-8")
        gb_io = io.StringIO(gb_string)
        gb_records = list(SeqIO.parse(gb_io, "genbank"))
        
        if not gb_records:
            st.error("❌ The uploaded reference file is NOT a valid GenBank format. Please ensure it has a .gb or .gbk extension and follows standard NCBI formatting.")
            st.stop()
            
        gb_record = gb_records[0]
        ref_id = gb_record.id
        ref_seq = str(gb_record.seq).upper()
        
    except Exception as e:
        st.error(f"Error reading GenBank file: {e}")
        st.stop()
        
    strain_ids = [ref_id]
    sequences = [ref_seq]
    summary_data = [{"Strain ID": f"{ref_id} (Reference - GenBank)", "Length (bp)": len(ref_seq)}]
    
    for f_file in fasta_files:
        try:
            fasta_string = f_file.getvalue().decode("utf-8")
            fasta_io = io.StringIO(fasta_string)
            for record in SeqIO.parse(fasta_io, "fasta"):
                strain_ids.append(record.id)
                sequences.append(str(record.seq).upper())
                summary_data.append({"Strain ID": record.id, "Length (bp)": len(record.seq)})
        except Exception as e:
            st.error(f"Error reading FASTA file {f_file.name}: {e}")
else:
    # ----------------------------------------------------
    # الحالة الافتراضية: عرض البيانات التجريبية فور فتح الموقع
    # ----------------------------------------------------
    is_sample = True
    st.info("💡 Showing sample data. Upload your own GenBank and FASTA files above to analyze your custom strains.")
    
    # قراءة البيانات التجريبية الثابتة وتحليلها
    fasta_io = io.StringIO(SAMPLE_FASTA)
    records = list(SeqIO.parse(fasta_io, "fasta"))
    
    strain_ids = [record.id for record in records]
    sequences = [str(record.seq).upper() for record in records]
    summary_data = [{"Strain ID": r.id, "Length (bp)": len(r.seq)} for r in records]

# طالما توفرت بيانات (سواء افتراضية أو مرفوعة) سيتم الحساب والعرض تلقائياً
num_strains = len(strain_ids)

st.subheader("📋 Loaded Strains Summary")
st.dataframe(pd.DataFrame(summary_data), use_container_width=True)

# حساب مصفوفة التشابه الجيني
similarity_matrix = np.zeros((num_strains, num_strains))
progress_bar = st.progress(0)

total_ops = (num_strains * (num_strains + 1)) // 2
current_op = 0

for i in range(num_strains):
    for j in range(i, num_strains):
        if i == j:
            similarity_matrix[i][j] = 100.0
        else:
            score = aligner.score(sequences[i], sequences[j])
            max_score_i = aligner.score(sequences[i], sequences[i])
            max_score_j = aligner.score(sequences[j], sequences[j])
            max_possible_score = max(max_score_i, max_score_j)
            
            if max_possible_score > 0:
                percentage = (score / max_possible_score) * 100
                percentage = max(0.0, min(100.0, percentage))
            else:
                percentage = 0.0

            similarity_matrix[i][j] = percentage
            similarity_matrix[j][i] = percentage
        
        current_op += 1
        progress_bar.progress(current_op / total_ops)

progress_bar.empty()

if is_sample:
    st.success("✅ Showing alignment matrix for built-in sample data!")
else:
    st.success("✅ Sequence Alignment Matrix Calculated Successfully!")

# تحويل المصفوفة إلى DataFrame وعرضها
df_similarity = pd.DataFrame(similarity_matrix, columns=strain_ids, index=strain_ids)

col1, col2 = st.columns([1, 1])

with col1:
    st.subheader("📊 Similarity Percentage Matrix")
    st.dataframe(df_similarity.style.format("{:.2f}%").background_gradient(cmap="viridis"), use_container_width=True)

with col2:
    st.subheader("🎨 Similarity Heatmap")
    fig, ax = plt.subplots(figsize=(6, 4.5))
    sns.heatmap(
        df_similarity, 
        annot=True, 
        fmt=".1f", 
        cmap="YlGnBu", 
        vmin=0, 
        vmax=100, 
        cbar_kws={'label': 'Similarity (%)'},
        ax=ax
    )
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()
    st.pyplot(fig)

# --- قسم تفسير النتائج الثابت في أسفل الصفحة ---
st.divider()
st.header("📘 Interpretation and Explanation of Results")

st.markdown("""
### 1. How to Read the Outputs
* **Diagonal Values (100.0%):** The diagonal elements (running from top-left to bottom-right) represent each strain aligned against itself (e.g., `Reference` vs `Reference`). This is your 100% identical control line.
* **Matrix Symmetry:** The matrix is perfectly mirrored across the diagonal. The similarity score of *Strain A vs Strain B* is identical to *Strain B vs Strain A*.
* **Color Depth:** In the heatmap, darker/warmer hues indicate high sequence preservation. Sudden shifts to lighter shades point to divergence, significant mutations, or evolutionary insertions/deletions (indels).

### 2. Biological Insights to Look For
* **High Similarity (> 95%):** Suggests these microorganism strains are closely related. They likely belong to the same species and share a very recent common ancestor. Observed differences might just be Single Nucleotide Polymorphisms (SNPs).
* **Moderate Similarity (70% - 95%):** Points to potential divergence. This level of variation often occurs when comparing different serotypes, variants of concern, or strains subjected to different environmental selective pressures.
* **Low Similarity (< 70%):** Indicates significant genomic rearrangements, large-scale gene deletions, or that the sequences may originate from entirely different species or genera altogether.

### 3. Understanding the Algorithm & Sidebar Parameters
The app uses Needleman-Wunsch-style global pairwise alignment (`Bio.Align.PairwiseAligner`). It matches the sequences from end-to-end based on the parameters set in the sidebar:
* **Matches (+1.0 by default):** Reward the score, pulling the final identity percentage closer to 100%.
* **Mismatches (-1.0 by default):** Penalize the score, representing **Point Mutations / Substitutions** where nucleotide letters vary between strains.
* **Gaps (-2.0 for Open / -0.5 for Extend):** Represent **Insertion or Deletion (Indel)** events during evolutionary drift.

The final percentage normalizes the maximum possible raw alignment score against the length and composition of the sequences being compared.
""")