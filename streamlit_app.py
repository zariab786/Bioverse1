
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from collections import Counter
import re
import numpy as np
from Bio import Entrez, SeqIO
import requests
import json

# Set NCBI email (required)
Entrez.email = "yyen90211@gmail.com"

st.set_page_config(
    page_title="Bioverse 2.0 - DNA to Protein Explorer",
    page_icon="🧬",
    layout="wide"
)

# Codon table
codon_table = {
    'ATA':'I','ATC':'I','ATT':'I','ATG':'M',
    'ACA':'T','ACC':'T','ACG':'T','ACT':'T',
    'AAC':'N','AAT':'N','AAA':'K','AAG':'K',
    'AGC':'S','AGT':'S','AGA':'R','AGG':'R',
    'CTA':'L','CTC':'L','CTG':'L','CTT':'L',
    'CCA':'P','CCC':'P','CCG':'P','CCT':'P',
    'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q',
    'CGA':'R','CGC':'R','CGG':'R','CGT':'R',
    'GTA':'V','GTC':'V','GTG':'V','GTT':'V',
    'GCA':'A','GCC':'A','GCG':'A','GCT':'A',
    'GAC':'D','GAT':'D','GAA':'E','GAG':'E',
    'GGA':'G','GGC':'G','GGG':'G','GGT':'G',
    'TCA':'S','TCC':'S','TCG':'S','TCT':'S',
    'TTC':'F','TTT':'F','TTA':'L','TTG':'L',
    'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*',
    'TGC':'C','TGT':'C','TGA':'*','TGG':'W'
}

aa_names = {
    'A':'Alanine','R':'Arginine','N':'Asparagine','D':'Aspartic Acid',
    'C':'Cysteine','E':'Glutamic Acid','Q':'Glutamine','G':'Glycine',
    'H':'Histidine','I':'Isoleucine','L':'Leucine','K':'Lysine',
    'M':'Methionine','F':'Phenylalanine','P':'Proline','S':'Serine',
    'T':'Threonine','W':'Tryptophan','Y':'Tyrosine','V':'Valine','*':'Stop'
}

aa_colors = {
    'A':'#4CAF50','R':'#2196F3','N':'#9C27B0','D':'#F44336',
    'C':'#FF9800','E':'#FF5722','Q':'#8BC34A','G':'#CDDC39',
    'H':'#3F51B5','I':'#009688','L':'#00BCD4','K':'#FF4081',
    'M':'#FFC107','F':'#E91E63','P':'#FF6F00','S':'#795548',
    'T':'#607D8B','W':'#9E9E9E','Y':'#FFD54F','V':'#26A69A','*':'#000000'
}

def translate_dna(dna):
    dna = dna.upper().strip()
    dna = re.sub(r"\s+", "", dna)
    amino_acids = []
    codons = []
    for i in range(0, len(dna) - 2, 3):
        codon = dna[i:i+3]
        codons.append(codon)
        amino_acids.append(codon_table.get(codon, '?'))
    return codons, amino_acids

def fetch_from_ncbi(accession):
    """Fetch DNA sequence from NCBI using accession number"""
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq), record.description
    except Exception as e:
        return None, str(e)

def get_gene_info(accession):
    """Fetch gene information from NCBI"""
    try:
        handle = Entrez.esummary(db="nucleotide", id=accession, retmode="json")
        data = handle.read()
        handle.close()
        return json.loads(data)
    except Exception as e:
        return None

def ai_analyze_protein(amino_acids):
    """Simple AI analysis using rules"""
    if not amino_acids:
        return {}
    
    # Count properties
    hydrophobic = ['A','V','L','I','F','W','M','P']
    hydrophilic = ['R','N','D','E','Q','K','H']
    charged = ['R','K','D','E','H']
    polar = ['S','T','Y','C','N','Q']
    
    counts = Counter(amino_acids)
    total = len(amino_acids)
    
    analysis = {
        'total_aa': total,
        'unique_aa': len(set(amino_acids)),
        'hydrophobic': sum(counts.get(aa, 0) for aa in hydrophobic),
        'hydrophilic': sum(counts.get(aa, 0) for aa in hydrophilic),
        'charged': sum(counts.get(aa, 0) for aa in charged),
        'polar': sum(counts.get(aa, 0) for aa in polar),
        'most_common': counts.most_common(3),
        'stop_codons': amino_acids.count('*')
    }
    
    # Generate insights
    insights = []
    if analysis['hydrophobic'] / total > 0.5:
        insights.append("🧬 This protein is hydrophobic - likely membrane-bound")
    if analysis['charged'] / total > 0.3:
        insights.append("⚡ High charge content - may interact with DNA/RNA")
    if analysis['polar'] / total > 0.4:
        insights.append("💧 Polar-rich protein - likely soluble in water")
    if analysis['unique_aa'] > 15:
        insights.append("🌈 High amino acid diversity - complex protein structure")
    if analysis['stop_codons'] > 0:
        insights.append("⏹️ Contains stop codons - may be incomplete sequence")
    
    analysis['insights'] = insights
    return analysis

# Custom CSS
st.markdown("""
<style>
    .stApp { background: linear-gradient(135deg, #0a0f0a 0%, #1a2a1a 50%, #0a1a0a 100%); }
    .main-header { text-align: center; padding: 2rem 0; }
    .main-header h1 { font-size: 4rem; font-weight: 800; background: linear-gradient(90deg, #4ade80, #22c55e, #16a34a); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
    .main-header p { color: #9ca3af; font-size: 1.2rem; }
    .stat-box { background: rgba(34, 197, 94, 0.1); border: 1px solid rgba(34, 197, 94, 0.2); border-radius: 10px; padding: 15px; text-align: center; }
    .stat-box .number { font-size: 2rem; font-weight: bold; color: #4ade80; }
    .stat-box .label { color: #9ca3af; font-size: 0.9rem; }
    .insight-box { background: rgba(34, 197, 94, 0.05); border-left: 3px solid #4ade80; padding: 10px 15px; margin: 5px 0; border-radius: 5px; }
</style>
""", unsafe_allow_html=True)

# Header
st.markdown("""
<div class="main-header">
    <h1>🧬 Bioverse 2.0</h1>
    <p>DNA to Protein Explorer with NCBI Integration & AI Analysis</p>
    <div style="display: flex; justify-content: center; gap: 1rem; margin-top: 1rem;">
        <span style="background: rgba(34,197,94,0.2); padding: 5px 15px; border-radius: 20px; color: #4ade80; font-size: 0.9rem;">🧬 CGPA: 3.73</span>
        <span style="background: rgba(168,85,247,0.2); padding: 5px 15px; border-radius: 20px; color: #a855f7; font-size: 0.9rem;">⚡ Biotech Graduate</span>
        <span style="background: rgba(59,130,246,0.2); padding: 5px 15px; border-radius: 20px; color: #3b82f6; font-size: 0.9rem;">🧪 NCBI + AI</span>
    </div>
</div>
""", unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.markdown("### 🔍 NCBI Accession Search")
    accession_input = st.text_input(
        "Enter Accession Number",
        placeholder="e.g., NM_000518, NP_000509, XM_005249",
        help="NCBI accession numbers for DNA sequences"
    )
    
    if st.button("🔍 Fetch from NCBI", use_container_width=True):
        if accession_input:
            with st.spinner("Fetching from NCBI..."):
                sequence, description = fetch_from_ncbi(accession_input.strip())
                if sequence:
                    st.session_state.dna_input = sequence
                    st.session_state.ncbi_description = description
                    st.session_state.accession = accession_input
                    st.success(f"✅ Fetched: {accession_input}")
                    st.rerun()
                else:
                    st.error(f"❌ Error: {description}")
    
    st.markdown("---")
    st.markdown("### 🧪 Quick Presets")
    examples = ['ATGGCGTAA', 'ATGGCCATTGTAATGGGCCGCTAA', 'ATGAAGTTTGGCACTTAA']
    for seq in examples:
        if st.button(seq, key=f"preset_{seq}", use_container_width=True):
            st.session_state.dna_input = seq
            st.rerun()
    
    st.markdown("---")
    st.markdown("### 👨‍🔬 About")
    st.info("""
    **Bioverse 2.0**
    DNA to Protein Explorer
    
    - 🧬 DNA Translation
    - 🔍 NCBI Integration
    - 🤖 AI Analysis
    - 🌈 3D Visualization
    
    CGPA: 3.73
    Biotech Graduate
    """)

# Main content
col_left, col_right = st.columns([1, 1], gap="large")

with col_left:
    st.markdown("### 📝 DNA Sequence Input")
    
    # Show NCBI info if available
    if 'accession' in st.session_state:
        st.success(f"📌 Accession: {st.session_state.accession}")
        if 'ncbi_description' in st.session_state:
            st.caption(st.session_state.ncbi_description[:100] + "...")
    
    dna_input = st.text_area(
        "Enter DNA sequence",
        value=st.session_state.get('dna_input', 'ATGGCGTAA'),
        height=120,
        key="dna_input"
    )
    
    col_btn1, col_btn2, col_btn3 = st.columns([2, 1, 1])
    with col_btn1:
        translate_clicked = st.button("🧬 Translate & Analyze", use_container_width=True)
    with col_btn2:
        if st.button("Clear", use_container_width=True):
            st.session_state.dna_input = ""
            if 'accession' in st.session_state:
                del st.session_state.accession
                del st.session_state.ncbi_description
            st.rerun()
    with col_btn3:
        if st.button("Example", use_container_width=True):
            st.session_state.dna_input = "ATGGCGTAA"
            st.rerun()
    
    if dna_input:
        valid = bool(re.match(r'^[ATCGatcg\s]+$', dna_input))
        if not valid:
            st.error("⚠️ Invalid sequence! Only A, T, C, G allowed.")
        elif len(dna_input.replace(' ', '')) < 3:
            st.warning("⚠️ Sequence must be at least 3 bases long.")
        else:
            st.success("✅ Valid DNA sequence")
    
    if translate_clicked and dna_input:
        clean_dna = re.sub(r"\s+", "", dna_input.upper())
        if len(clean_dna) >= 3 and re.match(r'^[ATCG]+$', clean_dna):
            codons, amino_acids = translate_dna(clean_dna)
            
            # AI Analysis
            analysis = ai_analyze_protein(amino_acids)
            
            st.markdown("### 🔬 Translation Results")
            
            # Display amino acids
            result = " → ".join([f"{codon}({aa})" for codon, aa in zip(codons, amino_acids)])
            st.code(result, language="text")
            
            st.session_state.amino_acids = amino_acids
            st.session_state.analysis = analysis

with col_right:
    st.markdown("### 🌈 3D Protein Structure")
    
    if 'amino_acids' in st.session_state and st.session_state.amino_acids:
        amino_acids = st.session_state.amino_acids
        
        n = len(amino_acids)
        if n > 0:
            t = np.linspace(0, 4*np.pi, n)
            x = np.cos(t) * 2.5 + np.random.normal(0, 0.1, n)
            y = np.sin(t) * 2.5 + np.random.normal(0, 0.1, n)
            z = np.linspace(-3, 3, n)
            
            colors = [aa_colors.get(aa, '#888888') for aa in amino_acids]
            
            fig = go.Figure()
            
            fig.add_trace(go.Scatter3d(
                x=x, y=y, z=z,
                mode='markers',
                marker=dict(
                    size=15,
                    color=colors,
                    symbol='circle',
                    line=dict(width=1, color='white')
                ),
                text=[f"{aa}: {aa_names.get(aa, 'Unknown')}" for aa in amino_acids],
                hoverinfo='text',
                name='Amino Acids'
            ))
            
            fig.add_trace(go.Scatter3d(
                x=x, y=y, z=z,
                mode='lines',
                line=dict(color='rgba(74, 222, 128, 0.5)', width=3),
                name='Peptide Bond',
                hoverinfo='skip'
            ))
            
            fig.update_layout(
                template='plotly_dark',
                scene=dict(
                    xaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
                    yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
                    zaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
                    camera=dict(eye=dict(x=1.5, y=1.5, z=1.5)),
                    bgcolor='rgba(0,0,0,0)'
                ),
                showlegend=False,
                margin=dict(l=0, r=0, t=0, b=0),
                height=400,
                paper_bgcolor='rgba(0,0,0,0)',
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Statistics
            st.markdown("### 📊 Statistics")
            col1, col2, col3 = st.columns(3)
            col1.metric("Amino Acids", len(amino_acids))
            col2.metric("Unique AA", len(set(amino_acids)))
            col3.metric("Stop Codons", amino_acids.count('*'))
            
            # AI Insights
            if 'analysis' in st.session_state:
                analysis = st.session_state.analysis
                st.markdown("### 🤖 AI Insights")
                for insight in analysis.get('insights', []):
                    st.markdown(f'<div class="insight-box">{insight}</div>', unsafe_allow_html=True)
            
            # Composition
            counts = Counter(amino_acids)
            df = pd.DataFrame({
                'Amino Acid': list(counts.keys()),
                'Count': list(counts.values())
            })
            df = df.sort_values('Count', ascending=False)
            
            fig2 = px.bar(df, x='Amino Acid', y='Count', title='Amino Acid Distribution')
            fig2.update_layout(
                template='plotly_dark',
                paper_bgcolor='rgba(0,0,0,0)',
                plot_bgcolor='rgba(0,0,0,0)',
                height=200,
            )
            st.plotly_chart(fig2, use_container_width=True)
    else:
        st.info("🧬 Enter a DNA sequence and click 'Translate & Analyze' to see the 3D structure")

# Codon Table
st.markdown("---")
st.markdown("### 📚 Codon Table")

cols = st.columns(4)
for idx, (codon, aa) in enumerate(codon_table.items()):
    col_idx = idx % 4
    with cols[col_idx]:
        color = aa_colors.get(aa, '#888888')
        st.markdown(f"""
        <div style="display: flex; justify-content: space-between; background: rgba(0,0,0,0.2); padding: 5px 10px; border-radius: 5px; margin: 2px 0;">
            <span style="font-family: monospace; color: #4ade80;">{codon}</span>
            <span style="background: {color}; padding: 0 8px; border-radius: 4px; color: white; font-weight: bold;">{aa}</span>
            <span style="font-size: 10px; color: #9ca3af;">{aa_names.get(aa, '')[:6]}</span>
        </div>
        """, unsafe_allow_html=True)

st.markdown("---")
st.markdown("Built with ❤️ by Zariab | Biotechnologist & Developer")
