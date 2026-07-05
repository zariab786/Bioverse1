
import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from collections import Counter
import re
import numpy as np

st.set_page_config(
    page_title="Bioverse 2.0",
    page_icon="🧬",
    layout="wide"
)

# Title
st.title("🧬 Bioverse 2.0")
st.subheader("DNA to Protein Explorer")

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

# Sidebar
with st.sidebar:
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
    
    CGPA: 3.73
    Biotech Graduate
    """)

# Main content
col_left, col_right = st.columns([1, 1], gap="large")

with col_left:
    st.markdown("### 📝 DNA Sequence Input")
    
    dna_input = st.text_area(
        "Enter DNA sequence",
        value=st.session_state.get('dna_input', 'ATGGCGTAA'),
        height=100,
        key="dna_input"
    )
    
    translate_clicked = st.button("🧬 Translate", use_container_width=True)
    
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
            
            st.markdown("### 🔬 Translation Results")
            
            # Display amino acids
            result = " → ".join([f"{codon}({aa})" for codon, aa in zip(codons, amino_acids)])
            st.code(result, language="text")
            
            st.session_state.amino_acids = amino_acids

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
        st.info("🧬 Enter a DNA sequence and click 'Translate' to see the 3D structure")

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
