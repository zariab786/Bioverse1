
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
import py3Dmol
from stmol import showmol
import time

# Set NCBI email (REQUIRED)
Entrez.email = "yyen90211@gmail.com"

# Try to get NCBI API key from secrets (optional but recommended)
try:
    NCBI_API_KEY = st.secrets["NCBI_API_KEY"]
except:
    NCBI_API_KEY = None

st.set_page_config(
    page_title="Bioverse 2.0 - REAL Protein Explorer",
    page_icon="🧬",
    layout="wide"
)

# REAL NCBI Fetch Function
def fetch_from_ncbi_real(accession):
    """Actually fetch from NCBI with proper error handling"""
    try:
        # Use API key if available
        if NCBI_API_KEY:
            Entrez.api_key = NCBI_API_KEY
        
        # First try to fetch as nucleotide
        try:
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return str(record.seq), record.description, "nucleotide"
        except:
            pass
        
        # Then try as protein
        try:
            handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
            record = SeqIO.read(handle, "fasta")
            handle.close()
            return str(record.seq), record.description, "protein"
        except:
            pass
        
        # Try UniProt API
        try:
            url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
            response = requests.get(url, timeout=10)
            if response.status_code == 200:
                lines = response.text.strip().split('\n')
                seq = ''.join(lines[1:])
                desc = lines[0] if lines else "Protein from UniProt"
                return seq, desc, "protein"
        except:
            pass
        
        return None, "Accession not found in NCBI or UniProt", None
        
    except Exception as e:
        return None, f"Error: {str(e)}", None

# Pre-loaded example sequences (for demo when NCBI fails)
EXAMPLE_SEQUENCES = {
    "P42212": {
        'sequence': "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK",
        'description': "Green Fluorescent Protein (GFP) - P42212",
        'type': 'protein'
    },
    "NM_000518": {
        'sequence': "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATGTGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGCACTGACTCTCTCTGCCTATTGGTCTATTTTCCCACCCTTAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTTTCTTGCTGTCCAATTTCTATTAAAGGTTCCTTTGTTCCCTAAGTCCAACTACTAAACTGGGGGATATTATGAAGGGCCTTGAGCATCTGGATTCTGCCTAATAAAAAACATTTATTTTCATTGC",
        'description': "Human beta-globin gene - NM_000518",
        'type': 'nucleotide'
    },
    "P01308": {
        'sequence': "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
        'description': "Insulin - P01308",
        'type': 'protein'
    }
}

def get_pdb_structure(pdb_id):
    """Fetch real PDB structure from RCSB"""
    try:
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.text
        return None
    except:
        return None

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

# Amino acid properties (same as before)
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

# Header
st.markdown("""
<style>
    .stApp { background: linear-gradient(135deg, #0a0f0a 0%, #1a2a1a 50%, #0a1a0a 100%); }
    .main-header { text-align: center; padding: 2rem 0; }
    .main-header h1 { font-size: 3.5rem; font-weight: 800; background: linear-gradient(90deg, #4ade80, #22c55e, #16a34a); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }
    .main-header p { color: #9ca3af; font-size: 1.1rem; }
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
    <p>REAL DNA to Protein Explorer with NCBI + PDB Integration</p>
    <div style="display: flex; justify-content: center; gap: 1rem; margin-top: 1rem;">
        <span style="background: rgba(34,197,94,0.2); padding: 5px 15px; border-radius: 20px; color: #4ade80; font-size: 0.9rem;">🧬 CGPA: 3.73</span>
        <span style="background: rgba(168,85,247,0.2); padding: 5px 15px; border-radius: 20px; color: #a855f7; font-size: 0.9rem;">⚡ Biotech Graduate</span>
        <span style="background: rgba(59,130,246,0.2); padding: 5px 15px; border-radius: 20px; color: #3b82f6; font-size: 0.9rem;">🔬 REAL NCBI</span>
        <span style="background: rgba(239,68,68,0.2); padding: 5px 15px; border-radius: 20px; color: #ef4444; font-size: 0.9rem;">🎨 PDB 3D</span>
    </div>
</div>
""", unsafe_allow_html=True)

# Sidebar
with st.sidebar:
    st.markdown("### 🔍 Search NCBI")
    st.caption("Try: P42212, NM_000518, P01308")
    
    accession_input = st.text_input(
        "Accession Number",
        placeholder="e.g., P42212"
    )
    
    col1, col2 = st.columns(2)
    with col1:
        if st.button("🔍 Fetch REAL", use_container_width=True):
            if accession_input:
                with st.spinner(f"Fetching {accession_input} from NCBI..."):
                    sequence, description, seq_type = fetch_from_ncbi_real(accession_input.strip())
                    if sequence:
                        st.session_state.dna_input = sequence
                        st.session_state.ncbi_description = description
                        st.session_state.accession = accession_input
                        st.session_state.seq_type = seq_type
                        st.session_state.is_real = True
                        st.success(f"✅ Fetched: {accession_input}")
                        time.sleep(0.5)
                        st.rerun()
                    else:
                        st.error(f"❌ {description}")
                        st.info("💡 Trying pre-loaded example...")
                        if accession_input in EXAMPLE_SEQUENCES:
                            seq_data = EXAMPLE_SEQUENCES[accession_input]
                            st.session_state.dna_input = seq_data['sequence']
                            st.session_state.ncbi_description = seq_data['description']
                            st.session_state.accession = accession_input
                            st.session_state.seq_type = seq_data['type']
                            st.session_state.is_real = False
                            st.success(f"✅ Loaded example: {accession_input}")
                            st.rerun()
    
    st.markdown("---")
    st.markdown("### 📌 Pre-loaded Examples")
    st.caption("Click to load instantly:")
    for acc in ['P42212', 'NM_000518', 'P01308']:
        if st.button(f"🧬 {acc}", key=f"example_{acc}", use_container_width=True):
            seq_data = EXAMPLE_SEQUENCES[acc]
            st.session_state.dna_input = seq_data['sequence']
            st.session_state.ncbi_description = seq_data['description']
            st.session_state.accession = acc
            st.session_state.seq_type = seq_data['type']
            st.session_state.is_real = False
            st.rerun()
    
    st.markdown("---")
    st.markdown("### 📊 Status")
    if 'accession' in st.session_state:
        st.success(f"Loaded: {st.session_state.accession}")
        if st.session_state.get('is_real', False):
            st.info("✅ From REAL NCBI")
        else:
            st.warning("📦 From pre-loaded data")

# Main content
st.markdown("### 🧬 Sequence Viewer")

col_left, col_right = st.columns([1, 2])

with col_left:
    dna_input = st.text_area(
        "DNA or Protein Sequence",
        value=st.session_state.get('dna_input', ''),
        height=200,
        key="dna_input"
    )
    
    if st.button("🧬 Translate & Analyze", use_container_width=True):
        if dna_input:
            clean_input = re.sub(r"\s+", "", dna_input.upper())
            is_dna = bool(re.match(r'^[ATCG]+$', clean_input))
            is_protein = bool(re.match(r'^[ARNDCEQGHILKMFPSTWYV]+$', clean_input))
            
            if is_dna and len(clean_input) >= 3:
                codons, amino_acids = translate_dna(clean_input)
                st.session_state.amino_acids = amino_acids
                st.session_state.seq_type = "DNA"
            elif is_protein:
                st.session_state.amino_acids = list(clean_input)
                st.session_state.seq_type = "Protein"
            else:
                st.error("Invalid input!")
    
    if 'amino_acids' in st.session_state:
        aa = st.session_state.amino_acids
        st.markdown("### 📊 Stats")
        col1, col2, col3 = st.columns(3)
        col1.metric("Length", len(aa))
        col2.metric("Unique", len(set(aa)))
        col3.metric("Stop", aa.count('*'))

with col_right:
    st.markdown("### 🎨 3D Structure Viewer")
    
    if 'amino_acids' in st.session_state and st.session_state.amino_acids:
        amino_acids = st.session_state.amino_acids
        
        # Show 3D viewer
        n = len(amino_acids)
        if n > 0:
            t = np.linspace(0, 4*np.pi, n)
            x = np.cos(t) * 2.5 + np.random.normal(0, 0.1, n)
            y = np.sin(t) * 2.5 + np.random.normal(0, 0.1, n)
            z = np.linspace(-4, 4, n)
            
            colors = [aa_colors.get(aa, '#888888') for aa in amino_acids]
            
            fig = go.Figure()
            
            fig.add_trace(go.Scatter3d(
                x=x, y=y, z=z,
                mode='markers+lines',
                marker=dict(
                    size=20,
                    color=colors,
                    symbol='circle',
                    line=dict(width=2, color='white')
                ),
                line=dict(color='rgba(74, 222, 128, 0.6)', width=4),
                text=[f"{aa}: {aa_names.get(aa, 'Unknown')}" for aa in amino_acids],
                hoverinfo='text'
            ))
            
            fig.update_layout(
                template='plotly_dark',
                scene=dict(
                    xaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
                    yaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
                    zaxis=dict(showticklabels=False, showgrid=False, zeroline=False),
                    camera=dict(
                        eye=dict(x=1.8, y=1.8, z=1.8),
                        center=dict(x=0, y=0, z=0),
                        up=dict(x=0, y=0, z=1)
                    ),
                    bgcolor='rgba(0,0,0,0)'
                ),
                showlegend=False,
                margin=dict(l=0, r=0, t=0, b=0),
                height=500,
                paper_bgcolor='rgba(0,0,0,0)',
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            # Show description if available
            if 'ncbi_description' in st.session_state:
                st.info(f"📌 {st.session_state.ncbi_description}")
            
            st.caption("🔄 Drag to rotate | Scroll to zoom | Color-coded by amino acid type")
    
    else:
        st.info("🧬 Enter a sequence or search an accession number")

st.markdown("---")
st.markdown("🔬 **Data Sources:** NCBI (nucleotide/protein) | UniProt | RCSB PDB")
