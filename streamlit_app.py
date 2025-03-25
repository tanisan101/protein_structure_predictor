# Importing necessary libraries
import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
import pandas as pd  # Added for DataFrame

# Function to parse FASTA file
def parse_fasta(uploaded_file):
    """Parses a FASTA file and returns a list of sequences."""
    sequences = []
    sequence = ""
    for line in uploaded_file:
        line = line.decode("utf-8").strip()  # Decode bytes to string
        if line.startswith(">"):  # Identifier line
            if sequence:  # If there's a sequence, add it to the list
                sequences.append(sequence)
            sequence = ""  # Reset for the next sequence
        else:
            sequence += line  # Append sequence data
    if sequence:  # Add the last sequence
        sequences.append(sequence)
    return sequences


# Streamlit app configuration
st.title("🎈 Predict Protein Structure")
st.sidebar.title("ESMFold")
st.sidebar.write(
    "[*ESMFold*](https://esmatlas.com/about) is an end-to-end single sequence protein structure predictor "
    "based on the ESM-2 language model. For more information, read the [research article]"
    "(https://www.biorxiv.org/content/10.1101/2022.07.20.500902v2) and the [news article]"
    "(https://www.nature.com/articles/d41586-022-03539-1) published in *Nature*."
)
# st.set_page_config(layout = 'wide')


# Function to render molecular structure
def render_mol(pdb_string):  # Expects pdb_string, not just pdb
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb_string, 'pdb')
    pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height=500, width=800)


# Protein sequence input and FASTA upload
Default_seq = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
txt = st.sidebar.text_area('Input Sequence', Default_seq, height=270)

# FASTA file upload
uploaded_file = st.sidebar.file_uploader("Upload FASTA File", type=["fasta"])

# ESMfold prediction function
def update(sequence):  # Expects sequence as argument
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }
    try:  # Added try-except block for error handling
        response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
        response.raise_for_status()  # Raise HTTPError for bad responses (4xx or 5xx)
        pdb_string = response.content.decode('utf-8')

        # Added temporary file handling using NamedTemporaryFile
        # Write the pdb_string to a temporary file
        with open('predicted.pdb', 'w') as f:
            f.write(pdb_string)  # Write string to file

        # Load the structure from the temporary file
        struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
        b_value = round(struct.b_factor.mean(), 4)

        # Output and Visualization
        st.subheader("Kafi interesting sequnece hai na!! ><")
        render_mol(pdb_string)  # Pass string to render_mol

        st.subheader('plDDT')
        st.write('plDDT is a per-residue estimate of the confidence in prediction on a scale from 0-100.')
        st.info(f'plDDT: {b_value}')

        st.download_button(
            label="Download PDB",
            data=pdb_string,
            file_name='predicted.pdb',
            mime='text/plain',
        )
    except requests.exceptions.RequestException as e:  # Catch network errors
        st.error(f"Request failed: {e}")
    except Exception as e:  # Catch other errors
        st.error(f"An error occurred: {e}")


# Predict button logic
predict = st.sidebar.button('Predict')

# Check if FASTA file is uploaded
if uploaded_file is not None:
    sequences = parse_fasta(uploaded_file)
    if sequences:
        st.success(f"Found {len(sequences)} sequences in the FASTA file.")
        for i, seq in enumerate(sequences):
            st.write(f"Processing sequence {i + 1}/{len(sequences)}")
            with st.spinner(f"Predicting structure for sequence {i + 1}..."):
                update(seq)  # Pass individual sequence to update()
    else:
        st.warning("The uploaded FASTA file does not contain valid sequences.")
elif predict:  # If no FASTA file, proceed with the sequence from the text area
    with st.spinner("Predicting structure..."):
        update(txt)  # Pass text area sequence to update()
else:
    st.warning('👈 Enter protein sequence data or upload a FASTA file!')
