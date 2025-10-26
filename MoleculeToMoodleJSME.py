# -*- coding: utf-8 -*-
"""
Created on Tue Sep 30 19:14:39 2025

@author: Carlos Fernandez
"""

import streamlit as st
import xml.etree.ElementTree as ET
import pandas as pd
import requests

# Configure the page to use the wide layout
st.set_page_config(layout="wide")


# ==============================================================================
# 1. CONSTANTS AND DEPENDENCY CHECKS
# ==============================================================================

# --- Dependency Check: requests (for NCI CIR API) ---
try:
    NCI_CIR_AVAILABLE = True 
except NameError:
    st.warning("Warning: The 'requests' module is not installed. Automatic SMILES lookup will not be available.")
    NCI_CIR_AVAILABLE = False
    requests = None


# --- Dependency Check: RDKit ---
try:
    from rdkit import Chem
except ImportError:
    st.warning("Warning: The 'rdkit' module is not installed. "
               "It will not be possible to validate molecule structures or canonicalize SMILES.")
    Chem = None


# Interface texts in different languages
TEXTS = {
    "es": {
        "title": "Generador de preguntas Moodle con SMILES (NCI CIR Compatible) 🧪",
        "section_individual": "1. Entrada Individual",
        "section_bulk": "2. Carga Masiva (CSV / Excel)",
        "expand_bulk": "Abrir Conversión Masiva",
        "section_export": "3. Exportar",
        "molecule_name": "Nombre de la molécula",
        "name_note": "La búsqueda se realiza a través de la API del Chemical Identifier Resolver (CIR) de NCI. El resultado se canonicaliza con RDKit para asegurar la máxima compatibilidad con el editor JSME de Moodle.",
        "search_smiles": "Buscar SMILES automáticamente (NCI CIR API + RDKit Canonicalización)",
        "smiles_input": "SMILES (si no usas la API o si falló la búsqueda)",
        "add_question": "Agregar pregunta",
        "upload_file": "Seleccionar archivo CSV/Excel (.csv, .xlsx, .xls)",
        "bulk_note": "Asegúrate de que el archivo contenga una columna llamada **'nombre'** con los nombres de las moléculas. Se recomienda usar nombres en inglés para la API de CIR.",
        "file_format_error": "Formato de archivo no soportado. Por favor, sube un archivo CSV, XLSX o XLS.",
        "column_error": "El archivo debe contener la columna requerida: {}. Si estás usando un archivo en inglés, puede ser 'name'.",
        "bulk_success": "Conversión masiva completada: **{}** éxitos, **{}** fallos.",
        "smiles_not_found": "Molécula no encontrada por NCI CIR o la estructura no se pudo generar/canonicalizar.", 
        "invalid_smiles": "SMILES inválido o no encontrado.",
        "question_added": "Pregunta añadida: {}",
        "questions_added_subtitle": "Preguntas añadidas",
        "clear_all": "Borrar todas las preguntas",
        "download_xml": "Descargar XML Moodle",
        "xml_error": "Error al generar el archivo XML: {}",
        "question_title": "Estructura de {}",
        "question_text": "Dibuja la estructura de <strong>{}</strong> usando el editor molecular.",
        "delete_tooltip": "Borrar"
    },
    "en": {
        "title": "Moodle Question Generator with SMILES (NCI CIR Compatible) 🧪",
        "section_individual": "1. Individual Entry",
        "section_bulk": "2. Bulk Upload (CSV / Excel)",
        "expand_bulk": "Open Bulk Conversion",
        "section_export": "3. Export",
        "molecule_name": "Molecule Name",
        "name_note": "The search is performed using the NCI Chemical Identifier Resolver (CIR) API. The result is canonicalized with RDKit to ensure maximum compatibility with Moodle's JSME editor.",
        "search_smiles": "Search for SMILES automatically (NCI CIR API + RDKit Canonicalization)",
        "smiles_input": "SMILES (if not using the API or if the search failed)",
        "add_question": "Add question",
        "upload_file": "Select CSV/Excel file (.csv, .xlsx, .xls)",
        "bulk_note": "Ensure the file contains a column named **'name'** with the molecule names. It is recommended to use English names for the CIR API.",
        "file_format_error": "Unsupported file format. Please upload a CSV, XLSX, or XLS file.",
        "column_error": "The file must contain the required column: {}. If you're using a Spanish file, it might be 'nombre'.",
        "bulk_success": "Bulk conversion complete: **{}** successful, **{}** failed.",
        "smiles_not_found": "Molecule not found by NCI CIR or the structure could not be generated/canonicalized.",
        "invalid_smiles": "Invalid or not found SMILES.",
        "question_added": "Question added: {}",
        "questions_added_subtitle": "Questions added",
        "clear_all": "Clear all questions",
        "download_xml": "Download Moodle XML",
        "xml_error": "Error generating XML file: {}",
        "question_title": "Structure of {}",
        "question_text": "Draw the structure of <strong>{}</strong> using the molecular editor.",
        "delete_tooltip": "Delete"
    }
}


# ==============================================================================
# 2. CORE UTILITY FUNCTIONS
# ==============================================================================

def _canonicalize_smiles(smiles):
    """Internal helper to canonicalize SMILES using RDKit."""
    if not Chem:
        return smiles  # Return as is if RDKit is not available
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Generate the canonical SMILES for Moodle/JSME standardization.
            return Chem.MolToSmiles(mol, canonical=True)
        return None
    except Exception:
        return None

def name_to_smiles(compound_name):
    """
    Looks up the SMILES of a compound using the NCI CIR API.
    Canonicalizes the SMILES using RDKit for Moodle compatibility.
    """
    if not NCI_CIR_AVAILABLE:
        return None
        
    try:
        encoded_name = requests.utils.quote(compound_name)
        url = f"http://cactus.nci.nih.gov/chemical/structure/{encoded_name}/smiles"
        response = requests.get(url, timeout=10)
        
        if response.status_code != 200:
            return None
        
        smiles = response.text.strip()
        
        # Handle multiple results or error strings returned by CIR
        if "\n" in smiles:
            smiles = smiles.split('\n')[0]
            
        if not smiles or any(err in smiles for err in ["Error", "Server Error", "NOT_FOUND"]):
            return None
                
        # Canonicalize if RDKit is available
        return _canonicalize_smiles(smiles)
        
    except (requests.exceptions.Timeout, requests.exceptions.RequestException):
        return None
    except Exception:
        return None


def escape_smiles(smiles):
    """Adds an escape symbol to Moodle's special characters in a SMILES string."""
    special_characters = "()[]\\"
    new_smiles = ""
    for char in smiles:
        if char in special_characters:
            new_smiles += "\\" + char
        else:
            new_smiles += char
    return new_smiles


def process_bulk_file(uploaded_file, texts):
    """Processes a bulk file (CSV/XLSX) to convert names to SMILES and update session state."""
    
    # 1. Read the file
    try:
        if uploaded_file.name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        elif uploaded_file.name.endswith(('.xlsx', '.xls')):
            df = pd.read_excel(uploaded_file)
        else:
            st.error(texts["file_format_error"])
            return
    except Exception as e:
        st.error(f"Error reading file: {e}")
        return

    # 2. Determine the required column
    possible_cols = ['name', 'nombre']
    required_col = next((col for col in possible_cols if col in df.columns), None)
    
    if not required_col:
        col_name = "'name'" if st.session_state.lang == 'en' else "'nombre'"
        st.error(texts["column_error"].format(col_name))
        return

    # 3. Process rows
    new_questions = []
    successes = 0
    failures = 0
    total_rows = len(df)
    
    with st.spinner(f'Processing {total_rows} molecules. This may take a moment...'):
        progress_bar = st.progress(0)
        
        for i, row in df.iterrows():
            name = str(row[required_col]).strip()
            progress_bar.progress((i + 1) / total_rows, text=f"Processing {i + 1} of {total_rows}: {name}...")

            if name and name.lower() != 'nan':
                smiles = name_to_smiles(name)
                if smiles:
                    new_questions.append((name, smiles))
                    successes += 1
                else:
                    failures += 1
                    # print(f"FAILURE looking up SMILES for: {name} (Row {i+2})") # Keep discrete in console
        
    progress_bar.empty()
    
    # 4. Integrate results
    st.session_state.questions.extend(new_questions)
    st.success(texts["bulk_success"].format(successes, failures))
    st.info(f"Omitted {failures} names that could not be resolved by NCI CIR.")


def generate_xml(questions_list, lang):
    """Generates an XML element of type <quiz> from a list of questions."""
    quiz = ET.Element("quiz")
    
    for name, smiles in questions_list:
        question = ET.Element("question", type="pmatchjme")
        
        q_name = ET.SubElement(question, "name")
        text = ET.SubElement(q_name, "text")
        text.text = TEXTS[lang]["question_title"].format(name)
        
        questiontext = ET.SubElement(question, "questiontext", format="html")
        qtext = ET.SubElement(questiontext, "text")
        qtext.text = TEXTS[lang]["question_text"].format(name)
        
        answer = ET.SubElement(question, "answer")
        answer.set("fraction", "100")
        answer.set("format", "moodle_auto_format")
        answer_text = ET.SubElement(answer, "text")
        
        modelanswer = ET.SubElement(question, "modelanswer")
        modelanswer.text = smiles 
        
        escaped_smiles = escape_smiles(smiles)
        answer_text.text = f"match({escaped_smiles})"
        
        quiz.append(question)
        
    tree = ET.ElementTree(quiz)
    ET.indent(tree, space="  ")
    # Return bytes for download button
    return ET.tostring(tree.getroot(), encoding="utf-8", xml_declaration=True)


# ==============================================================================
# 3. STREAMLIT APPLICATION LOGIC
# ==============================================================================

# --- Initialization and Language Selector ---

# Initialize the language and questions session state
if "lang" not in st.session_state:
    st.session_state.lang = "en" 
if "questions" not in st.session_state:
    st.session_state.questions = []

# Get texts for the current session language
texts = TEXTS[st.session_state.lang]

# Language selector at the top
st.markdown("###### Select language / Selecciona tu idioma")
flag_col1, flag_col2 = st.columns([1, 1])

with flag_col1:
    st.image("https://upload.wikimedia.org/wikipedia/commons/a/a5/Flag_of_the_United_Kingdom_%281-2%29.svg", width=60)
    en_btn = st.button("EN", key="en_btn")

with flag_col2:
    st.image("https://upload.wikimedia.org/wikipedia/commons/8/89/Bandera_de_Espa%C3%B1a.svg", width=50)
    es_btn = st.button("ES", key="es_btn")


# Use buttons to change the language
if en_btn:
    st.session_state.lang = "en"
    st.rerun()
if es_btn:
    st.session_state.lang = "es"
    st.rerun()

st.title(texts["title"])

# Main container with two columns for the interface
main_col, list_col = st.columns([1, 1])


# --- Individual Question Submission Handler ---
def handle_individual_submission(molecule_name, smiles_input, use_api, texts):
    """Handles the form submission logic for a single question."""
    
    if not molecule_name and not smiles_input:
        st.error("Please enter the molecule name or SMILES.")
        return

    final_smiles = smiles_input.strip() if smiles_input else ""
    
    if use_api and molecule_name:
        with st.spinner('Searching and canonicalizing SMILES...'):
            final_smiles = name_to_smiles(molecule_name.strip())
        
        if not final_smiles:
            st.error(texts["smiles_not_found"])
            return
    elif use_api and not molecule_name:
        st.error("Please enter a molecule name to search the API.")
        return
    
    if final_smiles:
        if Chem and not Chem.MolFromSmiles(final_smiles):
            st.error(texts["invalid_smiles"])
            return
        
        st.session_state.questions.append((molecule_name.strip(), final_smiles))
        st.success(texts["question_added"].format(molecule_name.strip()))
    else:
        st.error(texts["invalid_smiles"])
        

# --- Left Column: Controls ---
with main_col:
    # 1. Individual Entry Section
    st.subheader(texts["section_individual"]) 
    with st.form("question_form"):
        molecule_name = st.text_input(texts["molecule_name"])
        st.info(texts["name_note"])
        col_form1, col_form2 = st.columns([1, 1])
        with col_form1:
            use_api = st.checkbox(texts["search_smiles"], value=True) 
        with col_form2:
            smiles_input = st.text_input(texts["smiles_input"])
            
        submitted = st.form_submit_button(texts["add_question"], type= "primary", icon=":material/library_add:")

    if submitted:
        handle_individual_submission(molecule_name, smiles_input, use_api, texts)

    st.markdown("---") 

    # 2. Bulk Upload Section - Now automatically processes on file upload
    st.subheader(texts["section_bulk"]) 
    with st.expander(texts["expand_bulk"]): 
        st.info(texts["bulk_note"])
        # File uploader is outside the form, allowing immediate processing
        uploaded_file = st.file_uploader(texts["upload_file"], type=['csv', 'xlsx', 'xls'])
        
        # Check if a file was uploaded or changed
        if uploaded_file is not None:
            # Check if this file has already been processed in the session to avoid re-running logic on every interaction
            if "last_uploaded_file_name" not in st.session_state or st.session_state.last_uploaded_file_name != uploaded_file.name:
                
                # Update session state to mark this file as processed
                st.session_state.last_uploaded_file_name = uploaded_file.name
                
                uploaded_file.seek(0)
                # Call the bulk processing function
                process_bulk_file(uploaded_file, texts)


    st.markdown("---")

    # 3. Export Section - Combined XML generation and download
    if st.session_state.get("questions"):
        st.subheader(texts["section_export"])
        btn_col1, btn_col2 = st.columns(2)
        
        with btn_col1:
            if st.button(texts["clear_all"], icon=":material/delete:"):
                st.session_state.questions = []
                st.rerun() 
                
        with btn_col2:
            try:
                # Generate XML bytes dynamically when the download button is rendered
                xml_bytes = generate_xml(st.session_state.questions, st.session_state.lang)
                
                st.download_button(
                    label=texts["download_xml"],
                    data=xml_bytes,
                    file_name="moodle_questions.xml",
                    mime="application/xml",
                    type="primary", # Make the download button primary
                    icon= ":material/download:"
                )
            except Exception as e:
                st.error(texts["xml_error"].format(e))
    st.markdown("---")


# --- Right Column: Questions List ---
with list_col:
    if st.session_state.get("questions"):
        st.subheader(texts["questions_added_subtitle"])
        
        def delete_question(index):
            st.session_state.questions.pop(index)
            st.rerun()
            
        for i, (n, s) in enumerate(st.session_state.questions):
            cols = st.columns([5, 1])
            with cols[0]:
                st.write(f"**{i+1}.** {n} → `{s}`")
            with cols[1]:
                st.button("🗑️", help=texts["delete_tooltip"], key=f"del_{i}", on_click=delete_question, args=(i,))
            st.markdown("---")