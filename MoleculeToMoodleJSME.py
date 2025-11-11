# -*- coding: utf-8 -*-
"""
Moodle Question Generator for Chemistry (JSME + NCI CIR).
Unifies SMILES lookup with the JSME component for standardization.
"""

import streamlit as st
import xml.etree.ElementTree as ET
import pandas as pd
import requests
import time
import uuid
import json 

# Custom JSME component import (assumes it's in the Python path)
from my_component import jsme_editor 

# Configure the page to use the wide layout
st.set_page_config(layout="wide")


# ==============================================================================
# 1. DEPENDENCY CHECKS AND CONSTANTS
# ==============================================================================

# --- Dependency Check: requests (for NCI CIR API) ---
try:
    NCI_CIR_AVAILABLE = True 
except NameError:
    # This warning is intentionally left in English as it refers to a missing Python module installation.
    st.warning("Warning: The 'requests' module is not installed. Automatic SMILES lookup will not be available.")
    NCI_CIR_AVAILABLE = False
    requests = None


# --- Dependency Check: RDKit ---
try:
    from rdkit import Chem
except ImportError:
    # This warning is intentionally left in English as it refers to a missing Python module installation.
    st.warning("Warning: The 'rdkit' module is not installed. Structure validation and initial canonicalization will be limited.")
    Chem = None


# Interface texts in different languages (Spanish and English)
TEXTS = {
    "es": {
        "title": "Generador de preguntas Moodle con SMILES (JSME + NCI CIR) 🧪",
        "section_individual": "1. Entrada Individual",
        "section_bulk": "2. Carga Masiva (CSV / Excel)",
        "expand_bulk": "Abrir Conversión Masiva",
        "section_export": "3. Exportar",
        "molecule_name": "Nombre de la molécula",
        "custom_question_name": "Nombre de la pregunta (para el enunciado)",
        "name_note_individual": "Usará el Nombre para buscar el SMILES NCI si la casilla está marcada.",
        "search_smiles": "Modo: Usar Nombre (buscar en NCI CIR y estandarizar)",
        "smiles_input": "SMILES (manual)",
        "add_pending": "Buscar/Validar y Añadir a Pendientes",
        "start_standardization": "Iniciar Estandarización JSME ({} pendientes)",
        "pending_list_title": "Lista de preguntas pendientes de estandarizar",
        "upload_file": "Seleccionar archivo CSV/Excel (.csv, .xlsx, .xls)",
        "bulk_note": "Asegúrate de que el archivo contenga una columna llamada **'nombre'** con los nombres de las moléculas. Se recomienda usar nombres en inglés para la API de CIR.",
        "file_format_error": "Formato de archivo no soportado. Por favor, sube un archivo CSV, XLSX o XLS.",
        "column_error": "El archivo debe contener la columna requerida: {}. Si estás usando un archivo en inglés, puede ser 'name'.",
        "bulk_success": "Procesamiento masivo finalizado: **{}** éxitos, **{}** fallos.",
        "smiles_not_found": "Molécula no encontrada por NCI CIR o la estructura no se pudo generar/estandarizar.", 
        "invalid_smiles": "SMILES inválido o no encontrado.",
        "question_added_pending": "Añadido a pendientes: {}",
        "questions_added_subtitle": "Preguntas añadidas (Estandarizadas con JSME)",
        "clear_all": "Borrar todas las preguntas",
        "download_xml": "Descargar XML Moodle",
        "xml_error": "Error al generar el archivo XML: {}",
        "question_title": "Estructura de {}",
        "question_text": "Dibuja la estructura de <strong>{}</strong> usando el editor molecular.",
        "delete_tooltip": "Borrar",
        "jsme_status": "Procesador JSME: ",
        "processing_wait": "Procesando... Esperando respuesta del componente JSME.",
        "start_bulk": "Iniciar Procesamiento Masivo",
        "bulk_progress": "Progreso Masivo: {} de {} moléculas procesadas.",
        "bulk_api_lookup": "Buscando SMILES para **{}** ({} de {})",
        "bulk_jsme_process": "Estandarizando SMILES para **{}** ({} de {})",
        "bulk_finished": "Procesamiento Masivo Finalizado.",
        "smiles_hint": "Si el modo es **Usar Nombre**, se busca el SMILES en NCI. Si es **Manual**, debe proporcionar el SMILES y el nombre para el enunciado.",
        "missing_question_name": "Modo Manual: Debe introducir un Nombre de la pregunta para el enunciado. Por favor, escriba el nombre en el campo superior antes de volver a pulsar 'Añadir a Pendientes'.",
        "standardization_summary": "Estandarización de pendientes finalizada: **{}** éxitos, **{}** fallos (de {} elementos).",
        # --- NEW ERROR MESSAGES FOR LOCALIZATION ---
        "processing_running": "El procesador JSME está ocupado. Por favor, espere a que finalice el proceso actual.",
        "missing_molecule_name": "Modo Búsqueda: Debe introducir un Nombre de la molécula para realizar la consulta a la API.",
        "missing_smiles_manual": "Modo Manual: Debe introducir un SMILES válido.",
        "bulk_empty_file": "El archivo está vacío o no contiene nombres válidos.",
        "bulk_no_valid_smiles": "No se encontró ningún SMILES válido para procesar después de la consulta NCI CIR.",
        "generic_file_error": "Error durante la carga o procesamiento del archivo: {}",
        "jsme_error_prefix": "Error de procesamiento JSME para '{}': {}"
    },
    "en": {
        "title": "Moodle Question Generator with SMILES (JSME + NCI CIR) 🧪",
        "section_individual": "1. Individual Entry",
        "section_bulk": "2. Bulk Upload (CSV / Excel)",
        "expand_bulk": "Open Bulk Conversion",
        "section_export": "3. Export",
        "molecule_name": "Molecule name",
        "custom_question_name": "Question Name (for the statement)",
        "name_note_individual": "It will use the Name to search for the NCI SMILES if the checkbox is checked.",
        "search_smiles": "Mode: Use Name (search NCI CIR and standardize)",
        "smiles_input": "SMILES (manual)",
        "add_pending": "Search/Validate and Add to Pending",
        "start_standardization": "Start JSME Standardization ({} pending)",
        "pending_list_title": "List of pending questions for standardization",
        "upload_file": "Select CSV/Excel file (.csv, .xlsx, .xls)",
        "bulk_note": "Ensure the file contains a column called **'name'** with the molecule names. It is recommended to use English names for the CIR API.",
        "file_format_error": "Unsupported file format. Please upload a CSV, XLSX, or XLS file.",
        "column_error": "The file must contain the required column: {}. If you are using an English file, it can be 'name'.",
        "bulk_success": "Bulk processing finished: **{}** successes, **{}** failures.",
        "smiles_not_found": "Molecule not found by NCI CIR or the structure could not be generated/standardized.", 
        "invalid_smiles": "Invalid or not found SMILES.",
        "question_added_pending": "Added to pending: {}",
        "questions_added_subtitle": "Questions Added (Standardized with JSME)",
        "clear_all": "Clear all questions",
        "download_xml": "Download Moodle XML",
        "xml_error": "Error generating XML file: {}",
        "question_title": "Structure of {}",
        "question_text": "Draw the structure of <strong>{}</strong> using the molecular editor.",
        "delete_tooltip": "Delete",
        "jsme_status": "JSME Processor: ",
        "processing_wait": "Processing... Waiting for response from JSME component.",
        "start_bulk": "Start Bulk Processing",
        "bulk_progress": "Bulk Progress: {} of {} molecules processed.",
        "bulk_api_lookup": "Searching SMILES for **{}** ({} of {})",
        "bulk_jsme_process": "Standardizing SMILES for **{}** ({} of {})",
        "bulk_finished": "Bulk Processing Finished.",
        "smiles_hint": "If the mode is **Use Name**, the SMILES is searched in NCI. If it is **Manual**, you must provide the SMILES and the name for the statement.",
        "missing_question_name": "Manual Mode: You must enter a Question Name for the statement. Please write the name in the field above before pressing 'Add to Pending' again.",
        "standardization_summary": "Pending standardization finished: **{}** successes, **{}** failures (out of {} items).",
        # --- NEW ERROR MESSAGES FOR LOCALIZATION ---
        "processing_running": "The JSME processor is busy. Please wait for the current process to finish.",
        "missing_molecule_name": "Search Mode: You must enter a Molecule Name to query the API.",
        "missing_smiles_manual": "Manual Mode: You must enter a valid SMILES.",
        "bulk_empty_file": "The file is empty or does not contain valid names.",
        "bulk_no_valid_smiles": "No valid SMILES were found to process after the NCI CIR lookup.",
        "generic_file_error": "Error during file upload or processing: {}",
        "jsme_error_prefix": "JSME processing error for '{}': {}"
    }
}


# ==============================================================================
# 2. CORE UTILITY FUNCTIONS
# ==============================================================================

def _canonicalize_smiles(smiles):
    """Internal helper to canonicalize SMILES using RDKit."""
    if not Chem:
        return smiles 
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Returns the canonical SMILES.
            return Chem.MolToSmiles(mol, canonical=True)
        # If parsing fails, return None (indicating invalid SMILES)
        return None
    except Exception:
        return None


def name_to_smiles(compound_name):
    """
    Looks up the SMILES of a compound using the NCI CIR API.
    Canonicalizes the SMILES using RDKit for initial cleanup.
    """
    if not NCI_CIR_AVAILABLE:
        return None
        
    try:
        encoded_name = requests.utils.quote(compound_name)
        url = f"http://cactus.nci.nih.gov/chemical/structure/{encoded_name}/smiles"
        response = requests.get(url, timeout=15) 
        
        if response.status_code != 200:
            return None
        
        smiles = response.text.strip()
        
        if "\n" in smiles:
            smiles = smiles.split('\n')[0]
            
        if not smiles or any(err in smiles for err in ["Error", "Server Error", "NOT_FOUND"]):
            return None
            
        # Initial canonicalization with RDKit (if available)
        return _canonicalize_smiles(smiles)
            
    except (requests.exceptions.Timeout, requests.exceptions.RequestException):
        # Handle connection errors or timeouts
        return None
    except Exception:
        return None


def escape_smiles(smiles):
    """Escapes special SMILES characters for Moodle's pmatchjme answer format."""
    if smiles is None:
        return ""
    s = str(smiles)
    # Escape backslashes, parentheses, and brackets for Moodle's regex-like matching
    s = s.replace("\\", "\\\\")
    s = s.replace("(", "\\(").replace(")", "\\)").replace("[", "\\[").replace("]", "\\]")
    return s


def generate_xml_local(questions_list, lang):
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
    return ET.tostring(tree.getroot(), encoding="utf-8", xml_declaration=True)


# ==============================================================================
# 3. INITIALIZATION AND STATE MANAGEMENT
# ==============================================================================

# Initialize the language and questions session state
if "lang" not in st.session_state:
    st.session_state.lang = "en" 
if "questions" not in st.session_state:
    st.session_state.questions = []
    
# State variables for the asynchronous JSME processing flow
if "jsme_input" not in st.session_state:
    st.session_state.jsme_input = None 
if "jsme_processed_output" not in st.session_state:
    st.session_state.jsme_processed_output = None 

# State tracking the current molecule being processed by JSME (ID, name, bulk info)
if "pending_request" not in st.session_state:
    st.session_state.pending_request = None 
if "is_processing" not in st.session_state:
    st.session_state.is_processing = False

# List of molecules awaiting JSME standardization (after NCI lookup or manual entry)
if "pending_smiles_list" not in st.session_state:
    st.session_state.pending_smiles_list = [] # [{'name': 'MolX', 'smiles': 'C1=...'}]

# Bulk queue and results tracking for the asynchronous processing block
if "bulk_queue" not in st.session_state:
    st.session_state.bulk_queue = [] 
if "bulk_results" not in st.session_state:
    st.session_state.bulk_results = {'success': 0, 'failure': 0, 'total': 0}
    
# Checkbox state for input mode
if "use_api_search_form" not in st.session_state:
    # Initialize the default state to True (API Search)
    st.session_state.use_api_search_form = True

# Persistent error message for the individual entry form
if "individual_error_message" not in st.session_state:
    st.session_state.individual_error_message = None

# States for persistent job summaries
if "last_job_type" not in st.session_state:
    st.session_state.last_job_type = None # 'pending' or 'bulk_upload'

if "last_job_summary" not in st.session_state:
    st.session_state.last_job_summary = None # The persistent summary text to display

texts = TEXTS[st.session_state.lang]


# ==============================================================================
# 4. JSME ASYNCHRONOUS PROCESSOR BLOCK
# ==============================================================================

def start_jsme_processing(molecule_info):
    """Sets the state to start asynchronous processing for a molecule."""
    
    # 1. Generate unique ID for this request
    request_id = str(uuid.uuid4())
    
    # 2. Input JSON structure for React
    jsme_input_json = json.dumps({
        'smiles': molecule_info['smiles'],
        'id': request_id
    })
    
    # 3. Set the pending request marker (includes queue info)
    st.session_state.pending_request = {
        'request_id': request_id,
        'name': molecule_info['name'],
        'is_bulk': molecule_info.get('is_bulk', False),
        'index': molecule_info.get('index', 0),
        'total': molecule_info.get('total', 0)
    }
    
    st.session_state.jsme_input = jsme_input_json # <- This triggers the action in React
    st.session_state.jsme_processed_output = None # Clear previous output
    st.session_state.is_processing = True

# Renders the JSME component only once with a fixed key.
jsme_output = jsme_editor(
    smiles_json=st.session_state.jsme_input, # <-- Expects a JSON string
    key="jsme_global_processor"
)

# Update session state with the actual output
st.session_state.jsme_processed_output = jsme_output

# -----------------
# Completion and Transition Logic
# -----------------
if st.session_state.jsme_input and st.session_state.pending_request and st.session_state.jsme_processed_output:
    
    st.session_state.is_processing = True # Ensure processing state is active

    pending_req = st.session_state.pending_request
    pending_id = pending_req['request_id']
    pending_name = pending_req['name']
    is_bulk = pending_req.get('is_bulk', False)
    
    output_json = st.session_state.jsme_processed_output
    output_id = None
    output_smiles = None

    try:
        output_data = json.loads(output_json)
        output_smiles = output_data.get('smiles')
        output_id = output_data.get('id')
    except (json.JSONDecodeError, TypeError):
        # Handle transient value or error
        if output_json in ["", "null"]: pass
        elif output_json.startswith("JSME_ERROR:"):
            # This is a JSME error, process it as a failure
            st.error(texts["jsme_error_prefix"].format(pending_name, output_json)) # Localized error
            st.session_state.bulk_results['failure'] += 1
    
    
    # === ID VALIDATION ===
    if output_id == pending_id:
        
        # 1. Result Registration Process
        is_valid_smiles = output_smiles and not str(output_smiles).startswith("JSME_ERROR:")

        if is_valid_smiles:
            smiles_to_add = output_smiles
            
            # Final SMILES validation with RDKit if available
            if Chem and not Chem.MolFromSmiles(smiles_to_add):
                # If RDKit validation fails (invalid structure)
                st.error(texts["invalid_smiles"] + f" (Output SMILES: {smiles_to_add})") # Localized error
                if is_bulk: st.session_state.bulk_results['failure'] += 1
            else:
                # Success: Add to the final list of standardized questions
                st.session_state.questions.append((pending_name, smiles_to_add))
                if not is_bulk:
                    # Only show success for an immediate individual request (not in queue)
                    st.success(f"Question **{pending_name}** standardized and added.")
                if is_bulk: st.session_state.bulk_results['success'] += 1
                
        else:
             # Output is a specific JSME error (not standardizable)
             if not is_bulk and not str(output_json).startswith("JSME_ERROR:"):
                # Use a more descriptive error message here
                st.error(texts["jsme_error_prefix"].format(pending_name, output_smiles or 'Empty/invalid output')) # Localized error
             if is_bulk and not str(output_json).startswith("JSME_ERROR:"):
                 st.session_state.bulk_results['failure'] += 1


        # 2. Clear processing state
        st.session_state.jsme_input = None 
        st.session_state.jsme_processed_output = None 
        st.session_state.pending_request = None
        st.session_state.is_processing = False
        
        # 3. Queue Transition Logic
        if st.session_state.bulk_queue:
            # Elements in the queue, start the next one immediately.
            next_molecule = st.session_state.bulk_queue.pop(0)
            start_jsme_processing(next_molecule)
            st.rerun() # Force rerun to start the next process
        else:
            # The queue is empty, processing is finished.
            if is_bulk:
                # Calculate final results and store the persistent summary message
                t = st.session_state.bulk_results['success'] + st.session_state.bulk_results['failure']
                st.session_state.bulk_results['total'] = t
                
                s = st.session_state.bulk_results['success']
                f = st.session_state.bulk_results['failure']
                job_type = st.session_state.last_job_type

                if t > 0:
                    if job_type == 'pending':
                         # Pending Standardization Summary
                         st.session_state.last_job_summary = texts["standardization_summary"].format(s, f, t)
                    elif job_type == 'bulk_upload':
                         # Bulk Upload Summary
                         st.session_state.last_job_summary = texts["bulk_success"].format(s, f) + f" (Total processed by JSME: **{t}**)"
                
                # Clear temporary tracking state
                st.session_state.bulk_results = {'success': 0, 'failure': 0, 'total': 0}
                st.session_state.last_job_type = None 
            
            # Clear the processing status placeholder below
            # This will be done by the UI logic (Section 7) on the next rerun
            st.rerun() 
    
    # PENDING: Output has not arrived (ID does not match or is None/empty).
    if st.session_state.is_processing:
        # Halt execution to wait for JSME component response
        pass
else:
    # If the process ended (success/error) or hasn't started
    st.session_state.is_processing = False


# ==============================================================================
# 5. INDIVIDUAL ENTRY HANDLER
# ==============================================================================

def handle_add_to_pending(molecule_name, smiles_input, custom_name_input, use_api, texts):
    """
    Searches for SMILES (if necessary) and adds the unstandardized result 
    to the pending list (pending_smiles_list). Persists errors via session_state.
    """
    
    # 1. Clear previous state messages
    st.session_state.individual_error_message = None 
    st.session_state.bulk_results = {'success': 0, 'failure': 0, 'total': 0} 
    st.session_state.last_job_summary = None # Clear persistent summary
    st.session_state.last_job_type = None 

    if st.session_state.is_processing:
        st.session_state.individual_error_message = texts["processing_running"] # Localized error
        return

    molecule_name = molecule_name.strip()
    smiles_input = smiles_input.strip()
    custom_name_input = custom_name_input.strip() 
    
    final_smiles_to_store = None
    name_for_question = None

    # Use a specific placeholder for the "status while searching"
    local_status_placeholder = st.empty()
    local_status_placeholder.info(f"{texts['jsme_status']} Preparing SMILES...")

    # 2. CENTRALIZED DECISION LOGIC
    if use_api:
        # NAME SEARCH MODE (NCI CIR API)
        if not molecule_name:
            local_status_placeholder.empty()
            st.session_state.individual_error_message = texts["missing_molecule_name"] # Localized error
            return

        name_for_question = molecule_name
        with st.spinner('Searching NCI-CIR and canonicalizing SMILES...'):
            final_smiles_to_store = name_to_smiles(name_for_question)
        
        if not final_smiles_to_store:
            local_status_placeholder.empty()
            st.session_state.individual_error_message = texts["smiles_not_found"] # Localized error
            return
            
    else:
        # MANUAL SMILES ENTRY MODE
        if not smiles_input:
            local_status_placeholder.empty()
            st.session_state.individual_error_message = texts["missing_smiles_manual"] # Localized error
            return

        # VALIDATION: If SMILES is present but Question Name is missing
        if not custom_name_input:
            local_status_placeholder.empty()
            st.session_state.individual_error_message = "🚨 " + texts["missing_question_name"] # Localized error
            return

        name_for_question = custom_name_input
        
        with st.spinner('Canonicalizing SMILES with RDKit...'):
             # RDKit canonicalization for manual input
             final_smiles_to_store = _canonicalize_smiles(smiles_input)

        if not final_smiles_to_store:
            local_status_placeholder.empty()
            st.session_state.individual_error_message = texts["invalid_smiles"] # Localized error
            return
            
    # 3. Final Validation (Applies to API and Manual)
    if not final_smiles_to_store:
        local_status_placeholder.empty()
        st.session_state.individual_error_message = texts["invalid_smiles"] # Localized error
        return
        
    # 4. Add to the pending list
    st.session_state.pending_smiles_list.append({
        'name': name_for_question, 
        'smiles': final_smiles_to_store
    })
    
    # 5. Show success and clean up
    st.success(texts["question_added_pending"].format(name_for_question))
        
    local_status_placeholder.empty()
    st.session_state.individual_error_message = None # Clear error on success

def start_pending_standardization(texts):
    """Moves the pending list to the JSME processing queue and starts the cycle."""
    if not st.session_state.pending_smiles_list:
        st.warning("The list of pending questions is empty.") # Already localized through hardcoded string (Spanish is the default warning language)
        return
        
    if st.session_state.is_processing:
        st.warning(texts["processing_running"]) # Localized error
        return

    # Clear previous error and summary state before starting the new job
    st.session_state.individual_error_message = None 
    st.session_state.bulk_results = {'success': 0, 'failure': 0, 'total': 0}
    st.session_state.last_job_summary = None 

    # 1. Prepare the JSME queue
    total_mols = len(st.session_state.pending_smiles_list)
    
    # Format items for the JSME queue
    jsme_queue_items = [{
        'name': item['name'], 
        'smiles': item['smiles'],
        'is_bulk': True, # Treat as "bulk" because it's in a queue
        'total': total_mols,
        'index': i + 1 
    } for i, item in enumerate(st.session_state.pending_smiles_list)]
    
    # 2. Initialize queue and bulk state
    st.session_state.bulk_queue = jsme_queue_items
    st.session_state.pending_smiles_list = [] # Empty the pending list
    st.session_state.bulk_results = {'success': 0, 'failure': 0, 'total': len(jsme_queue_items)} # Set total for initial progress bar
    st.session_state.last_job_type = 'pending' # Set job type for summary
    
    # 3. Start the first element in the queue
    first_molecule = st.session_state.bulk_queue.pop(0)
    start_jsme_processing(first_molecule)
    st.rerun() # Starts the asynchronous cycle


# ==============================================================================
# 6. BULK PROCESSING LOGIC
# ==============================================================================

def handle_bulk_upload(uploaded_file, texts):
    """
    Reads the file, searches for SMILES, and fills the processing queue.
    """
    if st.session_state.is_processing or st.session_state.bulk_queue:
        st.warning(texts["processing_running"]) # Localized error
        return
        
    # Clear previous error and summary state before starting the new job
    st.session_state.individual_error_message = None 
    st.session_state.bulk_results = {'success': 0, 'failure': 0, 'total': 0}
    st.session_state.last_job_summary = None 

    try:
        # 1. Read file
        if uploaded_file.name.endswith('.csv'):
            df = pd.read_csv(uploaded_file)
        elif uploaded_file.name.endswith(('.xlsx', '.xls')):
            df = pd.read_excel(uploaded_file)
        else:
            st.error(texts["file_format_error"]) # Localized error
            return

        # 2. Validate column
        name_col = 'nombre'
        if name_col not in df.columns:
            # Try with 'name' if it's in English
            name_col = 'name'
            if name_col not in df.columns:
                st.error(texts["column_error"].format("'nombre' o 'name'")) # Localized error
                return
        
        names = df[name_col].astype(str).tolist()
        total_mols = len(names)
        
        if total_mols == 0:
             st.warning(texts["bulk_empty_file"]) # Localized error
             return
             
        # 3. Clean and Search SMILES (Initial Blocking)
        jsme_queue_temp = []
        bulk_lookup_success = 0
        
        progress_bar = st.progress(0, text="Starting NCI CIR lookup...")
        
        for i, name in enumerate(names):
            name = name.strip()
            progress_bar.progress((i + 1) / total_mols, text=texts['bulk_api_lookup'].format(name, i + 1, total_mols))
            
            # API Lookup (Blocking)
            smiles = name_to_smiles(name)
            
            if smiles:
                bulk_lookup_success += 1
                jsme_queue_temp.append({
                    'name': name, 
                    'smiles': smiles,
                    'is_bulk': True,
                    'total': total_mols,
                    'index': i + 1 
                })
        
        progress_bar.empty()
        # The result message is kept in English/Spanish based on texts:
        st.info(f"NCI CIR search finished: {bulk_lookup_success} SMILES found out of {total_mols} names.")

        if not jsme_queue_temp:
            st.error(texts["bulk_no_valid_smiles"]) # Localized error
            return

        # 4. Initialize JSME queue and bulk state
        st.session_state.bulk_queue = jsme_queue_temp
        st.session_state.pending_smiles_list = [] # Ensure the manual list is empty
        st.session_state.bulk_results = {'success': 0, 'failure': 0, 'total': len(jsme_queue_temp)}
        st.session_state.last_job_type = 'bulk_upload' # Set job type for summary
        
        # 5. Start the first element in the queue
        first_molecule = st.session_state.bulk_queue.pop(0)
        start_jsme_processing(first_molecule)
        st.rerun() # Starts the asynchronous cycle
        
    except Exception as e:
        st.error(texts["generic_file_error"].format(e)) # Localized error
        # Clean up states in case of failure
        st.session_state.is_processing = False
        st.session_state.bulk_queue = [] 
        st.session_state.pending_smiles_list = []
        st.session_state.bulk_results = {'success': 0, 'failure': 0, 'total': 0}
        st.session_state.last_job_summary = None
        st.session_state.last_job_type = None


# ==============================================================================
# 7. STREAMLIT UI LAYOUT
# ==============================================================================

# Language selector at the top
st.markdown("###### Select language / Selecciona tu idioma")
flag_col1, flag_col2 = st.columns([1, 1])

with flag_col1:
    st.image("https://upload.wikimedia.org/wikipedia/commons/a/a5/Flag_of_the_United_Kingdom_%281-2%29.svg", width=60)
    en_btn = st.button("EN", key="en_btn", disabled=st.session_state.is_processing)

with flag_col2:
    st.image("https://upload.wikimedia.org/wikipedia/commons/8/89/Bandera_de_Espa%C3%B1a.svg", width=50)
    es_btn = st.button("ES", key="es_btn", disabled=st.session_state.is_processing)


if en_btn and not st.session_state.is_processing:
    st.session_state.lang = "en"
    st.rerun()
if es_btn and not st.session_state.is_processing:
    st.session_state.lang = "es"
    st.rerun()

st.title(texts["title"])

# Placeholder for JSME processing status (in Section 4)
if st.session_state.is_processing and st.session_state.pending_request:
    pending_req = st.session_state.pending_request
    pending_name = pending_req['name']

    if pending_req.get('is_bulk', False):
        current_count = pending_req['total'] - len(st.session_state.bulk_queue)
        total = pending_req['total']
        # Use global queue progress
        st.info(texts['bulk_jsme_process'].format(pending_name, current_count, total))
        progress = current_count / total
        st.progress(progress)
    else:
        # Immediate individual process
        st.info(f"⏳ **{texts['processing_wait']}** ({pending_name})")

# Main container with two columns for the interface
main_col, list_col = st.columns([1, 1])


with main_col:
    # 1. Individual Entry Section
    st.subheader(texts["section_individual"]) 
    
    disabled_input = st.session_state.is_processing
    
    # Checkbox for API search mode
    st.checkbox(
        texts["search_smiles"], 
        key="use_api_search_form",  
        disabled=disabled_input,
    ) 
    
    st.info(f"ℹ️ {texts['smiles_hint']}") 
    
    # Dedicated placeholder for submission errors
    submission_error_placeholder = st.empty()
    
    # Display persistent error messages
    if st.session_state.individual_error_message:
        submission_error_placeholder.error(st.session_state.individual_error_message)
    else:
        submission_error_placeholder.empty()
    
    # Form to add to the pending list
    with st.form("pending_form", clear_on_submit=True): 
        
        molecule_name = ""
        smiles_input = ""
        custom_name_input = ""
        
        is_api_mode = st.session_state.use_api_search_form
        
        if is_api_mode:
            # --- API SEARCH MODE ---
            st.markdown("#### Búsqueda por Nombre (NCI CIR)")
            molecule_name = st.text_input(
                texts["molecule_name"], 
                disabled=disabled_input,
                key="api_molecule_name_input"
            )
            
        else:
            # --- MANUAL SMILES MODE ---
            st.markdown("#### Entrada Manual de SMILES")
            
            # 1. SMILES Field
            smiles_input = st.text_input(
                texts["smiles_input"], 
                disabled=disabled_input,
                key="manual_smiles_input"
            )
            
            # 2. Question Name Field
            custom_name_input = st.text_input(
                texts["custom_question_name"], 
                placeholder="Ej: Aspirina",
                disabled=disabled_input,
                key="manual_question_name_input"
            )
            
        st.markdown("---") # Visual separator

        # --- SUBMISSION BUTTON ---
        submitted = st.form_submit_button(
            texts["add_pending"], 
            type= "primary", 
            icon=":material/add_task:", 
            disabled=disabled_input
        )

    if submitted:
        use_api_value = st.session_state.use_api_search_form
        
        handle_add_to_pending(
            molecule_name, 
            smiles_input, 
            custom_name_input, 
            use_api_value, 
            texts, 
        )
        
        if not st.session_state.is_processing:
             st.rerun() 
    
    st.markdown("---")
    
    # Pending List and Standardization Button
    num_pending = len(st.session_state.pending_smiles_list)
    
    st.subheader(texts["pending_list_title"])
    
    # Display persistent job summary logic
    if st.session_state.last_job_summary:
         st.info(st.session_state.last_job_summary)
    
    st.button(
        texts["start_standardization"].format(num_pending),
        type="primary",
        icon=":material/rocket_launch:",
        disabled=st.session_state.is_processing or num_pending == 0,
        on_click=start_pending_standardization,
        args=(texts,)
    )
    
    # List of items pending
    if num_pending > 0:
        for i, item in enumerate(st.session_state.pending_smiles_list):
             st.markdown(f"**{i+1}.** {item['name']} → `{item['smiles']}` (SMILES Canónico)")
        st.markdown("---")
    else:
        st.markdown("_No hay moléculas pendientes de estandarización JSME._")


    # 2. Bulk Upload Section
    st.subheader(texts["section_bulk"]) 
    with st.expander(texts["expand_bulk"]): 
        st.info(texts["bulk_note"])
        
        uploaded_file = st.file_uploader(texts["upload_file"], type=['csv', 'xlsx', 'xls'], disabled=st.session_state.is_processing)
        
        start_bulk_btn = st.button(
            texts["start_bulk"], 
            type="secondary", 
            icon=":material/cloud_upload:",
            disabled=st.session_state.is_processing or not uploaded_file
        )
        
        if start_bulk_btn:
            handle_bulk_upload(uploaded_file, texts)

    st.markdown("---")
    
    # 3. Export Section
    if st.session_state.get("questions"):
        st.subheader(texts["section_export"])
        btn_col1, btn_col2 = st.columns(2)
        
        def clear_data_and_rerun():
            # Clears not only questions but all tracking states
            st.session_state.questions = []
            st.session_state.bulk_queue = []
            st.session_state.pending_smiles_list = []
            st.session_state.bulk_results = {'success': 0, 'failure': 0, 'total': 0}
            st.session_state.individual_error_message = None 
            st.session_state.last_job_summary = None # Clear persistent summary
            st.session_state.last_job_type = None
            st.rerun()
        
        with btn_col1:
            if st.button(texts["clear_all"], icon=":material/delete:", disabled=st.session_state.is_processing, on_click=clear_data_and_rerun):
                pass
                
        with btn_col2:
            try:
                xml_bytes = generate_xml_local(st.session_state.questions, st.session_state.lang)
                
                st.download_button(
                    label=texts["download_xml"],
                    data=xml_bytes,
                    file_name="moodle_questions.xml",
                    mime="application/xml",
                    type="primary", 
                    icon= ":material/download:",
                    disabled=st.session_state.is_processing
                )
            except Exception as e:
                st.error(texts["xml_error"].format(e)) # Localized error
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
                st.button("🗑️", help=texts["delete_tooltip"], key=f"del_{i}", on_click=delete_question, args=(i,), disabled=st.session_state.is_processing)
            st.markdown("---")