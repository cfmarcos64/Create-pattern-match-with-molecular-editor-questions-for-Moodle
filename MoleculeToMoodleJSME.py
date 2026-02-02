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
        "add_pending": "Añadir a la lista pendiente",
        "api_error_both_names": "⚠️ Por favor, introduce el nombre en inglés y el nombre en español.",
        "bulk_api_lookup": "Buscando en NCI CIR: **{0}** ({1}/{2})",
        "bulk_empty_file": "⚠️ El archivo no contiene nombres de moléculas.",
        "bulk_jsme_process": "Procesando: **{0}** ({1}/{2}). Esperando respuesta de JSME...",
        "bulk_no_valid_smiles": "🚨 No se pudo encontrar ningún SMILES válido en NCI CIR desde el archivo cargado.",
        "bulk_note": "Sube un archivo .csv o .xlsx con dos columnas: **'name'** (nombre en inglés para buscar el SMILES) y **'nombre'** (nombre en español para la pregunta).",
        "bulk_success": "✅ Carga masiva finalizada. Éxito: **{0}**, Fallo: **{1}**.",
        "clear_all": "Limpiar todas las preguntas",
        "column_error": "🚨 El archivo debe contener las columnas {0}.",
        "custom_question_name": "Nombre de la Pregunta (Ej: Ácido acetilsalicílico)",
        "delete_tooltip": "Eliminar esta pregunta",
        "download_xml": "Descargar XML de Moodle",
        "expand_bulk": "Instrucciones de Carga Masiva",
        "file_format_error": "🚨 Formato de archivo no soportado. Usa .csv, .xlsx o .xls.",
        "generic_file_error": "🚨 Error al procesar el archivo: {0}",
        "invalid_smiles": "🚨 SMILES inválido detectado (falló la validación con RDKit).",
        "jsme_error_prefix": "❌ Error al procesar '{0}': El editor JSME devolvió un error: {1}",
        "jsme_status": "Estado:",
        "missing_molecule_name": "🚨 Por favor, introduce el nombre de una molécula para buscar.",
        "missing_question_name": "Debes especificar un nombre para la pregunta cuando introduces el SMILES manualmente.",
        "missing_smiles_manual": "🚨 Por favor, introduce una cadena SMILES válida.",
        "molecule_name": "Nombre de la Molécula (Ej: Cafeína)",
        "name_search": "Búsqueda por Nombre (NCI CIR)",
        "no_pending_text": "No hay moléculas pendientes.",
        "pending_list_title": "2. Lista Pendiente de Estandarización",
        "processing_running": "⚠️ Ya hay un proceso de estandarización en curso. Espera a que termine.",
        "processing_wait": "Esperando respuesta del editor JSME para la estandarización...",
        "question_added_pending": "✅ Pregunta '{0}' añadida a la lista pendiente.",
        "question_text": "Utiliza el editor JSME para dibujar la estructura molecular del compuesto: {0}",
        "question_title": "Dibuja la estructura de {0}",
        "questions_added_subtitle": "Preguntas Estandarizadas Listas para Exportar",
        "search_smiles": "Buscar SMILES por Nombre (API NCI CIR)",
        "section_bulk": "3. Carga Masiva (Búsqueda por API)",
        "section_export": "4. Exportar Cuestionario Moodle",
        "section_individual": "1. Entrada Individual (Nombre o SMILES)",
        "smiles_hint": "Las cadenas SMILES de salida se envían a un componente JSME temporal para ser estandarizadas a la forma canónica preferida por el editor.",
        "smiles_input": "Cadena SMILES (Ej: O=C(C)Oc1ccccc1C(=O)O)",
        "smiles_not_found": "🚨 No se pudo encontrar el SMILES para ese nombre en NCI CIR.",
        "smiles_search": "Entrada Manual de SMILES",
        "start_bulk": "Iniciar Búsqueda y Estandarización Masiva",
        "start_NCI_lookup": "Comenzando búsqueda en NCI CIR...",
        "start_standardization": "Estandarizar {0} Pregunta(s) Pendiente(s) con JSME",
        "standardization_summary": "✅ Estandarización finalizada. Éxito: **{0}**, Fallo: **{1}** (Total: **{2}**).",
        "title": "Generador de Preguntas Químicas para Moodle (JSME + NCI CIR) 🧪",
        "xml_error": "🚨 Error al generar el XML: {0}"
    },
    "en": {
        "add_pending": "Add to Pending List",
        "api_error_both_names": "⚠️ Please enter both the English name and the Spanish name.",
        "bulk_api_lookup": "Searching NCI CIR: **{0}** ({1}/{2})",
        "bulk_empty_file": "⚠️ The file contains no molecule names.",
        "bulk_jsme_process": "Processing: **{0}** ({1}/{2}). Waiting for JSME response...",
        "bulk_no_valid_smiles": "🚨 No valid SMILES could be found in NCI CIR from the uploaded file.",
        "bulk_note": "Upload a .csv or .xlsx file with a column named **'name'**. The system will search for the SMILES for each name and then standardize it with JSME.",
        "bulk_success": "✅ Bulk upload finished. Success: **{0}**, Failure: **{1}**.",
        "clear_all": "Clear all questions",
        "column_error": "🚨 The file must contain a column named {0}.",
        "custom_question_name": "Question Name (Ex: Acetylsalicylic Acid)",
        "delete_tooltip": "Delete this question",
        "download_xml": "Download Moodle XML",
        "expand_bulk": "Bulk Upload Instructions",
        "file_format_error": "🚨 Unsupported file format. Please use .csv, .xlsx, or .xls.",
        "generic_file_error": "🚨 Error processing the file: {0}",
        "invalid_smiles": "🚨 Invalid SMILES detected (RDKit validation failed).",
        "jsme_error_prefix": "❌ Error processing '{0}': The JSME editor returned an error: {1}",
        "jsme_status": "Status:",
        "missing_molecule_name": "🚨 Please enter a molecule name to search.",
        "missing_question_name": "You must specify a name for the question when entering SMILES manually.",
        "missing_smiles_manual": "🚨 Please enter a valid SMILES string.",
        "molecule_name": "Molecule Name (Ex: Caffeine)",
        "name_search": "Molecule Search (NCI CIR)",
        "no_pending_text": "There are no pending molecules.",
        "pending_list_title": "2. Pending Standardization List",
        "processing_running": "⚠️ A standardization process is already running. Please wait for it to finish.",
        "processing_wait": "Waiting for response from the JSME editor for standardization...",
        "question_added_pending": "✅ Question '{0}' added to the pending list.",
        "question_text": "Use the JSME editor to draw the molecular structure of the compound: {0}",
        "question_title": "Draw the structure of {0}",
        "questions_added_subtitle": "Standardized Questions Ready for Export",
        "search_smiles": "Search SMILES by Name (NCI CIR API)",
        "section_bulk": "3. Bulk Upload (API Search)",
        "section_export": "4. Export Moodle Quiz",
        "section_individual": "1. Individual Entry (Name or SMILES)",
        "smiles_hint": "The output SMILES strings are sent to a temporary JSME component to be standardized to the editor's preferred canonical form.",
        "smiles_input": "SMILES String (Ex: O=C(C)Oc1ccccc1C(=O)O)",
        "smiles_not_found": "🚨 Could not find the SMILES for that name in NCI CIR.",
        "smiles_search": "SMILES Manual Input",
        "start_bulk": "Start Bulk Search and Standardization",
        "start_NCI_lookup": "Starting NCI CIR lookup...",
        "start_standardization": "Standardize {0} Pending Question(s) with JSME",
        "standardization_summary": "✅ Standardization finished. Success: **{0}**, Failure: **{1}** (Total: **{2}**).",
        "title": "Moodle Chemical Question Generator (JSME + NCI CIR) 🧪",
        "xml_error": "🚨 Error generating XML: {0}"
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

# IMPORTANT: Ensure TEXTS is defined before this line
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
    
    lang = st.session_state.lang
    
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
        if lang == "es":
            # REQUIRED: Both name fields filled
            if not molecule_name or not custom_name_input:
                st.session_state.individual_error_message = texts["api_error_both_names"]
                return
            lookup_name = molecule_name      # Name in English for the API
            name_for_question = custom_name_input # Name in Spanish for the question
        else:
            # REQUIRED: At least name in English
            if not molecule_name:
                st.session_state.individual_error_message = texts["missing_molecule_name"]
                return
            lookup_name = molecule_name
            name_for_question = molecule_name
        
        # Call the API
        local_status_placeholder.info(f"{texts['jsme_status']} Searching SMILES for '{lookup_name}'...")
        with st.spinner('Searching NCI-CIR...'):
            final_smiles_to_store = name_to_smiles(lookup_name)
        
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
        st.warning(texts["bulk_empty_file"]) # Using a relevant localized warning
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
    Bulk search conditional by language.
    """
    if st.session_state.is_processing or st.session_state.bulk_queue:
        st.warning(texts["processing_running"]) # Localized error
        return
        
    lang = st.session_state.lang  # Detects current language
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

        # 2. Validate column according to language
        if lang == "es":
            if 'name' not in df.columns or 'nombre' not in df.columns:
                st.error(texts["column_error"].format("'name' y 'nombre'"))
                return
            required_cols = ['name', 'nombre']
        else: # English
            if 'name' not in df.columns:
                st.error(texts["column_error"].format("'name'"))
                return
            required_cols = ['name']
        
        # Filter empty rows
        # Remove rows where required columns are void (NaN)
        df = df.dropna(subset=required_cols)
        
        # Remove rows where required columns are blank
        for col in required_cols:
            df = df[df[col].astype(str).str.strip() != ""]
        
        total_mols = len(df)
        if total_mols == 0:
              st.warning(texts["bulk_empty_file"])
              return
              
        # 3. Clean and Search SMILES (Initial Blocking)
        jsme_queue_temp = []
        bulk_lookup_success = 0
        
        progress_bar = st.progress(0, text=texts["start_NCI_lookup"])
        
        for i, row in df.iterrows():
            # LÓGICA DE NOMBRES
            if lang == "es":
                lookup_name = str(row['name']).strip()   # Name in Engllish for API
                display_name = str(row['nombre']).strip() # Name in Spanish for the question
            else:
                lookup_name = str(row['name']).strip()   # Name in Engllish for API
                display_name = lookup_name               # Name in Engllish for the question

            progress_bar.progress((i + 1) / total_mols, text=texts['bulk_api_lookup'].format(lookup_name, i + 1, total_mols))            
            
            # API Lookup (Blocking)
            smiles = name_to_smiles(lookup_name)
            
            if smiles:
                bulk_lookup_success += 1
                jsme_queue_temp.append({
                    'name': display_name, 
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
        lang = st.session_state.lang
        
        if is_api_mode:
            # --- API SEARCH MODE ---
            st.subheader(texts["name_search"])
            if lang == "es":
                molecule_name = st.text_input(
                    "Nombre en Inglés (para búsqueda API)", 
                    disabled=disabled_input,
                    key="api_mol_en",
                    placeholder="Ej: Caffeine"
                )
                custom_name_input = st.text_input(
                    "Nombre en Español (para la pregunta)", 
                    disabled=disabled_input,
                    key="api_mol_es",
                    placeholder="Ej: Cafeína"
                )
            else:
                # English
                molecule_name = st.text_input(
                    texts["molecule_name"], 
                    disabled=disabled_input,
                    key="api_molecule_name_input"
                )
            
        else:
            # --- MANUAL SMILES MODE ---
            st.subheader(texts["smiles_search"])
            
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
