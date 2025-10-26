# Moodle Question Generator with SMILES (NCI CIR Compatible)

A Streamlit-based web application to generate Moodle-compatible XML files for molecular structure questions using SMILES notation. The tool supports both individual molecule entry and bulk processing from CSV/Excel files, leveraging the NCI Chemical Identifier Resolver (CIR) API for SMILES lookup and RDKit for canonicalization.

## Features
- **Individual Entry**: Enter a molecule name to fetch its SMILES string via the NCI CIR API, with optional manual SMILES input.
- **Bulk Upload**: Upload CSV or Excel files containing molecule names to process multiple molecules at once.
- **SMILES Canonicalization**: Uses RDKit (optional) to standardize SMILES strings for compatibility with Moodle's JSME molecular editor.
- **Moodle XML Export**: Generates XML files formatted for Moodle's `pmatchjme` question type.
- **Multilingual Support**: Interface available in English and Spanish.

## Requirements
- Python 3.8 or higher
- Dependencies listed in `requirements.txt`:
  - `streamlit`
  - `pandas`
  - `requests`
  - `rdkit` (optional, for SMILES validation and canonicalization)
  - `lxml` (for Excel file support)

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/<your-repo-name>.git
   cd <your-repo-name>
   ```

2. Create a virtual environment (recommended):
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

4. (Optional) Install RDKit for SMILES canonicalization:
   ```bash
   pip install rdkit
   ```

## Usage
1. Run the Streamlit app:
   ```bash
   streamlit run app.py
   ```

2. Access the app in your browser (typically at `http://localhost:8501`).

3. **Individual Entry**:
   - Enter a molecule name (e.g., "benzene") and check "Search for SMILES automatically" to fetch SMILES via the NCI CIR API.
   - Alternatively, manually enter a SMILES string.
   - Click "Add question" to include the molecule in the question list.

4. **Bulk Upload**:
   - Upload a CSV or Excel file with a column named `name` (or `nombre` for Spanish files) containing molecule names.
   - The app will process the file and fetch SMILES strings for each molecule.

5. **Export**:
   - Once questions are added, download the generated Moodle XML file using the "Download Moodle XML" button.
   - Import the XML file into Moodle to create `pmatchjme` questions.

## Notes
- The NCI CIR API is used for automatic SMILES lookup. Ensure a stable internet connection for this feature.
- RDKit is optional but recommended for validating and canonicalizing SMILES strings to ensure compatibility with Moodle's JSME editor.
- For bulk uploads, ensure the input file has a column named `name` (or `nombre` for Spanish files) with molecule names, preferably in English for better API compatibility.

## Contributing
Contributions are welcome! Please open an issue or submit a pull request with improvements or bug fixes.

## License
This project is licensed under the CC BY_NC_SA 4.0 License. See the [LICENSE](LICENSE) file for details.