# 🧪 Moodle Chemistry Question Generator (SMILES/JSME Standardization)

This Streamlit application is designed to streamline the process of creating Moodle quiz questions of the pmatchjme type. These questions require students to draw a chemical structure using the integrated JSME molecular editor. The application automatically handles the heavy lifting: looking up the canonical SMILES for common compounds via the NCI Chemical Identifier Resolver (CIR) API and standardizing the structure using a hidden JSME component before exporting the final Moodle XML file.

## ✨ Features

- Multilingual Support: Interface available in English (EN) and Spanish (ES).

- SMILES Lookup: Automatically searches the NCI Chemical Identifier Resolver (CIR) API to find the SMILES string for a given common molecule name.

- JSME Standardization: Uses the integrated jsme_editor component to convert the initial SMILES string into a standardized, canonical format that works reliably as the Moodle answer key. This ensures robustness and consistency.

- Individual Entry: Easily look up and process one molecule at a time.

- Bulk Processing: Upload a CSV or Excel file containing a list of molecule names for fast, automated generation of multiple questions.

- Manual Entry: Bypass the API lookup by manually entering a SMILES string and a question name.

- Moodle XML Export: Generates a ready-to-import Moodle XML file containing all the standardized structure-drawing questions (pmatchjme).

## 🛠️ How It Works

The core functionality relies on a three-step chemical data pipeline:

1. Input & Lookup (NCI CIR): The application takes a molecule name and queries the NCI CIR API to retrieve the initial, potentially non-canonical, SMILES string.

2. RDKit Canonicalization: The retrieved or manually entered SMILES is first cleaned and canonicalized using RDKit.

3. Final Standardization (JSME Component): The application sends the RDKit-canonicalized SMILES to a non-visible JSME component. This component outputs the final, highly standardized SMILES that is proven to work correctly with Moodle's embedded JSME structure validation.

## 🚀 How to Run the Application

There are two primary ways to access and use this tool:

### Option 1: Use the Public Web Application (Recommended)

The application is deployed publicly and can be accessed directly through this link:

👉 https://github.com/jsme-editor/jsme-editor.github.io

### Option 2: Run Locally (Requires Python and Node.js)

To run the application in local development mode, you must run two processes simultaneously in separate terminals: the Streamlit server (Python) and the frontend component development server (Node/npm).

1. Clone the Repository and Install Python Dependencies:

git clone https://github.com/cfmarcos64/Create-pattern-match-with-molecular-editor-questions-for-Moodle

cd [repository-name]

// Install Python dependencies (skip if already done)
pip install -r requirements.txt


2. Run the Component Frontend (TERMINAL 1):

This step starts the component development server on http://localhost:3001. This is necessary for Streamlit to connect to the React component and see live changes.

// Navigate to the frontend directory
cd my_component/frontend
// Install JavaScript dependencies (only the first time)
npm install
// Start the component development server
npm run start


**Note:** Keep this terminal open and running while using the Streamlit application.

3. Run the Streamlit Application (TERMINAL 2):

Open a second terminal. Navigate back to the project root folder and run the application.

// Go back to the root directory
cd ../..
// Execute the main Streamlit application
streamlit run MoleculeToMoodleJSME.py


The Streamlit server will automatically connect to the component development server (Terminal 1).

## 📁 File Structure

The repository is organized into two main parts: the Streamlit Python application and the custom component frontend. The project directory should contain at least these files:

/ (repository root)

├── MoleculeToMoodleJSME.py  # The main Streamlit applicationThe main Streamlit application file. This script defines the user interface (UI), manages input/output data flow, handles the high-level logic for question generation, and calls the custom Python component function (jsme_processor_component) to trigger standardization.

├── requirements.txt         # Lists all necessary Python dependencies (e.g., streamlit, rdkit, etc.) required to run the main Streamlit application.

├──  __init__.py             # The Python Component Bridge. This file contains the Python wrapper function (jsme_processor_component) that links the Python application to the compiled React/TypeScript component files in the frontend/ directory.

├── frontend/                # The root directory for the Streamlit custom component's source code. This includes all necessary files for building the React/TypeScript application that hosts the JSME logic.

│    └── index.html         # The entry point HTML file for the custom component. It's the file Streamlit loads to display the component in the browser. It typically loads the bundled JavaScript assets.

│    └── vite.config.ts     # The Vite bundler configuration file. This dictates how the TypeScript and React files are compiled and optimized into the static JavaScript, HTML, and CSS assets used by the component.

│    └── src/               # Source directory for the core component files.

│    └──    └── MyComponent.tsx    # The Core React Component containing the primary client-side logic. This component is responsible for loading the hidden JSME applet, listening for SMILES input from Python, invoking the JSME standardization, and reliably returning the processed SMILES (with a unique request ID) back to the Streamlit backend.

│    └── node_modules/      # Directory containing all local JavaScript dependencies installed via npm. This directory is necessary for development and building but should generally be excluded from version control (e.g., via .gitignore).

└── README.md               # This file, providing an overview of the project and its structure.  




*This tool was created to assist educators and chemists in quickly generating high-quality Moodle quiz content.*

## License

This project is licensed under the CC BY_NC_SA 4.0 License. See the [LICENSE](LICENSE) file for details.
