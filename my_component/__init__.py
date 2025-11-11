import os
import streamlit.components.v1 as components

# Determines if we are in development or deployment mode
# _RELEASE = True activates the deployment mode.
_RELEASE = True

if not _RELEASE:
    # If we are in development, we aim the Vite local de server (we we do not use it here, but the structure is kept)
    _component_func = components.declare_component(
        "jsme_editor",
        url="http://localhost:3001",
    )

else:
    # Deployment mode:
    parent_dir = os.path.dirname(os.path.abspath(__file__)) 
    build_dir = os.path.join(parent_dir, "frontend/build")
    _component_func = components.declare_component("jsme_editor", path=build_dir)


def jsme_editor(smiles_json, key=None):
    """
    Creates an instance of JSME molecular editor.
    
    Parameters
    ----------
    smiles_json: str
        The JSON string containing the SMILES and ID of the requirement
        for the asyncronous processing.
    key: str or None
        Unique key to identify the component
    
    Returns
    -------
    str
        The resulting SMILES code returned by the component React.
    """
    # CRUCIAL: The argument 'smiles_json' is passed as keyword argument to the component
    component_value = _component_func(smiles_json=smiles_json, key=key, default="")
    
    return component_value
