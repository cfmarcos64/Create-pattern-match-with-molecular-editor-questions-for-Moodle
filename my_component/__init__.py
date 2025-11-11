import os
import streamlit.components.v1 as components

_RELEASE = True

if not _RELEASE:
    _component_func = components.declare_component(
        "jsme_editor",
        url="http://localhost:3001",
    )
else:
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
