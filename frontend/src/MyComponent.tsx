import {
  Streamlit,
  StreamlitComponentBase,
  withStreamlitConnection,
} from "streamlit-component-lib"
import React, { ReactNode } from "react"

// Declaración global para el applet JSME
declare global {
  interface Window {
    JSApplet: any;
  }
}

// Interfaz para la estructura JSON de entrada/salida
interface SmilesRequest {
    smiles: string;
    id: string; // Identificador único de la solicitud
}

interface State {
  jsmeLoaded: boolean;
  status: string; // Estado de procesamiento interno
  currentRequestId: string | null; // El ID de la última solicitud procesada/en curso
}

class JSMEComponent extends StreamlitComponentBase<State> {
  private jsmeApplet: any = null;
  // SMILES que Streamlit usa para comunicar la entrada es ahora 'smiles_json'
  private currentInputJson: string = ''; 

  constructor(props: any) {
    super(props);
    this.state = { jsmeLoaded: false, status: 'Initial', currentRequestId: null }; 
  }

  componentDidMount(): void {
    this.loadJSME();
  }

  componentDidUpdate(): void {
    // El input de Python es props.args.smiles_json
    const nextInputJson = this.props.args?.smiles_json || null;
    let nextRequest: SmilesRequest | null = null;
    
    // 1. Intentar parsear el JSON de entrada (si es diferente)
    if (nextInputJson && nextInputJson !== this.currentInputJson) {
        try {
            nextRequest = JSON.parse(nextInputJson);
        } catch (e) {
            console.error("Failed to parse input JSON:", e);
            this.currentInputJson = nextInputJson; // Bloquea el procesamiento
            return; 
        }
    }
    
    // 2. Lógica de Disparo de Procesamiento
    if (this.state.jsmeLoaded && nextRequest && nextRequest.smiles.trim() !== '') {
        // Disparar procesamiento solo si el SMILES es nuevo y el ID es nuevo
        if (nextRequest.id !== this.state.currentRequestId) {
            this.currentInputJson = nextInputJson; // Actualizar tracker de JSON completo
            this.setState({ 
                status: 'Processing...',
                currentRequestId: nextRequest.id 
            }); 
            // Iniciar el procesamiento seguro inmediatamente
            this.safeProcessSmiles(nextRequest.smiles, nextRequest.id); 
        }
    }
    
    // 3. Lógica de Limpieza (Python envía null para resetear)
    // El componente se reinicia solo si el input es null O si el estado de requestID no coincide
    // con el ID de la última respuesta enviada (para evitar reintentos de IDs ya completados).
    if (nextInputJson === null) {
        // Enviar JSON vacío a Python. Esto forzará un RERUN de limpieza en Python.
        Streamlit.setComponentValue(JSON.stringify({})); 
        this.currentInputJson = ''; // Reset tracker de input
        this.setState({ status: 'Idle', currentRequestId: null }); 
    }
  }

  loadJSME = (): void => {
    if (window.JSApplet) {
      this.initJSME();
      return;
    }

    const script = document.createElement('script');
    script.src = 'https://jsme-editor.github.io/dist/jsme/jsme.nocache.js';
    script.async = true;
    script.onload = () => this.waitForJSME();
    document.head.appendChild(script);
  }

  waitForJSME = (): void => {
    if (window.JSApplet) {
      this.initJSME();
    } else {
      setTimeout(() => this.waitForJSME(), 100);
    }
  }

  initJSME = (): void => {
    // CRUCIAL: Crear el JSME en un DIV oculto para que no interfiera con la UI
    const hiddenDiv = document.createElement('div');
    hiddenDiv.style.position = 'absolute';
    hiddenDiv.style.left = '-9999px';
    hiddenDiv.style.width = '200px'; 
    hiddenDiv.style.height = '200px';
    hiddenDiv.id = 'jsme-processor';
    document.body.appendChild(hiddenDiv);

    this.jsmeApplet = new window.JSApplet.JSME("jsme-processor", "200px", "200px");
    
    this.setState({ jsmeLoaded: true, status: 'Ready' }); 
    
    Streamlit.setFrameHeight(10);
  }
  
  isJSMEAppletReady = (): boolean => {
    return this.jsmeApplet && typeof this.jsmeApplet.smiles === 'function';
  }

  safeProcessSmiles = (smiles: string, requestId: string, attempt: number = 0): void => {
    const MAX_ATTEMPTS = 20;
    
    if (attempt >= MAX_ATTEMPTS) {
        console.error("JSME Processing failed: Max retry attempts reached.");
        // Devolver un valor de error específico con el ID para rastreo
        const errorOutput = JSON.stringify({ smiles: `JSME_ERROR:TIMEOUT`, id: requestId });
        Streamlit.setComponentValue(errorOutput); 
        this.setState({ status: 'Failed (Timeout)' });
        return;
    }

    if (this.isJSMEAppletReady()) {
      try {
        this.jsmeApplet.readGenericMolecularInput(smiles);
        
        // Esperamos un tick del evento JSME (50ms es seguro)
        setTimeout(() => {
             const processedSmiles = this.jsmeApplet.smiles() || "JSME_ERROR:EMPTY_OUTPUT"; 
        
             // Devolver el resultado ESTUCTURADO con el ID
             const finalOutput = JSON.stringify({ 
                 smiles: processedSmiles,
                 id: requestId
             });
             Streamlit.setComponentValue(finalOutput);
             this.setState({ status: 'Processed' });
        }, 50); 
        
      } catch (error) {
        console.error("Error procesando SMILES en JSME:", error);
        const errorOutput = JSON.stringify({ smiles: `JSME_ERROR:EXCEPTION`, id: requestId });
        Streamlit.setComponentValue(errorOutput);
        this.setState({ status: 'Failed (Exception)' });
      }
    } else {
      // Reintentar después de 100ms
      console.log(`JSME not ready, retrying... Attempt: ${attempt + 1}`);
      setTimeout(() => this.safeProcessSmiles(smiles, requestId, attempt + 1), 100);
    }
  }

  public render = (): ReactNode => {
    const { jsmeLoaded, status } = this.state;

    // Renderizamos un texto minimalista para indicar que el procesador está listo o trabajando.
    return (
      <div style={{
        padding: '5px 10px',
        fontSize: '10px',
        color: '#666',
        backgroundColor: jsmeLoaded ? '#e0ffe0' : '#ffeee0',
        borderRadius: '2px',
        display: 'inline-block',
        minWidth: '200px'
      }}>
        {jsmeLoaded ? `✓ Procesador JSME listo. Estado: ${status}` : '⏳ Cargando procesador JSME...'}
      </div>
    );
  }
}

export default withStreamlitConnection(JSMEComponent)
