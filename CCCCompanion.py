import streamlit as st
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from sklearn.preprocessing import StandardScaler, LabelEncoder
from rdkit import Chem
from rdkit.Chem import AllChem
import joblib
import base64
from datetime import datetime

class KDPredictor:
    def __init__(self, fingerprint_bits=2048, fingerprint_radius=2):
        self.fingerprint_bits = fingerprint_bits
        self.fingerprint_radius = fingerprint_radius
        self.solvent_encoder = LabelEncoder()
        self.composition_encoder = LabelEncoder()
        self.kd_scaler = StandardScaler()
        self.model = None
        self.is_trained = False
        self.valid_combinations = {}
        self.solvent_composition_map = {}
        
    def smiles_to_fingerprint(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            fingerprint = AllChem.GetMorganFingerprintAsBitVect(
                mol, self.fingerprint_radius, nBits=self.fingerprint_bits
            )
            return np.array(fingerprint)
        except Exception:
            return None
    
    def load_model(self, filepath):
        try:
            custom_objects = {
                'mse': tf.keras.losses.MeanSquaredError(),
                'mae': tf.keras.losses.MeanAbsoluteError(),
            }
            
            self.model = tf.keras.models.load_model(
                f'{filepath}_model.h5', 
                custom_objects=custom_objects,
                compile=False
            )
            self.model.compile(optimizer='adam', loss='mse', metrics=['mae'])
            
            preprocessors = joblib.load(f'{filepath}_preprocessors.pkl')
            self.solvent_encoder = preprocessors['solvent_encoder']
            self.composition_encoder = preprocessors['composition_encoder']
            self.kd_scaler = preprocessors['kd_scaler']
            self.fingerprint_bits = preprocessors['fingerprint_bits']
            self.fingerprint_radius = preprocessors['fingerprint_radius']
            self.solvent_composition_map = preprocessors['solvent_composition_map']
            
            self.valid_combinations = joblib.load(f'{filepath}_combinations.pkl')
            self.is_trained = True
            return True
            
        except Exception as e:
            st.error(f"Error loading model: {str(e)}")
            return False
    
    def get_available_solvents(self):
        return list(self.valid_combinations.keys())
    
    def get_available_compositions_for_solvent(self, solvent):
        if solvent in self.valid_combinations:
            return self.valid_combinations[solvent]
        return []
    
    def predict(self, smiles, solvent_system, composition):
        if not self.is_trained:
            return None
        
        smiles_fp = self.smiles_to_fingerprint(smiles)
        if smiles_fp is None:
            return None
        
        if solvent_system not in self.solvent_encoder.classes_:
            return None
        
        if composition not in self.valid_combinations.get(solvent_system, []):
            return None
        
        smiles_fp = smiles_fp.reshape(1, -1)
        solvent_encoded = self.solvent_encoder.transform([solvent_system]).reshape(1, -1)
        composition_encoded = self.solvent_composition_map[solvent_system]['mapping'][composition]
        composition_encoded = np.array([composition_encoded]).reshape(1, -1)
        
        try:
            prediction_scaled = self.model.predict({
                'smiles': smiles_fp,
                'solvent': solvent_encoded,
                'composition': composition_encoded
            }, verbose=0)
            
            prediction_original = self.kd_scaler.inverse_transform(prediction_scaled.reshape(-1, 1))
            return prediction_original[0][0]
            
        except Exception:
            return None

def main():
    st.set_page_config(
        page_title="🧪 KD Prediction Web App",
        page_icon="🧪",
        layout="wide",
        initial_sidebar_state="expanded"
    )
    
    # CSS personnalisé
    st.markdown("""
    <style>
    .main-header {
        font-size: 2.5rem;
        color: #2c3e50;
        text-align: center;
        margin-bottom: 2rem;
    }
    .result-box {
        background-color: #f8f9fa;
        border-radius: 10px;
        padding: 20px;
        border-left: 5px solid #3498db;
        margin: 10px 0;
    }
    .success-box {
        border-left: 5px solid #27ae60;
    }
    .warning-box {
        border-left: 5px solid #f39c12;
    }
    .error-box {
        border-left: 5px solid #e74c3c;
    }
    </style>
    """, unsafe_allow_html=True)
    
    st.markdown('<h1 class="main-header">🧪 Prediction of Partitioning Coefficient (log KD)</h1>', unsafe_allow_html=True)
    
    # Initialisation du prédicteur
    if 'predictor' not in st.session_state:
        st.session_state.predictor = KDPredictor()
        st.session_state.model_loaded = False
        st.session_state.load_model()
    
    def load_model():
        with st.spinner('🔄 Loading AI model...'):
            success = st.session_state.predictor.load_model('kd_predictor_model')
            st.session_state.model_loaded = success
    
    # Sidebar
    with st.sidebar:
        st.header("ℹ️ About")
        st.info("""
        This application predicts the partitioning coefficient (log KD) 
        of molecules in biphasic solvent systems using deep learning.
        
        **Features:**
        - Single system prediction
        - Optimal system search
        - Real-time SMILES validation
        """)
        
        st.header("📊 Model Status")
        if st.session_state.model_loaded:
            st.success("✅ Model loaded successfully")
            solvents = st.session_state.predictor.get_available_solvents()
            st.write(f"**Available systems:** {len(solvents)}")
        else:
            st.error("❌ Model not loaded")
            if st.button("Retry loading"):
                load_model()
    
    # Main content
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.header("🔬 System Configuration")
        
        # SMILES Input
        smiles_input = st.text_input(
            "**Enter SMILES string:**",
            placeholder="CCO for ethanol, C1=CC=CC=C1 for benzene...",
            help="Enter a valid SMILES string for your molecule"
        )
        
        # SMILES validation
        if smiles_input:
            mol = Chem.MolFromSmiles(smiles_input)
            if mol is None:
                st.error("❌ Invalid SMILES string")
            else:
                st.success("✅ Valid SMILES")
                
                # Display molecule info
                from rdkit.Chem import Draw
                from rdkit.Chem import Descriptors
                
                img = Draw.MolToImage(mol, size=(300, 200))
                st.image(img, caption="Molecule Structure")
                
                # Basic descriptors
                mol_weight = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                st.write(f"**Molecular Weight:** {mol_weight:.2f}")
                st.write(f"**LogP:** {logp:.2f}")
        
        # Solvent system selection
        if st.session_state.model_loaded:
            solvents = st.session_state.predictor.get_available_solvents()
            selected_solvent = st.selectbox(
                "**Select biphasic solvent system:**",
                solvents,
                index=0 if solvents else None
            )
            
            # Composition selection
            if selected_solvent:
                compositions = st.session_state.predictor.get_available_compositions_for_solvent(selected_solvent)
                selected_composition = st.selectbox(
                    "**Select composition:**",
                    compositions,
                    index=0 if compositions else None
                )
    
    with col2:
        st.header("🎯 Prediction Actions")
        
        # Single prediction
        if st.session_state.model_loaded and smiles_input and selected_solvent and selected_composition:
            if st.button("🎯 Predict KD for Selected System", type="primary", use_container_width=True):
                with st.spinner("Calculating prediction..."):
                    prediction = st.session_state.predictor.predict(smiles_input, selected_solvent, selected_composition)
                    
                    if prediction is not None:
                        # Display results
                        if prediction < -1:
                            interpretation = "Affinity with aqueous phase"
                            box_class = "warning-box"
                            color = "#e74c3c"
                        elif prediction < 1:
                            interpretation = "Optimal partitioning"
                            box_class = "success-box"
                            color = "#27ae60"
                        else:
                            interpretation = "Affinity with organic phase"
                            box_class = "error-box"
                            color = "#f39c12"
                        
                        st.markdown(f"""
                        <div class="result-box {box_class}">
                            <h3>📊 Prediction Results</h3>
                            <p><strong>SMILES:</strong> {smiles_input}</p>
                            <p><strong>System:</strong> {selected_solvent} + {selected_composition}</p>
                            <p><strong>Predicted log KD:</strong> <span style="color: {color}; font-weight: bold">{prediction:.4f}</span></p>
                            <p><strong>Interpretation:</strong> {interpretation}</p>
                        </div>
                        """, unsafe_allow_html=True)
                    else:
                        st.error("❌ Prediction failed. Please check your inputs.")
        
        # Complete scan
        if st.session_state.model_loaded and smiles_input:
            if st.button("🔍 Find Optimal Systems (-1 < KD < 1)", use_container_width=True):
                with st.spinner("Scanning all systems... This may take a while."):
                    results = []
                    total_combinations = 0
                    
                    solvents = st.session_state.predictor.get_available_solvents()
                    progress_bar = st.progress(0)
                    
                    for i, solvent in enumerate(solvents):
                        compositions = st.session_state.predictor.get_available_compositions_for_solvent(solvent)
                        total_combinations += len(compositions)
                        
                        for composition in compositions:
                            prediction = st.session_state.predictor.predict(smiles_input, solvent, composition)
                            if prediction is not None and -1 <= prediction <= 1:
                                results.append({
                                    'solvent': solvent,
                                    'composition': composition,
                                    'kd': prediction
                                })
                        
                        progress_bar.progress((i + 1) / len(solvents))
                    
                    # Display scan results
                    if results:
                        results.sort(key=lambda x: x['kd'])
                        
                        st.success(f"🎉 Found {len(results)} optimal systems!")
                        
                        # Create results dataframe
                        df_results = pd.DataFrame(results)
                        st.dataframe(df_results, use_container_width=True)
                        
                        # Download button
                        csv = df_results.to_csv(index=False)
                        b64 = base64.b64encode(csv.encode()).decode()
                        href = f'<a href="data:file/csv;base64,{b64}" download="optimal_systems.csv">📥 Download Results as CSV</a>'
                        st.markdown(href, unsafe_allow_html=True)
                        
                    else:
                        st.warning("❌ No systems found with -1 < log KD < 1")
        
        # Model information
        st.header("📈 About log KD")
        st.markdown("""
        **Interpretation guide:**
        - **log KD < -1**: Strong affinity with aqueous phase
        - **-1 < log KD < 1**: Good partitioning (optimal for CPC)
        - **log KD > 1**: Strong affinity with organic phase
        
        **Optimal range for CPC:** -1 to 1
        """)

if __name__ == "__main__":
    main()
