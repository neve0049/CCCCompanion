import streamlit as st
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from sklearn.preprocessing import StandardScaler, LabelEncoder
from rdkit import Chem
from rdkit.Chem import AllChem
import joblib
import requests
import io
import tempfile
import os
from datetime import datetime

# Configuration de la page
st.set_page_config(
    page_title="🧪 KD Prediction App",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

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
    
    def load_model_from_github(self, base_url):
        """Charge le modèle depuis les fichiers GitHub"""
        try:
            # URLs des fichiers sur GitHub
            model_url = f"{base_url}_model.h5"
            preprocessors_url = f"{base_url}_preprocessors.pkl"
            combinations_url = f"{base_url}_combinations.pkl"
            
            # Télécharger et charger le modèle
            with tempfile.NamedTemporaryFile(delete=False, suffix='.h5') as tmp_model:
                response = requests.get(model_url)
                response.raise_for_status()
                tmp_model.write(response.content)
                tmp_model.flush()
                
                custom_objects = {
                    'mse': tf.keras.losses.MeanSquaredError(),
                    'mae': tf.keras.losses.MeanAbsoluteError(),
                }
                
                self.model = tf.keras.models.load_model(
                    tmp_model.name, 
                    custom_objects=custom_objects,
                    compile=False
                )
                self.model.compile(optimizer='adam', loss='mse', metrics=['mae'])
            
            # Télécharger et charger les préprocesseurs
            response = requests.get(preprocessors_url)
            response.raise_for_status()
            preprocessors = joblib.load(io.BytesIO(response.content))
            
            self.solvent_encoder = preprocessors['solvent_encoder']
            self.composition_encoder = preprocessors['composition_encoder']
            self.kd_scaler = preprocessors['kd_scaler']
            self.fingerprint_bits = preprocessors['fingerprint_bits']
            self.fingerprint_radius = preprocessors['fingerprint_radius']
            self.solvent_composition_map = preprocessors['solvent_composition_map']
            
            # Télécharger et charger les combinaisons valides
            response = requests.get(combinations_url)
            response.raise_for_status()
            self.valid_combinations = joblib.load(io.BytesIO(response.content))
            
            self.is_trained = True
            return True
            
        except Exception as e:
            st.error(f"Erreur lors du chargement du modèle: {str(e)}")
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
    # Titre principal
    st.title("🧪 Prediction of Partitioning Coefficient (log KD)")
    st.markdown("---")
    
    # Initialisation du prédicteur
    if 'predictor' not in st.session_state:
        st.session_state.predictor = KDPredictor()
        st.session_state.model_loaded = False
    
    # Sidebar pour la configuration
    with st.sidebar:
        st.header("⚙️ Configuration")
        
        # URL GitHub pour les fichiers du modèle
        github_url = st.text_input(
            "URL de base GitHub pour les fichiers du modèle:",
            value="https://raw.githubusercontent.com/tonusername/tonrepo/main/models/kd_predictor_model",
            help="URL vers les fichiers .h5, .pkl sans l'extension"
        )
        
        if st.button("🔄 Charger le modèle depuis GitHub"):
            with st.spinner("Chargement du modèle..."):
                if st.session_state.predictor.load_model_from_github(github_url):
                    st.session_state.model_loaded = True
                    st.success("✅ Modèle chargé avec succès!")
                else:
                    st.session_state.model_loaded = False
                    st.error("❌ Erreur lors du chargement du modèle")
    
    # Section principale
    col1, col2 = st.columns([1, 1])
    
    with col1:
        st.header("🔬 Paramètres du système")
        
        # Sélection du solvant
        if st.session_state.model_loaded:
            solvents = st.session_state.predictor.get_available_solvents()
            selected_solvent = st.selectbox(
                "1. Sélectionnez un système de solvant biphasique:",
                options=solvents,
                index=0 if solvents else None
            )
            
            # Sélection de la composition
            if selected_solvent:
                compositions = st.session_state.predictor.get_available_compositions_for_solvent(selected_solvent)
                selected_composition = st.selectbox(
                    "2. Sélectionnez une composition:",
                    options=compositions,
                    index=0 if compositions else None
                )
            else:
                selected_composition = None
        else:
            selected_solvent = None
            selected_composition = None
            st.warning("Veuillez d'abord charger le modèle dans la sidebar")
        
        # Saisie du SMILES
        smiles_input = st.text_input(
            "3. Entrez le SMILES:",
            placeholder="Ex: CCO pour l'éthanol",
            help="Entrez une structure chimique au format SMILES"
        )
        
        # Validation du SMILES
        if smiles_input:
            mol = Chem.MolFromSmiles(smiles_input)
            if mol is None:
                st.error("❌ SMILES invalide")
                smiles_valid = False
            else:
                st.success("✅ SMILES valide")
                smiles_valid = True
        else:
            smiles_valid = False
    
    with col2:
        st.header("🎯 Actions")
        
        # Bouton de prédiction simple
        if st.button(
            "🎯 Prédire KD pour le système sélectionné",
            disabled=not (st.session_state.model_loaded and selected_solvent and selected_composition and smiles_valid),
            use_container_width=True
        ):
            with st.spinner("Calcul en cours..."):
                prediction = st.session_state.predictor.predict(
                    smiles_input, selected_solvent, selected_composition
                )
                
                if prediction is not None:
                    # Affichage du résultat
                    st.success("✅ Prédiction terminée!")
                    
                    # Interprétation
                    if prediction < -1:
                        interpretation = "Affinité avec la phase aqueuse"
                        color = "red"
                    elif prediction < 1:
                        interpretation = "Partitionnement optimal"
                        color = "orange"
                    else:
                        interpretation = "Affinité avec la phase organique"
                        color = "green"
                    
                    # Résultats
                    st.subheader("📊 Résultats de la prédiction")
                    st.metric("log KD prédit", f"{prediction:.4f}")
                    
                    col_a, col_b = st.columns(2)
                    with col_a:
                        st.info(f"**Interprétation:** {interpretation}")
                    with col_b:
                        st.info(f"**Couleur:** :{color}[{color}]")
                    
                    # Détails
                    with st.expander("📋 Détails de la prédiction"):
                        st.write(f"**SMILES:** {smiles_input}")
                        st.write(f"**Système de solvant:** {selected_solvent}")
                        st.write(f"**Composition:** {selected_composition}")
                        st.write(f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
                        
                        st.markdown("""
                        **📈 Plage de KD:**
                        - < -1.0 : Affinité avec la phase aqueuse
                        - -1.0 - 1.0 : Bon partitionnement
                        - > 1.0 : Affinité avec la phase organique
                        """)
                else:
                    st.error("❌ Erreur lors de la prédiction")
        
        # Bouton de scan complet
        if st.button(
            "🔍 Rechercher les systèmes optimaux (-1 < KD < 1)",
            disabled=not (st.session_state.model_loaded and smiles_valid),
            use_container_width=True
        ):
            with st.spinner("Recherche en cours... Cela peut prendre quelques minutes"):
                results = []
                total_combinations = 0
                valid_combinations = 0
                
                solvents = st.session_state.predictor.get_available_solvents()
                progress_bar = st.progress(0)
                
                for i, solvent in enumerate(solvents):
                    compositions = st.session_state.predictor.get_available_compositions_for_solvent(solvent)
                    total_combinations += len(compositions)
                    
                    for composition in compositions:
                        prediction = st.session_state.predictor.predict(
                            smiles_input, solvent, composition
                        )
                        
                        if prediction is not None and -1 <= prediction <= 1:
                            results.append({
                                'solvent': solvent,
                                'composition': composition,
                                'kd': prediction
                            })
                            valid_combinations += 1
                    
                    progress_bar.progress((i + 1) / len(solvents))
                
                # Affichage des résultats du scan
                st.subheader("🔍 Résultats de la recherche")
                st.write(f"**Compositions testées:** {total_combinations}")
                st.write(f"**Compositions avec -1 < KD < 1:** {valid_combinations}")
                
                if results:
                    # Trier les résultats
                    results.sort(key=lambda x: x['kd'])
                    
                    # Créer un DataFrame pour l'affichage
                    df_results = pd.DataFrame(results)
                    df_results['KD'] = df_results['kd'].round(4)
                    df_results = df_results[['solvent', 'composition', 'KD']]
                    
                    st.success(f"✅ {len(results)} systèmes optimaux trouvés!")
                    
                    # Affichage sous forme de tableau
                    st.dataframe(
                        df_results,
                        use_container_width=True,
                        hide_index=True
                    )
                    
                    # Téléchargement des résultats
                    csv = df_results.to_csv(index=False)
                    st.download_button(
                        label="📥 Télécharger les résultats en CSV",
                        data=csv,
                        file_name=f"kd_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                        mime="text/csv"
                    )
                    
                    # Recommandations
                    st.subheader("🎯 Recommandations")
                    col_rec1, col_rec2 = st.columns(2)
                    
                    with col_rec1:
                        st.info(
                            f"**Système le plus bas:**\n"
                            f"{results[0]['solvent']} + {results[0]['composition']}\n"
                            f"KD = {results[0]['kd']:.4f}"
                        )
                    
                    with col_rec2:
                        st.info(
                            f"**Système le plus haut:**\n"
                            f"{results[-1]['solvent']} + {results[-1]['composition']}\n"
                            f"KD = {results[-1]['kd']:.4f}"
                        )
                    
                else:
                    st.warning("⚠️ Aucun système optimal trouvé")
                    st.info("""
                    **Suggestions:**
                    - Essayez un autre composé
                    - Sélectionnez manuellement des compositions
                    - Vérifiez le SMILES
                    """)
        
        # Bouton de réinitialisation
        if st.button("🔄 Réinitialiser", use_container_width=True):
            st.rerun()
    
    # Section d'information
    st.markdown("---")
    st.header("📚 Informations")
    
    col_info1, col_info2, col_info3 = st.columns(3)
    
    with col_info1:
        st.subheader("ℹ️ À propos")
        st.markdown("""
        Cette application prédit le coefficient de partage (log KD) 
        pour des systèmes de solvants biphasiques utilisés en CPC.
        """)
    
    with col_info2:
        st.subheader("🎯 Plage optimale")
        st.markdown("""
        **-1 < log KD < 1**
        
        Cette plage indique un bon partitionnement pour la
        Chromatographie par Centrifugation Partition (CPC).
        """)
    
    with col_info3:
        st.subheader("🔧 Utilisation")
        st.markdown("""
        1. Chargez le modèle depuis GitHub
        2. Sélectionnez le système de solvant
        3. Entrez un SMILES valide
        4. Lancez la prédiction ou la recherche
        """)

if __name__ == "__main__":
    main()
