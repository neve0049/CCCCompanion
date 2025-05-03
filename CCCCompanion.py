import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from io import BytesIO
import json

# Configuration de l'application
st.set_page_config(
    page_title="CCC Companion",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Chemins des fichiers
EXCEL_PATH = "KDDB.xlsx"
DBDT_PATH = "DBDT.xlsx"

# Configuration des colonnes par système
SYSTEM_CONFIG = {
    "CPME Ethanol Water": {
        "required_columns": {
            "number": "Number",
            "vol1_up": "%Vol1 - UP",
            "vol2_up": "%Vol2 - UP",
            "vol3_up": "%Vol3 - UP",
            "vol1_lp": "%Vol1 - LP",
            "vol2_lp": "%Vol2 - LP",
            "vol3_lp": "%Vol3 - LP"
        },
        "labels": {
            "vol1": "CPME",
            "vol2": "Ethanol",
            "vol3": "Water"
        }
    },
    # Ajoutez d'autres systèmes ici au même format
}

# =============================================
# Fonctions utilitaires améliorées
# =============================================

@st.cache_data
def load_excel_sheets(file_path):
    """Charge les noms de feuilles d'un fichier Excel avec gestion d'erreur"""
    try:
        return pd.ExcelFile(file_path).sheet_names
    except Exception as e:
        st.error(f"Erreur de chargement du fichier {file_path}: {str(e)}")
        return []

@st.cache_data
def load_system_data(system_name):
    """Charge les données d'un système avec validation des colonnes"""
    try:
        df = pd.read_excel(DBDT_PATH, sheet_name=system_name)
        
        # Vérification des colonnes requises
        if system_name in SYSTEM_CONFIG:
            required = SYSTEM_CONFIG[system_name]["required_columns"]
            missing = [col for col in required.values() if col not in df.columns]
            
            if missing:
                st.error(f"Colonnes manquantes dans {system_name}: {', '.join(missing)}")
                available = [col for col in df.columns if not col.startswith('Unnamed')]
                st.info(f"Colonnes disponibles: {', '.join(available)}")
                return None
                
        return df.dropna(subset=[col for col in df.columns if '%Vol' in col])
    except Exception as e:
        st.error(f"Erreur lors du chargement de {system_name}: {str(e)}")
        return None

def create_ternary_plot(df_system, selected_row, system_name):
    """Crée le diagramme ternaire avec gestion des configurations système"""
    if system_name not in SYSTEM_CONFIG:
        st.error(f"Aucune configuration trouvée pour le système {system_name}")
        return None
    
    config = SYSTEM_CONFIG[system_name]
    cols = config["required_columns"]
    labels = config["labels"]
    
    try:
        fig = go.Figure()

        # Traces pour les phases UP et LP
        for phase, color in [('UP', 'red'), ('LP', 'blue')]:
            fig.add_trace(go.Scatter(
                x=df_system[cols[f'vol3_{phase.lower()}']],
                y=df_system[cols[f'vol2_{phase.lower()}']],
                mode='markers',
                name=f'Phase {phase}',
                marker=dict(color=color, size=10),
                customdata=np.stack((
                    df_system[cols['number']],
                    df_system[cols[f'vol1_{phase.lower()}']],
                    df_system[cols[f'vol2_{phase.lower()}']],
                    df_system[cols[f'vol3_{phase.lower()}']]
                ), axis=-1),
                hovertemplate=(
                    f"<b>Number</b>: %{{customdata[0]}}<br>"
                    f"{labels['vol1']}: %{{customdata[1]:.2f}}%<br>"
                    f"{labels['vol2']}: %{{customdata[2]:.2f}}%<br>"
                    f"{labels['vol3']}: %{{customdata[3]:.2f}}%<extra></extra>"
                )
            ))

        # Lignes de connexion
        for _, row in df_system.iterrows():
            fig.add_trace(go.Scatter(
                x=[row[cols['vol3_up']], row[cols['vol3_lp']]],
                y=[row[cols['vol2_up']], row[cols['vol2_lp']]],
                mode='lines',
                line=dict(color='gray', width=1, dash='dot'),
                showlegend=False,
                hoverinfo='none'
            ))

        # Mise en page
        fig.update_layout(
            title=f"Diagramme de Phase Ternaire - {system_name}",
            xaxis_title=f"% {labels['vol3']}",
            yaxis_title=f"% {labels['vol2']}",
            showlegend=True,
            height=600,
            margin=dict(l=60, r=60, t=80, b=60)
        )
        
        return fig
        
    except Exception as e:
        st.error(f"Erreur lors de la création du graphique: {str(e)}")
        return None

# =============================================
# Pages de l'application
# =============================================

def show_home_page():
    st.title("CCC Companion")
    col1, col2 = st.columns(2)
    with col1:
        if st.button("KD Database", use_container_width=True):
            st.session_state.current_page = "kddb"
    with col2:
        if st.button("Ternary Phase Diagrams", use_container_width=True):
            st.session_state.current_page = "dbdt"

def show_kddb_page():
    st.title("KD Database Explorer")
    
    # Chargement des données KDDB
    kddb_sheets = load_excel_sheets(EXCEL_PATH)
    if not kddb_sheets:
        return
        
    selected_sheet = st.selectbox("Sélectionnez une feuille", kddb_sheets)
    
    try:
        df = pd.read_excel(EXCEL_PATH, sheet_name=selected_sheet)
        
        # Vérification des colonnes requises
        required_cols = ['Compound', 'SMILES', 'Number', 'System']
        missing = [col for col in required_cols if col not in df.columns]
        
        if missing:
            st.error(f"Colonnes requises manquantes: {', '.join(missing)}")
            return
            
        # Affichage des données
        st.dataframe(df[required_cols], use_container_width=True)
        
        # Sélection d'un composé
        selected_number = st.selectbox("Sélectionnez un numéro", df['Number'].unique())
        selected_row = df[df['Number'] == selected_number].iloc[0]
        system_name = selected_row['System']
        
        # Chargement des données du système
        df_system = load_system_data(system_name)
        if df_system is None:
            return
            
        # Création de la visualisation
        col1, col2 = st.columns([0.7, 0.3])
        
        with col1:
            fig = create_ternary_plot(df_system, selected_row, system_name)
            if fig:
                st.plotly_chart(fig, use_container_width=True)
        
        with col2:
            if system_name in SYSTEM_CONFIG:
                config = SYSTEM_CONFIG[system_name]
                cols = config["required_columns"]
                labels = config["labels"]
                
                st.subheader("Compositions")
                for phase in ['UP', 'LP']:
                    with st.expander(f"Phase {phase}", expanded=True):
                        st.markdown(f"""
                        **{labels['vol1']}:** {selected_row[cols[f'vol1_{phase.lower()}']]:.2f}%  
                        **{labels['vol2']}:** {selected_row[cols[f'vol2_{phase.lower()}']]:.2f}%  
                        **{labels['vol3']}:** {selected_row[cols[f'vol3_{phase.lower()}']]:.2f}%
                        """)
    
    except Exception as e:
        st.error(f"Erreur lors du traitement des données: {str(e)}")

def show_dbdt_page():
    st.title("Ternary Phase Diagram Database")
    
    # Chargement des systèmes disponibles
    dbdt_sheets = load_excel_sheets(DBDT_PATH)
    if not dbdt_sheets:
        return
        
    selected_system = st.selectbox("Sélectionnez un système", dbdt_sheets)
    
    # Chargement des données
    df_system = load_system_data(selected_system)
    if df_system is None:
        return
        
    # Configuration du système
    system_config = SYSTEM_CONFIG.get(selected_system, {})
    if not system_config:
        st.warning(f"Aucune configuration trouvée pour {selected_system}")
        return
        
    # Visualisation
    col1, col2 = st.columns([0.7, 0.3])
    
    with col1:
        fig = create_ternary_plot(df_system, None, selected_system)
        if fig:
            st.plotly_chart(fig, use_container_width=True)
            
    with col2:
        st.subheader("Sélection de point")
        selected_number = st.selectbox(
            "Choisir un numéro", 
            df_system[system_config["required_columns"]["number"]].unique()
        )
        
        selected_row = df_system[
            df_system[system_config["required_columns"]["number"]] == selected_number
        ].iloc[0]
        
        st.subheader("Compositions")
        for phase in ['UP', 'LP']:
            with st.expander(f"Phase {phase}", expanded=True):
                st.markdown(f"""
                **{system_config['labels']['vol1']}:** {selected_row[system_config['required_columns'][f'vol1_{phase.lower()}']]:.2f}%  
                **{system_config['labels']['vol2']}:** {selected_row[system_config['required_columns'][f'vol2_{phase.lower()}']]:.2f}%  
                **{system_config['labels']['vol3']}:** {selected_row[system_config['required_columns'][f'vol3_{phase.lower()}']]:.2f}%
                """)

# =============================================
# Application principale
# =============================================

def main():
    # Initialisation de l'état
    if 'current_page' not in st.session_state:
        st.session_state.current_page = "home"
    
    # Navigation
    with st.sidebar:
        st.title("Navigation")
        if st.button("Accueil"):
            st.session_state.current_page = "home"
        if st.button("KD Database"):
            st.session_state.current_page = "kddb"
        if st.button("Ternary Diagrams"):
            st.session_state.current_page = "dbdt"
    
    # Router vers la page active
    if st.session_state.current_page == "home":
        show_home_page()
    elif st.session_state.current_page == "kddb":
        show_kddb_page()
    elif st.session_state.current_page == "dbdt":
        show_dbdt_page()

if __name__ == "__main__":
    main()