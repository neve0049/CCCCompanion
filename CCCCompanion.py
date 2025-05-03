import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from st_aggrid import AgGrid, GridOptionsBuilder

# Configuration des chemins des fichiers
EXCEL_PATH = "KDDB.xlsx"
DBDT_PATH = "DBDT.xlsx"

# Configuration de l'application
st.set_page_config(
    page_title="CCC Companion",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# =============================================
# Fonctions utilitaires
# =============================================

@st.cache_data
def load_excel_sheets(file_path):
    """Charge les noms de feuilles d'un fichier Excel avec gestion d'erreur"""
    try:
        return pd.ExcelFile(file_path).sheet_names
    except Exception as e:
        st.error(f"Erreur de chargement du fichier {file_path}: {str(e)}")
        return []

def configure_aggrid(df):
    """Configure AgGrid pour la sélection de lignes"""
    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_selection(
        selection_mode='single',
        use_checkbox=True,
        header_checkbox=False
    )
    gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
    return gb.build()

def create_ternary_plot(df_system, selected_number=None):
    """Crée le diagramme ternaire avec les points sélectionnés"""
    fig = go.Figure()
    
    # Ajout des traces pour les phases UP et LP
    for phase, color in [('UP', 'red'), ('LP', 'blue')]:
        fig.add_trace(go.Scatter(
            x=df_system[f'%Vol3 - {phase}'],
            y=df_system[f'%Vol2 - {phase}'],
            mode='markers',
            name=f'Phase {phase}',
            marker=dict(color=color, size=10),
            customdata=np.stack((
                df_system['Number'],
                df_system[f'%Vol1 - {phase}'],
                df_system[f'%Vol2 - {phase}'],
                df_system[f'%Vol3 - {phase}']
            ), axis=-1),
            hovertemplate=(
                "<b>Number</b>: %{customdata[0]}<br>"
                "<b>%Vol1</b>: %{customdata[1]:.2f}%<br>"
                "<b>%Vol2</b>: %{customdata[2]:.2f}%<br>"
                "<b>%Vol3</b>: %{customdata[3]:.2f}%<extra></extra>"
            )
        ))
    
    # Lignes de connexion
    for _, row in df_system.iterrows():
        fig.add_trace(go.Scatter(
            x=[row['%Vol3 - UP'], row['%Vol3 - LP']],
            y=[row['%Vol2 - UP'], row['%Vol2 - LP']],
            mode='lines',
            line=dict(color='gray', width=1, dash='dot'),
            showlegend=False,
            hoverinfo='none'
        ))
    
    # Ligne x + y = 1
    x_vals = np.linspace(0, 1, 100)
    fig.add_trace(go.Scatter(
        x=x_vals,
        y=1 - x_vals,
        mode='lines',
        name='x + y = 1',
        line=dict(color='green', width=2),
        hoverinfo='none'
    ))
    
    # Mise en évidence des points sélectionnés si un numéro est spécifié
    if selected_number is not None:
        df_selected = df_system[df_system['Number'] == selected_number]
        if not df_selected.empty:
            for phase in ['UP', 'LP']:
                fig.add_trace(go.Scatter(
                    x=[df_selected[f'%Vol3 - {phase}'].values[0]],
                    y=[df_selected[f'%Vol2 - {phase}'].values[0]],
                    mode='markers',
                    name=f'Selected {phase}',
                    marker=dict(
                        color='black',
                        size=16,
                        symbol='circle-open',
                        line=dict(width=2)
                    ),
                    hoverinfo='none'
                ))
    
    # Mise en forme du graphique
    fig.update_layout(
        title="Ternary Phase Diagram",
        xaxis_title='%Vol3',
        yaxis_title='%Vol2',
        showlegend=True,
        height=600,
        margin=dict(l=60, r=60, t=80, b=60)
    )
    
    return fig

# =============================================
# Pages de l'application
# =============================================

def show_home_page():
    st.title("CCC Companion")
    col1, col2 = st.columns(2)
    with col1:
        if st.button("KD Database Explorer", use_container_width=True):
            st.session_state.current_page = "kddb"
    with col2:
        if st.button("Ternary Phase Diagrams", use_container_width=True):
            st.session_state.current_page = "dbdt"

def show_kddb_page():
    st.title("KD Database Explorer")
    
    # Chargement des données KDDB
    sheet_names = load_excel_sheets(EXCEL_PATH)
    if not sheet_names:
        return
        
    selected_sheet = st.selectbox("Sélectionnez une feuille", sheet_names)
    
    try:
        df = pd.read_excel(EXCEL_PATH, sheet_name=selected_sheet)
        
        # Vérification des colonnes requises
        required_cols = ['Compound', 'SMILES', 'Number', 'System']
        if not all(col in df.columns for col in required_cols):
            st.error("Colonnes requises manquantes")
            return
            
        # Configuration AgGrid
        grid_options = configure_aggrid(df[required_cols])
        grid_response = AgGrid(
            df[required_cols],
            gridOptions=grid_options,
            height=400,
            update_mode='MODEL_CHANGED'
        )
        
        selected_rows = grid_response['selected_rows']
        
        if selected_rows:
            selected_row = selected_rows[0]
            selected_number = selected_row['Number']
            system_name = selected_row['System']
            
            # Chargement des données du système
            try:
                df_system = pd.read_excel(DBDT_PATH, sheet_name=system_name)
                
                # Création du graphique
                fig = create_ternary_plot(df_system, selected_number)
                
                # Affichage
                col1, col2 = st.columns([0.7, 0.3])
                
                with col1:
                    st.plotly_chart(fig, use_container_width=True)
                
                with col2:
                    st.subheader("Compositions")
                    df_selected = df_system[df_system['Number'] == selected_number]
                    if not df_selected.empty:
                        for phase in ['UP', 'LP']:
                            with st.expander(f"Phase {phase}"):
                                st.write(f"%Vol1: {df_selected[f'%Vol1 - {phase}'].values[0]:.2f}%")
                                st.write(f"%Vol2: {df_selected[f'%Vol2 - {phase}'].values[0]:.2f}%")
                                st.write(f"%Vol3: {df_selected[f'%Vol3 - {phase}'].values[0]:.2f}%")
                    else:
                        st.warning("Données non trouvées")
            
            except Exception as e:
                st.error(f"Erreur système: {str(e)}")
    
    except Exception as e:
        st.error(f"Erreur données: {str(e)}")
    
    if st.button("Retour"):
        st.session_state.current_page = "home"

def show_dbdt_page():
    st.title("Ternary Phase Diagrams")
    
    # Chargement des systèmes disponibles
    sheet_names = load_excel_sheets(DBDT_PATH)
    if not sheet_names:
        return
        
    selected_system = st.selectbox("Sélectionnez un système", sheet_names)
    
    try:
        df = pd.read_excel(DBDT_PATH, sheet_name=selected_system)
        lib_row = pd.read_excel(DBDT_PATH, sheet_name=selected_system, header=None).iloc[1]
        labels = {
            'vol1': lib_row[0],
            'vol2': lib_row[1],
            'vol3': lib_row[2]
        }
        
        # Nettoyage des données
        df = df.dropna(subset=['%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', 
                             '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP'])
        
        # Création du graphique interactif
        fig = go.Figure()
        
        # Ajout des traces
        for phase, color in [('UP', 'red'), ('LP', 'blue')]:
            fig.add_trace(go.Scatter(
                x=df[f'%Vol3 - {phase}'],
                y=df[f'%Vol2 - {phase}'],
                mode='markers',
                name=phase,
                marker=dict(color=color, size=10),
                customdata=df[['Number', 
                             f'%Vol1 - {phase}', 
                             f'%Vol2 - {phase}', 
                             f'%Vol3 - {phase}']].values,
                hovertemplate=(
                    f"<b>Number</b>: %{{customdata[0]}}<br>"
                    f"{labels['vol1']}: %{{customdata[1]:.2f}}%<br>"
                    f"{labels['vol2']}: %{{customdata[2]:.2f}}%<br>"
                    f"{labels['vol3']}: %{{customdata[3]:.2f}}%<extra></extra>"
                )
            ))
        
        # Lignes de connexion
        for _, row in df.iterrows():
            fig.add_trace(go.Scatter(
                x=[row['%Vol3 - UP'], row['%Vol3 - LP']],
                y=[row['%Vol2 - UP'], row['%Vol2 - LP']],
                mode='lines',
                line=dict(color='gray', width=1),
                showlegend=False,
                hoverinfo='none'
            ))
        
        # Mise en forme
        fig.update_layout(
            title=f"{labels['vol1']}/{labels['vol2']}/{labels['vol3']}",
            xaxis_title=f"% {labels['vol3']}",
            yaxis_title=f"% {labels['vol2']}",
            height=700
        )
        
        # Affichage
        col1, col2 = st.columns([0.7, 0.3])
        
        with col1:
            selected_point = st.plotly_chart(
                fig, 
                use_container_width=True,
                on_select="rerun",
                key="ternary_chart"
            )
            
            st.dataframe(df, use_container_width=True)
        
        with col2:
            st.subheader("Compositions")
            
            if st.session_state.get("ternary_chart_select"):
                point_data = st.session_state.ternary_chart_select["points"][0]
                point_number = point_data["customdata"][0]
                
                selected_row = df[df['Number'] == point_number].iloc[0]
                
                for phase in ['UP', 'LP']:
                    with st.expander(f"Phase {phase}"):
                        st.write(f"{labels['vol1']}: {selected_row[f'%Vol1 - {phase}']:.2f}%")
                        st.write(f"{labels['vol2']}: {selected_row[f'%Vol2 - {phase}']:.2f}%")
                        st.write(f"{labels['vol3']}: {selected_row[f'%Vol3 - {phase}']:.2f}%")
            else:
                st.info("Cliquez sur un point pour voir ses compositions")
    
    except Exception as e:
        st.error(f"Erreur: {str(e)}")
    
    if st.button("Retour"):
        st.session_state.current_page = "home"

# =============================================
# Application principale
# =============================================

def main():
    if 'current_page' not in st.session_state:
        st.session_state.current_page = "home"
    
    # Navigation
    if st.session_state.current_page == "home":
        show_home_page()
    elif st.session_state.current_page == "kddb":
        show_kddb_page()
    elif st.session_state.current_page == "dbdt":
        show_dbdt_page()

if __name__ == "__main__":
    main()