import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from io import BytesIO

# Configuration de l'application
st.set_page_config(
    page_title="CCC Companion",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Chemins des fichiers (adaptés pour Streamlit Cloud)
EXCEL_PATH = "KDDB.xlsx"
DBDT_PATH = "DBDT.xlsx"

# =============================================
# Fonctions utilitaires
# =============================================

@st.cache_data
def load_kddb_sheets():
    """Charge les noms de feuilles du fichier KDDB"""
    try:
        return pd.ExcelFile(EXCEL_PATH).sheet_names
    except Exception as e:
        st.error(f"Erreur lors du chargement de KDDB.xlsx: {str(e)}")
        return []

@st.cache_data
def load_dbdt_sheets():
    """Charge les noms de feuilles du fichier DBDT"""
    try:
        return pd.ExcelFile(DBDT_PATH).sheet_names
    except Exception as e:
        st.error(f"Erreur lors du chargement de DBDT.xlsx: {str(e)}")
        return []

def create_ternary_plot(df_system, selected_row, labels):
    """Crée le diagramme ternaire avec Plotly"""
    fig = go.Figure()

    # Ajout des traces UP et LP
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

    # Points sélectionnés
    for phase in ['UP', 'LP']:
        fig.add_trace(go.Scatter(
            x=[selected_row[f'%Vol3 - {phase}']],
            y=[selected_row[f'%Vol2 - {phase}']],
            mode='markers',
            name=f'Sélection {phase}',
            marker=dict(
                color='black',
                size=16,
                symbol='circle-open',
                line=dict(width=2)
            ),
            hoverinfo='none'
        ))

    fig.update_layout(
        title=dict(
            text=f"Diagramme de Phase Ternaire - Système {selected_row['System']}",
            x=0.5,
            font=dict(size=16)
        ),
        xaxis=dict(
            title='%Vol3',
            range=[0, 1],
            constrain='domain'
        ),
        yaxis=dict(
            title='%Vol2',
            range=[0, 1],
            scaleanchor='x',
            scaleratio=1
        ),
        showlegend=True,
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=1.02,
            xanchor='right',
            x=1
        ),
        margin=dict(l=60, r=60, t=80, b=60, pad=20),
        height=600
    )
    
    return fig

# =============================================
# Interface principale
# =============================================

def main():
    # Initialisation de l'état de session
    if 'current_page' not in st.session_state:
        st.session_state.current_page = "home"
    if 'selected_sheet' not in st.session_state:
        st.session_state.selected_sheet = None
    if 'selected_row' not in st.session_state:
        st.session_state.selected_row = None

    # Barre latérale de navigation
    with st.sidebar:
        st.title("Navigation")
        if st.button("Accueil"):
            st.session_state.current_page = "home"
            st.rerun()
        if st.button("KD Database"):
            st.session_state.current_page = "kddb"
            st.rerun()
        if st.button("Ternary Phase Diagrams"):
            st.session_state.current_page = "dbdt"
            st.rerun()

    # Gestion des pages
    if st.session_state.current_page == "home":
        show_home_page()
    elif st.session_state.current_page == "kddb":
        show_kddb_page()
    elif st.session_state.current_page == "dbdt":
        show_dbdt_page()

def show_home_page():
    st.title("CCC Companion")
    st.markdown("""
    Bienvenue dans l'application CCC Companion. Sélectionnez une base de données à explorer:
    """)
    
    col1, col2 = st.columns(2)
    with col1:
        if st.button("KD Database", use_container_width=True, key="kddb_home"):
            st.session_state.current_page = "kddb"
            st.rerun()
    with col2:
        if st.button("Ternary Phase Diagrams", use_container_width=True, key="dbdt_home"):
            st.session_state.current_page = "dbdt"
            st.rerun()

def show_kddb_page():
    st.title("KD Database Explorer")
    
    # Chargement des données
    sheet_names = load_kddb_sheets()
    if not sheet_names:
        st.error("Aucune feuille trouvée dans le fichier KDDB.xlsx")
        return
    
    # Recherche et sélection de feuille
    search_query = st.text_input("Rechercher une feuille par nom")
    filtered_sheets = [s for s in sheet_names if search_query.lower() in s.lower()] if search_query else sheet_names
    
    if not filtered_sheets:
        st.warning("Aucune feuille ne correspond à votre recherche")
        return
    
    selected_sheet = st.selectbox("Sélectionnez une feuille", filtered_sheets, key="kddb_sheet_select")
    
    # Chargement des données de la feuille sélectionnée
    try:
        df = pd.read_excel(EXCEL_PATH, sheet_name=selected_sheet)
        required_cols = ['Compound', 'SMILES', 'Number', 'System']
        
        if not all(col in df.columns for col in required_cols):
            st.error("La feuille sélectionnée ne contient pas toutes les colonnes requises")
            return
            
        # Affichage des données
        st.data_editor(
            df,
            use_container_width=True,
            hide_index=True,
            column_config={
                "Number": st.column_config.NumberColumn(format="%d"),
            },
            disabled=df.columns
        )
        
        # Sélection d'une ligne pour visualisation
        selected_number = st.selectbox(
            "Sélectionnez un numéro pour visualisation",
            df['Number'].unique(),
            key="kddb_number_select"
        )
        
        selected_row = df[df['Number'] == selected_number].iloc[0]
        st.session_state.selected_row = selected_row
        
        # Chargement des données du système correspondant
        system_name = selected_row['System']
        try:
            df_system = pd.read_excel(DBDT_PATH, sheet_name=system_name)
            
            # Création du layout de visualisation
            col1, col2 = st.columns([0.7, 0.3])
            
            with col1:
                # Diagramme ternaire
                fig = create_ternary_plot(df_system, selected_row, {})
                st.plotly_chart(fig, use_container_width=True)
            
            with col2:
                # Informations sur les compositions
                st.subheader("Compositions des phases")
                
                for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                    with st.expander(f"Phase {phase}", expanded=True):
                        st.markdown(f"""
                        **%Vol1:** {selected_row[f'%Vol1 - {phase}']:.2f}%  
                        **%Vol2:** {selected_row[f'%Vol2 - {phase}']:.2f}%  
                        **%Vol3:** {selected_row[f'%Vol3 - {phase}']:.2f}%
                        """)
                        
        except Exception as e:
            st.error(f"Erreur lors du chargement du système {system_name}: {str(e)}")
            
    except Exception as e:
        st.error(f"Erreur lors du chargement des données: {str(e)}")

def show_dbdt_page():
    st.title("Ternary Phase Diagram Database")
    
    # Chargement des données
    sheet_names = load_dbdt_sheets()
    if not sheet_names:
        st.error("Aucune feuille trouvée dans le fichier DBDT.xlsx")
        return
    
    # Sélection de la feuille
    selected_sheet = st.selectbox("Sélectionnez un système", sheet_names, key="dbdt_sheet_select")
    
    try:
        # Lecture des données
        df = pd.read_excel(DBDT_PATH, sheet_name=selected_sheet)
        lib_row = pd.read_excel(DBDT_PATH, sheet_name=selected_sheet, header=None).iloc[1]
        labels = {
            'vol1': lib_row[0],
            'vol2': lib_row[1],
            'vol3': lib_row[2],
            'sheet': selected_sheet
        }
        
        # Nettoyage des données
        df = df.dropna(subset=['%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP'])
        df = df.sort_values('%Vol3 - UP')
        
        # Création du graphique
        fig = go.Figure()

        # UP Phase
        fig.add_trace(go.Scatter(
            x=df['%Vol3 - UP'],
            y=df['%Vol2 - UP'],
            mode='lines+markers',
            name='UP',
            line=dict(color='red'),
            marker=dict(size=8),
            customdata=df[['Number', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP']].values,
            hovertemplate=(
                f'Phase: UP<br>{labels["vol1"]}: %{{customdata[1]:.2f}}%<br>'
                f'{labels["vol2"]}: %{{customdata[2]:.2f}}%<br>'
                f'{labels["vol3"]}: %{{customdata[3]:.2f}}%<extra></extra>'
            )
        ))

        # LP Phase
        fig.add_trace(go.Scatter(
            x=df['%Vol3 - LP'],
            y=df['%Vol2 - LP'],
            mode='lines+markers',
            name='LP',
            line=dict(color='blue'),
            marker=dict(size=8),
            customdata=df[['Number', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP']].values,
            hovertemplate=(
                f'Phase: LP<br>{labels["vol1"]}: %{{customdata[1]:.2f}}%<br>'
                f'{labels["vol2"]}: %{{customdata[2]:.2f}}%<br>'
                f'{labels["vol3"]}: %{{customdata[3]:.2f}}%<extra></extra>'
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

        # Ligne x + y = 1
        fig.add_trace(go.Scatter(
            x=[0, 1], y=[1, 0],
            mode='lines',
            line=dict(color='black', dash='dash'),
            name='x + y = 1',
            hoverinfo='none'
        ))

        fig.update_layout(
            title=f"Ternary Phase Diagram of {labels['vol1']} / {labels['vol2']} / {labels['vol3']} - {labels['sheet']}",
            xaxis_title=f'% {labels["vol3"]}',
            yaxis_title=f'% {labels["vol2"]}',
            hovermode='closest',
            height=700
        )

        # Affichage
        col1, col2 = st.columns([0.7, 0.3])
        
        with col1:
            st.plotly_chart(fig, use_container_width=True)
            
            # Affichage des données brutes
            st.subheader("Données brutes")
            st.dataframe(df, use_container_width=True)
        
        with col2:
            # Sélection interactive
            st.subheader("Sélection de point")
            selected_number = st.selectbox(
                "Choisir un numéro", 
                df['Number'].unique(),
                key="dbdt_number_select"
            )
            
            selected_row = df[df['Number'] == selected_number].iloc[0]
            
            # Affichage des compositions
            st.subheader("Compositions des phases")
            
            for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                with st.expander(f"Phase {phase}", expanded=True):
                    st.markdown(f"""
                    **{labels['vol1']}:** {selected_row[f'%Vol1 - {phase}']:.2f}%  
                    **{labels['vol2']}:** {selected_row[f'%Vol2 - {phase}']:.2f}%  
                    **{labels['vol3']}:** {selected_row[f'%Vol3 - {phase}']:.2f}%
                    """)
    
    except Exception as e:
        st.error(f"Erreur lors du traitement des données: {str(e)}")

if __name__ == "__main__":
    main()