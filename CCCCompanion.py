import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import sys

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

def create_phase_display(row, labels):
    """Crée l'affichage des compositions des phases"""
    return {
        "UP": {
            "vol1": row['%Vol1 - UP'],
            "vol2": row['%Vol2 - UP'],
            "vol3": row['%Vol3 - UP']
        },
        "LP": {
            "vol1": row['%Vol1 - LP'],
            "vol2": row['%Vol2 - LP'],
            "vol3": row['%Vol3 - LP']
        },
        "labels": labels
    }

# =============================================
# Pages de l'application
# =============================================

def show_home_page():
    st.title("CCC Companion")
    st.markdown("""
    Bienvenue dans l'application CCC Companion. Sélectionnez une base de données à explorer:
    """)
    
    col1, col2 = st.columns(2)
    with col1:
        if st.button("KD Database Explorer", use_container_width=True, key="kddb_home"):
            st.session_state.current_page = "kddb"
            st.rerun()
    with col2:
        if st.button("Ternary Phase Diagrams", use_container_width=True, key="dbdt_home"):
            st.session_state.current_page = "dbdt"
            st.rerun()

def show_kddb_page():
    """Page KD Database Explorer - Version complète avec recherche"""
    st.title("KD Database Explorer")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(EXCEL_PATH)
    if not sheet_names:
        return
        
    # Zone de recherche
    col1, col2 = st.columns([0.7, 0.3])
    with col1:
        search_query = st.text_input(
            "Entrez un nom de molécule ou système...",
            key="search_input",
            placeholder="Rechercher"
        )
    with col2:
        st.write("")  # Pour l'alignement
        if st.button("Rechercher", key="search_button"):
            st.session_state.search_triggered = True

    # Gestion de la recherche
    if 'search_triggered' not in st.session_state:
        st.session_state.search_triggered = False

    if st.session_state.search_triggered or search_query:
        search_value = search_query.strip().lower()
        matching_sheets = [sheet for sheet in sheet_names if search_value in sheet.lower()]
        
        if not matching_sheets:
            st.warning("Aucune correspondance trouvée.")
        else:
            # Affichage des résultats
            selected_sheet = st.radio(
                "Feuilles correspondantes:",
                matching_sheets,
                key="sheet_selection"
            )
            
            # Chargement des données de la feuille sélectionnée
            try:
                df = pd.read_excel(EXCEL_PATH, sheet_name=selected_sheet)
                
                # Colonnes requises et optionnelles
                required_cols = ['Compound', 'SMILES', 'Number', 'System']
                additional_cols = ['Log Pow', 'Log KD']
                
                # Vérification des colonnes disponibles
                available_cols = [col for col in required_cols + additional_cols if col in df.columns]
                
                if len([col for col in required_cols if col in available_cols]) < len(required_cols):
                    st.error("Les colonnes requises ne sont pas toutes présentes dans la feuille.")
                else:
                    # Affichage du tableau de données
                    st.dataframe(
                        df[available_cols],
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            "Number": st.column_config.NumberColumn(format="%d"),
                        }
                    )
                    
                    # Sélection d'une ligne
                    selected_number = st.selectbox(
                        "Sélectionnez un numéro",
                        df['Number'].unique(),
                        key="number_selection"
                    )
                    
                    selected_row = df[df['Number'] == selected_number].iloc[0]
                    
                    # Chargement des données du système correspondant
                    system_name = selected_row['System']
                    try:
                        df_system = pd.read_excel(DBDT_PATH, sheet_name=system_name)
                        df_filtered = df_system[df_system['Number'] == selected_number]
                        
                        if df_filtered.empty:
                            st.error(f"Aucune donnée trouvée pour le numéro {selected_number} dans le système {system_name}")
                        else:
                            # Création du graphique ternaire
                            fig = go.Figure()
                            
                            # Traces pour les phases UP et LP
                            for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                                fig.add_trace(go.Scatter(
                                    x=df_system[f'%Vol3 - {phase}'],
                                    y=df_system[f'%Vol2 - {phase}'],
                                    mode='markers',
                                    name=f'Phase {phase}',
                                    marker=dict(color=color, size=10, line=dict(width=1, color='DarkSlateGrey')),
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
                                    x=[df_filtered[f'%Vol3 - {phase}'].values[0]],
                                    y=[df_filtered[f'%Vol2 - {phase}'].values[0]],
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
                            
                            # Mise en forme du graphique
                            fig.update_layout(
                                title=dict(
                                    text=f"Diagramme de Phase Ternaire - Système {system_name}",
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
                            
                            # Affichage en deux colonnes (graphique + compositions)
                            col1, col2 = st.columns([0.7, 0.3])
                            
                            with col1:
                                st.plotly_chart(fig, use_container_width=True)
                            
                            with col2:
                                st.subheader("Compositions des phases")
                                
                                for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                                    with st.expander(f"Phase {phase}", expanded=True):
                                        st.markdown(f"""
                                        **%Vol1:** {df_filtered[f'%Vol1 - {phase}'].values[0]:.2f}%  
                                        **%Vol2:** {df_filtered[f'%Vol2 - {phase}'].values[0]:.2f}%  
                                        **%Vol3:** {df_filtered[f'%Vol3 - {phase}'].values[0]:.2f}%
                                        """)
                    
                    except Exception as e:
                        st.error(f"Erreur lors du chargement du système {system_name}: {str(e)}")
            
            except Exception as e:
                st.error(f"Erreur lors du chargement des données: {str(e)}")

    # Bouton pour effacer la recherche
    if st.session_state.search_triggered:
        if st.button("Effacer la recherche"):
            st.session_state.search_triggered = False
            st.rerun()
    
    # Bouton de retour
    if st.button("Retour à l'accueil", key="kddb_back"):
        st.session_state.current_page = "home"
        st.rerun()

def show_dbdt_page():
    """Page Ternary Phase Diagrams - Version complète"""
    st.title("Ternary Phase Diagrams")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(DBDT_PATH)
    if not sheet_names:
        return
    
    # Gestion des arguments passés
    initial_sheet = None
    selected_number = None
    if len(sys.argv) > 2:
        initial_sheet = sys.argv[1]
        selected_number = int(sys.argv[2])
    
    # Sélection de la feuille
    selected_sheet = st.selectbox(
        "Sélectionnez un système",
        sheet_names,
        index=sheet_names.index(initial_sheet) if initial_sheet in sheet_names else 0
    )
    
    # Chargement des données
    try:
        df = pd.read_excel(DBDT_PATH, sheet_name=selected_sheet)
        lib_row = pd.read_excel(DBDT_PATH, sheet_name=selected_sheet, header=None).iloc[1]
        labels = {
            'vol1': lib_row[0],
            'vol2': lib_row[1],
            'vol3': lib_row[2],
            'sheet': selected_sheet
        }
        
        # Nettoyage et tri
        df = df.dropna(subset=['%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP'])
        df = df.sort_values('%Vol3 - UP')
        
        # Création du graphique
        fig = go.Figure()

        # Styles
        line_style = dict(color='gray', width=1)
        up_marker = dict(color='red', size=8, symbol='circle')
        lp_marker = dict(color='blue', size=8, symbol='circle')

        # Courbe UP
        fig.add_trace(go.Scatter(
            x=df['%Vol3 - UP'],
            y=df['%Vol2 - UP'],
            mode='lines+markers',
            name='UP',
            line=line_style,
            marker=up_marker,
            customdata=df[['Number', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP']].values,
            hovertemplate=(
                f'Phase: UP<br>{labels["vol1"]}: %{{customdata[1]:.2f}}%<br>'
                f'{labels["vol2"]}: %{{customdata[2]:.2f}}%<br>'
                f'{labels["vol3"]}: %{{customdata[3]:.2f}}%<extra></extra>'
            )
        ))

        # Courbe LP
        fig.add_trace(go.Scatter(
            x=df['%Vol3 - LP'],
            y=df['%Vol2 - LP'],
            mode='lines+markers',
            name='LP',
            line=line_style,
            marker=lp_marker,
            customdata=df[['Number', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP']].values,
            hovertemplate=(
                f'Phase: LP<br>{labels["vol1"]}: %{{customdata[4]:.2f}}%<br>'
                f'{labels["vol2"]}: %{{customdata[5]:.2f}}%<br>'
                f'{labels["vol3"]}: %{{customdata[6]:.2f}}%<extra></extra>'
            )
        ))

        # Lignes de connexion
        for _, row in df.iterrows():
            fig.add_trace(go.Scatter(
                x=[row['%Vol3 - UP'], row['%Vol3 - LP']],
                y=[row['%Vol2 - UP'], row['%Vol2 - LP']],
                mode='lines',
                line=line_style,
                showlegend=False,
                hoverinfo='none'
            ))

        # Triangle de référence
        fig.add_trace(go.Scatter(
            x=[0, 1], y=[1, 0],
            mode='lines',
            line=dict(color='black', dash='dash'),
            name='x + y = 1',
            hoverinfo='none'
        ))

        # Configuration des axes avec grille améliorée
        fig.update_layout(
            xaxis=dict(
                title=f'% {labels["vol3"]} (%Vol3)',
                range=[0, 1],
                dtick=0.05,
                tick0=0,
                tickmode='linear',
                tickvals=np.arange(0, 1.05, 0.05).round(2),
                showgrid=True,
                gridwidth=1.5,
                gridcolor='#666666',
                griddash='solid',
                minor=dict(
                    ticklen=4,
                    dtick=0.01,
                    gridcolor='LightGrey',
                    gridwidth=0.5
                )
            ),
            yaxis=dict(
                title=f'% {labels["vol2"]} (%Vol2)',
                range=[0, 1],
                dtick=0.05,
                tick0=0,
                tickmode='linear',
                tickvals=np.arange(0, 1.05, 0.05).round(2),
                showgrid=True,
                gridwidth=1.5,
                gridcolor='#666666',
                griddash='solid',
                minor=dict(
                    ticklen=4,
                    dtick=0.01,
                    gridcolor='LightGrey',
                    gridwidth=0.5
                )
            ),
            title=f"Ternary Phase Diagram of {labels['vol1']} / {labels['vol2']} / {labels['vol3']} - {labels['sheet']}",
            width=1400,
            height=700,
            hovermode='closest'
        )

        # Affichage en deux colonnes
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
                index=list(df['Number']).index(selected_number) if selected_number in df['Number'].values else 0
            )
            
            selected_row = df[df['Number'] == selected_number].iloc[0]
            phase_data = create_phase_display(selected_row, labels)
            
            # Affichage des compositions
            st.subheader("Compositions des phases")
            
            for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                with st.expander(f"Phase {phase}", expanded=True):
                    st.markdown(f"""
                    **{labels['vol1']}:** {phase_data[phase]['vol1']:.2f}%  
                    **{labels['vol2']}:** {phase_data[phase]['vol2']:.2f}%  
                    **{labels['vol3']}:** {phase_data[phase]['vol3']:.2f}%
                    """)
    
    except Exception as e:
        st.error(f"Erreur lors du traitement des données: {str(e)}")
    
    # Bouton de retour
    if st.button("Retour à l'accueil", key="dbdt_back"):
        st.session_state.current_page = "home"
        st.rerun()

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
            st.rerun()
        if st.button("KD Database Explorer"):
            st.session_state.current_page = "kddb"
            st.rerun()
        if st.button("Ternary Diagrams"):
            st.session_state.current_page = "dbdt"
            st.rerun()
    
    # Router vers la page active
    if st.session_state.current_page == "home":
        show_home_page()
    elif st.session_state.current_page == "kddb":
        show_kddb_page()
    elif st.session_state.current_page == "dbdt":
        show_dbdt_page()

if __name__ == "__main__":
    main()