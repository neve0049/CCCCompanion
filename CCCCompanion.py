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

def configure_aggrid(df, selection_mode='single'):
    """Configure AgGrid pour la sélection de lignes"""
    gb = GridOptionsBuilder.from_dataframe(df)
    gb.configure_selection(
        selection_mode=selection_mode,
        use_checkbox=True,
        header_checkbox=True if selection_mode == 'multiple' else False,
        header_checkbox_filtered_only=True
    )
    gb.configure_pagination(paginationAutoPageSize=False, paginationPageSize=10)
    return gb.build()

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
    """Page KD Database Explorer avec sélection via cases à cocher"""
    st.title("KD Database Explorer")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(EXCEL_PATH)
    if not sheet_names:
        return
        
    # Sélection de la feuille
    selected_sheet = st.selectbox("Sélectionnez une feuille", sheet_names)
    
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
            # Configuration d'AgGrid avec cases à cocher
            grid_options = configure_aggrid(df[available_cols])
            
            # Affichage du tableau avec sélection
            grid_response = AgGrid(
                df[available_cols],
                gridOptions=grid_options,
                height=400,
                width='100%',
                data_return_mode='FILTERED_AND_SORTED',
                update_mode='MODEL_CHANGED'
            )
            
            selected_rows = grid_response['selected_rows']
            
            if selected_rows:
                selected_row = selected_rows[0]
                selected_number = selected_row['Number']
                system_name = selected_row['System']
                
                # Chargement des données du système correspondant
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
    
    # Bouton de retour
    if st.button("Retour à l'accueil", key="kddb_back"):
        st.session_state.current_page = "home"
        st.rerun()

def show_dbdt_page():
    """Page Ternary Phase Diagrams avec sélection interactive"""
    st.title("Ternary Phase Diagrams")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(DBDT_PATH)
    if not sheet_names:
        return
    
    # Sélection de la feuille
    selected_sheet = st.selectbox("Sélectionnez un système", sheet_names)
    
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
            customdata=df[['Number', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP']].values,
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
            # Affichage du graphique
            st.plotly_chart(fig, use_container_width=True, key="ternary_chart")
            
            # Affichage des données brutes
            st.subheader("Données brutes")
            st.dataframe(df, use_container_width=True)
        
        with col2:
            st.subheader("Compositions des phases")
            
            # Sélection manuelle du point
            selected_number = st.selectbox(
                "Sélectionnez un numéro", 
                df['Number'].unique(),
                key="dbdt_number_select"
            )
            
            selected_row = df[df['Number'] == selected_number].iloc[0]
            
            for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                with st.expander(f"Phase {phase}", expanded=True):
                    st.markdown(f"""
                    **{labels['vol1']}:** {selected_row[f'%Vol1 - {phase}']:.2f}%  
                    **{labels['vol2']}:** {selected_row[f'%Vol2 - {phase}']:.2f}%  
                    **{labels['vol3']}:** {selected_row[f'%Vol3 - {phase}']:.2f}%
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