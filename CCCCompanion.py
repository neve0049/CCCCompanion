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
    initial_sidebar_state="collapsed"
)

# =============================================
# Fonctions utilitaires
# =============================================

def is_mobile():
    """Détecte si l'utilisateur est sur un appareil mobile"""
    user_agent = st.experimental_get_query_params().get("user_agent", [""])[0]
    mobile_keywords = ['mobile', 'android', 'iphone', 'ipad', 'windows phone']
    return any(keyword in user_agent.lower() for keyword in mobile_keywords)

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
    
    if is_mobile():
        # Version mobile
        if st.button("KD Database Explorer", use_container_width=True, key="kddb_home"):
            st.session_state.current_page = "kddb"
            st.rerun()
        if st.button("Ternary Phase Diagrams", use_container_width=True, key="dbdt_home"):
            st.session_state.current_page = "dbdt"
            st.rerun()
    else:
        # Version desktop
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
    """Page KD Database Explorer - Version adaptée mobile"""
    st.title("KD Database Explorer")
    
    if is_mobile():
        st.markdown("ℹ️ Utilisez le mode paysage pour une meilleure expérience")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(EXCEL_PATH)
    if not sheet_names:
        return
        
    # Zone de recherche adaptée mobile
    if is_mobile():
        search_query = st.text_input(
            "Recherche...",
            key="search_input",
            placeholder="Molécule/système"
        )
        if st.button("Rechercher", use_container_width=True, key="search_button"):
            st.session_state.search_triggered = True
    else:
        col1, col2 = st.columns([0.7, 0.3])
        with col1:
            search_query = st.text_input(
                "Entrez un nom de molécule ou système...",
                key="search_input",
                placeholder="Rechercher"
            )
        with col2:
            st.write("")
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
            
            # Chargement des données
            try:
                df = pd.read_excel(EXCEL_PATH, sheet_name=selected_sheet)
                
                # Colonnes requises et optionnelles
                required_cols = ['Compound', 'SMILES', 'Number', 'System']
                additional_cols = ['Log P (Pubchem)', 'Log P (COSMO-RS)', 'Log KD']
                
                # Vérification des colonnes
                available_cols = [col for col in required_cols + additional_cols if col in df.columns]
                
                if len([col for col in required_cols if col in available_cols]) < len(required_cols):
                    st.error("Colonnes requises manquantes.")
                else:
                    # Sélection interactive
                    st.subheader("Sélectionnez une entrée")
                    
                    selection_df = df[available_cols].copy()
                    selection_df.insert(0, 'Select', False)
                    
                    column_config = {
                        "Select": st.column_config.CheckboxColumn(
                            "Sélection",
                            help="Sélectionnez une ligne",
                            default=False,
                            required=True
                        ),
                        "Number": st.column_config.NumberColumn(format="%d"),
                    }
                    
                    edited_df = st.data_editor(
                        selection_df,
                        column_config=column_config,
                        use_container_width=True,
                        hide_index=True,
                        disabled=[col for col in available_cols],
                        key="data_editor"
                    )
                    
                    selected_rows = edited_df[edited_df['Select']]
                    
                    if not selected_rows.empty:
                        selected_row = selected_rows.iloc[0]
                        system_name = selected_row['System']
                        selected_number = selected_row['Number']
                        
                        try:
                            all_sheets = pd.read_excel(DBDT_PATH, sheet_name=None)
                            
                            if system_name not in all_sheets:
                                st.error(f"Système {system_name} non trouvé")
                            else:
                                df_system = all_sheets[system_name]
                                df_filtered = df_system[df_system['Number'] == selected_number]
                                
                                if df_filtered.empty:
                                    st.error(f"Numéro {selected_number} non trouvé")
                                else:
                                    # Création du graphique
                                    fig = go.Figure()
                                    
                                    # Traces pour les phases
                                    for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                                        fig.add_trace(go.Scatter(
                                            x=df_system[f'%Vol3 - {phase}'],
                                            y=df_system[f'%Vol2 - {phase}'],
                                            mode='markers',
                                            name=f'Phase {phase}',
                                            marker=dict(color=color, size=8 if is_mobile() else 10),
                                            customdata=np.stack((
                                                df_system['Number'],
                                                df_system[f'%Vol1 - {phase}'],
                                                df_system[f'%Vol2 - {phase}'],
                                                df_system[f'%Vol3 - {phase}']
                                            ), axis=-1),
                                            hovertemplate=(
                                                "<b>Number</b>: %{customdata[0]}<br>"
                                                "<b>Vol1</b>: %{customdata[1]:.2f}<br>"
                                                "<b>Vol2</b>: %{customdata[2]:.2f}<br>"
                                                "<b>Vol3</b>: %{customdata[3]:.2f}<extra></extra>"
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
                                        line=dict(color='green', width=1.5 if is_mobile() else 2),
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
                                                size=12 if is_mobile() else 16,
                                                symbol='circle-open',
                                                line=dict(width=1.5 if is_mobile() else 2)
                                            ),
                                            hoverinfo='none'
                                        ))
                                    
                                    # Mise en forme
                                    fig.update_layout(
                                        title=dict(
                                            text=f"Diagramme de Phase - {system_name}",
                                            x=0.5,
                                            font=dict(size=14 if is_mobile() else 16)
                                        ),
                                        xaxis=dict(
                                            title='Vol3',
                                            range=[0, 1]
                                        ),
                                        yaxis=dict(
                                            title='Vol2',
                                            range=[0, 1]
                                        ),
                                        showlegend=True,
                                        legend=dict(
                                            orientation='h',
                                            yanchor='bottom',
                                            y=1.02,
                                            xanchor='right',
                                            x=1
                                        ),
                                        margin=dict(
                                            l=40 if is_mobile() else 60,
                                            r=40 if is_mobile() else 60,
                                            t=60 if is_mobile() else 80,
                                            b=40 if is_mobile() else 60
                                        ),
                                        height=400 if is_mobile() else 600
                                    )
                                    
                                    # Affichage adapté mobile
                                    if is_mobile():
                                        st.plotly_chart(fig, use_container_width=True)
                                        st.subheader("Compositions")
                                        
                                        for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                                            with st.expander(f"Phase {phase}", expanded=False):
                                                st.markdown(f"""
                                                **Vol1:** {df_filtered[f'%Vol1 - {phase}'].values[0]:.2f}  
                                                **Vol2:** {df_filtered[f'%Vol2 - {phase}'].values[0]:.2f}  
                                                **Vol3:** {df_filtered[f'%Vol3 - {phase}'].values[0]:.2f}
                                                """)
                                    else:
                                        col1, col2 = st.columns([0.7, 0.3])
                                        with col1:
                                            st.plotly_chart(fig, use_container_width=True)
                                        with col2:
                                            st.subheader("Compositions des phases")
                                            for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                                                with st.expander(f"Phase {phase}", expanded=True):
                                                    st.markdown(f"""
                                                    **Vol1:** {df_filtered[f'%Vol1 - {phase}'].values[0]:.2f}  
                                                    **Vol2:** {df_filtered[f'%Vol2 - {phase}'].values[0]:.2f}  
                                                    **Vol3:** {df_filtered[f'%Vol3 - {phase}'].values[0]:.2f}
                                                    """)
                        except Exception as e:
                            st.error(f"Erreur: {str(e)}")
                    else:
                        st.info("Sélectionnez une ligne pour voir les détails")
            
            except Exception as e:
                st.error(f"Erreur: {str(e)}")

    # Boutons
    if st.session_state.search_triggered:
        if st.button("Effacer la recherche", use_container_width=is_mobile()):
            st.session_state.search_triggered = False
            st.rerun()
    
    if st.button("Retour à l'accueil", key="kddb_back", use_container_width=is_mobile()):
        st.session_state.current_page = "home"
        st.rerun()

def show_dbdt_page():
    """Page Ternary Phase Diagrams - Version adaptée mobile"""
    st.title("Ternary Phase Diagrams")
    
    if is_mobile():
        st.markdown("ℹ️ Utilisez le mode paysage pour une meilleure visualisation")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(DBDT_PATH)
    if not sheet_names:
        return
    
    # Gestion des arguments
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
        
        # Nettoyage
        df = df.dropna(subset=['%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP'])
        df = df.sort_values('%Vol3 - UP')
        
        # Création du graphique
        fig = go.Figure()

        # Courbe UP
        fig.add_trace(go.Scatter(
            x=df['%Vol3 - UP'],
            y=df['%Vol2 - UP'],
            mode='lines+markers',
            name='UP',
            line=dict(color='black', width=1.5 if is_mobile() else 2),
            marker=dict(color='red', size=6 if is_mobile() else 8),
            customdata=df[['Number', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP']].values,
            hovertemplate=(
                f'Phase: UP<br>{labels["vol1"]}: %{{customdata[1]:.2f}}<br>'
                f'{labels["vol2"]}: %{{customdata[2]:.2f}}<br>'
                f'{labels["vol3"]}: %{{customdata[3]:.2f}}<extra></extra>'
            )
        ))

        # Courbe LP
        fig.add_trace(go.Scatter(
            x=df['%Vol3 - LP'],
            y=df['%Vol2 - LP'],
            mode='lines+markers',
            name='LP',
            line=dict(color='black', width=1.5 if is_mobile() else 2),
            marker=dict(color='blue', size=6 if is_mobile() else 8),
            customdata=df[['Number', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP']].values,
            hovertemplate=(
                f'Phase: LP<br>{labels["vol1"]}: %{{customdata[1]:.2f}}<br>'
                f'{labels["vol2"]}: %{{customdata[2]:.2f}}<br>'
                f'{labels["vol3"]}: %{{customdata[3]:.2f}}<extra></extra>'
            )
        ))

        # Lignes de connexion
        for _, row in df.iterrows():
            fig.add_trace(go.Scatter(
                x=[row['%Vol3 - UP'], row['%Vol3 - LP']],
                y=[row['%Vol2 - UP'], row['%Vol2 - LP']],
                mode='lines',
                line=dict(color='blue', width=1.5 if is_mobile() else 2),
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

        # Configuration
        fig.update_layout(
            xaxis=dict(
                title=f'{labels["vol3"]} (Vol3)',
                range=[0, 1],
                dtick=0.1,
                showgrid=True
            ),
            yaxis=dict(
                title=f'{labels["vol2"]} (Vol2)',
                range=[0, 1],
                dtick=0.05,
                showgrid=True
            ),
            title=f"Diagramme de {labels['vol1']} / {labels['vol2']} / {labels['vol3']}",
            height=500 if is_mobile() else 700,
            hovermode='closest',
            margin=dict(
                l=40 if is_mobile() else 60,
                r=40 if is_mobile() else 60,
                t=80 if is_mobile() else 100,
                b=40 if is_mobile() else 60
            )
        )

        # Affichage adapté mobile
        if is_mobile():
            st.plotly_chart(fig, use_container_width=True)
            
            # Sélection
            selected_number = st.selectbox(
                "Choisir un numéro", 
                df['Number'].unique(),
                index=list(df['Number']).index(selected_number) if selected_number in df['Number'].values else 0
            )
            
            selected_row = df[df['Number'] == selected_number].iloc[0]
            phase_data = create_phase_display(selected_row, labels)
            
            # Compositions
            st.subheader("Compositions")
            for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                with st.expander(f"Phase {phase}", expanded=False):
                    st.markdown(f"""
                    **{labels['vol1']}:** {phase_data[phase]['vol1']:.2f}  
                    **{labels['vol2']}:** {phase_data[phase]['vol2']:.2f}  
                    **{labels['vol3']}:** {phase_data[phase]['vol3']:.2f}
                    """)
            
            # Données brutes
            st.subheader("Données")
            with st.container(height=300, border=True):
                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True,
                    column_config={
                        "Number": st.column_config.NumberColumn("Numéro", format="%d"),
                        "%Vol1 - UP": st.column_config.NumberColumn("UP Vol1", format="%.2f"),
                        "%Vol2 - UP": st.column_config.NumberColumn("UP Vol2", format="%.2f"),
                        "%Vol3 - UP": st.column_config.NumberColumn("UP Vol3", format="%.2f"),
                        "%Vol1 - LP": st.column_config.NumberColumn("LP Vol1", format="%.2f"),
                        "%Vol2 - LP": st.column_config.NumberColumn("LP Vol2", format="%.2f"),
                        "%Vol3 - LP": st.column_config.NumberColumn("LP Vol3", format="%.2f")
                    }
                )
        else:
            # Version desktop
            col1, col2 = st.columns([0.8, 0.2])
            with col1:
                st.plotly_chart(fig, use_container_width=True)
                
                # Données brutes
                st.subheader("Données brutes")
                with st.container(height=400, border=True):
                    st.dataframe(
                        df,
                        use_container_width=True,
                        hide_index=True,
                        column_config={
                            "Number": st.column_config.NumberColumn("Numéro", format="%d"),
                            "%Vol1 - UP": st.column_config.NumberColumn("UP Vol1", format="%.2f"),
                            "%Vol2 - UP": st.column_config.NumberColumn("UP Vol2", format="%.2f"),
                            "%Vol3 - UP": st.column_config.NumberColumn("UP Vol3", format="%.2f"),
                            "%Vol1 - LP": st.column_config.NumberColumn("LP Vol1", format="%.2f"),
                            "%Vol2 - LP": st.column_config.NumberColumn("LP Vol2", format="%.2f"),
                            "%Vol3 - LP": st.column_config.NumberColumn("LP Vol3", format="%.2f")
                        }
                    )
            with col2:
                # Sélection
                st.subheader("Sélection")
                selected_number = st.selectbox(
                    "Choisir un numéro", 
                    df['Number'].unique(),
                    index=list(df['Number']).index(selected_number) if selected_number in df['Number'].values else 0
                )
                
                selected_row = df[df['Number'] == selected_number].iloc[0]
                phase_data = create_phase_display(selected_row, labels)
                
                # Compositions
                st.subheader("Compositions")
                for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                    with st.expander(f"Phase {phase}", expanded=True):
                        st.markdown(f"""
                        **{labels['vol1']}:** {phase_data[phase]['vol1']:.2f}  
                        **{labels['vol2']}:** {phase_data[phase]['vol2']:.2f}  
                        **{labels['vol3']}:** {phase_data[phase]['vol3']:.2f}
                        """)
    
    except Exception as e:
        st.error(f"Erreur: {str(e)}")
    
    # Bouton de retour
    if st.button("Retour à l'accueil", key="dbdt_back", use_container_width=is_mobile()):
        st.session_state.current_page = "home"
        st.rerun()

# =============================================
# Application principale
# =============================================

def main():
    # Initialisation de l'état
    if 'current_page' not in st.session_state:
        st.session_state.current_page = "home"
    
    # Navigation adaptée mobile
    if is_mobile():
        # Menu en haut pour mobile
        menu = ["Accueil", "KD Database", "Ternary Diagrams"]
        choice = st.selectbox("Menu", menu, label_visibility="collapsed")
        
        if choice == "Accueil":
            st.session_state.current_page = "home"
        elif choice == "KD Database":
            st.session_state.current_page = "kddb"
        elif choice == "Ternary Diagrams":
            st.session_state.current_page = "dbdt"
    else:
        # Sidebar pour desktop
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
