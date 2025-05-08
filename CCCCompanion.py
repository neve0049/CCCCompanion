import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import sys

# Configuration des chemins des fichiers
EXCEL_PATH = "KDDB.xlsx"
DBDT_PATH = "DBDT.xlsx"
DBDQ_PATH = "DBDQ.xlsx"

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

def create_phase_display(row, labels, is_quaternary=False):
    """Crée l'affichage des compositions des phases"""
    if is_quaternary:
        return {
            "UP": {
                "vol1": row['%Vol1 - UP'],
                "vol2": row['%Vol2 - UP'],
                "vol3": row['%Vol3 - UP'],
                "vol4": row['%Vol4 - UP']
            },
            "LP": {
                "vol1": row['%Vol1 - LP'],
                "vol2": row['%Vol2 - LP'],
                "vol3": row['%Vol3 - LP'],
                "vol4": row['%Vol4 - LP']
            },
            "labels": labels
        }
    else:
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

def quaternary_to_3d(v1, v2, v3, v4):
    """Convertit des coordonnées quaternaires (v1+v2+v3+v4=100) en coordonnées 3D"""
    # Normalisation des pourcentages (de 0-100 à 0-1)
    total = v1 + v2 + v3 + v4
    v1_norm = v1 / total
    v2_norm = v2 / total
    v3_norm = v3 / total
    v4_norm = v4 / total
    
    # Configuration des axes
    x = v4_norm  # %Vol4 sur l'axe X
    y = v3_norm  # %Vol3 sur l'axe Y
    z = v2_norm  # %Vol2 sur l'axe Z
    
    return x, y, z

# =============================================
# Pages de l'application
# =============================================

def show_home_page():
    st.title("CCCCompanion ")
    st.markdown("""
    Welcome to the CCCCompanion app. Please select a database to explore:
    """)
    
    col1, col2, col3 = st.columns(3)
    with col1:
        if st.button("KD Database Explorer", use_container_width=True, key="kddb_home"):
            st.session_state.current_page = "kddb"
            st.rerun()
    with col2:
        if st.button("Ternary Phase Diagrams", use_container_width=True, key="dbdt_home"):
            st.session_state.current_page = "dbdt"
            st.rerun()
    with col3:
        if st.button("Quaternary Phase Diagrams", use_container_width=True, key="dbdq_home"):
            st.session_state.current_page = "dbdq"
            st.rerun()

def show_kddb_page():
    """Page KD Database Explorer - Version corrigée"""
    st.title("KD Database Explorer")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(EXCEL_PATH)
    if not sheet_names:
        return
        
    # Zone de recherche améliorée
    search_query = st.text_input(
        "Enter a compound or system name",
        value=st.session_state.get("search_query", ""),
        key="search_input",
        placeholder="Search by compound or system"
    )
    
    # Bouton de recherche avec état
    if st.button("Search", key="search_button") or search_query:
        st.session_state.search_triggered = True
        st.session_state.search_query = search_query

    # Si recherche déclenchée
    if st.session_state.get("search_triggered", False):
        search_value = st.session_state.search_query.strip().lower()
        matching_sheets = [sheet for sheet in sheet_names if search_value in sheet.lower()]
        
        if not matching_sheets:
            st.warning("No matching sheets found.")
        else:
            # Sélection de la feuille avec valeur par défaut
            selected_sheet = st.selectbox(
                "Select a sheet:",
                matching_sheets,
                key="sheet_selection"
            )
            
            try:
                # Chargement des données avec gestion d'erreur améliorée
                df = pd.read_excel(EXCEL_PATH, sheet_name=selected_sheet)
                
                # Colonnes requises avec alternatives possibles
                required_cols = {
                    'compound': ['Compound', 'Composé', 'Molecule'],
                    'system': ['System', 'Système', 'Phase System'],
                    'number': ['Number', 'Numéro', 'ID']
                }
                
                # Trouver les noms de colonnes réels
                actual_cols = {}
                for col_type, alternatives in required_cols.items():
                    for alt in alternatives:
                        if alt in df.columns:
                            actual_cols[col_type] = alt
                            break
                    if col_type not in actual_cols:
                        st.error(f"Required column not found (tried: {', '.join(alternatives)})")
                        return
                
                # Colonnes supplémentaires optionnelles
                optional_cols = ['Log KD', 'Log Pow', 'Log P (Pubchem)', 'Log P (COSMO-RS)']
                display_cols = [actual_cols['compound'], actual_cols['number'], actual_cols['system']] + \
                             [col for col in optional_cols if col in df.columns]
                
                # Affichage du tableau de sélection
                st.subheader("Select an entry")
                selection_df = df[display_cols].copy()
                selection_df.insert(0, 'Select', False)
                
                # Éditeur de données avec sélection
                edited_df = st.data_editor(
                    selection_df,
                    column_config={
                        "Select": st.column_config.CheckboxColumn(
                            "Select",
                            default=False,
                            disabled=[col for col in display_cols]
                        ),
                        actual_cols['number']: st.column_config.NumberColumn(format="%d"),
                    },
                    use_container_width=True,
                    hide_index=True,
                    key="data_editor"
                )
                
                # Récupération de la sélection
                selected_rows = edited_df[edited_df['Select']]
                
                if not selected_rows.empty:
                    selected_row = df.iloc[selected_rows.index[0]]
                    system_name = selected_row[actual_cols['system']]
                    selected_number = selected_row[actual_cols['number']]
                    
                    # Vérification du type de système
                    is_quaternary = st.checkbox(
                        "Show as quaternary diagram", 
                        key="quaternary_check",
                        help="Check if this is a 4-component system"
                    )
                    
                    # Chargement des données du diagramme
                    target_file = DBDQ_PATH if is_quaternary else DBDT_PATH
                    try:
                        all_sheets = pd.read_excel(target_file, sheet_name=None)
                        
                        # Nettoyage du nom du système pour la correspondance
                        clean_system_name = str(system_name).strip()
                        
                        if clean_system_name not in all_sheets:
                            st.error(f"No diagram data found for system: {clean_system_name}")
                            st.info(f"Available systems: {', '.join(all_sheets.keys())}")
                        else:
                            df_system = all_sheets[clean_system_name]
                            
                            # Conversion cohérente du numéro (en cas de types différents)
                            df_system['Number'] = pd.to_numeric(df_system['Number'], errors='coerce')
                            selected_number = pd.to_numeric(selected_number, errors='coerce')
                            
                            df_filtered = df_system[df_system['Number'] == selected_number]
                            
                            if df_filtered.empty:
                                st.error(f"No data found for number {selected_number} in system {clean_system_name}")
                                st.info(f"Available numbers: {', '.join(map(str, df_system['Number'].unique()))}")
                            else:
                                if is_quaternary:
                                    show_quaternary_diagram(df_system, df_filtered, clean_system_name, selected_number)
                                else:
                                    show_ternary_diagram(df_system, df_filtered, clean_system_name, selected_number)
                    
                    except Exception as e:
                        st.error(f"Error loading diagram data: {str(e)}")
                
            except Exception as e:
                st.error(f"Error processing data: {str(e)}")
                
        # Bouton pour effacer la recherche
        if st.button("Clear search"):
            st.session_state.search_triggered = False
            st.session_state.search_query = ""
            st.rerun()
    
    # Bouton de retour
    if st.button("Back to Home", key="kddb_back"):
        st.session_state.current_page = "home"
        st.rerun()

def show_ternary_diagram(df_system, df_filtered, system_name, selected_number):
    """Affiche un diagramme ternaire"""
    lib_row = pd.read_excel(DBDT_PATH, sheet_name=system_name, header=None).iloc[1]
    labels = {
        'vol1': lib_row[0],
        'vol2': lib_row[1],
        'vol3': lib_row[2],
        'sheet': system_name
    }
    
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
                "<b>%Vol1</b>: %{customdata[1]:.2f}<br>"
                "<b>%Vol2</b>: %{customdata[2]:.2f}<br>"
                "<b>%Vol3</b>: %{customdata[3]:.2f}<extra></extra>"
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
        
        phase_data = create_phase_display(df_filtered.iloc[0], labels)
        
        for phase, color in [('UP', 'red'), ('LP', 'blue')]:
            with st.expander(f"Phase {phase}", expanded=True):
                st.markdown(f"""
                **{labels['vol1']}:** {phase_data[phase]['vol1']:.2f}  
                **{labels['vol2']}:** {phase_data[phase]['vol2']:.2f}  
                **{labels['vol3']}:** {phase_data[phase]['vol3']:.2f}
                """)

def show_quaternary_diagram(df_system, df_filtered, system_name, selected_number):
    """Affiche un diagramme quaternaire"""
    lib_row = pd.read_excel(DBDQ_PATH, sheet_name=system_name, header=None).iloc[1]
    labels = {
        'vol1': lib_row[0],
        'vol2': lib_row[1],
        'vol3': lib_row[2],
        'vol4': lib_row[3],
        'sheet': system_name
    }
    
    # Conversion des coordonnées
    def convert_row(row, phase):
        v1 = row[f'%Vol1 - {phase}']
        v2 = row[f'%Vol2 - {phase}']
        v3 = row[f'%Vol3 - {phase}']
        v4 = row[f'%Vol4 - {phase}']
        return quaternary_to_3d(v1, v2, v3, v4)
    
    # Appliquer la conversion
    df_system[['x_up', 'y_up', 'z_up']] = df_system.apply(lambda x: convert_row(x, 'UP'), axis=1, result_type='expand')
    df_system[['x_lp', 'y_lp', 'z_lp']] = df_system.apply(lambda x: convert_row(x, 'LP'), axis=1, result_type='expand')
    df_filtered[['x_up', 'y_up', 'z_up']] = df_filtered.apply(lambda x: convert_row(x, 'UP'), axis=1, result_type='expand')
    df_filtered[['x_lp', 'y_lp', 'z_lp']] = df_filtered.apply(lambda x: convert_row(x, 'LP'), axis=1, result_type='expand')
    
    # Création de la figure 3D
    fig = go.Figure()
    
    # Ajout des points UP (rouge)
    fig.add_trace(go.Scatter3d(
        x=df_system['x_up'],
        y=df_system['y_up'],
        z=df_system['z_up'],
        mode='markers',
        name='UP',
        marker=dict(color='red', size=5),
        customdata=df_system[['Number', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', '%Vol4 - UP']].values,
        hovertemplate=(
            f"<b>Phase UP</b><br>"
            f"{labels['vol1']}: %{{customdata[1]:.2f}}%<br>"
            f"{labels['vol2']}: %{{customdata[2]:.2f}}%<br>"
            f"{labels['vol3']}: %{{customdata[3]:.2f}}%<br>"
            f"{labels['vol4']}: %{{customdata[4]:.2f}}%<extra></extra>"
        )
    ))
    
    # Ajout des points LP (bleu)
    fig.add_trace(go.Scatter3d(
        x=df_system['x_lp'],
        y=df_system['y_lp'],
        z=df_system['z_lp'],
        mode='markers',
        name='LP',
        marker=dict(color='blue', size=5),
        customdata=df_system[['Number', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP', '%Vol4 - LP']].values,
        hovertemplate=(
            f"<b>Phase LP</b><br>"
            f"{labels['vol1']}: %{{customdata[1]:.2f}}%<br>"
            f"{labels['vol2']}: %{{customdata[2]:.2f}}%<br>"
            f"{labels['vol3']}: %{{customdata[3]:.2f}}%<br>"
            f"{labels['vol4']}: %{{customdata[4]:.2f}}%<extra></extra>"
        )
    ))
    
    # Ajout des lignes de connexion
    for _, row in df_system.iterrows():
        fig.add_trace(go.Scatter3d(
            x=[row['x_up'], row['x_lp']],
            y=[row['y_up'], row['y_lp']],
            z=[row['z_up'], row['z_lp']],
            mode='lines',
            line=dict(color='gray', width=1),
            showlegend=False,
            hoverinfo='none'
        ))
    
    # Points sélectionnés
    for phase, color in [('UP', 'red'), ('LP', 'blue')]:
        fig.add_trace(go.Scatter3d(
            x=[df_filtered[f'x_{phase.lower()}'].values[0]],
            y=[df_filtered[f'y_{phase.lower()}'].values[0]],
            z=[df_filtered[f'z_{phase.lower()}'].values[0]],
            mode='markers',
            name=f'Sélection {phase}',
            marker=dict(
                color='black',
                size=10,
                symbol='circle-open',
                line=dict(width=2)
            ),
            hoverinfo='none'
        ))
    
    # Configuration du layout avec les nouveaux axes
    fig.update_layout(
        scene=dict(
            xaxis=dict(title=f"% {labels['vol4']} (X)"),
            yaxis=dict(title=f"% {labels['vol3']} (Y)"),
            zaxis=dict(title=f"% {labels['vol2']} (Z)"),
            aspectmode='cube',
            camera=dict(
                up=dict(x=0, y=0, z=1),
                center=dict(x=0, y=0, z=0),
                eye=dict(x=1.25, y=1.25, z=1.25)
            )
        ),
        title=f"Quaternary Phase Diagram: {labels['vol1']} / {labels['vol2']} / {labels['vol3']} / {labels['vol4']}",
        width=1200,
        height=800,
        hovermode='closest'
    )
    
    # Affichage en deux colonnes (graphique + compositions)
    col1, col2 = st.columns([0.7, 0.3])
    
    with col1:
        st.plotly_chart(fig, use_container_width=True)
    
    with col2:
        st.subheader("Compositions des phases")
        
        phase_data = create_phase_display(df_filtered.iloc[0], labels, is_quaternary=True)
        
        for phase, color in [('UP', 'red'), ('LP', 'blue')]:
            with st.expander(f"Phase {phase}", expanded=True):
                st.markdown(f"""
                **{labels['vol1']}:** {phase_data[phase]['vol1']:.2f}%  
                **{labels['vol2']}:** {phase_data[phase]['vol2']:.2f}%  
                **{labels['vol3']}:** {phase_data[phase]['vol3']:.2f}%  
                **{labels['vol4']}:** {phase_data[phase]['vol4']:.2f}%
                """)

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
            line=dict(
            color='black',  # Ligne noire
            width=2        # Épaisseur de 1 pixel
        ),
            marker=up_marker,
            customdata=df[['Number', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP']].values,
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
            line=dict(
            color='black',  # Ligne noire
            width=2        # Épaisseur de 1 pixel
        ),
            marker=lp_marker,
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
                line=dict(
                color='blue',  # Ligne bleue
                width=2        # Épaisseur de 1 pixel
            ),
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
                dtick=0.1,
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

        # Dans la partie "Affichage en deux colonnes"
        col1, col2 = st.columns([0.8, 0.2])

        # Dans la fonction show_dbdt_page(), remplacez la partie "Données brutes" par :

        with col1:
            st.plotly_chart(fig, use_container_width=True)
    
            # Section Données brutes avec utilisation optimale de l'espace
            st.subheader("Data of selected system")
    
            # Définir la hauteur en fonction du nombre de lignes (min 200px, max 600px)
            table_height = min(200 + len(df) * 35, 600)
    
            # Container avec bordure et défilement si nécessaire
            with st.container(height=table_height, border=True):
                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True,
                    height=table_height - 10,  # Compenser la bordure
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
            # Sélection interactive
            st.subheader("Select a system")
            selected_number = st.selectbox(
                "Select a number", 
                df['Number'].unique(),
                index=list(df['Number']).index(selected_number) if selected_number in df['Number'].values else 0
            )
            
            selected_row = df[df['Number'] == selected_number].iloc[0]
            phase_data = create_phase_display(selected_row, labels)
            
            # Affichage des compositions
            st.subheader("Phase composition")
            
            for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                with st.expander(f"Phase {phase}", expanded=True):
                    st.markdown(f"""
                    **{labels['vol1']}:** {phase_data[phase]['vol1']:.2f}  
                    **{labels['vol2']}:** {phase_data[phase]['vol2']:.2f}  
                    **{labels['vol3']}:** {phase_data[phase]['vol3']:.2f}
                    """)
    
    except Exception as e:
        st.error(f"Erreur lors du traitement des données: {str(e)}")
    
    # Bouton de retour
    if st.button("Retour à l'accueil", key="dbdt_back"):
        st.session_state.current_page = "home"
        st.rerun()

def show_dbdq_page():
    """Page Quaternary Phase Diagrams"""
    st.title("Quaternary Phase Diagrams")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(DBDQ_PATH)
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
        df = pd.read_excel(DBDQ_PATH, sheet_name=selected_sheet)
        lib_row = pd.read_excel(DBDQ_PATH, sheet_name=selected_sheet, header=None).iloc[1]
        labels = {
            'vol1': lib_row[0],
            'vol2': lib_row[1],
            'vol3': lib_row[2],
            'vol4': lib_row[3],
            'sheet': selected_sheet
        }
        
        # Nettoyage des données
        required_cols = ['%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', '%Vol4 - UP',
                        '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP', '%Vol4 - LP']
        df = df.dropna(subset=required_cols)
        
        # Conversion des coordonnées
        def convert_row(row, phase):
            v1 = row[f'%Vol1 - {phase}']
            v2 = row[f'%Vol2 - {phase}']
            v3 = row[f'%Vol3 - {phase}']
            v4 = row[f'%Vol4 - {phase}']
            return quaternary_to_3d(v1, v2, v3, v4)
        
        # Appliquer la conversion
        df[['x_up', 'y_up', 'z_up']] = df.apply(lambda x: convert_row(x, 'UP'), axis=1, result_type='expand')
        df[['x_lp', 'y_lp', 'z_lp']] = df.apply(lambda x: convert_row(x, 'LP'), axis=1, result_type='expand')
        
        # Création de la figure 3D
        fig = go.Figure()
        
        # Ajout des points UP (rouge)
        fig.add_trace(go.Scatter3d(
            x=df['x_up'],
            y=df['y_up'],
            z=df['z_up'],
            mode='markers',
            name='UP',
            marker=dict(color='red', size=5),
            customdata=df[['Number', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', '%Vol4 - UP']].values,
            hovertemplate=(
                f"<b>Phase UP</b><br>"
                f"{labels['vol1']}: %{{customdata[1]:.2f}}%<br>"
                f"{labels['vol2']}: %{{customdata[2]:.2f}}%<br>"
                f"{labels['vol3']}: %{{customdata[3]:.2f}}%<br>"
                f"{labels['vol4']}: %{{customdata[4]:.2f}}%<extra></extra>"
            )
        ))
        
        # Ajout des points LP (bleu)
        fig.add_trace(go.Scatter3d(
            x=df['x_lp'],
            y=df['y_lp'],
            z=df['z_lp'],
            mode='markers',
            name='LP',
            marker=dict(color='blue', size=5),
            customdata=df[['Number', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP', '%Vol4 - LP']].values,
            hovertemplate=(
                f"<b>Phase LP</b><br>"
                f"{labels['vol1']}: %{{customdata[1]:.2f}}%<br>"
                f"{labels['vol2']}: %{{customdata[2]:.2f}}%<br>"
                f"{labels['vol3']}: %{{customdata[3]:.2f}}%<br>"
                f"{labels['vol4']}: %{{customdata[4]:.2f}}%<extra></extra>"
            )
        ))
        
        # Ajout des lignes de connexion
        for _, row in df.iterrows():
            fig.add_trace(go.Scatter3d(
                x=[row['x_up'], row['x_lp']],
                y=[row['y_up'], row['y_lp']],
                z=[row['z_up'], row['z_lp']],
                mode='lines',
                line=dict(color='gray', width=1),
                showlegend=False,
                hoverinfo='none'
            ))
        
        # Configuration du layout avec les nouveaux axes
        fig.update_layout(
            scene=dict(
                xaxis=dict(title=f"% {labels['vol4']} (X)"),
                yaxis=dict(title=f"% {labels['vol3']} (Y)"),
                zaxis=dict(title=f"% {labels['vol2']} (Z)"),
                aspectmode='cube',
                camera=dict(
                    up=dict(x=0, y=0, z=1),
                    center=dict(x=0, y=0, z=0),
                    eye=dict(x=1.25, y=1.25, z=1.25)
                )
            ),
            title=f"Quaternary Phase Diagram: {labels['vol1']} / {labels['vol2']} / {labels['vol3']} / {labels['vol4']}",
            width=1200,
            height=800,
            hovermode='closest',
            margin=dict(l=60, r=60, t=80, b=60)
        )

        # Affichage en deux colonnes
        col1, col2 = st.columns([0.7, 0.3])

        with col1:
            st.plotly_chart(fig, use_container_width=True)

            # Affichage des données brutes
            st.subheader("Raw Data")
            with st.expander("View complete dataset", expanded=False):
                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True,
                    height=400,
                    column_config={
                        "Number": st.column_config.NumberColumn("Numéro", format="%d"),
                        "%Vol1 - UP": st.column_config.NumberColumn("UP Vol1", format="%.2f"),
                        "%Vol2 - UP": st.column_config.NumberColumn("UP Vol2", format="%.2f"),
                        "%Vol3 - UP": st.column_config.NumberColumn("UP Vol3", format="%.2f"),
                        "%Vol4 - UP": st.column_config.NumberColumn("UP Vol4", format="%.2f"),
                        "%Vol1 - LP": st.column_config.NumberColumn("LP Vol1", format="%.2f"),
                        "%Vol2 - LP": st.column_config.NumberColumn("LP Vol2", format="%.2f"),
                        "%Vol3 - LP": st.column_config.NumberColumn("LP Vol3", format="%.2f"),
                        "%Vol4 - LP": st.column_config.NumberColumn("LP Vol4", format="%.2f")
                    }
                )

        with col2:
            # Sélection interactive
            st.subheader("Select a point")
            selected_number = st.selectbox(
                "Select by number", 
                df['Number'].unique(),
                index=list(df['Number']).index(selected_number) if selected_number in df['Number'].values else 0
            )
            
            selected_row = df[df['Number'] == selected_number].iloc[0]
            phase_data = create_phase_display(selected_row, labels, is_quaternary=True)
            
            # Affichage des compositions
            st.subheader("Phase Compositions")
            
            for phase, color in [('UP', 'red'), ('LP', 'blue')]:
                with st.expander(f"Phase {phase}", expanded=True):
                    st.markdown(f"""
                    **{labels['vol1']}:** {phase_data[phase]['vol1']:.2f}%  
                    **{labels['vol2']}:** {phase_data[phase]['vol2']:.2f}%  
                    **{labels['vol3']}:** {phase_data[phase]['vol3']:.2f}%  
                    **{labels['vol4']}:** {phase_data[phase]['vol4']:.2f}%
                    """)

            # Bouton pour explorer dans KD Database
            if st.button("Find in KD Database"):
                # Vérifier si le système existe dans KDDB
                try:
                    kddb_sheets = load_excel_sheets(EXCEL_PATH)
                    system_name = selected_sheet
                    
                    if system_name in kddb_sheets:
                        st.session_state.current_page = "kddb"
                        st.session_state.search_query = system_name
                        st.session_state.search_triggered = True
                        st.rerun()
                    else:
                        st.warning(f"No matching system found in KD Database for {system_name}")
                except Exception as e:
                    st.error(f"Error accessing KD Database: {str(e)}")

    except Exception as e:
        st.error(f"Error processing data: {str(e)}")
    
    # Bouton de retour
    if st.button("Back to Home", key="dbdq_back"):
        st.session_state.current_page = "home"
        st.rerun()

# =============================================
# Application principale
# =============================================

def main():
    # Initialisation de l'état
    if 'current_page' not in st.session_state:
        st.session_state.current_page = "home"
    if 'search_triggered' not in st.session_state:
        st.session_state.search_triggered = False
    if 'search_query' not in st.session_state:
        st.session_state.search_query = ""
    
    # Navigation
    with st.sidebar:
        st.title("Navigation")
        if st.button("Home"):
            st.session_state.current_page = "home"
            st.rerun()
        if st.button("KD Database Explorer"):
            st.session_state.current_page = "kddb"
            st.rerun()
        if st.button("Ternary Diagrams"):
            st.session_state.current_page = "dbdt"
            st.rerun()
        if st.button("Quaternary Diagrams"):
            st.session_state.current_page = "dbdq"
            st.rerun()
    
    # Router vers la page active
    if st.session_state.current_page == "home":
        show_home_page()
    elif st.session_state.current_page == "kddb":
        show_kddb_page()
    elif st.session_state.current_page == "dbdt":
        show_dbdt_page()
    elif st.session_state.current_page == "dbdq":
        show_dbdq_page()

if __name__ == "__main__":
    main()
