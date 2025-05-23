import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import sys
from rdkit import Chem
from rdkit.Chem import Draw

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

def show_kddb_page():
    """Page KD Database Explorer - Version avec sélection par ligne"""
    st.title("KD Database Explorer")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(EXCEL_PATH)
    if not sheet_names:
        return
        
    # Zone de recherche
    col1, col2 = st.columns([0.7, 0.3])
    with col1:
        search_query = st.text_input(
            "Enter a compound name",
            key="search_input",
            placeholder="Search"
        )
    with col2:
        st.write("")  # Pour l'alignement
        if st.button("Search", key="search_button"):
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
                required_cols = ['Compound', 'Log KD', 'System', 'Composition', 'SMILES']
                additional_cols = ['Log P (Pubchem)', 'Log P (COSMO-RS)']
                
                # Vérification des colonnes disponibles
                available_cols = [col for col in required_cols + additional_cols if col in df.columns]
                
                if len([col for col in required_cols if col in available_cols]) < len(required_cols):
                    st.error("Les colonnes requises ne sont pas toutes présentes dans la feuille.")
                else:
                    # Création d'un dataframe pour la sélection
                    selection_df = df[available_cols].copy()
                    selection_df.insert(0, 'Select', False)
                    
                    # Configuration des colonnes pour st.data_editor
                    column_config = {
                        "Select": st.column_config.CheckboxColumn(
                            "Select",
                            help="Select a line by checking a box",
                            default=False,
                            required=True
                        ),
                        "Composition": st.column_config.TextColumn("Composition"),
                    }
                    
                    # Affichage du tableau avec case à cocher
                    edited_df = st.data_editor(
                        selection_df,
                        column_config=column_config,
                        use_container_width=True,
                        hide_index=True,
                        disabled=[col for col in available_cols],
                        key="data_editor"
                    )
                    
                    # Récupération de la ligne sélectionnée
                    selected_rows = edited_df[edited_df['Select']]
                    
                    # Afficher la structure au-dessus si selectionnée
                    if not selected_rows.empty:
                        selected_row = selected_rows.iloc[0]
                        selected_composition = selected_row['Composition']
                        smiles = selected_row['SMILES']

                        st.markdown("### Structure du composé sélectionné")
                        if pd.notna(smiles) and smiles.strip() != "":
                            mol = Chem.MolFromSmiles(smiles)
                            if mol is not None:
                                img = Draw.MolToImage(mol)
                                st.image(img, use_column_width=False, width=250)
                            else:
                                st.warning("Impossible de générer la structure à partir du SMILES.")
                        else:
                            st.info("SMILES non disponible pour ce composé.")
                    else:
                        st.info("Veuillez sélectionner une ligne pour afficher la structure du composé.")
                    
                    # Affichage du tableau sous l'image
                    st.subheader("Données des composés")
                    st.dataframe(edited_df[available_cols], use_container_width=True, hide_index=True)

                    # Si une ligne est sélectionnée, charger et afficher les diagrammes
                    if not selected_rows.empty:
                        selected_row = selected_rows.iloc[0]
                        system_name = selected_row['System']
                        selected_composition = selected_row['Composition']
                        is_quaternary = st.checkbox("Afficher en diagramme quaternaire", key="quaternary_check")
                        
                        try:
                            target_file = DBDQ_PATH if is_quaternary else DBDT_PATH
                            all_sheets = pd.read_excel(target_file, sheet_name=None)
                            
                            if system_name not in all_sheets:
                                st.error(f"Aucune donnée trouvée pour le système {system_name}")
                            else:
                                df_system = all_sheets[system_name]
                                df_system['Composition'] = df_system['Composition'].astype(str)
                                selected_composition_str = str(selected_composition)
                                df_filtered = df_system[df_system['Composition'] == selected_composition_str]
                                
                                if df_filtered.empty:
                                    st.error(f"Aucune donnée trouvée pour la composition {selected_composition} dans le système {system_name}")
                                    st.write("Compositions disponibles:", df_system['Composition'].unique())
                                else:
                                    if is_quaternary:
                                        show_quaternary_diagram(df_system, df_filtered, system_name, selected_composition)
                                    else:
                                        show_ternary_diagram(df_system, df_filtered, system_name, selected_composition)
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
    
    try:
        # Chargement des noms de feuilles
        sheet_names = load_excel_sheets(DBDT_PATH)
        if not sheet_names:
            st.warning("Aucune feuille trouvée dans le fichier DBDT.xlsx")
            return
        
        # Sélection de la feuille
        selected_sheet = st.selectbox(
            "Sélectionnez un système",
            sheet_names,
            key="dbdt_sheet_select"
        )
        
        # Chargement des données
        try:
            df = pd.read_excel(DBDT_PATH, sheet_name=selected_sheet)
            
            # Vérifier que les colonnes nécessaires existent
            required_cols = ['Composition', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', 
                            '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP']
            
            missing_cols = [col for col in required_cols if col not in df.columns]
            if missing_cols:
                st.error(f"Colonnes manquantes: {', '.join(missing_cols)}")
                return
            
            # Nettoyage des données
            df = df.dropna(subset=['%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', 
                                 '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP'])
            
            # Récupérer les noms des solvants (première ligne après l'en-tête)
            lib_row = pd.read_excel(DBDT_PATH, sheet_name=selected_sheet, header=None).iloc[1]
            labels = {
                'vol1': lib_row[0],
                'vol2': lib_row[1],
                'vol3': lib_row[2],
                'sheet': selected_sheet
            }
            
            # Création du graphique
            fig = go.Figure()

            # Courbe UP
            fig.add_trace(go.Scatter(
                x=df['%Vol3 - UP'],
                y=df['%Vol2 - UP'],
                mode='lines+markers',
                name='UP',
                line=dict(color='black', width=2),
                marker=dict(color='red', size=8, symbol='circle'),
                customdata=df[['Composition', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP']].values,
                hovertemplate=(
                    f'<b>Composition</b>: %{{customdata[0]}}<br>'
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
                line=dict(color='black', width=2),
                marker=dict(color='blue', size=8, symbol='circle'),
                customdata=df[['Composition', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP']].values,
                hovertemplate=(
                    f'<b>Composition</b>: %{{customdata[0]}}<br>'
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
                    line=dict(color='blue', width=2),
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

            # Configuration des axes
            fig.update_layout(
                xaxis=dict(
                    title=f'% {labels["vol3"]} (%Vol3)',
                    range=[0, 1],
                    dtick=0.1,
                    showgrid=True,
                    gridwidth=1.5,
                    gridcolor='#666666'
                ),
                yaxis=dict(
                    title=f'% {labels["vol2"]} (%Vol2)',
                    range=[0, 1],
                    dtick=0.1,
                    showgrid=True,
                    gridwidth=1.5,
                    gridcolor='#666666'
                ),
                title=f"Ternary Phase Diagram: {labels['vol1']} / {labels['vol2']} / {labels['vol3']}",
                width=1400,
                height=700,
                hovermode='closest'
            )

            # Affichage en deux colonnes
            col1, col2 = st.columns([0.8, 0.2])

            with col1:
                st.plotly_chart(fig, use_container_width=True)

                # Affichage des données brutes
                st.subheader("Raw Data")
                st.dataframe(
                    df,
                    use_container_width=True,
                    hide_index=True,
                    height=400
                )

            with col2:
                # Sélection interactive
                st.subheader("Select a composition")
                selected_composition = st.selectbox(
                    "Select by composition", 
                    df['Composition'].unique(),
                    key="dbdt_composition_select"
                )
                
                selected_row = df[df['Composition'] == selected_composition].iloc[0]
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
            st.error(f"Erreur lors du chargement des données: {str(e)}")
            st.error("Veuillez vérifier que le fichier DBDT.xlsx est correctement formaté.")

    except Exception as e:
        st.error(f"Erreur critique: {str(e)}")
    
    # Bouton de retour
    if st.button("Retour à l'accueil", key="dbdt_back"):
        st.session_state.current_page = "home"
        st.rerun()
        
def show_ternary_diagram(df_system, df_filtered, system_name, selected_composition):
    """Affiche un diagramme ternaire"""
    # Récupérer les labels des solvants depuis la première ligne
    lib_row = pd.read_excel(DBDT_PATH, sheet_name=system_name, header=None).iloc[1]
    labels = {
        'vol1': lib_row[0],
        'vol2': lib_row[1],
        'vol3': lib_row[2],
        'sheet': system_name
    }
    
    # Vérifier que les colonnes nécessaires existent
    required_cols = [f'%Vol1 - UP', f'%Vol2 - UP', f'%Vol3 - UP',
                    f'%Vol1 - LP', f'%Vol2 - LP', f'%Vol3 - LP']
    
    for col in required_cols:
        if col not in df_system.columns:
            st.error(f"Colonne manquante: {col}")
            return
    
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
                df_system['Composition'].astype(str),
                df_system[f'%Vol1 - {phase}'],
                df_system[f'%Vol2 - {phase}'],
                df_system[f'%Vol3 - {phase}']
            ), axis=-1),
            hovertemplate=(
                "<b>Composition</b>: %{customdata[0]}<br>"
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

def show_quaternary_diagram(df_system, df_filtered, system_name, selected_composition):
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
        customdata=df_system[['Composition', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', '%Vol4 - UP']].values,
        hovertemplate=(
            f"<b>Phase UP</b><br>"
            f"<b>Composition</b>: %{{customdata[0]}}<br>"
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
        customdata=df_system[['Composition', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP', '%Vol4 - LP']].values,
        hovertemplate=(
            f"<b>Phase LP</b><br>"
            f"<b>Composition</b>: %{{customdata[0]}}<br>"
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

def show_kddb_page():
    """Page KD Database Explorer - Version avec sélection par ligne"""
    st.title("KD Database Explorer")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(EXCEL_PATH)
    if not sheet_names:
        return
        
    # Zone de recherche
    col1, col2 = st.columns([0.7, 0.3])
    with col1:
        search_query = st.text_input(
            "Enter a compound name",
            key="search_input",
            placeholder="Search"
        )
    with col2:
        st.write("")  # Pour l'alignement
        if st.button("Search", key="search_button"):
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
                required_cols = ['Compound', 'Log KD', 'System', 'Composition']
                additional_cols = ['Log P (Pubchem)', 'Log P (COSMO-RS)']
                
                # Vérification des colonnes disponibles
                available_cols = [col for col in required_cols + additional_cols if col in df.columns]
                
                if len([col for col in required_cols if col in available_cols]) < len(required_cols):
                    st.error("Les colonnes requises ne sont pas toutes présentes dans la feuille.")
                else:
                    # Sélection interactive par ligne
                    st.subheader("Sélectionnez une entrée")
                    
                    # Création d'un dataframe pour la sélection
                    selection_df = df[available_cols].copy()
                    selection_df.insert(0, 'Select', False)
                    
                    # Configuration des colonnes pour st.data_editor
                    column_config = {
                        "Select": st.column_config.CheckboxColumn(
                            "Select",
                            help="Select a line by checking a box",
                            default=False,
                            required=True
                        ),
                        "Composition": st.column_config.TextColumn("Composition"),
                    }
                    
                    # Affichage du tableau avec case à cocher
                    edited_df = st.data_editor(
                        selection_df,
                        column_config=column_config,
                        use_container_width=True,
                        hide_index=True,
                        disabled=[col for col in available_cols],
                        key="data_editor"
                    )
                    
                    # Récupération de la ligne sélectionnée
                    selected_rows = edited_df[edited_df['Select']]
                    
                    if not selected_rows.empty:
                        # Ne garder que la première sélection si plusieurs cases cochées
                        selected_row = selected_rows.iloc[0]
                        system_name = selected_row['System']
                        selected_composition = str(selected_row['Composition'])  # Convertir en string pour la comparaison
                        
                        # Déterminer si c'est un système ternaire ou quaternaire
                        is_quaternary = st.checkbox("Afficher en diagramme quaternaire", key="quaternary_check")
                        
                        # Chargement des données du système correspondant
                        try:
                            # Charger toutes les feuilles du fichier approprié
                            target_file = DBDQ_PATH if is_quaternary else DBDT_PATH
                            all_sheets = pd.read_excel(target_file, sheet_name=None)
                            
                            if system_name not in all_sheets:
                                st.error(f"Aucune donnée trouvée pour le système {system_name}")
                            else:
                                df_system = all_sheets[system_name]
                                # Convertir la colonne Composition en string pour la comparaison
                                df_system['Composition'] = df_system['Composition'].astype(str)
                                df_filtered = df_system[df_system['Composition'] == selected_composition]
                                
                                if df_filtered.empty:
                                    st.error(f"Aucune donnée trouvée pour la composition {selected_composition} dans le système {system_name}")
                                    st.write("Compositions disponibles:", df_system['Composition'].unique())
                                else:
                                    if is_quaternary:
                                        show_quaternary_diagram(df_system, df_filtered, system_name, selected_composition)
                                    else:
                                        show_ternary_diagram(df_system, df_filtered, system_name, selected_composition)
                        
                        except Exception as e:
                            st.error(f"Erreur lors du chargement du système {system_name}: {str(e)}")
                    else:
                        st.info("Veuillez sélectionner une ligne dans le tableau pour afficher les détails")
            
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

def show_dbdq_page():
    """Page Quaternary Phase Diagrams"""
    st.title("Quaternary Phase Diagrams")
    
    # Chargement des noms de feuilles
    sheet_names = load_excel_sheets(DBDQ_PATH)
    if not sheet_names:
        return
    
    # Gestion des arguments passés
    initial_sheet = None
    selected_composition = None
    if len(sys.argv) > 2:
        initial_sheet = sys.argv[1]
        selected_composition = sys.argv[2]
    
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
            customdata=df[['Composition', '%Vol1 - UP', '%Vol2 - UP', '%Vol3 - UP', '%Vol4 - UP']].values,
            hovertemplate=(
                f"<b>Phase UP</b><br>"
                f"<b>Composition</b>: %{{customdata[0]}}<br>"
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
            customdata=df[['Composition', '%Vol1 - LP', '%Vol2 - LP', '%Vol3 - LP', '%Vol4 - LP']].values,
            hovertemplate=(
                f"<b>Phase LP</b><br>"
                f"<b>Composition</b>: %{{customdata[0]}}<br>"
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
                        "Composition": st.column_config.TextColumn("Composition"),
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
            selected_composition = st.selectbox(
                "Select by composition", 
                df['Composition'].unique(),
                index=list(df['Composition']).index(selected_composition) if selected_composition in df['Composition'].values else 0
            )
            
            selected_row = df[df['Composition'] == selected_composition].iloc[0]
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
    
    # Navigation
    with st.sidebar:
        st.title("Navigation")
        if st.button("Home"):
            st.session_state.current_page = "home"
        if st.button("KD Database Explorer"):
            st.session_state.current_page = "kddb"
        if st.button("Ternary Diagrams"):
            st.session_state.current_page = "dbdt"
        if st.button("Quaternary Diagrams"):
            st.session_state.current_page = "dbdq"
    
    # Router vers la page active
    try:
        if st.session_state.current_page == "home":
            show_home_page()
        elif st.session_state.current_page == "kddb":
            show_kddb_page()
        elif st.session_state.current_page == "dbdt":
            show_dbdt_page()
        elif st.session_state.current_page == "dbdq":
            show_dbdq_page()
    except Exception as e:
        st.error(f"Une erreur est survenue: {str(e)}")
        st.session_state.current_page = "home"
        st.rerun()

if __name__ == "__main__":
    main()
