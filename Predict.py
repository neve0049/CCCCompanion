import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from sklearn.preprocessing import StandardScaler, LabelEncoder
from rdkit import Chem
from rdkit.Chem import AllChem
import joblib
import threading
from datetime import datetime

class KDPredictorGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("🧪 KD Prediction - UI")
        self.root.geometry("900x750")
        self.root.configure(bg='#f0f0f0')
        
        # Initialisation du prédicteur
        self.predictor = KDPredictor()
        self.model_loaded = False
        
        # Variables pour les menus déroulants
        self.solvent_var = tk.StringVar()
        self.composition_var = tk.StringVar()
        self.smiles_var = tk.StringVar()
        
        self.setup_ui()
        self.load_model_async()
    
    def setup_ui(self):
        """Configurating UI"""
        # Titre principal
        title_frame = tk.Frame(self.root, bg='#2c3e50', height=80)
        title_frame.pack(fill='x', padx=10, pady=10)
        title_frame.pack_propagate(False)
        
        title_label = tk.Label(
            title_frame,
            text="🧪 Prediction of partitioning coefficient (log KD)",
            font=('Arial', 16, 'bold'),
            fg='white',
            bg='#2c3e50'
        )
        title_label.pack(expand=True)
        
        # Frame principale
        main_frame = tk.Frame(self.root, bg='#f0f0f0')
        main_frame.pack(fill='both', expand=True, padx=20, pady=10)
        
        # Status du modèle
        self.status_frame = tk.Frame(main_frame, bg='#f0f0f0')
        self.status_frame.pack(fill='x', pady=(0, 20))
        
        self.status_label = tk.Label(
            self.status_frame,
            text="🔄 Loading model...",
            font=('Arial', 10),
            fg='orange',
            bg='#f0f0f0'
        )
        self.status_label.pack(side='left')
        
        # Frame de configuration
        config_frame = tk.LabelFrame(
            main_frame,
            text=" System configuration ",
            font=('Arial', 12, 'bold'),
            bg='#f0f0f0',
            padx=15,
            pady=15
        )
        config_frame.pack(fill='x', pady=(0, 20))
        
        # Sélection du solvant
        solvent_frame = tk.Frame(config_frame, bg='#f0f0f0')
        solvent_frame.pack(fill='x', pady=5)
        
        tk.Label(
            solvent_frame,
            text="1. Select a biphasic solvent system:",
            font=('Arial', 10, 'bold'),
            bg='#f0f0f0'
        ).pack(anchor='w')
        
        self.solvent_combo = ttk.Combobox(
            solvent_frame,
            textvariable=self.solvent_var,
            state="readonly",
            width=50,
            font=('Arial', 10)
        )
        self.solvent_combo.pack(fill='x', pady=5)
        self.solvent_combo.bind('<<ComboboxSelected>>', self.on_solvent_selected)
        
        # Sélection de la composition
        composition_frame = tk.Frame(config_frame, bg='#f0f0f0')
        composition_frame.pack(fill='x', pady=5)
        
        tk.Label(
            composition_frame,
            text="2. Select a composition:",
            font=('Arial', 10, 'bold'),
            bg='#f0f0f0'
        ).pack(anchor='w')
        
        self.composition_combo = ttk.Combobox(
            composition_frame,
            textvariable=self.composition_var,
            state="disabled",
            width=50,
            font=('Arial', 10)
        )
        self.composition_combo.pack(fill='x', pady=5)
        
        # Saisie du SMILES
        smiles_frame = tk.Frame(config_frame, bg='#f0f0f0')
        smiles_frame.pack(fill='x', pady=5)
        
        tk.Label(
            smiles_frame,
            text="3. Enter SMILES:",
            font=('Arial', 10, 'bold'),
            bg='#f0f0f0'
        ).pack(anchor='w')
        
        self.smiles_entry = tk.Entry(
            smiles_frame,
            textvariable=self.smiles_var,
            width=50,
            font=('Arial', 10),
            relief='solid',
            bd=1
        )
        self.smiles_entry.pack(fill='x', pady=5)
        self.smiles_entry.bind('<KeyRelease>', self.validate_smiles)
        
        self.smiles_status = tk.Label(
            smiles_frame,
            text="",
            font=('Arial', 9),
            bg='#f0f0f0'
        )
        self.smiles_status.pack(anchor='w')
        
        # Boutons d'action
        button_frame = tk.Frame(config_frame, bg='#f0f0f0')
        button_frame.pack(fill='x', pady=10)
        
        # Bouton de prédiction simple
        self.predict_button = tk.Button(
            button_frame,
            text="🎯 Predict KD for selected system",
            command=self.launch_prediction,
            font=('Arial', 10, 'bold'),
            bg='#3498db',
            fg='white',
            relief='raised',
            bd=2,
            state='disabled',
            width=40,
            height=1
        )
        self.predict_button.pack(side='left', padx=(0, 10))
        
        # Bouton de scan complet
        self.scan_button = tk.Button(
            button_frame,
            text="🔍 Find suitable Log KD for CPC (-1 < KD < 1)",
            command=self.launch_complete_scan,
            font=('Arial', 10, 'bold'),
            bg='#9b59b6',
            fg='white',
            relief='raised',
            bd=2,
            state='disabled',
            width=40,
            height=1
        )
        self.scan_button.pack(side='left', padx=(0, 10))
        
        # Bouton de réinitialisation
        self.reset_button = tk.Button(
            button_frame,
            text="🔄 Reinitialize",
            command=self.reset_interface,
            font=('Arial', 10, 'bold'),
            bg='#e67e22',
            fg='white',
            relief='raised',
            bd=2,
            width=12,
            height=1
        )
        self.reset_button.pack(side='left')
        
        # Zone de résultats
        results_frame = tk.LabelFrame(
            main_frame,
            text=" Prediction results ",
            font=('Arial', 12, 'bold'),
            bg='#f0f0f0',
            padx=15,
            pady=15
        )
        results_frame.pack(fill='both', expand=True)
        
        self.results_text = scrolledtext.ScrolledText(
            results_frame,
            wrap=tk.WORD,
            width=80,
            height=12,
            font=('Consolas', 9),
            relief='solid',
            bd=1
        )
        self.results_text.pack(fill='both', expand=True)
        self.results_text.config(state='disabled')
        
        # Barre de statut en bas
        self.status_bar = tk.Label(
            self.root,
            text="Ready - waiting for the loading of the model",
            relief='sunken',
            anchor='w',
            font=('Arial', 9)
        )
        self.status_bar.pack(side='bottom', fill='x')
    
    def load_model_async(self):
        """Charge le modèle en arrière-plan"""
        def load_task():
            try:
                success = self.predictor.load_model('kd_predictor_model')
                if success:
                    self.root.after(0, self.on_model_loaded)
                else:
                    self.root.after(0, self.on_model_error)
            except Exception as e:
                self.root.after(0, lambda: self.on_model_error(str(e)))
        
        thread = threading.Thread(target=load_task)
        thread.daemon = True
        thread.start()
    
    def on_model_loaded(self):
        """Callback quand le modèle est chargé"""
        self.model_loaded = True
        self.status_label.config(text="✅ Model successfully loaded", fg='green')
        
        # Mettre à jour les combobox
        solvents = self.predictor.get_available_solvents()
        self.solvent_combo['values'] = solvents
        
        if solvents:
            self.solvent_combo.set(solvents[0])
            self.on_solvent_selected()
        
        self.update_status_bar(f"Model loaded - {len(solvents)} systems available")
        self.predict_button.config(state='normal', bg='#27ae60')
        self.scan_button.config(state='normal', bg='#8e44ad')
    
    def on_model_error(self, error_msg=""):
        """Callback en cas d'erreur de chargement"""
        self.status_label.config(text="❌ Error when loading the model", fg='red')
        messagebox.showerror(
            "Error",
            f"Impossible to load the model.\n{error_msg}\n\nVerify the presence of the files:\n"
            "- kd_predictor_model_model.h5\n"
            "- kd_predictor_model_preprocessors.pkl\n" 
            "- kd_predictor_model_combinations.pkl"
        )
    
    def on_solvent_selected(self, event=None):
        """Quand un solvant est sélectionné"""
        solvent = self.solvent_var.get()
        if solvent and self.model_loaded:
            compositions = self.predictor.get_available_compositions_for_solvent(solvent)
            self.composition_combo['values'] = compositions
            self.composition_combo['state'] = 'readonly'
            
            if compositions:
                self.composition_combo.set(compositions[0])
            else:
                self.composition_combo.set('')
            
            self.update_status_bar(f"Solvent selected: {solvent} - {len(compositions)} compositions available")
    
    def validate_smiles(self, event=None):
        """Valide le SMILES en temps réel"""
        smiles = self.smiles_var.get().strip()
        
        if not smiles:
            self.smiles_status.config(text="", fg='black')
            return
        
        # Validation basique du SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            self.smiles_status.config(text="❌ Invalid SMILES", fg='red')
        else:
            self.smiles_status.config(text="✅ Valid SMILES", fg='green')
    
    def launch_prediction(self):
        """Lance la prédiction simple"""
        if not self.model_loaded:
            messagebox.showerror("Error", "The model did not load")
            return
        
        # Récupération des valeurs
        solvent = self.solvent_var.get()
        composition = self.composition_var.get()
        smiles = self.smiles_var.get().strip()
        
        # Validation
        if not all([solvent, composition, smiles]):
            messagebox.showwarning("Warning", "Please fill the parameters")
            return
        
        # Validation du SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            messagebox.showerror("Error", "SMILES is invalid")
            return
        
        # Désactiver les boutons pendant la prédiction
        self.predict_button.config(state='disabled', text="🔄 Calculating...")
        self.scan_button.config(state='disabled')
        self.update_status_bar("Prediction in process...")
        
        # Lancer la prédiction en arrière-plan
        threading.Thread(target=self.run_prediction, args=(smiles, solvent, composition), daemon=True).start()
    
    def launch_complete_scan(self):
        """Lance le scan complet de tous les systèmes"""
        if not self.model_loaded:
            messagebox.showerror("Error", "The model did not load")
            return
        
        smiles = self.smiles_var.get().strip()
        
        # Validation du SMILES
        if not smiles:
            messagebox.showwarning("Warning", "Please enter a SMILES")
            return
        
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            messagebox.showerror("Error", "SMILES is not valid")
            return
        
        # Confirmation (le scan peut prendre du temps)
        response = messagebox.askyesno(
            "Search for optimal system", 
            f"Do you wish to search for an optimal system for:\n{smiles}\n\n"
            f"Your compound will evaluated for all systems in every composition "
            f"and will only display those for -1 < KD < 1.\n\n"
            f"This operation might take a while."
        )
        
        if not response:
            return
        
        # Désactiver les boutons pendant le scan
        self.predict_button.config(state='disabled')
        self.scan_button.config(state='disabled', text="🔍 Scanning...")
        self.update_status_bar("Searching for optimal systems...")
        
        # Lancer le scan en arrière-plan
        threading.Thread(target=self.run_complete_scan, args=(smiles,), daemon=True).start()
    
    def run_prediction(self, smiles, solvent, composition):
        """Exécute la prédiction simple dans un thread séparé"""
        try:
            prediction = self.predictor.predict(smiles, solvent, composition)
            
            if prediction is not None:
                self.root.after(0, lambda: self.display_prediction_result(
                    smiles, solvent, composition, prediction
                ))
            else:
                self.root.after(0, lambda: self.display_prediction_error())
                
        except Exception as e:
            self.root.after(0, lambda: self.display_prediction_error(str(e)))
        
        finally:
            self.root.after(0, self.reset_buttons)
    
    def run_complete_scan(self, smiles):
        """Exécute le scan complet dans un thread séparé"""
        try:
            results = []
            total_combinations = 0
            valid_combinations = 0
            
            # Récupérer tous les solvants et compositions
            solvents = self.predictor.get_available_solvents()
            
            for solvent in solvents:
                compositions = self.predictor.get_available_compositions_for_solvent(solvent)
                total_combinations += len(compositions)
                
                for composition in compositions:
                    # Faire la prédiction pour chaque combinaison
                    prediction = self.predictor.predict(smiles, solvent, composition)
                    
                    if prediction is not None and -1 <= prediction <= 1:
                        results.append({
                            'solvent': solvent,
                            'composition': composition,
                            'kd': prediction
                        })
                        valid_combinations += 1
            
            self.root.after(0, lambda: self.display_scan_results(
                smiles, results, total_combinations, valid_combinations
            ))
            
        except Exception as e:
            self.root.after(0, lambda: self.display_scan_error(str(e)))
        
        finally:
            self.root.after(0, self.reset_buttons)
    
    def display_prediction_result(self, smiles, solvent, composition, prediction):
        """Affiche le résultat de la prédiction simple"""
        # Interprétation
        if prediction < -1:
            interpretation = "Affinity with aqueous phase"
            color = "#e74c3c"
        elif prediction < 1:
            interpretation = "Optimal partitioning"
            color = "#f39c12"
        else:
            interpretation = "Affinity with organic phase" 
            color = "#27ae60"
        
        result_text = f"""
{'='*60}
🧪 PRÉDICTION TERMINÉE - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
{'='*60}

🔬 MOLÉCULE:
   SMILES: {smiles}

💧 SYSTÈME:
   Solvant: {solvent}
   Composition: {composition}

📊 Results:
   Predicted log KD: {prediction:.4f}

💡 Interpretation:
   {interpretation}

📈 PLAGE DE KD:
   < -1.0 : Affinity with aqueous phase
   -1.0 - 1.0 : Good with partitioning
   > 1.0 : Affinity with organic phase

{'='*60}
"""
        
        self.display_in_results(result_text, interpretation, color)
        self.update_status_bar(f"Prédiction terminée - KD: {prediction:.4f}")
    
    def display_scan_results(self, smiles, results, total_combinations, valid_combinations):
        """Affiche les résultats du scan complet"""
        if not results:
            result_text = f"""
{'='*60}
🔍 Search for optimal system ended - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
{'='*60}

🔬 MOLÉCULE:
   SMILES: {smiles}

📊 Results of search:
   Compositions tested: {total_combinations}
   Compositions with -1 < KD < 1: {valid_combinations}

❌ No system found
   No composition for any system gave a satisfactory result

💡 SUGGESTIONS:
   • Try another compound
   • Manually select compositions and try to find a system near -1 or 1
   • Verify SMILES

{'='*60}
"""
            self.display_in_results(result_text)
            self.update_status_bar(f"Search finished - No result for -1 < Log KD < 1")
            return
        
        # Trier les résultats par KD (croissant)
        results.sort(key=lambda x: x['kd'])
        
        result_text = f"""
{'='*60}
🔍 Search completed - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
{'='*60}

🔬 MOLECULE:
   SMILES: {smiles}

📊 Results of search:
   Compositions tested: {total_combinations}
   Compositions with -1 < log KD < 1: {valid_combinations}

🎯 Optimal systems (sorted by log KD):
"""
        
        # Ajouter chaque résultat
        for i, result in enumerate(results, 1):
            kd = result['kd']
            # Color code basé sur la valeur de log KD
            if kd < -0.5:
                kd_color = "🟢"  # Très bas
            elif kd < 0:
                kd_color = "🟡"  # Bas
            elif kd < 0.5:
                kd_color = "🟠"  # Modéré
            else:
                kd_color = "🔴"  # Élevé mais dans la plage
            
            result_text += f"\n{kd_color} {i:2d}. {result['solvent']} + {result['composition']}"
            result_text += f"\n     KD = {kd:.4f}\n"
        
        result_text += f"""
{'='*60}
💡 INTERPRÉTATION:
   • Log KD between -1 et 1 indicates a good partitioning
   • Negative values indicate a preference for the aqueous phase
   • Positive values indicate a preference for the organic phase

🎯 RECOMMANDATIONS:
   • Lower system: {results[0]['solvent']} + {results[0]['composition']} (KD = {results[0]['kd']:.4f})
   • Système le plus haut: {results[-1]['solvent']} + {results[-1]['composition']} (KD = {results[-1]['kd']:.4f})

{'='*60}
"""
        
        self.display_in_results(result_text)
        self.update_status_bar(f"Search finished - {valid_combinations} systems found with -1 < log KD < 1")
    
    def display_scan_error(self, error_msg=""):
        """Affiche une erreur de scan"""
        error_text = f"""
{'='*60}
❌ ERREUR DE SCAN - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
{'='*60}

Can't compute the complete search.

{error_msg}

Verify that:
• SMILES is valid
• Model is correctly loaded
• The file are accessible
{'='*60}
"""
        self.display_in_results(error_text)
        messagebox.showerror("Error", "Cannot perform the search for optimal compositions")
    
    def display_in_results(self, text, highlight_text="", highlight_color="#000000"):
        """Affiche du texte dans la zone de résultats"""
        self.results_text.config(state='normal')
        self.results_text.delete(1.0, tk.END)
        self.results_text.insert(1.0, text)
        
        if highlight_text:
            start_idx = text.find(highlight_text)
            if start_idx != -1:
                end_idx = start_idx + len(highlight_text)
                self.results_text.tag_add("highlight", f"1.0+{start_idx}c", f"1.0+{end_idx}c")
                self.results_text.tag_config("highlight", foreground=highlight_color, font=('Consolas', 9, 'bold'))
        
        self.results_text.config(state='disabled')
        self.results_text.see(1.0)  # Scroll to top
    
    def reset_buttons(self):
        """Réactive les boutons"""
        self.predict_button.config(state='normal', text="🎯 Prediction for selected system and composition", bg='#27ae60')
        self.scan_button.config(state='normal', text="🔍 Search optimal compositions for CPC (-1 < log KD < 1)", bg='#8e44ad')
    
    def reset_interface(self):
        """Réinitialise l'interface"""
        self.smiles_var.set("")
        self.smiles_status.config(text="", fg='black')
        self.results_text.config(state='normal')
        self.results_text.delete(1.0, tk.END)
        self.results_text.config(state='disabled')
        self.update_status_bar("Interface réinitialisée")
    
    def update_status_bar(self, message):
        """Met à jour la barre de statut"""
        self.status_bar.config(text=f" {message}")

# Classe KDPredictor (identique)
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
            
        except Exception:
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

if __name__ == "__main__":
    root = tk.Tk()
    app = KDPredictorGUI(root)
    root.mainloop()