import pandas as pd
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense, Embedding, Concatenate, Dropout, Flatten
from sklearn.preprocessing import StandardScaler, LabelEncoder
from sklearn.model_selection import train_test_split
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import joblib

class KDPredictor:
    def __init__(self, fingerprint_bits=2048, fingerprint_radius=2):
        self.fingerprint_bits = fingerprint_bits
        self.fingerprint_radius = fingerprint_radius
        self.solvent_encoder = LabelEncoder()
        self.composition_encoder = LabelEncoder()  # Encodage global principal
        self.kd_scaler = StandardScaler()
        self.model = None
        self.is_trained = False
        self.valid_combinations = {}
        self.solvent_composition_map = {}  # Pour la validation seulement
        
    def smiles_to_fingerprint(self, smiles):
        """Convertit un SMILES en fingerprint moléculaire"""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Attention: SMILES invalide: {smiles}")
                return np.zeros(self.fingerprint_bits)
            fingerprint = AllChem.GetMorganFingerprintAsBitVect(
                mol, self.fingerprint_radius, nBits=self.fingerprint_bits
            )
            return np.array(fingerprint)
        except Exception as e:
            print(f"Erreur avec SMILES {smiles}: {e}")
            return np.zeros(self.fingerprint_bits)
    
    def extract_valid_combinations(self, data):
        """Extrait les combinaisons solvant-composition valides du dataset"""
        valid_combinations = {}
        
        for solvent in data['Système de solvant'].unique():
            compositions = data[data['Système de solvant'] == solvent]['Composition'].unique()
            valid_combinations[solvent] = list(compositions)
        
        return valid_combinations
    
    def load_and_preprocess_data(self, csv_path, sep=';'):
        """Charge et prétraite les données depuis le fichier CSV"""
        # Chargement des données
        data = pd.read_csv(csv_path, sep=sep)
        print(f"Données chargées: {len(data)} entrées")
        print(f"Colonnes: {data.columns.tolist()}")
        
        # Vérification des données manquantes
        print("\nVérification des données manquantes:")
        print(data.isnull().sum())
        
        # Nettoyage des données
        data = data.dropna()
        print(f"Données après nettoyage: {len(data)} entrées")
        
        # Extraction des combinaisons valides (pour la validation seulement)
        self.valid_combinations = self.extract_valid_combinations(data)
        print(f"\n🔍 Combinaisons valides trouvées:")
        for solvent, compositions in self.valid_combinations.items():
            print(f"   {solvent}: {len(compositions)} compositions")
        
        # Conversion des SMILES en fingerprints
        print("\nConversion des SMILES en fingerprints...")
        X_smiles = np.array([self.smiles_to_fingerprint(smiles) for smiles in data['Smiles']])
        
        # Encodage des systèmes de solvant (GLOBAL)
        print("Encodage des systèmes de solvant...")
        X_solvent = self.solvent_encoder.fit_transform(data['Système de solvant'])
        
        # Encodage des compositions (GLOBAL - pour l'entraînement)
        print("Encodage GLOBAL des compositions...")
        X_composition = self.composition_encoder.fit_transform(data['Composition'])
        
        # Création du mapping pour la validation (séparé de l'entraînement)
        print("Création du mapping de validation...")
        self.solvent_composition_map = {}
        for solvent in data['Système de solvant'].unique():
            compositions = data[data['Système de solvant'] == solvent]['Composition'].unique()
            # On utilise l'encodage GLOBAL pour mapper
            encoded_compositions = self.composition_encoder.transform(compositions)
            self.solvent_composition_map[solvent] = {
                'compositions': list(compositions),
                'encoded_compositions': list(encoded_compositions),
                'mapping': dict(zip(compositions, encoded_compositions))
            }
        
        # Préparation de la target (KD)
        y_kd = data['KD'].values.reshape(-1, 1)
        y_kd_scaled = self.kd_scaler.fit_transform(y_kd).flatten()
        
        print(f"\nRésumé du prétraitement:")
        print(f"- Fingerprints: {X_smiles.shape}")
        print(f"- Solvants uniques: {len(self.solvent_encoder.classes_)}")
        print(f"- Compositions globales uniques: {len(self.composition_encoder.classes_)}")
        print(f"- Plage des KD: {data['KD'].min():.3f} à {data['KD'].max():.3f}")
        
        # Afficher le mapping des compositions par solvant
        print(f"\n📊 Mapping des compositions par solvant (pour validation):")
        for solvent, mapping_info in self.solvent_composition_map.items():
            print(f"   {solvent}: {mapping_info['mapping']}")
        
        return {
            'smiles': X_smiles,
            'solvent': X_solvent,
            'composition': X_composition  # Utilise l'encodage GLOBAL
        }, y_kd_scaled, data['KD'].values, data
    
    def build_model(self, n_solvents, n_compositions):
        """Construit le modèle de réseau de neurones"""
        # Dimensions
        fingerprint_dim = self.fingerprint_bits
        
        # Inputs
        smiles_input = Input(shape=(fingerprint_dim,), name='smiles')
        solvent_input = Input(shape=(1,), name='solvent')
        composition_input = Input(shape=(1,), name='composition')
        
        # Embedding pour les solvants
        solvent_embed = Embedding(n_solvents, 10)(solvent_input)
        solvent_embed = Flatten()(solvent_embed)
        
        # Embedding pour les compositions (GLOBAL)
        composition_embed = Embedding(n_compositions, 10)(composition_input)
        composition_embed = Flatten()(composition_embed)
        
        # Concatenation de toutes les features
        concatenated = Concatenate()([smiles_input, solvent_embed, composition_embed])
        
        # Couches cachées (adaptées à la taille de votre dataset)
        x = Dense(512, activation='relu')(concatenated)
        x = Dropout(0.4)(x)
        x = Dense(256, activation='relu')(x)
        x = Dropout(0.3)(x)
        x = Dense(128, activation='relu')(x)
        x = Dropout(0.2)(x)
        x = Dense(64, activation='relu')(x)
        
        # Couche de sortie
        output = Dense(1, activation='linear', name='kd_prediction')(x)
        
        # Construction du modèle
        model = Model(
            inputs=[smiles_input, solvent_input, composition_input],
            outputs=output
        )
        
        # Compilation
        model.compile(
            optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
            loss='mse',
            metrics=['mae', 'mse']
        )
        
        return model

    def prepare_data_for_training(self, X, y, test_size=0.2, random_state=42):
        """Prépare les données pour l'entraînement avec un split cohérent"""
        n_samples = len(y)
        indices = np.arange(n_samples)
        
        # Split des indices
        train_idx, test_idx = train_test_split(
            indices, test_size=test_size, random_state=random_state
        )
        
        # Split des données en utilisant les mêmes indices pour toutes les features
        X_train = {
            'smiles': X['smiles'][train_idx],
            'solvent': X['solvent'][train_idx],
            'composition': X['composition'][train_idx]
        }
        X_test = {
            'smiles': X['smiles'][test_idx],
            'solvent': X['solvent'][test_idx],
            'composition': X['composition'][test_idx]
        }
        
        y_train = y[train_idx]
        y_test = y[test_idx]
        
        return X_train, X_test, y_train, y_test
    
    def train(self, csv_path, test_size=0.2, random_state=42, epochs=100, batch_size=32):
        """Entraîne le modèle sur les données"""
        try:
            # Chargement et prétraitement des données
            X, y_scaled, y_original, original_data = self.load_and_preprocess_data(csv_path)
            
            # Vérification du nombre d'échantillons
            n_samples = len(y_scaled)
            print(f"\nNombre total d'échantillons: {n_samples}")
            
            # Ajustement de la taille du test set
            if n_samples < 100:
                test_size = 0.1  # Plus de données pour l'entraînement
                print(f"Dataset de taille moyenne, test_size ajusté à: {test_size}")
            
            # Utilisation de la méthode pour préparer les données
            X_train, X_test, y_train, y_test = self.prepare_data_for_training(
                X, y_scaled, test_size=test_size, random_state=random_state
            )
            
            print(f"Split des données:")
            print(f"- Entraînement: {len(y_train)} échantillons")
            print(f"- Test: {len(y_test)} échantillons")
            
            # Construction du modèle
            n_solvents = len(self.solvent_encoder.classes_)
            n_compositions = len(self.composition_encoder.classes_)
            
            print(f"📊 Dimensions du modèle:")
            print(f"- Solvants: {n_solvents}")
            print(f"- Compositions globales: {n_compositions}")
            
            self.model = self.build_model(n_solvents, n_compositions)
            
            print("\nArchitecture du modèle:")
            self.model.summary()
            
            # Ajustement de la batch_size
            actual_batch_size = min(batch_size, len(y_train))
            if actual_batch_size != batch_size:
                print(f"Batch_size ajusté à: {actual_batch_size}")
            
            # Callbacks
            callbacks = [
                tf.keras.callbacks.EarlyStopping(
                    monitor='val_loss',
                    patience=15,
                    restore_best_weights=True
                ),
                tf.keras.callbacks.ReduceLROnPlateau(
                    monitor='val_loss',
                    factor=0.5,
                    patience=8,
                    min_lr=1e-7
                )
            ]
            
            # Entraînement
            print("\nDébut de l'entraînement...")
            history = self.model.fit(
                X_train, y_train,
                epochs=epochs,
                batch_size=actual_batch_size,
                validation_data=(X_test, y_test),
                callbacks=callbacks,
                verbose=1
            )
            
            # Évaluation
            print("\nÉvaluation sur le jeu de test:")
            test_loss = self.model.evaluate(X_test, y_test, verbose=0)
            print(f"Loss: {test_loss[0]:.4f}, MAE: {test_loss[1]:.4f}")
            
            # Calcul de la MAE sur les données originales
            y_pred_scaled = self.model.predict(X_test, verbose=0)
            y_pred_original = self.kd_scaler.inverse_transform(y_pred_scaled.reshape(-1, 1)).flatten()
            y_test_original = self.kd_scaler.inverse_transform(y_test.reshape(-1, 1)).flatten()
            mae_original = np.mean(np.abs(y_pred_original - y_test_original))
            print(f"MAE sur KD original: {mae_original:.4f}")
            
            self.is_trained = True
            
            return history, X_test, y_test, y_original, original_data
            
        except Exception as e:
            print(f"❌ Erreur lors de l'entraînement: {e}")
            import traceback
            traceback.print_exc()
            raise
    
    def save_model(self, filepath, data=None):
        """Sauvegarde le modèle, les preprocesseurs et les combinaisons valides"""
        if not self.is_trained:
            raise ValueError("Le modèle doit être entraîné avant de être sauvegardé")
        
        # Sauvegarde du modèle Keras
        self.model.save(f'{filepath}_model.h5')
        
        # Sauvegarde des preprocesseurs
        joblib.dump({
            'solvent_encoder': self.solvent_encoder,
            'composition_encoder': self.composition_encoder,  # Encodage GLOBAL
            'kd_scaler': self.kd_scaler,
            'fingerprint_bits': self.fingerprint_bits,
            'fingerprint_radius': self.fingerprint_radius,
            'solvent_composition_map': self.solvent_composition_map  # Pour validation
        }, f'{filepath}_preprocessors.pkl')
        
        # Sauvegarde des combinaisons valides
        joblib.dump(self.valid_combinations, f'{filepath}_combinations.pkl')
        
        print(f"💾 Modèle sauvegardé: {filepath}_model.h5")
        print(f"💾 Preprocesseurs sauvegardés: {filepath}_preprocessors.pkl")
        print(f"💾 Combinaisons valides sauvegardées: {filepath}_combinations.pkl")
        print(f"📊 Solvants: {len(self.solvent_encoder.classes_)}")
        print(f"📊 Compositions globales: {len(self.composition_encoder.classes_)}")

def plot_training_history(history):
    """Visualise l'historique d'entraînement"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    
    # Plot loss
    ax1.plot(history.history['loss'], label='Train Loss')
    ax1.plot(history.history['val_loss'], label='Val Loss')
    ax1.set_title('Model Loss')
    ax1.set_xlabel('Epoch')
    ax1.set_ylabel('Loss')
    ax1.legend()
    
    # Plot MAE
    ax2.plot(history.history['mae'], label='Train MAE')
    ax2.plot(history.history['val_mae'], label='Val MAE')
    ax2.set_title('Model MAE')
    ax2.set_xlabel('Epoch')
    ax2.set_ylabel('MAE')
    ax2.legend()
    
    plt.tight_layout()
    plt.show()

# Exemple d'utilisation
if __name__ == "__main__":
    # Initialisation du prédicteur
    predictor = KDPredictor()
    
    try:
        # Entraînement du modèle
        print("🚀 Démarrage de l'entraînement...")
        history, X_test, y_test, y_original, original_data = predictor.train(
            csv_path='data.csv',
            test_size=0.2,
            epochs=100,
            batch_size=32
        )
        
        # Visualisation de l'entraînement
        plot_training_history(history)
        
        # Sauvegarde du modèle
        predictor.save_model('kd_predictor_model', original_data)
        
        print("\n✅ Entraînement terminé avec succès!")
        
    except Exception as e:
        print(f"\n❌ Erreur lors de l'entraînement: {e}")