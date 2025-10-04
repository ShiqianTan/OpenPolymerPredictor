# ============= models/predictor.py =============
"""
Molecular property prediction model
"""
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
import numpy as np
from sklearn.ensemble import RandomForestRegressor
import pickle
import os

class MolecularPropertyPredictor:
    """Random Forest based molecular property predictor"""
    
    def __init__(self, n_estimators=50, random_state=42, fingerprint_bits=2048):
        self.n_estimators = n_estimators
        self.random_state = random_state
        self.fingerprint_bits = fingerprint_bits
        
        self.models = {
            'molecular_weight': RandomForestRegressor(
                n_estimators=n_estimators, 
                random_state=random_state
            ),
            'logP': RandomForestRegressor(
                n_estimators=n_estimators, 
                random_state=random_state
            ),
            'tpsa': RandomForestRegressor(
                n_estimators=n_estimators, 
                random_state=random_state
            )
        }
        self.trained = False
        
    def smiles_to_fingerprint(self, smiles):
        """Convert SMILES to Morgan fingerprint"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, 2, nBits=self.fingerprint_bits
        )
        return np.array(fp)
    
    def calculate_properties(self, smiles):
        """Calculate actual molecular properties"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return {
            'molecular_weight': Descriptors.MolWt(mol),
            'logP': Descriptors.MolLogP(mol),
            'tpsa': Descriptors.TPSA(mol)
        }
    
    def get_sample_data(self):
        """Get sample SMILES for training"""
        return [
            'CCO', 'CC(C)O', 'CCCC', 'c1ccccc1', 'CC(=O)O',
            'CCN', 'CCC(=O)O', 'CCCCCC', 'c1ccc(O)cc1', 'CC(C)CC',
            'CCCCCCCC', 'CC(C)(C)O', 'c1ccccc1O', 'CCO', 'CCCCO',
            'c1ccccc1C', 'CC(=O)OC', 'CCCCCCO', 'c1ccc(N)cc1', 'CCCCCO',
            'CC(C)C', 'CCCCC', 'c1ccc(Cl)cc1', 'CCN(CC)CC', 'CCOC'
        ]
    
    def train_with_sample_data(self):
        """Train models with sample molecules"""
        sample_smiles = self.get_sample_data()
        
        X = []
        y_mw = []
        y_logp = []
        y_tpsa = []
        
        for smiles in sample_smiles:
            fp = self.smiles_to_fingerprint(smiles)
            props = self.calculate_properties(smiles)
            if fp is not None and props is not None:
                X.append(fp)
                y_mw.append(props['molecular_weight'])
                y_logp.append(props['logP'])
                y_tpsa.append(props['tpsa'])
        
        X = np.array(X)
        
        # Train models
        self.models['molecular_weight'].fit(X, y_mw)
        self.models['logP'].fit(X, y_logp)
        self.models['tpsa'].fit(X, y_tpsa)
        self.trained = True
        
        return len(X)
    
    def predict(self, smiles):
        """Predict properties for a given SMILES string"""
        if not self.trained:
            self.train_with_sample_data()
        
        fp = self.smiles_to_fingerprint(smiles)
        if fp is None:
            return None
        
        predictions = {
            'molecular_weight': float(self.models['molecular_weight'].predict([fp])[0]),
            'logP': float(self.models['logP'].predict([fp])[0]),
            'tpsa': float(self.models['tpsa'].predict([fp])[0])
        }
        return predictions
    
    def save_models(self, filepath='models/trained_models.pkl'):
        """Save trained models to disk"""
        os.makedirs(os.path.dirname(filepath), exist_ok=True)
        with open(filepath, 'wb') as f:
            pickle.dump(self.models, f)
    
    def load_models(self, filepath='models/trained_models.pkl'):
        """Load trained models from disk"""
        if os.path.exists(filepath):
            with open(filepath, 'rb') as f:
                self.models = pickle.load(f)
            self.trained = True
            return True
        return False