"""
Main Flask application entry point
"""
from flask import Flask, request, jsonify, render_template
from models.predictor import MolecularPropertyPredictor
from utils.validators import validate_smiles, generate_molecule_image
import config
import base64
from rdkit import Chem
from rdkit.Chem import Draw

app = Flask(__name__)
app.config.from_object(config)

# Initialize predictor
predictor = MolecularPropertyPredictor()

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/predict', methods=['POST'])
def predict():
    try:
        data = request.get_json()
        smiles = data.get('smiles', '')
        
        if not smiles:
            return jsonify({'error': 'No SMILES string provided'}), 400
        
        # Validate SMILES with detailed error reporting
        is_valid, mol, error_msg = validate_smiles(smiles)
        
        if not is_valid:
            return jsonify({'error': error_msg}), 400
        
        # Generate molecular structure image (heavy atoms only)
        img_base64 = generate_molecule_image(mol, remove_hydrogens=True)
        
        if img_base64 is None:
            return jsonify({'error': 'Could not generate molecular structure'}), 500
        
        # Predict properties
        predictions = predictor.predict(smiles)
        
        if predictions is None:
            return jsonify({'error': 'Could not generate predictions'}), 500
        
        return jsonify({
            'smiles': smiles,
            'canonical_smiles': Chem.MolToSmiles(mol),
            'predicted_properties': predictions,
            'molecule_image': img_base64
        })
    
    except Exception as e:
        return jsonify({'error': f'Unexpected error: {str(e)}'}), 500

@app.route('/health', methods=['GET'])
def health_check():
    """Health check endpoint for monitoring"""
    return jsonify({
        'status': 'healthy',
        'model_trained': predictor.trained
    })

if __name__ == '__main__':
    print("Training model with sample data...")
    predictor.train_with_sample_data()
    print("Model trained successfully!")
    
    app.run(
        debug=app.config['DEBUG'],
        host=app.config['HOST'],
        port=app.config['PORT']
    )