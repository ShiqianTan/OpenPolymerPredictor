# ============= app.py =============
"""
Main Flask application entry point
"""
from flask import Flask, request, jsonify, render_template
from models.predictor import MolecularPropertyPredictor
from rdkit import Chem
import config

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
        
        # Validate SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return jsonify({'error': 'Invalid SMILES string'}), 400
        
        # Predict properties
        predictions = predictor.predict(smiles)
        
        if predictions is None:
            return jsonify({'error': 'Could not generate predictions'}), 500
        
        return jsonify({
            'smiles': smiles,
            'predicted_properties': predictions
        })
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

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