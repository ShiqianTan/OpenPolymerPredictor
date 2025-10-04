# Molecular Property Predictor

A scalable Flask application for predicting molecular properties from SMILES strings using Random Forest regression with robust SMILES validation and molecular structure visualization.

## Features

- **Robust SMILES Validation**: Multi-level validation using RDKit built-in functions
  - Parse validation
  - Atom count check
  - Valence validation
  - Radical electron detection
  
- **Molecular Structure Visualization**: 2D structure rendering with heavy atoms only (hydrogens removed)

- **Property Prediction**: Random Forest models for:
  - Molecular Weight
  - LogP (lipophilicity)
  - TPSA (Topological Polar Surface Area)

## Project Structure

```
molecular_predictor/
├── app.py                  # Main Flask application
├── config.py              # Configuration settings
├── models/
│   └── predictor.py       # ML model implementation
├── utils/
│   └── validators.py      # SMILES validation & image generation
├── static/
│   ├── css/
│   │   └── style.css     # Stylesheets
│   └── js/
│       └── main.js       # Frontend JavaScript
├── templates/
│   └── index.html        # HTML template
├── requirements.txt       # Python dependencies
└── README.md             # This file
```

## Installation

1. Create a virtual environment:
```bash
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

Run the application:
```bash
python app.py
```

Visit `http://localhost:5000` in your browser.

## SMILES Validation

The app performs comprehensive validation:
1. **Format Check**: Ensures SMILES is a non-empty string
2. **Parsing**: Validates molecular structure can be parsed
3. **Sanitization**: Checks for valid valence states
4. **Quality Checks**: Detects radicals and unusual structures

## API Endpoints

- `GET /` - Main web interface
- `POST /predict` - Predict molecular properties
  - Request: `{"smiles": "CCO"}`
  - Response: `{"smiles": "CCO", "canonical_smiles": "CCO", "predicted_properties": {...}, "molecule_image": "base64..."}`
- `GET /health` - Health check endpoint

## Example SMILES to Try

Valid:
- `CCO` - ethanol
- `c1ccccc1` - benzene
- `CC(=O)O` - acetic acid
- `CC(C)(C)C` - neopentane

Invalid (will show specific error):
- `XYZ` - Invalid characters
- `C[` - Incomplete structure
- `C=C=C=C` - Valence issues

## Configuration

Edit `config.py` to modify:
- Model parameters (n_estimators, fingerprint_bits)
- Server settings (host, port, debug mode)
- Environment-specific configurations

## Future Enhancements

- Add model training endpoint
- Implement batch predictions
- Add user authentication
- Deploy with Docker
- Add database for storing predictions
- Support for 3D structure visualization