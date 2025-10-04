# Molecular Property Predictor

A scalable Flask application for predicting molecular properties from SMILES strings using Random Forest regression.

## Project Structure

```
molecular_predictor/
├── app.py                  # Main Flask application
├── config.py              # Configuration settings
├── models/
│   └── predictor.py       # ML model implementation
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

## API Endpoints

- `GET /` - Main web interface
- `POST /predict` - Predict molecular properties
- `GET /health` - Health check endpoint

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