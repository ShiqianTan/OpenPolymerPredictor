async function predict() {
    const smiles = document.getElementById('smiles').value;
    const resultDiv = document.getElementById('result');
    const moleculeDiv = document.getElementById('molecule-structure');
    const button = document.querySelector('button');
    
    if (!smiles) {
        showError('Please enter a SMILES string');
        moleculeDiv.style.display = 'none';
        return;
    }
    
    // Show loading state
    button.disabled = true;
    button.innerHTML = '<span class="loading"></span> Validating & Predicting...';
    
    try {
        const response = await fetch('/predict', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json',
            },
            body: JSON.stringify({ smiles: smiles })
        });
        
        const data = await response.json();
        
        if (data.error) {
            showError(data.error);
            moleculeDiv.style.display = 'none';
        } else {
            showResults(data);
            showMoleculeStructure(data);
        }
    } catch (error) {
        showError('Network error: ' + error.message);
        moleculeDiv.style.display = 'none';
    } finally {
        button.disabled = false;
        button.innerHTML = 'Predict Properties';
    }
}

function showResults(data) {
    const resultDiv = document.getElementById('result');
    resultDiv.className = 'result';
    resultDiv.innerHTML = `
        <h3>✅ Predicted Properties:</h3>
        <p><strong>Molecular Weight:</strong> ${data.predicted_properties.molecular_weight.toFixed(2)} g/mol</p>
        <p><strong>LogP:</strong> ${data.predicted_properties.logP.toFixed(2)}</p>
        <p><strong>TPSA:</strong> ${data.predicted_properties.tpsa.toFixed(2)} Ų</p>
        <hr>
        <p style="font-size: 12px; color: #666;"><strong>Input SMILES:</strong> ${data.smiles}</p>
        <p style="font-size: 12px; color: #666;"><strong>Canonical SMILES:</strong> ${data.canonical_smiles}</p>
    `;
    resultDiv.style.display = 'block';
}

function showMoleculeStructure(data) {
    const moleculeDiv = document.getElementById('molecule-structure');
    
    if (data.molecule_image) {
        moleculeDiv.innerHTML = `
            <h3>Molecular Structure (Heavy Atoms Only)</h3>
            <img src="data:image/png;base64,${data.molecule_image}" alt="Molecular Structure">
            <div class="molecule-info">2D structure generated with RDKit (hydrogens removed)</div>
        `;
        moleculeDiv.style.display = 'block';
    } else {
        moleculeDiv.style.display = 'none';
    }
}

function showError(message) {
    const resultDiv = document.getElementById('result');
    resultDiv.className = 'result error';
    resultDiv.innerHTML = `<h3>❌ Error</h3><p>${message}</p>`;
    resultDiv.style.display = 'block';
}

// Allow Enter key to submit
document.getElementById('smiles').addEventListener('keypress', function(event) {
    if (event.key === 'Enter') {
        predict();
    }
});