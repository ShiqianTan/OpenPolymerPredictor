async function predict() {
    const smiles = document.getElementById('smiles').value;
    const resultDiv = document.getElementById('result');
    const button = document.querySelector('button');
    
    if (!smiles) {
        showError('Please enter a SMILES string');
        return;
    }
    
    // Show loading state
    button.disabled = true;
    button.innerHTML = '<span class="loading"></span> Predicting...';
    
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
        } else {
            showResults(data);
        }
    } catch (error) {
        showError('Network error: ' + error.message);
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
        <p style="font-size: 12px; color: #666;">SMILES: ${data.smiles}</p>
    `;
    resultDiv.style.display = 'block';
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
