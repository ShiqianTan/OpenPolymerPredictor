"""
SMILES validation and molecule image generation utilities
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from io import BytesIO
import base64

def validate_smiles(smiles):
    """
    Validate SMILES string using RDKit built-in functions
    
    Args:
        smiles (str): SMILES string to validate
        
    Returns:
        tuple: (is_valid, mol_object, error_message)
    """
    if not smiles or not isinstance(smiles, str):
        return False, None, "SMILES string is empty or not a string"
    
    # Remove leading/trailing whitespace
    smiles = smiles.strip()
    
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        
        if mol is None:
            return False, None, "Invalid SMILES: Cannot parse molecular structure"
        
        # Additional validation checks
        num_atoms = mol.GetNumAtoms()
        if num_atoms == 0:
            return False, None, "Invalid SMILES: Molecule has no atoms"
        
        # Check for valid valence
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            return False, None, f"Invalid SMILES: Valence error - {str(e)}"
        
        # Check for radicals or unusual structures
        for atom in mol.GetAtoms():
            if atom.GetNumRadicalElectrons() > 0:
                return False, None, "Warning: Molecule contains radical electrons"
        
        # All checks passed
        return True, mol, None
        
    except Exception as e:
        return False, None, f"Invalid SMILES: {str(e)}"


def generate_molecule_image(mol, remove_hydrogens=True, img_size=(400, 300)):
    """
    Generate 2D molecular structure image
    
    Args:
        mol: RDKit molecule object
        remove_hydrogens (bool): If True, show only heavy atoms
        img_size (tuple): Image dimensions (width, height)
        
    Returns:
        str: Base64 encoded PNG image or None on error
    """
    try:
        # Create a copy to avoid modifying original
        mol_copy = Chem.Mol(mol)
        
        # Remove explicit hydrogens if requested
        if remove_hydrogens:
            mol_copy = Chem.RemoveHs(mol_copy)
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol_copy)
        
        # Draw molecule
        drawer = Draw.MolDraw2DCairo(img_size[0], img_size[1])
        drawer.DrawMolecule(mol_copy)
        drawer.FinishDrawing()
        
        # Get PNG data
        png_data = drawer.GetDrawingText()
        
        # Encode to base64
        img_base64 = base64.b64encode(png_data).decode('utf-8')
        
        return img_base64
        
    except Exception as e:
        print(f"Error generating molecule image: {e}")
        return None