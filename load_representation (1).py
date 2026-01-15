"""
load_representation.py

Loads Artin representation data from LMFDB-style JSON or text files.
Extracts character tables, Frobenius data, local factors, and Swan conductors.

Usage:
    python3 load_representation.py --file data/A5_Conductor_Data.txt
    python3 load_representation.py --lmfdb 3.2083.12t76.a
"""

import json
import argparse
from pathlib import Path

def parse_lmfdb_label(label):
    """
    Parse LMFDB label format: dimension.conductor.galois_group.variant
    Example: 3.2083.12t76.a → (dim=3, N=2083, group=12t76, var=a)
    """
    parts = label.split('.')
    if len(parts) != 4:
        raise ValueError(f"Invalid LMFDB label format: {label}")
    
    return {
        'dimension': int(parts[0]),
        'conductor': int(parts[1]),
        'galois_label': parts[2],
        'variant': parts[3],
        'label': label
    }

def load_from_json(filepath):
    """Load representation data from JSON file"""
    with open(filepath, 'r') as f:
        data = json.load(f)
    
    return normalize_representation_data(data)

def load_from_text(filepath):
    """Load representation data from text file (legacy format)"""
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Parse text format (assumes JSON embedded or key-value pairs)
    try:
        data = json.loads(content)
        return normalize_representation_data(data)
    except json.JSONDecodeError:
        # Parse key-value format
        return parse_text_format(content)

def parse_text_format(content):
    """Parse legacy text format with key: value pairs"""
    data = {}
    for line in content.strip().split('\n'):
        if ':' in line:
            key, value = line.split(':', 1)
            key = key.strip().lower().replace(' ', '_')
            value = value.strip()
            
            # Try to parse as number
            try:
                if '.' in value:
                    data[key] = float(value)
                else:
                    data[key] = int(value)
            except ValueError:
                data[key] = value
    
    return normalize_representation_data(data)

def normalize_representation_data(raw_data):
    """
    Normalize representation data to standard format.
    
    Standard format:
    {
        'label': str,
        'conductor': int,
        'dimension': int,
        'galois_group': str,
        'character_table': {
            'classes': [str],
            'traces': [float],
            'norms': [float]
        },
        'frobenius_data': [
            {'prime': int, 'conjugacy_class': str, 'trace': float}
        ],
        'local_factors': [
            {'prime': int, 'polynomial': [float], 'type': str}
        ],
        'ramification': {
            'ramified_primes': [
                {
                    'prime': int,
                    'tame': bool,
                    'wild': bool,
                    'swan_conductor': int,
                    'inertia_quotient_dim': int
                }
            ],
            'conductor_exponents': {int: int}
        },
        'kappa_inf': float,
        'kappa_local': float
    }
    """
    normalized = {}
    
    # Basic metadata
    normalized['label'] = raw_data.get('label', 'unknown')
    normalized['conductor'] = int(raw_data.get('conductor', 0))
    normalized['dimension'] = int(raw_data.get('dimension', 0))
    normalized['galois_group'] = raw_data.get('galois_group', raw_data.get('type', 'unknown'))
    
    # Character table
    if 'character_table' in raw_data:
        ct = raw_data['character_table']
        normalized['character_table'] = {
            'classes': ct.get('classes', []),
            'traces': [float(x) for x in ct.get('traces', [])],
            'norms': [float(x) for x in ct.get('norms', [])]
        }
    else:
        normalized['character_table'] = {
            'classes': [],
            'traces': [],
            'norms': []
        }
    
    # Frobenius data (if available)
    normalized['frobenius_data'] = raw_data.get('frobenius_data', [])
    
    # Local factors
    normalized['local_factors'] = raw_data.get('local_factors', [])
    
    # Ramification data
    if 'ramification' in raw_data:
        normalized['ramification'] = raw_data['ramification']
    else:
        # Construct from other fields
        ram_primes = []
        if 'ramified_primes' in raw_data:
            for p_data in raw_data['ramified_primes']:
                ram_primes.append({
                    'prime': p_data.get('prime', 0),
                    'tame': p_data.get('tame', True),
                    'wild': p_data.get('wild', False),
                    'swan_conductor': p_data.get('swan_conductor', 0),
                    'inertia_quotient_dim': p_data.get('inertia_quotient_dim', 1)
                })
        
        normalized['ramification'] = {
            'ramified_primes': ram_primes,
            'conductor_exponents': {}
        }
    
    # Curvature constants
    normalized['kappa_inf'] = float(raw_data.get('kappa_inf', raw_data.get('kappa_local', 0)))
    normalized['kappa_local'] = float(raw_data.get('kappa_local', 0))
    
    return normalized

def verify_schur_orthogonality(character_table, group_order=None):
    """
    Verify Schur orthogonality: (1/|G|) Σ |χ(g)|² = 1
    
    This is the fundamental algebraic requirement for deficiency indices (0,0).
    """
    if not character_table.get('norms') or not character_table.get('classes'):
        print("Warning: Insufficient character table data for orthogonality check")
        return None
    
    norms = character_table['norms']
    
    # If class sizes are provided, use them; otherwise assume from LMFDB format
    if 'class_sizes' in character_table:
        class_sizes = character_table['class_sizes']
    else:
        # Try to infer from group structure
        if 'classes' in character_table:
            # For known groups, use standard class sizes
            # This is placeholder - real implementation would parse group structure
            class_sizes = [1] * len(norms)  # Simplification
    
    if group_order is None:
        # Try to determine from Galois group label or sum of class sizes
        group_order = sum(class_sizes)
    
    # Compute Schur orthogonality value
    total = sum(size * norm for size, norm in zip(class_sizes, norms))
    ortho_value = total / group_order
    
    error = abs(1.0 - ortho_value)
    
    return {
        'value': ortho_value,
        'error': error,
        'passes': error < 1e-10
    }

def print_representation_summary(data):
    """Print human-readable summary of representation data"""
    print("=" * 70)
    print("ARTIN REPRESENTATION SUMMARY")
    print("=" * 70)
    print(f"\nLabel:          {data['label']}")
    print(f"Conductor:      {data['conductor']}")
    print(f"Dimension:      {data['dimension']}")
    print(f"Galois Group:   {data['galois_group']}")
    
    if data['character_table']['traces']:
        print(f"\nCharacter Table:")
        print(f"  Classes:  {data['character_table']['classes']}")
        print(f"  Traces:   {[f'{x:.3f}' for x in data['character_table']['traces']]}")
        print(f"  Norms:    {[f'{x:.3f}' for x in data['character_table']['norms']]}")
    
    if data['ramification']['ramified_primes']:
        print(f"\nRamification:")
        for p_data in data['ramification']['ramified_primes']:
            prime = p_data['prime']
            wild = " (WILD)" if p_data['wild'] else " (tame)"
            swan = f", Swan={p_data['swan_conductor']}" if p_data['swan_conductor'] > 0 else ""
            print(f"  p = {prime}{wild}{swan}")
    
    print(f"\nCurvature Constants:")
    print(f"  κ_∞ (Archimedean):  {data['kappa_inf']:.6f}")
    print(f"  κ_local (Ramified): {data['kappa_local']:.6f}")
    print(f"  κ(L) (Total):       {data['kappa_inf'] + data['kappa_local']:.6f}")
    
    print("=" * 70)

def export_to_json(data, output_path):
    """Export normalized data to JSON file"""
    with open(output_path, 'w') as f:
        json.dump(data, f, indent=2)
    print(f"\nExported to: {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Load and normalize Artin representation data"
    )
    parser.add_argument('--file', type=str,
                       help='Path to representation data file')
    parser.add_argument('--lmfdb', type=str,
                       help='LMFDB label (e.g., 3.2083.12t76.a)')
    parser.add_argument('--verify', action='store_true',
                       help='Verify Schur orthogonality')
    parser.add_argument('--export', type=str,
                       help='Export normalized data to JSON file')
    
    args = parser.parse_args()
    
    if args.file:
        # Load from file
        filepath = Path(args.file)
        if filepath.suffix == '.json':
            data = load_from_json(filepath)
        else:
            data = load_from_text(filepath)
        
        print_representation_summary(data)
        
        if args.verify:
            ortho = verify_schur_orthogonality(data['character_table'])
            if ortho:
                print(f"\nSchur Orthogonality Check:")
                print(f"  Value: {ortho['value']:.15f}")
                print(f"  Error: {ortho['error']:.2e}")
                print(f"  Status: {'PASS' if ortho['passes'] else 'FAIL'}")
        
        if args.export:
            export_to_json(data, args.export)
    
    elif args.lmfdb:
        # Parse LMFDB label
        label_data = parse_lmfdb_label(args.lmfdb)
        print(f"LMFDB Label: {args.lmfdb}")
        print(f"  Dimension: {label_data['dimension']}")
        print(f"  Conductor: {label_data['conductor']}")
        print(f"  Galois Group: {label_data['galois_label']}")
        print("\nNote: Fetching from LMFDB API not yet implemented")
        print("Please download data manually and use --file option")
    
    else:
        parser.print_help()
