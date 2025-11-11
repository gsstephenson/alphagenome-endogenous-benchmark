#!/usr/bin/env python3
"""Quick test to verify AlphaGenome API key is working."""

import os
from pathlib import Path
from dotenv import load_dotenv

# Load .env file
ENV_FILE = Path("/mnt/work_1/gest9386/CU_Boulder/rotations/LAYER/Alpha_genome_quickstart_notebook/.env")
if ENV_FILE.exists():
    load_dotenv(ENV_FILE)
    print(f"✓ Loaded .env file from {ENV_FILE}")
else:
    print(f"✗ .env file not found at {ENV_FILE}")

# Check if API key is set
API_KEY = os.getenv("ALPHA_GENOME_KEY")
if API_KEY:
    print(f"✓ API key found: {API_KEY[:20]}...{API_KEY[-4:]}")
else:
    print("✗ API key not found in environment")
    exit(1)

# Try to import AlphaGenome
try:
    from alphagenome.models.dna_client import create
    from alphagenome.models.dna_output import OutputType
    print("✓ AlphaGenome SDK imported successfully")
except ImportError as e:
    print(f"✗ Failed to import AlphaGenome: {e}")
    exit(1)

# Try to create client
try:
    client = create(api_key=API_KEY)
    print("✓ AlphaGenome client created successfully")
except Exception as e:
    print(f"✗ Failed to create client: {e}")
    exit(1)

# Try a simple prediction with a short test sequence
try:
    test_seq = "ATCG" * 512  # 2048 bp test sequence (exactly)
    print(f"→ Testing prediction on {len(test_seq)} bp sequence...")
    
    output = client.predict_sequence(
        sequence=test_seq,
        requested_outputs=[OutputType.DNASE],
        ontology_terms=None
    )
    
    print(f"✓ Prediction successful!")
    print(f"  - Output shape: {output.dnase.values.shape}")
    print(f"  - Mean DNase signal: {output.dnase.values.mean():.6f}")
    print("\n🎉 API key is working correctly!")
    
except Exception as e:
    print(f"✗ Prediction failed: {e}")
    exit(1)
