"""IDT Complexity Score API Integration"""

import requests
from typing import Dict, List, Any, Optional, Tuple
import logging

ENDPOINT_URL = 'https://www.idtdna.com/api/complexities/screengBlockSequences'

def evaluate_idt_gblock_complexity(credential: str, sequences: Dict[str, str], timeout: int = 60) -> List[Dict[str, Any]]:
    """
    Evaluate IDT gBlock complexity for given sequences.
    
    Args:
        credential: IDT API access token
        sequences: Dictionary mapping sequence names to sequences
        timeout: Request timeout in seconds
    
    Returns:
        List of complexity evaluation results from IDT API
    """
    sequences_enc = [
        {'Name': name, 'Sequence': seq}
        for name, seq in sequences.items()
    ]

    request_headers = {
        'Content-Type': 'application/json',
        'Authorization': 'Bearer ' + credential
    }

    resp = requests.post(
        ENDPOINT_URL,
        json=sequences_enc,
        timeout=timeout,
        headers=request_headers
    )

    resp.raise_for_status()
    response = resp.json()

    return response

def score_idt_complexity(seq: str, idt_token: Optional[str] = None, **kwargs) -> Tuple[float, Dict[str, Any]]:
    """
    Score function for IDT complexity that integrates with VaxLab scoring system.
    
    Args:
        seq: DNA/RNA sequence to evaluate
        idt_token: IDT API access token
        **kwargs: Additional parameters
    
    Returns:
        Tuple of (score, metrics_dict)
    """
    log = logging.getLogger(__name__)
    
    if not idt_token:
        log.warning("IDT API token not provided. Skipping IDT complexity evaluation.")
        return 0.0, {"idt_complexity_error": "No API token provided"}
    
    # Convert RNA to DNA for IDT API
    dna_seq = seq.replace('U', 'T')
    
    try:
        # Call IDT API
        results = evaluate_idt_gblock_complexity(idt_token, {'sequence': dna_seq})
        
        if results and len(results) > 0:
            result = results[0]  # Get first result
            
            # Extract all complexity data
            complexity_items = []
            if isinstance(result, list):
                complexity_items = result
            elif isinstance(result, dict) and 'ComplexityItems' in result:
                complexity_items = result['ComplexityItems']
            elif isinstance(result, dict):
                # If result is directly the complexity item
                complexity_items = [result]
            
            # Process complexity items
            idt_metrics = {
                "idt_complexity_items": complexity_items,
                "idt_complexity_count": len(complexity_items)
            }
            
            # Calculate aggregate score (higher score = more complex = worse)
            total_score = 0.0
            for item in complexity_items:
                if isinstance(item, dict) and 'Score' in item:
                    total_score += float(item['Score'])
            
            idt_metrics["idt_complexity_total_score"] = total_score
            
            # Return negative score (so minimization = better)
            return -total_score, idt_metrics
            
        else:
            return 0.0, {"idt_complexity_error": "No results from API"}
            
    except requests.exceptions.RequestException as e:
        log.error(f"IDT API request error: {e}")
        return 0.0, {"idt_complexity_error": f"API request failed: {str(e)}"}
    except Exception as e:
        log.error(f"Error evaluating IDT complexity: {e}")
        return 0.0, {"idt_complexity_error": f"Evaluation error: {str(e)}"}

# VaxLab scoring function interface
name = "idt_complexity"
description = "IDT gBlock complexity score"