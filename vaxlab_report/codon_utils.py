# Minimal codon table utilities for scoring functions

from Bio.Data import CodonTable

STOP = '*'

class MinimalCodonTable:
    """Minimal codon table for scoring functions that only need aa2codons mapping."""
    
    def __init__(self, codon_table='standard'):
        self.initialize_codon_table(codon_table)
    
    def initialize_codon_table(self, codon_table: str) -> None:
        table_var_name = f'{codon_table}_rna_table'
        if not hasattr(CodonTable, table_var_name):
            raise ValueError(f'Invalid codon table name: {codon_table}')

        self.codon_table = getattr(CodonTable, table_var_name)

        # Build aa2codons mapping without pandas
        self.aa2codons = {}
        
        # Add regular codons
        for codon, aa in self.codon_table.forward_table.items():
            if aa not in self.aa2codons:
                self.aa2codons[aa] = set()
            self.aa2codons[aa].add(codon)
        
        # Add stop codons
        if STOP not in self.aa2codons:
            self.aa2codons[STOP] = set()
        for stopcodon in self.codon_table.stop_codons:
            self.aa2codons[STOP].add(stopcodon)