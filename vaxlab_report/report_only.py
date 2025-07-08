import argparse
import os
# import pandas as pd # Not strictly needed if not using df_checkpoints
# import csv # Not strictly needed if not using checkpoints_list from tsv
from Bio import SeqIO
from Bio.Seq import Seq # Not directly used, but often comes with SeqRecord
from Bio.SeqRecord import SeqRecord # Not directly used for ReportGenerator input
import time
import json
import plotly.graph_objects as go
import plotly.offline as pyo
import numpy as np
from typing import Dict, Any, List, Optional
import traceback
import shlex # For command line string reconstruction
from vaxlab_report.log import initialize_logging, log

try:
    from vaxlab_report.reporting import ReportGenerator
    from vaxlab_report.presets import load_preset
    try:
        from vaxlab_report.evolution_chamber import ExecutionOptions
    except ImportError:
        # Fallback definition if the main one is not found (e.g. simpler environment)
        class ExecutionOptions:
            def __init__(self, n_iterations=0, n_population=1, n_survivors=1,
                         initial_mutation_rate=0.1, winddown_trigger=15, winddown_rate=0.9,
                         output='.', command_line='', overwrite=False, seed=0, processes=1,
                         random_initialization=False, conservative_start=None,
                         boost_loop_mutations="1.5:15", full_scan_interval=300,
                         species="Homo sapiens", codon_table="standard", protein=False, quiet=True,
                         seq_description=None, print_top_mutants=0, addons=None,
                         lineardesign_dir=None, lineardesign_lambda=None, lineardesign_omit_start=5,
                         folding_engine="vienna", **kwargs):
                self.n_iterations=n_iterations; self.n_population=n_population; self.n_survivors=n_survivors
                self.initial_mutation_rate=initial_mutation_rate; self.winddown_trigger=winddown_trigger; self.winddown_rate=winddown_rate
                self.output=output; self.command_line=command_line; self.overwrite=overwrite; self.seed=seed; self.processes=processes
                self.random_initialization=random_initialization; self.conservative_start=conservative_start
                self.boost_loop_mutations=boost_loop_mutations; self.full_scan_interval=full_scan_interval
                self.species=species; self.codon_table=codon_table; self.protein=protein; self.quiet=quiet
                self.seq_description=seq_description; self.print_top_mutants=print_top_mutants
                self.addons=addons if addons is not None else []; self.lineardesign_dir=lineardesign_dir
                self.lineardesign_lambda=lineardesign_lambda; self.lineardesign_omit_start=lineardesign_omit_start
                self.folding_engine=folding_engine; self._kwargs=kwargs
            def to_dict(self) -> Dict[str, Any]:
                data = self.__dict__.copy(); data.pop('_kwargs', None); return data
except ImportError as e:
    print(f"Error importing vaxlab_report components: {e}"); exit(1)

METRIC_LABELS = {
    "gc": "GC ratio", "ucount": "Uridine ratio", "bicodon": "Bicodon score",
    "repeat": "Repeat count", "longstem": "Long stem count",
    "start_str": "Start codon structure", "loop": "Loop Length", "mfe": "Minimum Free Energy",
    "cai": "CAI", "structure": "Secondary Structure", "degscore": "DegScore"
}
METRIC_DESCRIPTIONS = {
    "gc": "Ratio of G/C in the sequence", "degscore": "Predicted degradation score (Eterna DegScore)",
    "ucount": "Ratio of U in the sequence", "bicodon": "Codon pair usage bias",
    "repeat": "Total length of all tandem repeats found" , "longstem": "Count of long stem-loop structures (≥ 27bp)",
    "start_str": "Base-paired nucleotides near start codon (-14 to +14)",
    "loop": "Total length of unpaired regions (loops ≥ 2nt)", "mfe": "Minimum free energy of predicted structure",
    "cai": "Codon usage optimality based on species-specific frequencies (0-1 scale).",
    "structure": "Predicted secondary structure (dot-bracket)"
}
METRIC_RANGES = {
    "gc": "0.4–0.6", "ucount": "Lowest possible", "bicodon": "Maximize (≥ 0)",
    "repeat": "< 20 nt", "longstem": "0–1", "start_str": "≤ 3 bp paired",
    "loop": "≤ 100 per kb", "mfe": "≤ –30 kcal/mol per kb", "cai": "Maximize (closer to 1)",
    "degscore": "Lowest possible", "structure": "N/A"
}

metainfo = {"vaxpress_version": "unknown", "command_line": " ".join(shlex.quote(x) for x in os.sys.argv), "start_time": time.time(), "forna_url": None }
scoring_functions = {} 
scoring_options = {}

def generate_idt_complexity_table(idt_data: List[Dict[str, Any]]) -> str:
    """
    Generate HTML table for IDT complexity results with emoji visualization.
    
    Args:
        idt_data: List of IDT complexity items
    
    Returns:
        HTML string for the table
    """
    if not idt_data:
        return ""
    
    # Calculate total IDT complexity score
    total_score = sum(float(item.get('Score', 0)) for item in idt_data)
    
    # Determine overall status based on total score
    if total_score < 7:
        total_emoji = "😊"
        total_status_class = "text-success"
        total_status = "Low"
    elif 7 <= total_score < 20:
        total_emoji = "⚠️"
        total_status_class = "text-warning"
        total_status = "Moderate"
    else:
        total_emoji = "😞"
        total_status_class = "text-danger"
        total_status = "High"
    
    html = f"""
    <div class="idt-complexity-section" style="margin-top: 30px;">
        <h3>IDT Complexity Analysis</h3>
        <p style="margin-bottom: 15px;"><strong>IDT Complexity Score: <span class="{total_status_class}">{total_score:.1f} {total_emoji} ({total_status})</span></strong></p>
        <table class="table table-striped table-hover" style="width: 100%;">
            <thead>
                <tr>
                    <th>Name</th>
                    <th>Score</th>
                    <th>Status</th>
                    <th>Actual Value</th>
                    <th>Description</th>
                    <th>Solution</th>
                </tr>
            </thead>
            <tbody>
    """
    
    for item in idt_data:
        name = item.get('Name', 'Unknown')
        score = float(item.get('Score', 0))
        actual_value = item.get('ActualValue', 'N/A')
        display_text = item.get('DisplayText', 'N/A')
        
        # Split description and solution
        if 'Solution:' in display_text:
            description_part, solution_part = display_text.split('Solution:', 1)
            description = description_part.strip()
            solution = solution_part.strip()
        else:
            description = display_text
            solution = 'N/A'
        
        # Determine emoji based on score
        if score < 6:
            emoji = "😊"
            status_class = "text-success"
        elif 6 <= score < 15:
            emoji = "⚠️"
            status_class = "text-warning"
        else:
            emoji = "😞"
            status_class = "text-danger"
        
        html += f"""
                <tr>
                    <td>{name}</td>
                    <td>{score:.1f}</td>
                    <td class="{status_class}">{emoji}</td>
                    <td>{actual_value}</td>
                    <td>{description}</td>
                    <td>{solution}</td>
                </tr>
        """
    
    html += """
            </tbody>
        </table>
        <p class="text-muted">
            <small>
                Individual score interpretation: &lt;6 = Good 😊 | 6-15 = Caution ⚠️ | ≥15 = Poor 😞<br>
                Total complexity interpretation: &lt;7 = Low 😊 | 7-20 = Moderate ⚠️ | ≥20 = High 😞
            </small>
        </p>
    </div>
    """
    
    return html


def generate_positional_plot(positions: list, values: dict, title: str) -> Optional[str]:
    if not positions or not values: log.warning("No data for positional plot."); return None
    fig = go.Figure(); valid_trace_added = False
    for name, data in values.items():
        if data is not None and len(data) == len(positions):
             y_values = [v if not (v is None or np.isnan(v)) else None for v in data]
             fig.add_trace(go.Scatter(x=positions, y=y_values, mode='lines', name=name, connectgaps=False)); valid_trace_added = True
        elif data is not None: log.warning(f"Skipping trace '{name}' (length mismatch {len(data)} vs {len(positions)}).")
        else: log.warning(f"Skipping trace '{name}' (data is None).")
    if not valid_trace_added: log.warning("No valid traces added to plot."); return None
    fig.update_layout(title=title, xaxis_title="Position (approx. window center, nt)", yaxis_title="Metric Value", hovermode="x unified")
    try: return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    except Exception as e: log.error(f"Error converting Plotly figure to div: {e}"); traceback.print_exc(); return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file (containing 5UTR, CDS, 3UTR)")
    parser.add_argument("-o", "--output", required=True, help="Output directory (must contain evaluation_result_cds.json and evaluation_result_mrna.json)")
    parser.add_argument("--preset", help=argparse.SUPPRESS)
    parser.add_argument('--iters', type=int, default=0, help=argparse.SUPPRESS)
    parser.add_argument('--population', type=int, default=0, help=argparse.SUPPRESS)
    parser.add_argument('--cpu-count', type=int, default=1, help=argparse.SUPPRESS)
    parser.add_argument('--random-seed', type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument('--species', type=str, default='Homo sapiens', help=argparse.SUPPRESS)
    parser.add_argument('--forna', type=str, choices=['qbio', 'tbi'], default=None, help='Forna server to use for structure visualization (qbio or tbi)')
    args = parser.parse_args()

    log_file = os.path.join(args.output, "report_only_log.txt")
    initialize_logging(log_file, quiet=False)
    log.debug(f"Parsed args: {args}")

    execution_options_dict_from_preset = {}
    if args.preset:
        try:
            preset_data = load_preset(args.preset)
            scoring_options.update(preset_data.get("scoring", {})) # global scoring_options
            execution_options_dict_from_preset.update(preset_data.get("execution", {}))
            log.debug(f"Loaded preset '{args.preset}'")
        except Exception as e: print(f"Error loading preset '{args.preset}': {e}"); exit(1)

    utr5_seq_str = ""
    cds_seq_str = ""
    utr3_seq_str = ""
    cds_record_id = "CDS_Unknown"
    cds_record_desc = "CDS sequence"
    
    try:
        num_potential_cds = 0
        temp_cds_seq = None
        temp_cds_id = None
        temp_cds_desc = None
        for record in SeqIO.parse(args.input, "fasta"):
            seq_upper_u = str(record.seq).upper().replace("T", "U")
            if record.id.upper() == "5UTR":
                utr5_seq_str = seq_upper_u
                log.info(f"Loaded 5UTR: {record.id}, Length: {len(utr5_seq_str)}")
            elif record.id.upper() == "3UTR":
                utr3_seq_str = seq_upper_u
                log.info(f"Loaded 3UTR: {record.id}, Length: {len(utr3_seq_str)}")
            else: # Assumed CDS
                if num_potential_cds == 0: # Take the first non-UTR as CDS
                    temp_cds_seq = seq_upper_u
                    temp_cds_id = record.id
                    temp_cds_desc = record.description
                num_potential_cds +=1
        
        if temp_cds_seq:
            cds_seq_str = temp_cds_seq
            cds_record_id = temp_cds_id
            cds_record_desc = temp_cds_desc
            log.info(f"Identified CDS: {cds_record_id}, Length: {len(cds_seq_str)}")
            if num_potential_cds > 1:
                 log.warning(f"Multiple non-UTR sequences found. Using '{cds_record_id}' as CDS.")
        else:
            log.critical(f"CDS sequence not found in {args.input}"); exit(1)
    except FileNotFoundError: print(f"Error: Input FASTA file not found: {args.input}"); exit(1)
    except Exception as e: print(f"Error reading input FASTA file: {e}"); exit(1)

    entire_mrna_sequence = utr5_seq_str + cds_seq_str + utr3_seq_str
    outputseq_dict = { # This is for the entire mRNA, which is the primary output of the report
        'id': f"{cds_record_id}_mRNA_Report",
        'seq': entire_mrna_sequence,
        'description': f"{cds_record_desc} (Full mRNA: 5UTR-CDS-3UTR)"
    }
    inputseq_dict = {'id': cds_record_id, 'seq': cds_seq_str, 'description': cds_record_desc}
    log.info(f"Report target mRNA: ID={outputseq_dict['id']}, Length={len(outputseq_dict['seq'])}")

    # Initialize ExecutionOptions with all required arguments
    init_args = {
        'n_iterations': args.iters, 'n_population': args.population, 'n_survivors': 1,
        'output': args.output, 'command_line': "report_only execution", 'overwrite': True,
        'seed': args.random_seed if args.random_seed is not None else int(time.time()), # Ensure seed is int
        'processes': args.cpu_count, 'species': args.species, 'codon_table': "standard",
        'quiet': True, 'seq_description': outputseq_dict['description'], # Use entire mRNA description
        'initial_mutation_rate': 0.1, 'winddown_trigger': 15, 'winddown_rate': 0.9,
        'random_initialization': False, 'conservative_start': None,
        'boost_loop_mutations': "1.5:15", 'full_scan_interval': 0, 
        'protein': False, 'print_top_mutants': 0, 'addons': [], 
        'lineardesign_dir': None, 'lineardesign_lambda': None, 'lineardesign_omit_start': 5,
        'folding_engine': "vienna", 
    }
    init_args.update(execution_options_dict_from_preset)
    # Override with command-line args where applicable
    init_args.update({k: v for k, v in vars(args).items() if v is not None and k in init_args})


    try: execution_options = ExecutionOptions(**init_args); log.debug("Execution Options Initialized.")
    except Exception as e: print(f"Error: Failed to initialize ExecutionOptions: {e}"); exit(1)

    # --- Load Evaluation Data from JSON files ---
    cds_eval_path = os.path.join(args.output, "evaluation_result_cds.json")
    mrna_eval_path = os.path.join(args.output, "evaluation_result_mrna.json")

    if not os.path.exists(cds_eval_path):
        print(f"Error: CDS evaluation result file not found: {cds_eval_path}"); exit(1)
    if not os.path.exists(mrna_eval_path):
        print(f"Error: mRNA evaluation result file not found: {mrna_eval_path}"); exit(1)

    log.info(f"Loading CDS evaluation results from: {cds_eval_path}")
    with open(cds_eval_path, 'r') as f:
        cds_evaluation_result = json.load(f)
    log.info(f"Loading mRNA evaluation results from: {mrna_eval_path}")
    with open(mrna_eval_path, 'r') as f:
        mrna_evaluation_result = json.load(f)
    
    # Extract IDT complexity data if available
    idt_complexity_data = cds_evaluation_result.get("idt_complexity", None)

    # --- Prepare Global Metrics for the Report ---
    # Start with mRNA global metrics as the base
    final_global_metrics = mrna_evaluation_result.get("global_metrics", {}).copy()
    cds_global_metrics = cds_evaluation_result.get("global_metrics", {})
    # log.info(f"Base global metrics from mRNA: {final_global_metrics}")
    # log.info(f"Global metrics from CDS: {cds_global_metrics}")
    
    # Define which metrics to take from CDS results
    # Key: key in cds_global_metrics (from cds_evaluation_result.json)
    # Value: key to be used in final_global_metrics and METRIC_LABELS
    metrics_to_source_from_cds = {
        "cai": "cai",
        "bicodon": "bicodon",
        "start_str": "start_str"
    }

    for cds_key, target_key in metrics_to_source_from_cds.items():
        if cds_key in cds_global_metrics:
            value_from_cds = cds_global_metrics[cds_key]
            if target_key == "cai":
                try:
                    log_score_cai = float(value_from_cds)
                    # Assuming the value in JSON is a log-score, convert to 0-1 scale CAI
                    # This assumption needs to be verified based on how evaluate_only.py saves CAI for CDS.
                    # If evaluate_only.py already saves the exp() transformed value, this np.exp() is not needed.
                    # Given -0.0346..., it's likely a log score.
                    if -10 < log_score_cai < 10 : # Heuristic to check if it's a log score rather than already 0-1
                         standard_cai_value = np.exp(log_score_cai)
                         final_global_metrics[target_key] = standard_cai_value
                         log.info(f"Global metric '{target_key}' (CAI) from CDS: {standard_cai_value:.4f} (converted from log-score: {log_score_cai:.4f})")
                    else: # Assume it's already in the desired scale or not a typical log-score for exp transform
                         final_global_metrics[target_key] = float(value_from_cds)
                         log.info(f"Global metric '{target_key}' (CAI) from CDS: {value_from_cds} (used as is or float converted)")
                except (ValueError, TypeError) as e:
                    final_global_metrics[target_key] = value_from_cds # Fallback
                    log.warning(f"Could not convert/process CAI value '{value_from_cds}' from CDS. Using as is. Error: {e}")
            else:
                final_global_metrics[target_key] = value_from_cds
                log.info(f"Global metric '{target_key}' from CDS: {value_from_cds}")
        else:
            log.warning(f"Metric key '{cds_key}' for target '{target_key}' not found in CDS global metrics.")
            if target_key not in final_global_metrics: # If not even in mRNA metrics
                 final_global_metrics[target_key] = None
    
    # Ensure 'structure' comes from mRNA results for Forna
    if "structure" in mrna_evaluation_result.get("global_metrics", {}):
        final_global_metrics["structure"] = mrna_evaluation_result["global_metrics"]["structure"]
    elif "structure" not in final_global_metrics : # If it wasn't set by CDS override (it shouldn't be)
        log.warning("Structure (dot-bracket) string not found in mRNA results.")
        final_global_metrics["structure"] = "" # Default to empty string for Forna

    log.debug(f"Final combined global metrics for report: {final_global_metrics}")
    # log.info(f"Structure in final_global_metrics: '{final_global_metrics.get('structure', 'NOT FOUND')}'")

    # Local metrics for plot are from entire mRNA evaluation
    final_local_metrics = mrna_evaluation_result.get("local_metrics", {})

    checkpoints_list = [] # --evaluations option removed, so this remains empty or is not used.
    log.info("Checkpoints.tsv loading skipped as --evaluations option is not used for primary data.")


    # --- Prepare Positional Plot for Entire mRNA (excluding CAI) ---
    positional_plot_div = None
    plot_data = {}
    all_positions = []
    
    log.info("--- Processing Positional Metrics for Plot (Entire mRNA, CAI excluded) ---")
    metrics_to_process_for_plot: Dict[str, List[List[Any]]] = {}

    if isinstance(final_local_metrics, dict):
        for key, metric_data in final_local_metrics.items():
            # Exclude CAI from the plot as it's not evaluated for the entire mRNA
            if key == 'cai' or METRIC_LABELS.get(key, key).lower() == 'cai':
                log.info(f"Excluding '{key}' from positional plot for entire mRNA.")
                continue
            if key.endswith("_error"): continue # Skip error entries

            if isinstance(metric_data, (list, tuple)) and len(metric_data) == 2 and \
               all(isinstance(el, (list, np.ndarray, tuple)) for el in metric_data): # Check for [[positions], [values]]
                positions, values = metric_data
                # Ensure positions and values are lists/arrays of numbers
                if all(isinstance(p, (int, float)) for p in positions) and \
                   all(isinstance(v, (int, float, type(None))) or np.isnan(v) for v in values if isinstance(v, (int, float, type(None))) or hasattr(v, '__iter__') is False ): # Check for numbers or None/NaN
                    if len(positions) == len(values):
                        num_positions = [p for p in positions]
                        num_values = [float(v) if v is not None and not (isinstance(v, float) and np.isnan(v)) else np.nan for v in values]
                        metrics_to_process_for_plot[key] = [num_positions, num_values]
                    else: log.warning(f"Length mismatch for local metric '{key}' for plot.")
                else: log.warning(f"Non-numeric data in positions/values for local metric '{key}' for plot.")
            elif metric_data is not None:
                 log.warning(f"Unexpected data format for local metric '{key}' for plot: {type(metric_data)}.")
    else:
        log.warning("final_local_metrics (mRNA) is not a dictionary. Cannot process for positional plot.")
    
    # Determine base x-axis for the plot
    priority_keys_for_x_axis = ["gc", "degscore"] # Common metrics with good positional data
    base_metric_key_for_x_axis = None
    for key in priority_keys_for_x_axis:
        if key in metrics_to_process_for_plot:
            base_metric_key_for_x_axis = key; break
    if not base_metric_key_for_x_axis and metrics_to_process_for_plot: # Fallback to first available
        base_metric_key_for_x_axis = list(metrics_to_process_for_plot.keys())[0]

    if base_metric_key_for_x_axis:
        all_positions = metrics_to_process_for_plot[base_metric_key_for_x_axis][0]
        log.info(f"Using positions from '{base_metric_key_for_x_axis}' as x-axis for plot ({len(all_positions)} points).")
        
        plot_order_preference = ["gc", "degscore"] # Order of traces in plot
        plot_keys_ordered = [k for k in plot_order_preference if k in metrics_to_process_for_plot] + \
                            [k for k in metrics_to_process_for_plot if k not in plot_order_preference and k != base_metric_key_for_x_axis]
        
        # Add base metric first if not in preference, to ensure it's plotted
        if base_metric_key_for_x_axis not in plot_keys_ordered:
            plot_keys_ordered.insert(0, base_metric_key_for_x_axis)
        
        unique_plot_keys = []
        for k in plot_keys_ordered: # Ensure no duplicates if base_metric was also in preference
            if k not in unique_plot_keys:
                unique_plot_keys.append(k)

        for key in unique_plot_keys:
            positions, values = metrics_to_process_for_plot[key]
            plot_label = METRIC_LABELS.get(key, key)
            
            if list(positions) == list(all_positions): # No interpolation needed
                plot_values = [v if v is not None and not np.isnan(v) else None for v in values]
            else: # Interpolate to common x-axis
                log.debug(f"Interpolating values for '{plot_label}' (original key: {key})...")
                try:
                    pos_arr = np.array(positions)
                    val_arr = np.array([float(v) if v is not None and not np.isnan(v) else np.nan for v in values])
                    valid_mask = ~np.isnan(val_arr)
                    if np.any(valid_mask):
                         interp_positions = pos_arr[valid_mask]
                         interp_values = val_arr[valid_mask]
                         plot_values = np.interp(all_positions, interp_positions, interp_values, left=np.nan, right=np.nan).tolist()
                    else: plot_values = [None] * len(all_positions) # All NaNs
                except Exception as interp_e:
                    log.error(f"Error during interpolation for '{plot_label}': {interp_e}", exc_info=True)
                    plot_values = [None] * len(all_positions)
            
            if plot_values is not None: plot_data[plot_label] = plot_values
            else: log.warning(f"Could not determine plot values for '{plot_label}'.")

    if all_positions and plot_data:
        plot_title = f"Positional Metrics (Entire mRNA)"
        positional_plot_div = generate_positional_plot(all_positions, plot_data, plot_title)
        if positional_plot_div: log.info("Positional plot for entire mRNA generated successfully.")
        else: log.warning("Positional plot generation for entire mRNA failed.")
    else:
        log.warning("No valid data or common x-axis to generate positional plot for entire mRNA.")

    # --- Generate IDT Complexity Table HTML ---
    idt_complexity_html = ""
    idt_token_error = False
    if idt_complexity_data:
        log.info("Generating IDT complexity table...")
        # Check if it's a list of items or a single result containing items
        if isinstance(idt_complexity_data, list):
            idt_complexity_html = generate_idt_complexity_table(idt_complexity_data)
        elif isinstance(idt_complexity_data, dict):
            # Check for 401 Unauthorized error
            if 'error' in idt_complexity_data and '401' in str(idt_complexity_data.get('error', '')):
                idt_token_error = True
                idt_complexity_html = """
                <div class="alert alert-danger" style="margin-top: 30px; padding: 20px; background-color: #f8d7da; border: 1px solid #f5c6cb; border-radius: 5px; color: #721c24;">
                    <h4 style="margin-top: 0;">⚠️ IDT Complexity Score Error</h4>
                    <p><strong>Invalid or expired IDT API token.</strong></p>
                    <p>Obtain a new valid IDT API token and re-run VaxLab with the correct token value.</p>
                </div>
                """
            else:
                # Try to extract items from different possible structures
                items = idt_complexity_data.get('ComplexityItems', idt_complexity_data.get('items', [idt_complexity_data]))
                if items:
                    idt_complexity_html = generate_idt_complexity_table(items if isinstance(items, list) else [items])
    
    # --- Construct Status Dictionary for ReportGenerator ---
    status: Dict[str, Any] = {
        "checkpoints": checkpoints_list, # Will be empty
        "evaluations": {"optimized": {"global_metrics": final_global_metrics, "local_metrics": final_local_metrics}},
        "positional_plot_div": positional_plot_div,
        "idt_complexity_html": idt_complexity_html,
        "idt_token_error": idt_token_error
    }

    # --- Metainfo Update for Forna ---
    structure_mrna = final_global_metrics.get("structure", "") # Should be from mRNA
    seq_for_forna = outputseq_dict['seq'] # Entire mRNA sequence
    
    metainfo['structure'] = structure_mrna
    
    # Debug logging
    log.info(f"Forna option: {args.forna}")
    log.info(f"Structure available: {bool(structure_mrna)} (length: {len(structure_mrna) if structure_mrna else 0})")
    log.info(f"Sequence available: {bool(seq_for_forna)} (length: {len(seq_for_forna) if seq_for_forna else 0})")
    
    # Generate Forna URL based on server selection
    if args.forna is None:
        # No --forna option provided: don't show structure visualization
        metainfo['forna_url'] = ""
        log.info("No --forna option provided: structure visualization disabled")
    else:
        # --forna option provided: generate URL with structure and sequence
        if args.forna == 'qbio':
            metainfo['forna_url'] = f"https://pub-forna.qbio.io/?id=url/vaxpress&sequence={seq_for_forna}&structure={structure_mrna}"
            log.info("Generated qbio Forna URL")
        elif args.forna == 'tbi':
            metainfo['forna_url'] = f"http://nibiru.tbi.univie.ac.at/forna/forna.html?id=url/name&sequence={seq_for_forna}&structure={structure_mrna}"
            log.info("Generated TBI Forna URL")
    
    metainfo['end_time'] = time.time()

    # --- Instantiate ReportGenerator ---
    # log.info(f"Final metainfo forna_url: {metainfo.get('forna_url', 'NOT SET')}")
    
    try:
        log.debug("Initializing ReportGenerator...")
        generator = ReportGenerator(
            status=status,
            args=args, # Pass parsed args
            metainfo=metainfo,
            scoring_options=scoring_options, # Loaded from preset
            execution_options=execution_options, # Initialized
            inputseq=inputseq_dict,   # CDS info
            outputseq=outputseq_dict, # Entire mRNA info
            scoring_functions=scoring_functions, # Empty, not used by report_only directly
            metric_labels=METRIC_LABELS,
            metric_descriptions=METRIC_DESCRIPTIONS,
            metric_ranges=METRIC_RANGES
        )
        log.debug("ReportGenerator initialized.")
    except Exception as e: print(f"Error initializing ReportGenerator: {e}"); traceback.print_exc(); exit(1)

    # --- Generate Report ---
    try:
        log.debug("Generating report HTML..."); report_html = generator.generate()
        output_html_path = os.path.join(args.output, "report.html")
        with open(output_html_path, "w", encoding="utf-8") as f: f.write(report_html)
        print(f"✅ Report successfully generated: {output_html_path}")
    except Exception as e: print(f"\nError generating report: {e}"); traceback.print_exc(); exit(1)

if __name__ == "__main__":
    main()
