# report_only.py

import argparse
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time
import json
import plotly.graph_objects as go
import plotly.offline as pyo
import numpy as np
from typing import Dict, Any, List, Optional
import traceback
import shlex
from vaxlab_report.log import initialize_logging, log

# vaxlab_report imports
try:
    from vaxlab_report.reporting import ReportGenerator
    from vaxlab_report.presets import load_preset
    try:
        from vaxlab_report.evolution_chamber import ExecutionOptions
    except ImportError:
        # Keep local definition as fallback, remove print
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
    print(f"Error importing vaxlab_report components: {e}"); exit(1) # Keep critical exit error

# VaxPress data/code imports
try: from vaxlab_report.data import codon_usage_data
except ImportError:
    try: import codon_usage_data
    except ImportError: print("Error: Could not import codon_usage_data."); codon_usage_data = None # Keep error print
try: from Bio.Data import CodonTable; standard_table = CodonTable.unambiguous_rna_by_name["Standard"]
except ImportError: print("Warning: Biopython CodonTable not found."); standard_table = None # Keep warning print


# --- Hardcoded Metric Metadata ---
METRIC_LABELS = {
    "gc": "GC ratio", "ucount": "Uridine ratio", "bicodon": "Bicodon score",
    "repeat": "Repeat count", "longstem": "Long stem count",
    "start_str": "Start codon structure", "loop": "Loop Length", "mfe": "Minimum Free Energy",
    "cai": "CAI", "structure": "Secondary Structure", "degscore": "DegScore"
}
METRIC_DESCRIPTIONS = {
    "gc": "Ratio of G/C in the sequence", "degscore": "Predicted degradation score (Eterna DegScore)",
    "ucount": "Ratio of U in the sequence", "bicodon": "Codon pair usage bias",
    "repeat": "Length of longest tandem repeat", "longstem": "Count of long stem-loop structures (≥ 27bp)",
    "start_str": "Base-paired nucleotides near start codon (-14 to +14)",
    "loop": "Total length of unpaired regions (loops ≥ 2nt)", "mfe": "Minimum free energy of predicted structure",
    "cai": "Codon usage optimality based on species-specific frequencies (0-1 scale).",
    "structure": "Predicted secondary structure (dot-bracket)"
}
METRIC_RANGES = {
    "gc": "0.4–0.6", "ucount": "Lowest possible", "bicodon": "Maximize (≥ 0)",
    "repeat": "< 15-20 nt", "longstem": "0–1", "start_str": "≤ 3 bp paired",
    "loop": "≤ 100 per kb", "mfe": "≤ –30 kcal/mol per kb", "cai": "Maximize (closer to 1)",
    "degscore": "Lowest possible", "structure": "N/A"
}

# --- Metainfo, Globals ---
metainfo = {"vaxpress_version": "unknown", "command_line": " ".join(os.sys.argv), "start_time": time.time(), "forna_url": None }
scoring_functions = {}
scoring_options = {}

# --- Helper Functions ---
def calculate_sliding_window_gc(sequence: str, window: int, step: int) -> (list, list):
    """Calculates GC content in sliding windows."""
    positions = []; gc_values = []
    if len(sequence) < window: log.warning(f"Seq length {len(sequence)} < GC window {window}."); return positions, gc_values
    for i in range(0, len(sequence) - window + 1, step):
        sub_seq = sequence[i:i+window].upper(); gc_count = sub_seq.count('G') + sub_seq.count('C')
        gc_ratio = gc_count / window if window > 0 else 0
        positions.append(i + window // 2); gc_values.append(gc_ratio)
    return positions, gc_values

def get_vaxpress_codon_scores(species: str) -> Optional[Dict[str, float]]:
    """Calculates VaxPress-style log-scores for codons based on species usage data."""
    if codon_usage_data is None or standard_table is None: log.error("Cannot get VaxPress codon scores (missing data/table)."); return None
    if not hasattr(codon_usage_data, 'codon_usage'): log.error("codon_usage_data missing 'codon_usage' attr."); return None
    if species not in codon_usage_data.codon_usage: log.error(f'Species "{species}" not in codon_usage_data.'); return None
    log.debug(f"Loading VaxPress codon usage data for {species}..."); species_codon_usage = codon_usage_data.codon_usage[species]
    codon_scores = {}; aa_to_codons = {}
    for codon, aa in standard_table.forward_table.items(): codon_rna = codon.replace('T', 'U'); aa_to_codons.setdefault(aa, []).append(codon_rna)
    stop_codons_rna = [c.replace('T', 'U') for c in standard_table.stop_codons]
    log.debug("Calculating VaxPress codon log-scores...");
    for aa, codons in aa_to_codons.items():
        valid_codons = [c for c in codons if c in species_codon_usage];
        if not valid_codons: continue
        freqs = np.array([species_codon_usage[c] for c in valid_codons]); max_freq = freqs.max()
        if max_freq <= 0: log_rel_scores = np.zeros(len(valid_codons))
        else: epsilon = 1e-9; relative_freqs = (freqs + epsilon) / max_freq; log_rel_scores = np.log(relative_freqs); log_rel_scores[np.isneginf(log_rel_scores)] = np.log(epsilon / max_freq if max_freq > 0 else epsilon)
        codon_scores.update(dict(zip(valid_codons, log_rel_scores)))
    all_codons = [c.replace('T','U') for c in standard_table.forward_table.keys()] + stop_codons_rna
    for codon in all_codons: codon_scores.setdefault(codon, np.nan)
    log.debug(f"Finished calculating {len(codon_scores)} VaxPress codon scores.")
    return codon_scores

def calculate_sliding_window_vaxpress_cai(sequence: str, window_codons: int, step_codons: int, codon_scores: Dict[str, float]) -> (list, list):
    """Calculates mean VaxPress log-score in sliding windows."""
    positions = []; cai_scores = []
    seq_len = len(sequence); cds_seq = sequence[:seq_len - (seq_len % 3)].upper().replace('T', 'U'); cds_len = len(cds_seq)
    if not codon_scores: log.error("Empty codon scores for VaxPress CAI calc."); return positions, cai_scores
    if cds_len < window_codons * 3: log.warning(f"CDS length {cds_len} < CAI window {window_codons} codons."); return positions, cai_scores
    window_nt = window_codons * 3; step_nt = step_codons * 3
    log.debug(f"Starting VaxPress CAI sliding window calc ({window_codons}c window, {step_codons}c step)...")
    for i in range(0, cds_len - window_nt + 1, step_nt):
        sub_seq = cds_seq[i:i+window_nt]
        window_log_scores = [codon_scores.get(sub_seq[j:j+3], np.nan) for j in range(0, len(sub_seq), 3)]
        mean_log_score = np.nanmean(window_log_scores)
        positions.append(i + window_nt // 2); cai_scores.append(mean_log_score)
    log.debug(f"Finished VaxPress CAI calculation ({len(positions)} points).")
    return positions, cai_scores

def generate_positional_plot(positions: list, values: dict, title: str) -> Optional[str]:
    """Generates a Plotly line plot for positional metrics."""
    if not positions or not values: log.warning("No data for positional plot."); return None
    fig = go.Figure(); valid_trace_added = False
    for name, data in values.items():
        if data is not None and len(data) == len(positions):
             y_values = [v if not (v is None or np.isnan(v)) else None for v in data]
             fig.add_trace(go.Scatter(x=positions, y=y_values, mode='lines', name=name, connectgaps=False)); valid_trace_added = True; log.debug(f"Added trace '{name}' to plot.")
        elif data is not None: log.warning(f"Skipping trace '{name}' (length mismatch {len(data)} vs {len(positions)}).")
        else: log.warning(f"Skipping trace '{name}' (data is None).")
    if not valid_trace_added: log.warning("No valid traces added to plot."); return None
    fig.update_layout(title=title, xaxis_title="Position (approx. window center, nt)", yaxis_title="Metric Value", hovermode="x unified")
    try: return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    except Exception as e: log.error(f"Error converting Plotly figure to div: {e}"); traceback.print_exc(); return None

def calculate_global_vaxpress_cai(sequence: str, codon_scores: Dict[str, float]) -> Optional[float]:
    """Calculates the overall mean VaxPress log-score (Global CAI) for the sequence."""
    seq_len = len(sequence)
    cds_seq = sequence[:seq_len - (seq_len % 3)].upper().replace('T', 'U')
    cds_len = len(cds_seq)

    if not codon_scores:
        log.error("Cannot calculate global CAI: Empty codon scores provided.")
        return None
    if cds_len == 0:
        log.warning("Cannot calculate global CAI: CDS length is zero.")
        return None

    log.info(f"Calculating Global VaxPress CAI for sequence length {cds_len}...")
    all_log_scores = [codon_scores.get(cds_seq[j:j+3], np.nan) for j in range(0, cds_len, 3)]
    mean_log_score = np.nanmean(all_log_scores)

    if np.isnan(mean_log_score):
        log.warning("Global CAI calculation resulted in NaN (possibly no valid codons found).")
        return None
    else:
        log.info(f"Global VaxPress CAI calculated: {mean_log_score:.4f}")
        return float(mean_log_score)

# --- Main Function ---
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    parser.add_argument("--evaluations", required=True, help="TSV file with evaluation metrics")
    parser.add_argument("--preset", help="Preset JSON file for options")
    parser.add_argument('--iters', type=int, default=0, help='Placeholder (n_iterations)')
    parser.add_argument('--population', type=int, default=0, help='Placeholder (n_population)')
    parser.add_argument('--cpu-count', type=int, default=1, help='Placeholder (processes)')
    parser.add_argument('--random-seed', type=int, default=None, help='Placeholder (seed)')
    parser.add_argument('--species', type=str, default='Homo sapiens', help='Species for CAI calculation')
    args = parser.parse_args()
    log.debug(f"Parsed args: {args}")

    # --- Load Preset and Options ---
    execution_options_dict_from_preset = {}
    if args.preset:
        try: preset_data = load_preset(args.preset); scoring_options.update(preset_data.get("scoring", {})); execution_options_dict_from_preset.update(preset_data.get("execution", {})); log.debug(f"Loaded preset '{args.preset}'")
        except FileNotFoundError: print(f"Error: Preset file not found: {args.preset}"); exit(1)
        except json.JSONDecodeError: print(f"Error: Could not decode JSON from preset file: {args.preset}"); exit(1)
        except Exception as e: print(f"Error loading preset '{args.preset}': {e}"); exit(1)

    # --- Load Input Sequence ---
    try:
        inputseq = SeqIO.read(args.input, "fasta")
        input_seq_description = inputseq.description
        outputseq = SeqRecord(inputseq.seq, id=inputseq.id + "_reported", description="Reported sequence");
        inputseq_dict = {'id': inputseq.id, 'seq': str(inputseq.seq)};
        outputseq_dict = {'id': outputseq.id, 'seq': str(outputseq.seq), 'description': outputseq.description};
        log.debug(f"Loaded sequence '{inputseq.id}', length {len(inputseq.seq)}")
    except FileNotFoundError: print(f"Error: Input FASTA file not found: {args.input}"); exit(1)
    except Exception as e: print(f"Error reading input FASTA file: {e}"); exit(1)

    # --- Create ExecutionOptions instance ---
    init_args = {
        'n_iterations': args.iters, 'n_population': args.population, 'n_survivors': 1,
        'initial_mutation_rate': 0.1, 'winddown_trigger': 15, 'winddown_rate': 0.9,
        'output': args.output, 'command_line': "report_only execution", 'overwrite': True,
        'seed': args.random_seed if args.random_seed is not None else 0,
        'processes': args.cpu_count, 'random_initialization': False, 'conservative_start': None,
        'boost_loop_mutations': "1.5:15", 'full_scan_interval': 300, 'species': args.species,
        'codon_table': "standard", 'protein': False, 'quiet': True,
        'seq_description': input_seq_description,
        'print_top_mutants': 0, 'addons': [], 'lineardesign_dir': None,
        'lineardesign_lambda': None, 'lineardesign_omit_start': 5, 'folding_engine': "vienna"
    }
    init_args.update(execution_options_dict_from_preset)
    init_args['output'] = args.output; init_args['n_iterations'] = args.iters; init_args['n_population'] = args.population; init_args['processes'] = args.cpu_count
    if args.random_seed is not None: init_args['seed'] = args.random_seed
    init_args['species'] = args.species

    try: execution_options = ExecutionOptions(**init_args); log.debug(f"Execution Options Initialized.")
    except TypeError as e: print(f"Error: Failed to initialize ExecutionOptions: {e}"); exit(1)

    log.debug("Using hardcoded metric metadata.")

    # --- Load Evaluation Data ---
    checkpoints_list = []; final_global_metrics = {}; final_local_metrics = {}
    try:
        log.debug(f"Loading metrics from TSV file: {args.evaluations}")
        df = pd.read_csv(args.evaluations, sep='\t')
        if not df.empty:
            last_row = df.iloc[-1]; final_global_metrics = {k.replace('metric:', ''): v for k, v in last_row.items() if k.startswith('metric:')}
            for k, v in final_global_metrics.items():
                 if isinstance(v, str):
                      try: final_global_metrics[k] = int(v)
                      except ValueError:
                           try: final_global_metrics[k] = float(v)
                           except ValueError: pass
            checkpoints_list = df.to_dict(orient='records')
            final_sequence = last_row.get("sequence", str(inputseq.seq)); outputseq_dict['seq'] = final_sequence; log.debug(f"Loaded metrics and sequence from TSV.")
        else: log.warning(f"TSV file '{args.evaluations}' is empty. Trying JSON.")
    except FileNotFoundError: print(f"Error: Evaluations TSV file not found: {args.evaluations}"); exit(1)
    except Exception as e: print(f"Error reading or processing evaluations TSV file: {e}"); exit(1)

    # Load local metrics from JSON
    json_eval_path = os.path.join(args.output, "evaluation_result.json")
    if os.path.exists(json_eval_path):
        try:
            with open(json_eval_path, 'r') as f: evaluation_result = json.load(f)
            final_local_metrics = evaluation_result.get("local_metrics", {})
            log.debug(f"Loaded local metrics from {json_eval_path}")
            if not final_global_metrics:
                 final_global_metrics = evaluation_result.get("global_metrics", {})
                 log.debug("Loaded global metrics from JSON as TSV was empty.")
                 for k, v in final_global_metrics.items():
                      if isinstance(v, str):
                           try: final_global_metrics[k] = int(v)
                           except ValueError:
                                try: final_global_metrics[k] = float(v)
                                except ValueError: pass
        except Exception as e: log.warning(f"Could not load or parse {json_eval_path}. Error: {e}")

    if not final_global_metrics: print("Error: Could not load any global evaluation metrics."); exit(1)

    sequence_str = str(outputseq_dict.get('seq', ''))
    if not sequence_str:
        log.critical("Error: Sequence string is empty after loading data.")
        exit(1)

    # --- Calculate Global CAI ---
    log.info(f"Attempting to calculate standard CAI (0-1 scale) for species: {args.species}")
    vaxpress_codon_scores = get_vaxpress_codon_scores(args.species)
    if vaxpress_codon_scores:
        mean_log_score = calculate_global_vaxpress_cai(sequence_str, vaxpress_codon_scores)
        if mean_log_score is not None:
            standard_cai_value = np.exp(mean_log_score)
            final_global_metrics['cai'] = float(standard_cai_value)
            log.info(f"Added calculated standard CAI ({standard_cai_value:.4f}) to final_global_metrics.")
        else:
            log.warning("Mean log score for CAI could not be calculated. Standard CAI will be missing.")
            final_global_metrics.pop('cai', None)
    else:
        log.error("Could not get VaxPress codon scores. CAI calculation skipped.")
        final_global_metrics.pop('cai', None)

    # --- Calculate Positional Metrics & Generate Plot ---
    positional_plot_div = None
    plot_data = {}
    all_positions = []
    sequence_str = str(outputseq_dict.get('seq', ''))

    log.info("--- Processing Positional Metrics for Plot ---")

    # 1. Load and validate local metrics data
    metrics_to_process: Dict[str, List[List[Any]]] = {}
    if isinstance(final_local_metrics, dict):
        for key, metric_data in final_local_metrics.items():
            if key.endswith("_error"):
                 log.error(f"Found error key in local metrics: {key} = {metric_data}")
                 continue
            if isinstance(metric_data, (list, tuple)) and len(metric_data) == 2 and \
               all(isinstance(el, (list, tuple)) for el in metric_data):
                positions, values = metric_data
                if positions and values and len(positions) == len(values):
                    try:
                        num_positions = [p for p in positions]
                        num_values = [float(v) if v is not None else float('nan') for v in values]
                        metrics_to_process[key] = [num_positions, num_values]
                        log.debug(f"Successfully loaded and validated data for '{key}' ({len(positions)} points).")
                    except (ValueError, TypeError) as e:
                         log.warning(f"Could not convert values to numeric for local metric '{key}'. Skipping. Error: {e}")
                else:
                    log.warning(f"Incomplete data (positions or values empty/mismatched) for local metric '{key}'. Skipping.")
            elif metric_data is not None:
                 log.warning(f"Unexpected data format for local metric '{key}': {type(metric_data)}. Skipping.")
    else:
        log.warning("final_local_metrics is not a dictionary. Cannot process positional metrics.")

    # 2. Set base positions (x-axis)
    priority_keys = ["gc", "cai", "degscore"]
    base_metric_key = None
    for key in priority_keys:
        if key in metrics_to_process:
            base_metric_key = key
            break
    if not base_metric_key and metrics_to_process:
        base_metric_key = list(metrics_to_process.keys())[0]

    if base_metric_key:
        all_positions = metrics_to_process[base_metric_key][0]
        log.info(f"Using positions from '{base_metric_key}' as the primary x-axis ({len(all_positions)} points).")
    else:
        log.warning("No valid positional metric data found to establish x-axis for the plot.")

    # 3. Prepare plot data (including interpolation and CAI transformation)
    if all_positions:
        processed_keys_for_plot = set()
        plot_order_preference = ["gc", "cai", "degscore"]
        plot_keys_ordered = [k for k in plot_order_preference if k in metrics_to_process] + \
                            [k for k in metrics_to_process if k not in plot_order_preference]

        for key in plot_keys_ordered:
            if key not in metrics_to_process or key in processed_keys_for_plot: continue
            if key == 'ucount':
                log.debug("Skipping 'ucount' for positional plotting as it's not needed/available in positional format.")
                continue
            positions, values = metrics_to_process[key]
            plot_label = METRIC_LABELS.get(key, key)
            plot_values = None
            current_values_to_process = values

            if key == 'cai':
                log.debug(f"Exponentiating local log scores for '{plot_label}' (original key: {key}) for plotting.")
                try:
                    current_values_to_process = [np.exp(v) if v is not None and not np.isnan(v) else None for v in values]
                except Exception as exp_e:
                    log.error(f"Error exponentiating CAI values for plot: {exp_e}", exc_info=True)
                    current_values_to_process = [None] * len(values)

            if list(positions) != list(all_positions):
                log.debug(f"Interpolating values for '{plot_label}' (original key: {key})...")
                try:
                    pos_arr = np.array(positions)
                    val_arr = np.array([float(v) if v is not None else np.nan for v in current_values_to_process])
                    valid_mask = ~np.isnan(val_arr)

                    if np.any(valid_mask):
                         interp_positions = pos_arr[valid_mask]
                         interp_values = val_arr[valid_mask]
                         plot_values = np.interp(all_positions, interp_positions, interp_values).tolist()
                         log.debug(f"Interpolation successful for '{plot_label}'.")
                    else:
                         log.warning(f"Cannot interpolate '{plot_label}' (all NaN or invalid values).")
                         plot_values = [None] * len(all_positions)
                except Exception as interp_e:
                    log.error(f"Error during interpolation for '{plot_label}': {interp_e}", exc_info=True)
                    plot_values = None
            else:
                plot_values = [v if v is not None and not np.isnan(v) else None for v in current_values_to_process]

            if plot_values is not None:
                 plot_data[plot_label] = plot_values
                 processed_keys_for_plot.add(key)
                 log.debug(f"Added '{plot_label}' to plot_data.")
            else:
                 log.warning(f"Could not determine final plot values for '{plot_label}'. Skipping trace.")

    # --- Generate Plotly graph ---
    if all_positions and plot_data:
        plot_title = f"Positional Metrics"
        positional_plot_div = generate_positional_plot(all_positions, plot_data, plot_title)
        if positional_plot_div:
            log.info("Positional plot generated successfully.")
        else:
            log.warning("Positional plot generation failed (generate_positional_plot returned None).")
    else:
        log.warning("No valid data available to generate positional plot.")

    log.info("--- Finished Positional Metric Calculation & Plot Generation ---")

    # --- Construct Status Dictionary ---
    status: Dict[str, Any] = {
        "checkpoints": checkpoints_list,
        "evaluations": {"optimized": {"global_metrics": final_global_metrics, "local_metrics": final_local_metrics}},
        "positional_plot_div": positional_plot_div
    }

    # --- Metainfo Update ---
    structure = final_global_metrics.get("structure", "")
    seq_for_forna = sequence_str
    metainfo['structure'] = structure
    metainfo['forna_url'] = f"https://pub-forna.qbio.io/?id=url/vaxpress&sequence={seq_for_forna}&structure={structure}" if structure and seq_for_forna else ""
    metainfo['end_time'] = time.time()

    # --- Instantiate ReportGenerator ---
    try:
        log.debug("Initializing ReportGenerator...")
        generator = ReportGenerator(
            status=status,
            args=args,
            metainfo=metainfo,
            scoring_options=scoring_options,
            execution_options=execution_options,
            inputseq=inputseq_dict,
            outputseq=outputseq_dict,
            scoring_functions=scoring_functions,
            metric_labels=METRIC_LABELS,
            metric_descriptions=METRIC_DESCRIPTIONS,
            metric_ranges=METRIC_RANGES
        )
        log.debug("ReportGenerator initialized.")
    except Exception as e: print(f"Error initializing ReportGenerator: {e}"); traceback.print_exc(); exit(1)

    # --- Generate Report ---
    try:
        log.debug("Generating report HTML..."); report_html = generator.generate()
        output_html_path = os.path.join(args.output, "report.html"); os.makedirs(args.output, exist_ok=True)
        with open(output_html_path, "w", encoding="utf-8") as f: f.write(report_html)
        print(f"✅ Report successfully generated: {output_html_path}") # Keep success message
    except Exception as e: print(f"\nError generating report: {e}"); traceback.print_exc(); exit(1)

if __name__ == "__main__":
    main()