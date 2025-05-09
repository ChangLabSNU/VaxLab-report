import argparse
import os
# import pandas as pd # checkpoints_list 로딩에 pandas를 사용하지 않는다면 제거 가능
# import csv # checkpoints_list 로딩에 csv 모듈을 사용하지 않는다면 제거 가능 (현재는 사용 안함)
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

try:
    from vaxlab_report.reporting import ReportGenerator
    from vaxlab_report.presets import load_preset
    try:
        from vaxlab_report.evolution_chamber import ExecutionOptions
    except ImportError:
        class ExecutionOptions: # Fallback definition
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
    "repeat": "Length of longest tandem repeat", "longstem": "Count of long stem-loop structures (≥ 27bp)",
    "start_str": "Base-paired nucleotides near start codon (0 to +14)",
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

metainfo = {"vaxpress_version": "unknown", "command_line": " ".join(os.sys.argv), "start_time": time.time(), "forna_url": None }
scoring_functions = {} 
scoring_options = {}

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
    parser.add_argument("--preset", help="Preset JSON file for options")
    parser.add_argument('--iters', type=int, default=0, help='Placeholder (n_iterations)')
    parser.add_argument('--population', type=int, default=0, help='Placeholder (n_population)')
    parser.add_argument('--cpu-count', type=int, default=1, help='Placeholder (processes)')
    parser.add_argument('--random-seed', type=int, default=None, help='Placeholder (seed)')
    parser.add_argument('--species', type=str, default='Homo sapiens', help='Species for CAI calculation')
    args = parser.parse_args()

    log_file = os.path.join(args.output, "report_only_log.txt")
    initialize_logging(log_file, quiet=False)
    log.debug(f"Parsed args: {args}")

    execution_options_dict_from_preset = {}
    if args.preset:
        try:
            preset_data = load_preset(args.preset)
            scoring_options.update(preset_data.get("scoring", {}))
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
            else:
                if num_potential_cds == 0:
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
    outputseq_dict = {
        'id': f"{cds_record_id}",
        'seq': entire_mrna_sequence,
        'description': f"{cds_record_desc}"
    }
    inputseq_dict = {'id': cds_record_id, 'seq': cds_seq_str, 'description': cds_record_desc}
    log.info(f"Constructed entire mRNA: ID={outputseq_dict['id']}, Length={len(outputseq_dict['seq'])}")

    init_args = {
        'n_iterations': args.iters, 
        'n_population': args.population, 
        'n_survivors': 1,
        'output': args.output, 
        'command_line': "report_only execution", 
        'overwrite': True,
        'seed': args.random_seed if args.random_seed is not None else 0,
        'processes': args.cpu_count, 
        'species': args.species, 
        'codon_table': "standard",
        'quiet': True, 
        'seq_description': outputseq_dict['description'],
        'initial_mutation_rate': 0.1,
        'winddown_trigger': 15,
        'winddown_rate': 0.9,
        'random_initialization': False,
        'conservative_start': None,
        'boost_loop_mutations': "1.5:15",
        'full_scan_interval': 0, 
        'protein': False,
        'print_top_mutants': 0,
        'addons': [], 
        'lineardesign_dir': None,
        'lineardesign_lambda': None,
        'lineardesign_omit_start': 5,
        'folding_engine': "vienna", 
    }
    init_args.update(execution_options_dict_from_preset) # Preset 값으로 덮어쓰기
    # 명령줄 인자로 다시 덮어쓰기 (우선순위 가장 높음)
    init_args['output'] = args.output; init_args['n_iterations'] = args.iters;
    init_args['n_population'] = args.population; init_args['processes'] = args.cpu_count
    if args.random_seed is not None: init_args['seed'] = args.random_seed
    init_args['species'] = args.species

    try: execution_options = ExecutionOptions(**init_args); log.debug("Execution Options Initialized.")
    except Exception as e: print(f"Error: Failed to initialize ExecutionOptions: {e}"); exit(1)

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

    final_global_metrics = mrna_evaluation_result.get("global_metrics", {}).copy()
    cds_global_metrics = cds_evaluation_result.get("global_metrics", {})
    
    metrics_to_override_from_cds = {
        "cai": "cai", 
        "bicodon": "bicodon",
        "start_str": "start_str"
    }

    for cds_key, final_key in metrics_to_override_from_cds.items():
        if cds_key in cds_global_metrics:
            final_global_metrics[final_key] = cds_global_metrics[cds_key]
            log.info(f"Overridden global metric '{final_key}' with value from CDS results: {cds_global_metrics[cds_key]}")
        else:
            log.warning(f"Metric '{cds_key}' not found in CDS global metrics. Cannot use for overriding.")
            if final_key not in final_global_metrics:
                 final_global_metrics[final_key] = None

    final_local_metrics = mrna_evaluation_result.get("local_metrics", {})
    log.debug(f"Final global metrics for report: {final_global_metrics}")

    checkpoints_list = [] # --evaluations 옵션 제거로 인해 항상 빈 리스트
    log.info("Skipping checkpoints.tsv loading as the --evaluations option is removed.")

    positional_plot_div = None
    plot_data = {}
    all_positions = []
    
    log.info("--- Processing Positional Metrics for Plot (mRNA based, CAI excluded) ---")
    metrics_to_process_for_plot: Dict[str, List[List[Any]]] = {}
    if isinstance(final_local_metrics, dict):
        for key, metric_data in final_local_metrics.items():
            if key == 'cai' or METRIC_LABELS.get(key, key).lower() == 'cai':
                log.info("Excluding 'CAI' from positional plot as it's not available or requested for mRNA context.")
                continue
            if key.endswith("_error"): continue
            if isinstance(metric_data, (list, tuple)) and len(metric_data) == 2 and \
               all(isinstance(el, (list, tuple)) for el in metric_data):
                positions, values = metric_data
                if positions and values and len(positions) == len(values):
                    try:
                        num_positions = [p for p in positions]
                        num_values = [float(v) if v is not None else float('nan') for v in values]
                        metrics_to_process_for_plot[key] = [num_positions, num_values]
                    except (ValueError, TypeError) as e:
                         log.warning(f"Could not convert values for local metric '{key}' for plot. Error: {e}")
                else: log.warning(f"Incomplete data for local metric '{key}' for plot.")
            elif metric_data is not None:
                 log.warning(f"Unexpected data format for local metric '{key}' for plot: {type(metric_data)}.")
    else:
        log.warning("final_local_metrics (mRNA) is not a dictionary. Cannot process for positional plot.")

    priority_keys_for_x_axis = ["gc", "degscore"]
    base_metric_key_for_x_axis = None
    for key in priority_keys_for_x_axis:
        if key in metrics_to_process_for_plot:
            base_metric_key_for_x_axis = key; break
    if not base_metric_key_for_x_axis and metrics_to_process_for_plot:
        base_metric_key_for_x_axis = list(metrics_to_process_for_plot.keys())[0]

    if base_metric_key_for_x_axis:
        all_positions = metrics_to_process_for_plot[base_metric_key_for_x_axis][0]
        log.info(f"Using positions from '{base_metric_key_for_x_axis}' as x-axis for plot ({len(all_positions)} points).")
    else:
        log.warning("No valid positional metric data from mRNA results to establish x-axis for the plot.")

    if all_positions:
        processed_keys_for_plot_data = set()
        plot_order_preference = ["gc", "degscore"]
        plot_keys_ordered = [k for k in plot_order_preference if k in metrics_to_process_for_plot] + \
                            [k for k in metrics_to_process_for_plot if k not in plot_order_preference]

        for key in plot_keys_ordered:
            if key in processed_keys_for_plot_data: continue
            
            positions, values = metrics_to_process_for_plot[key]
            plot_label = METRIC_LABELS.get(key, key)
            plot_values = None
            current_values_to_process = values

            if list(positions) != list(all_positions):
                log.debug(f"Interpolating values for '{plot_label}' (original key: {key})...")
                try:
                    pos_arr = np.array(positions)
                    val_arr = np.array([float(v) if v is not None else np.nan for v in current_values_to_process])
                    valid_mask = ~np.isnan(val_arr)
                    if np.any(valid_mask):
                         interp_positions = pos_arr[valid_mask]
                         interp_values = val_arr[valid_mask]
                         plot_values = np.interp(all_positions, interp_positions, interp_values, left=np.nan, right=np.nan).tolist()
                    else:
                         plot_values = [None] * len(all_positions)
                except Exception as interp_e:
                    log.error(f"Error during interpolation for '{plot_label}': {interp_e}", exc_info=True)
                    plot_values = [None] * len(all_positions)
            else:
                plot_values = [v if v is not None and not np.isnan(v) else None for v in current_values_to_process]

            if plot_values is not None: plot_data[plot_label] = plot_values
            processed_keys_for_plot_data.add(key)
    
    if all_positions and plot_data:
        plot_title = f"Positional Metrics"
        positional_plot_div = generate_positional_plot(all_positions, plot_data, plot_title)
        if positional_plot_div: log.info("Positional plot for entire mRNA generated.")
        else: log.warning("Positional plot generation failed.")
    else:
        log.warning("No valid data to generate positional plot for entire mRNA.")

    status: Dict[str, Any] = {
        "checkpoints": checkpoints_list,
        "evaluations": {"optimized": {"global_metrics": final_global_metrics, "local_metrics": final_local_metrics}},
        "positional_plot_div": positional_plot_div
    }

    structure_mrna = final_global_metrics.get("structure", "")
    seq_for_forna = outputseq_dict['seq']
    metainfo['structure'] = structure_mrna
    metainfo['forna_url'] = f"https://pub-forna.qbio.io/?id=url/vaxpress&sequence={seq_for_forna}&structure={structure_mrna}" if structure_mrna and seq_for_forna else ""
    metainfo['end_time'] = time.time()

    try:
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

    try:
        report_html = generator.generate()
        output_html_path = os.path.join(args.output, "report.html")
        with open(output_html_path, "w", encoding="utf-8") as f: f.write(report_html)
        print(f"✅ Report successfully generated: {output_html_path}")
    except Exception as e: print(f"\nError generating report: {e}"); traceback.print_exc(); exit(1)

if __name__ == "__main__":
    main()