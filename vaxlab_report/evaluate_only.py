# evaluate_only.py

import argparse
import os
import csv
import random
import shlex
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
# import subprocess # subprocess는 현재 사용되지 않으므로 주석 처리 또는 제거 가능
import json
import numpy as np
from typing import Dict, List, Tuple, Any, Optional
import copy

from vaxlab_report import config, scoring, __version__
from vaxlab_report.sequence_evaluator import SequenceEvaluator
try:
    from vaxlab_report.evolution_chamber import ExecutionOptions
except ImportError as e:
     print(f"FATAL ERROR: Failed to import ExecutionOptions from vaxlab_report.evolution_chamber: {e}")
     print("        Please ensure this module and class exist and the Python path is correct.")
     exit(1)
from vaxlab_report.mutant_generator import MutantGenerator
from vaxlab_report.presets import load_preset
from vaxlab_report.log import initialize_logging, log
from vaxlab_report.scoring.idt_complexity import evaluate_idt_gblock_complexity


def evaluate_idt_complexity_for_sequence(sequence: str, sequence_id: str, idt_token: str, log) -> Optional[Dict[str, Any]]:
    """
    Evaluate IDT complexity for a given sequence.
    
    Args:
        sequence: RNA sequence to evaluate
        sequence_id: ID for the sequence
        idt_token: IDT API access token
        log: Logger instance
    
    Returns:
        Dictionary with IDT complexity results or None if evaluation fails
    """
    if not idt_token:
        log.info("IDT API token not provided. Skipping IDT complexity evaluation.")
        return None
    
    # Convert RNA to DNA for IDT API
    dna_sequence = sequence.replace('U', 'T')
    
    try:
        log.info(f"Evaluating IDT complexity for {sequence_id}...")
        results = evaluate_idt_gblock_complexity(idt_token, {sequence_id: dna_sequence})
        
        if results and len(results) > 0:
            # Process the first result
            result = results[0]
            log.info(f"IDT complexity evaluation successful for {sequence_id}")
            return result
        else:
            log.warning(f"No IDT complexity results returned for {sequence_id}")
            return None
            
    except Exception as e:
        log.error(f"Error evaluating IDT complexity for {sequence_id}: {e}", exc_info=True)
        # Return error information instead of None
        return {"error": str(e)}


def convert_numpy_to_list(data: Any) -> Any:
    """Recursively converts NumPy types in data structures to JSON serializable types."""
    if isinstance(data, dict):
        return {k: convert_numpy_to_list(v) for k, v in data.items()}
    elif isinstance(data, (list, tuple)):
        if len(data) == 2 and all(isinstance(el, (list, tuple, np.ndarray)) for el in data):
             idx_list = data[0].tolist() if isinstance(data[0], np.ndarray) else list(data[0])
             val_list = data[1].tolist() if isinstance(data[1], np.ndarray) else list(data[1])
             return [idx_list, val_list]
        else:
             return [convert_numpy_to_list(item) for item in data]
    elif isinstance(data, np.ndarray):
        return data.tolist()
    elif isinstance(data, (np.int8, np.int16, np.int32, np.int64, np.uint8, np.uint16, np.uint32, np.uint64)):
        return int(data)
    elif isinstance(data, (np.float16, np.float32, np.float64)):
        if np.isnan(data): return None
        if np.isinf(data): return "Infinity" if data > 0 else "-Infinity"
        return float(data)
    elif isinstance(data, (np.complex64, np.complex128)):
        return {'real': data.real, 'imag': data.imag}
    elif isinstance(data, np.bool_):
         return bool(data)
    elif isinstance(data, np.void):
         return None
    else:
        return data

def run_evaluation(
    sequence_to_evaluate: str,      # 실제 평가 대상 시퀀스 (CDS 또는 전체 mRNA)
    original_cds_sequence_for_mutgen: str, # MutantGenerator 초기화용 원본 CDS
    sequence_id: str,
    sequence_description: str,
    output_suffix: str,
    execopts: ExecutionOptions,
    base_scoring_funcs: Dict[str, callable],
    base_scoring_options: Dict,
    metrics_to_skip_user_names: List[str],
    args_namespace: argparse.Namespace,
    random_state_obj: random.Random,
    original_cds_len: int, # SequenceEvaluator에 전달될 원본 CDS 길이
    executor: ThreadPoolExecutor,
    idt_token: Optional[str] = None
):
    """
    Runs a single evaluation pipeline for a given sequence and configuration.
    """
    log.info(f"--- Starting evaluation for {sequence_id} ({output_suffix.strip('_').upper()}) ---")
    log.info(f"Sequence to evaluate length: {len(sequence_to_evaluate)}")
    log.info(f"Original CDS sequence for MutantGenerator length: {len(original_cds_sequence_for_mutgen)}")
    log.info(f"Original CDS length for SequenceEvaluator: {original_cds_len}")


    metric_name_to_key_map = {
        "cai": "cai",
        "bicodon score": "bicodon",
        "start codon structure": "start_str"
    }

    excluded_actual_keys = []
    if metrics_to_skip_user_names:
        log.info(f"Attempting to exclude metrics for {output_suffix}: {metrics_to_skip_user_names}")
        for user_metric_name in metrics_to_skip_user_names:
            key = metric_name_to_key_map.get(user_metric_name.lower())
            if key and key in base_scoring_funcs:
                excluded_actual_keys.append(key)
                log.info(f"Will exclude metric with key: '{key}' (from '{user_metric_name}')")
            elif key:
                log.warning(f"Metric key '{key}' (from '{user_metric_name}') not found in discovered scoring functions. Cannot exclude.")
            else:
                log.warning(f"User metric name '{user_metric_name}' not mapped to a known key. Cannot exclude.")

    current_scoring_funcs = {
        k: v for k, v in base_scoring_funcs.items() if k not in excluded_actual_keys
    }
    if excluded_actual_keys:
        log.info(f"Effective scoring functions for this run: {list(current_scoring_funcs.keys())}")

    current_scoring_options = copy.deepcopy(base_scoring_options)
    for key_to_remove in excluded_actual_keys:
        if "weights" in current_scoring_options and key_to_remove in current_scoring_options["weights"]:
            current_scoring_options["weights"][key_to_remove] = 0.0
            log.info(f"Set weight for '{key_to_remove}' to 0.0 in current_scoring_options.")
        if "active_scores" in current_scoring_options and key_to_remove in current_scoring_options["active_scores"]:
            current_scoring_options["active_scores"][key_to_remove] = False
            log.info(f"Set active_score for '{key_to_remove}' to False in current_scoring_options.")

    log.info("Initializing MutantGenerator for this run...")
    try:
        # MutantGenerator는 항상 원본 CDS 서열로 초기화
        mutantgen = MutantGenerator(original_cds_sequence_for_mutgen, random_state_obj)
        log.info("MutantGenerator initialized successfully.")
    except ValueError as ve:
        log.critical(f"ValueError during MutantGenerator initialization for {output_suffix} with sequence length {len(original_cds_sequence_for_mutgen)}: {ve}", exc_info=True)
        log.critical("This might indicate the original CDS sequence is not valid (e.g., not multiple of 3).")
        return
    except Exception as e:
        log.critical(f"Failed to initialize MutantGenerator for {output_suffix}: {e}", exc_info=True)
        return

    log.info("Initializing SequenceEvaluator for this run...")
    try:
        # SequenceEvaluator 생성 시 len_cds 키워드 인자 대신 위치 인자로 original_cds_len 전달
        evaluator = SequenceEvaluator(
            current_scoring_funcs,
            current_scoring_options,
            execopts,
            mutantgen, # 원본 CDS로 초기화된 mutantgen
            args_namespace.species,
            original_cds_len, # 6번째 위치 인자로 원본 CDS 길이 전달
            quiet=False
        )
        log.info("SequenceEvaluator initialized successfully for this run.")
    except Exception as e:
         log.critical(f"Failed to initialize SequenceEvaluator for {output_suffix}: {e}", exc_info=True)
         return

    log.info(f"Starting sequence evaluation for {sequence_id}...")
    global_metrics: Dict[str, Any] = {}
    local_metrics_output: Dict[str, Any] = {}
    evaluation_successful = False

    try:
        results = evaluator.evaluate([sequence_to_evaluate], executor) # 평가 대상 시퀀스 전달

        if results is not None and len(results) == 5:
            total_scores, scores_list, metrics_list, foldings_list, local_metrics_list = results
            if metrics_list: global_metrics = metrics_list[0]
            else: log.warning("SequenceEvaluator did not return global metrics.")
            if local_metrics_list: local_metrics_output = local_metrics_list[0]
            else: log.warning("SequenceEvaluator did not return local metrics.")
            
            # Add structure information from foldings to global_metrics
            if foldings_list and len(foldings_list) > 0 and foldings_list[0]:
                folding_data = foldings_list[0]
                if isinstance(folding_data, dict) and 'folding' in folding_data:
                    global_metrics['structure'] = folding_data['folding']
                    log.info(f"Added structure information to global_metrics (length: {len(folding_data['folding'])})")
                else:
                    log.warning("Folding data does not contain 'folding' key")
            else:
                log.warning("No folding information available to add to global_metrics")
            
            evaluation_successful = True
            log.info("Evaluation process completed for this run.")
        else:
            log.error("evaluator.evaluate() returned None or unexpected data. Evaluation failed for this run.")
            global_metrics = {'error': 'Evaluation failed or returned None'}
            local_metrics_output = {'error': 'Evaluation failed or returned None'}
    except Exception as e:
        log.critical(f"An critical exception occurred during evaluator.evaluate() for {output_suffix}: {e}", exc_info=True)
        global_metrics = {'error': f'Evaluation exception: {e}'}
        local_metrics_output = {'error': f'Evaluation exception: {e}'}

    # Evaluate IDT complexity if token is provided
    idt_complexity_result = None
    if idt_token and output_suffix == "_cds":  # Only evaluate IDT complexity for CDS
        idt_complexity_result = evaluate_idt_complexity_for_sequence(
            sequence_to_evaluate,
            sequence_id,
            idt_token,
            log
        )
    
    log.debug("Converting results for JSON serialization...")
    try:
        serializable_global_metrics = convert_numpy_to_list(global_metrics)
        serializable_local_metrics = convert_numpy_to_list(local_metrics_output)
        serializable_idt_complexity = convert_numpy_to_list(idt_complexity_result) if idt_complexity_result else None
    except Exception as e:
         log.error(f"Error converting evaluation results for JSON ({output_suffix}): {e}", exc_info=True)
         serializable_global_metrics = {'error': f'Failed to convert global metrics: {e}'}
         serializable_local_metrics = {'error': f'Failed to convert local metrics: {e}'}
         serializable_idt_complexity = None

    evaluation_result_path: str = os.path.join(args_namespace.output, f"evaluation_result{output_suffix}.json")
    log.info(f"Preparing to save evaluation results to {evaluation_result_path}")
    evaluation_result: Dict[str, Any] = {
        "input_sequence_id": sequence_id,
        "input_sequence_description": sequence_description,
        "evaluated_sequence_length": len(sequence_to_evaluate),
        "original_cds_length_reference": original_cds_len,
        "evaluation_status": "Success" if evaluation_successful and 'error' not in serializable_global_metrics and 'error' not in serializable_local_metrics else "Failure",
        "metrics_excluded": metrics_to_skip_user_names if metrics_to_skip_user_names else "None",
        "global_metrics": serializable_global_metrics,
        "local_metrics": serializable_local_metrics,
    }
    
    # Add IDT complexity results if available
    if serializable_idt_complexity is not None:
        evaluation_result["idt_complexity"] = serializable_idt_complexity
    if output_suffix == "_mrna":
        evaluation_result["notes"] = "Evaluation for entire mRNA (5UTR-CDS-3UTR)"

    try:
        with open(evaluation_result_path, "w") as f:
            json.dump(evaluation_result, f, indent=4)
        log.info(f"Evaluation results successfully saved to {evaluation_result_path}")
    except Exception as e:
        log.error(f"Failed to save evaluation results to JSON ({output_suffix}): {e}", exc_info=True)

    checkpoints_path: str = os.path.join(args_namespace.output, f"checkpoints{output_suffix}.tsv")
    log.info(f"Saving evaluation summary to {checkpoints_path}")
    try:
        with open(checkpoints_path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            gm_dict = serializable_global_metrics
            if isinstance(gm_dict, dict) and gm_dict and 'error' not in gm_dict:
                 gm_keys = [k for k in gm_dict.keys() if '_error' not in k]
                 gm_values = [str(gm_dict[k]) if isinstance(gm_dict[k], (list, dict)) else gm_dict[k] for k in gm_keys]
            elif isinstance(gm_dict, dict) and 'error' in gm_dict:
                 gm_keys = ['error']
                 gm_values = [gm_dict['error']]
            else:
                 gm_keys = ['status']
                 gm_values = ['No valid global metrics']
            header: List[str] = ["description", "sequence_id", "sequence_preview"] + [f"metric:{k}" for k in gm_keys]
            writer.writerow(header)
            seq_preview = (sequence_to_evaluate[:75] + '...') if len(sequence_to_evaluate) > 75 else sequence_to_evaluate
            row: List[Any] = [sequence_description, sequence_id, seq_preview] + gm_values
            writer.writerow(row)
        log.info(f"✅ Evaluation summary successfully saved to TSV: {checkpoints_path}")
    except Exception as e:
        log.error(f"Failed to save evaluation summary to TSV ({output_suffix}): {e}", exc_info=True)
    log.info(f"--- Finished evaluation for {sequence_id} ({output_suffix.strip('_').upper()}) ---")


def main():
    parser = argparse.ArgumentParser(description="Vaxlab Report - Evaluate Only Script")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file (can contain CDS, 5UTR, 3UTR)")
    parser.add_argument("-o", "--output", required=True, help="Output directory for results")
    parser.add_argument("--preset", required=False, help="Preset JSON file path")
    parser.add_argument("--folding-engine", default="vienna", choices=["vienna", "linearfold"], help=argparse.SUPPRESS)
    parser.add_argument("--overwrite", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--seed", type=int, default=None, help=argparse.SUPPRESS)
    parser.add_argument("--codon-table", default="standard", help=argparse.SUPPRESS)
    parser.add_argument("--species", default="Homo sapiens", help=argparse.SUPPRESS)
    parser.add_argument("--token", type=str, default=None, help="IDT API access token for complexity score")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    log_file = os.path.join(args.output, "evaluate_only_log.txt")
    initialize_logging(log_file, quiet=False)
    log.info(f"Starting evaluate_only script. Vaxlab Report Version: {__version__}")
    log.info(f"Arguments: {args}")

    sequences: Dict[str, str] = {}
    seq_meta: Dict[str, Dict[str, str]] = {"CDS": {}, "5UTR": {}, "3UTR": {}}

    try:
        fasta_records = list(SeqIO.parse(args.input, "fasta"))
        if not fasta_records:
            log.critical(f"No sequences found in FASTA file: {args.input}")
            exit(1)
        found_cds = False
        for record in fasta_records:
            seq_upper_u = str(record.seq).upper().replace("T", "U")
            if record.id.upper() == "5UTR":
                sequences["5UTR"] = seq_upper_u
                seq_meta["5UTR"] = {"id": record.id, "description": record.description}
                log.info(f"Loaded 5UTR: {record.id}, Length: {len(sequences['5UTR'])}")
            elif record.id.upper() == "3UTR":
                sequences["3UTR"] = seq_upper_u
                seq_meta["3UTR"] = {"id": record.id, "description": record.description}
                log.info(f"Loaded 3UTR: {record.id}, Length: {len(sequences['3UTR'])}")
            elif not found_cds:
                sequences["CDS"] = seq_upper_u
                seq_meta["CDS"] = {"id": record.id, "description": record.description}
                found_cds = True
                log.info(f"Identified CDS: {record.id}, Length: {len(sequences['CDS'])}")
            else:
                log.warning(f"Multiple potential CDS sequences found. Using '{seq_meta['CDS']['id']}' as CDS. Sequence '{record.id}' will be ignored as CDS.")
        if "CDS" not in sequences:
            log.critical(f"CDS sequence not identified in {args.input}. A sequence not labeled '5UTR' or '3UTR' must be present.")
            exit(1)
    except FileNotFoundError:
        log.critical(f"Input FASTA file not found: {args.input}")
        exit(1)
    except Exception as e:
         log.critical(f"Error reading FASTA file {args.input}: {e}", exc_info=True)
         exit(1)

    cds_sequence_str = sequences["CDS"]
    cds_id = seq_meta["CDS"]["id"]
    cds_description = seq_meta["CDS"]["description"]
    original_cds_length = len(cds_sequence_str)
    
    # CDS 길이 유효성 검사 (MutantGenerator를 위해)
    if original_cds_length % 3 != 0:
        log.warning(f"Identified CDS sequence ('{cds_id}') length ({original_cds_length}) is not a multiple of 3. This might cause issues with MutantGenerator or CDS-specific metrics.")


    current_seed = args.seed if args.seed is not None else random.randint(0, 2**32 - 1)
    random_state = random.Random(current_seed)
    log.info(f"Using random seed: {current_seed}")

    preset_config = config.load_config()
    if args.preset:
        try:
            with open(args.preset) as f:
                 preset_data = load_preset(f.read())
                 preset_config.update(preset_data)
                 log.info(f"Loaded and applied preset from {args.preset}")
        except Exception as e:
             log.error(f"Error loading preset file {args.preset}: {e}. Using default or current settings.", exc_info=True)
    
    addon_paths = preset_config.get("addons", [])
    base_scoring_options = preset_config.get("fitness", {})
    execution_options_from_preset = preset_config.get("execution", {})
    if not base_scoring_options: log.warning("No 'fitness' section found in preset.")

    try:
        base_scoring_funcs = scoring.discover_scoring_functions(addon_paths)
        log.info(f"Discovered scoring functions: {list(base_scoring_funcs.keys())}")
        if not base_scoring_funcs: log.warning("No scoring functions were discovered.")
    except Exception as e:
        log.error(f"Error discovering scoring functions: {e}", exc_info=True)
        base_scoring_funcs = {}

    log.info("Initializing ExecutionOptions...")
    try:
        execopts_dict = {
            'n_iterations': 0, 'n_population': 1, 'n_survivors': 1,
            'output': args.output,
            'command_line': " ".join(shlex.quote(x) for x in os.sys.argv),
            'overwrite': args.overwrite, 'seed': current_seed, 
            'processes': execution_options_from_preset.get('processes', 1),
            'species': args.species, 'codon_table': args.codon_table,
            'quiet': False, 'seq_description': cds_description,
            'print_top_mutants': 0, 'folding_engine': args.folding_engine,
            'addons': addon_paths,
            # ... (rest of execopts_dict as before) ...
            'initial_mutation_rate': execution_options_from_preset.get('initial_mutation_rate', 0.1),
            'winddown_trigger': execution_options_from_preset.get('winddown_trigger', 15),
            'winddown_rate': execution_options_from_preset.get('winddown_rate', 0.9),
            'random_initialization': execution_options_from_preset.get('random_initialization', False),
            'conservative_start': execution_options_from_preset.get('conservative_start', None),
            'boost_loop_mutations': execution_options_from_preset.get('boost_loop_mutations', "1.5:15"),
            'full_scan_interval': execution_options_from_preset.get('full_scan_interval', 0),
            'protein': execution_options_from_preset.get('protein', False),
            'lineardesign_dir': execution_options_from_preset.get('lineardesign_dir', None),
            'lineardesign_lambda': execution_options_from_preset.get('lineardesign_lambda', None),
            'lineardesign_omit_start': execution_options_from_preset.get('lineardesign_omit_start', 5),
        }
        execopts = ExecutionOptions(**execopts_dict)
        log.info("ExecutionOptions initialized successfully.")
    except Exception as e:
         log.critical(f"Error during ExecutionOptions initialization: {e}", exc_info=True)
         exit(1)

    max_workers = execopts.processes if execopts.processes > 0 else os.cpu_count() or 1
    log.info(f"Using ThreadPoolExecutor with max_workers={max_workers}")
    with ThreadPoolExecutor(max_workers=max_workers) as executor:

        log.info("=== Starting CDS Only Evaluation ===")
        run_evaluation(
            sequence_to_evaluate=cds_sequence_str,
            original_cds_sequence_for_mutgen=cds_sequence_str, # CDS 평가이므로 동일
            sequence_id=cds_id,
            sequence_description=cds_description,
            output_suffix="_cds",
            execopts=execopts,
            base_scoring_funcs=base_scoring_funcs,
            base_scoring_options=base_scoring_options,
            metrics_to_skip_user_names=[],
            args_namespace=args,
            random_state_obj=random_state,
            original_cds_len=original_cds_length,
            executor=executor,
            idt_token=args.token
        )

        if "5UTR" in sequences and "3UTR" in sequences:
            log.info("=== Starting Entire mRNA Evaluation ===")
            utr5_seq = sequences["5UTR"]
            utr3_seq = sequences["3UTR"]
            entire_mrna_seq = utr5_seq + cds_sequence_str + utr3_seq
            
            mrna_id = f"{cds_id}_mRNA"
            mrna_description = f"{cds_description} (Entire mRNA: 5UTR + CDS + 3UTR)"
            metrics_to_skip_for_mrna = ["CAI", "Bicodon score", "start codon structure"]
            
            run_evaluation(
                sequence_to_evaluate=entire_mrna_seq, # 평가 대상은 전체 mRNA
                original_cds_sequence_for_mutgen=cds_sequence_str, # MutantGenerator는 원본 CDS 사용
                sequence_id=mrna_id,
                sequence_description=mrna_description,
                output_suffix="_mrna",
                execopts=execopts,
                base_scoring_funcs=base_scoring_funcs,
                base_scoring_options=base_scoring_options,
                metrics_to_skip_user_names=metrics_to_skip_for_mrna,
                args_namespace=args,
                random_state_obj=random_state,
                original_cds_len=original_cds_length, # SequenceEvaluator에는 원본 CDS 길이 전달
                executor=executor,
                idt_token=args.token
            )
        elif "5UTR" in sequences or "3UTR" in sequences:
            log.warning("One UTR found but not both (5UTR and 3UTR needed). Skipping entire mRNA evaluation.")
        else:
            log.info("No 5UTR and/or 3UTR sequences found. Skipping entire mRNA evaluation.")

    log.info("evaluate_only script finished.")

if __name__ == "__main__":
    main()