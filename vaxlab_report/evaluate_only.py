# evaluate_only.py

import argparse
import os
import csv
import random
import shlex
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
import subprocess
import json
import numpy as np
from typing import Dict, List, Tuple, Any, Optional

from vaxlab_report import config, scoring, __version__
from vaxlab_report.sequence_evaluator import SequenceEvaluator
try:
    from vaxlab_report.evolution_chamber import ExecutionOptions
except ImportError as e:
     # Keep critical error handling, but remove debug prints
     print(f"FATAL ERROR: Failed to import ExecutionOptions from vaxlab_report.evolution_chamber: {e}")
     print("        Please ensure this module and class exist and the Python path is correct.")
     exit(1)
from vaxlab_report.mutant_generator import MutantGenerator
from vaxlab_report.presets import load_preset
from vaxlab_report.log import initialize_logging, log


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

def main():
    # --- 인수 파싱 ---
    parser = argparse.ArgumentParser(description="Vaxlab Report - Evaluate Only Script")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output directory for results")
    parser.add_argument("--preset", required=False, help="Preset JSON file path")
    parser.add_argument("--folding-engine", default="vienna", choices=["vienna", "linearfold"], help="RNA folding engine")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")
    parser.add_argument("--seed", type=int, default=None, help="Random seed (optional)")
    parser.add_argument("--codon-table", default="standard", help="Codon table name")
    parser.add_argument("--species", default="Homo sapiens", help="Species for calculations")
    args = parser.parse_args()

    # --- 초기 설정 ---
    os.makedirs(args.output, exist_ok=True)
    log_file = os.path.join(args.output, "evaluate_only_log.txt")
    initialize_logging(log_file, quiet=False)
    log.info(f"Starting evaluate_only script. Vaxlab Report Version: {__version__}")
    log.info(f"Arguments: {args}")

    # 시퀀스 로드
    try:
        input_seq = next(SeqIO.parse(args.input, "fasta"))
        cds_sequence = str(input_seq.seq).upper().replace("T", "U")
        log.info(f"Loaded sequence '{input_seq.id}' from {args.input}, Length: {len(cds_sequence)}")
    except FileNotFoundError:
        log.critical(f"Input FASTA file not found: {args.input}")
        exit(1)
    except StopIteration:
         log.critical(f"No sequences found in FASTA file: {args.input}")
         exit(1)
    except Exception as e:
         log.critical(f"Error reading FASTA file {args.input}: {e}", exc_info=True)
         exit(1)

    # 랜덤 시드 설정
    current_seed = args.seed if args.seed is not None else random.randint(0, 2**32 - 1)
    random_state = random.Random(current_seed)
    log.info(f"Using random seed: {current_seed}")

    # --- 설정 로드 (Preset) ---
    preset = config.load_config()
    if args.preset:
        try:
            with open(args.preset) as f:
                 preset_data = load_preset(f.read())
                 preset.update(preset_data)
                 log.info(f"Loaded and applied preset from {args.preset}")
        except FileNotFoundError:
             log.error(f"Preset file not found: {args.preset}. Using default settings.")
        except json.JSONDecodeError as e:
             log.error(f"Error decoding JSON from preset file {args.preset}: {e}. Using default settings.")
        except Exception as e:
             log.error(f"Error loading preset file {args.preset}: {e}. Using default settings.", exc_info=True)

    addon_paths = preset.get("addons", [])
    scoring_options = preset.get("fitness", {})
    execution_options_from_preset = preset.get("execution", {})
    if not scoring_options: log.warning("No 'fitness' section found in preset.")
    if not execution_options_from_preset: log.warning("No 'execution' section found in preset.")


    # --- 스코어링 함수 로드 ---
    try:
        scoring_funcs = scoring.discover_scoring_functions(addon_paths)
        log.info(f"Discovered scoring functions: {list(scoring_funcs.keys())}")
        if not scoring_funcs: log.warning("No scoring functions were discovered.")
    except Exception as e:
        log.error(f"Error discovering scoring functions: {e}", exc_info=True)
        scoring_funcs = {}


    # --- ExecutionOptions 인스턴스 생성 ---
    log.info("Initializing ExecutionOptions...")
    try:
        execopts_dict = {
            'n_iterations': 0, 'n_population': 1, 'n_survivors': 1,
            'output': args.output,
            'command_line': " ".join(shlex.quote(x) for x in os.sys.argv),
            'overwrite': args.overwrite, 'seed': current_seed, 'processes': 1,
            'species': args.species, 'codon_table': args.codon_table,
            'quiet': False, 'seq_description': input_seq.description,
            'print_top_mutants': 0, 'folding_engine': args.folding_engine,
            'addons': addon_paths,
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
    except TypeError as e:
         log.critical(f"TypeError during ExecutionOptions initialization: {e}", exc_info=True)
         log.critical("        Check 'vaxlab_report.evolution_chamber.ExecutionOptions' definition and import.")
         exit(1)
    except Exception as e:
         log.critical(f"Unexpected error during ExecutionOptions initialization: {e}", exc_info=True)
         exit(1)


    # --- SequenceEvaluator 인스턴스 생성 ---
    log.info("Initializing SequenceEvaluator...")
    mutantgen = MutantGenerator(cds_sequence, random_state)
    try:
        evaluator = SequenceEvaluator(
            scoring_funcs, scoring_options, execopts, mutantgen,
            args.species, len(cds_sequence), quiet=False
        )
        log.info("SequenceEvaluator initialized successfully.")
    except Exception as e:
         log.critical(f"Failed to initialize SequenceEvaluator: {e}", exc_info=True)
         exit(1)


    # --- 평가 실행 ---
    log.info("Starting sequence evaluation...")
    global_metrics: Dict[str, Any] = {}
    local_metrics_output: Dict[str, Any] = {}
    evaluation_successful = False

    try:
        max_workers = execopts.processes
        log.info(f"Using ThreadPoolExecutor with max_workers={max_workers}")
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = evaluator.evaluate([cds_sequence], executor)

        if results is not None and len(results) == 5:
            total_scores, scores_list, metrics_list, foldings_list, local_metrics_list = results

            if metrics_list: global_metrics = metrics_list[0]
            else: log.warning("SequenceEvaluator did not return global metrics.")

            if local_metrics_list: local_metrics_output = local_metrics_list[0]
            else: log.warning("SequenceEvaluator did not return local metrics.")

            evaluation_successful = True
            log.info("Evaluation process completed.")
            log.debug(f"Received Global Metrics Keys: {list(global_metrics.keys())}")
            log.debug(f"Received Local Metrics Keys: {list(local_metrics_output.keys())}")

        else:
            log.error("evaluator.evaluate() returned None or unexpected data. Evaluation failed.")
            global_metrics = {'error': 'Evaluation failed or returned None'}
            local_metrics_output = {'error': 'Evaluation failed or returned None'}

    except Exception as e:
        log.critical(f"An critical exception occurred during evaluator.evaluate(): {e}", exc_info=True)
        global_metrics = {'error': f'Evaluation exception: {e}'}
        local_metrics_output = {'error': f'Evaluation exception: {e}'}


    # --- NumPy 타입을 JSON 직렬화 가능하게 변환 ---
    log.debug("Converting results for JSON serialization...")
    try:
        serializable_global_metrics = convert_numpy_to_list(global_metrics)
        serializable_local_metrics = convert_numpy_to_list(local_metrics_output)
        log.debug("Conversion complete.")
    except Exception as e:
         log.error(f"Error converting evaluation results for JSON: {e}", exc_info=True)
         serializable_global_metrics = {'error': f'Failed to convert global metrics: {e}'}
         serializable_local_metrics = {'error': f'Failed to convert local metrics: {e}'}


    # --- 평가 결과 저장 (JSON 파일) ---
    evaluation_result_path: str = os.path.join(args.output, "evaluation_result.json")
    log.info(f"Preparing to save evaluation results to {evaluation_result_path}")
    evaluation_result: Dict[str, Any] = {
        "input_sequence_id": input_seq.id,
        "input_sequence_description": input_seq.description,
        "cds_length": len(cds_sequence),
        "evaluation_status": "Success" if evaluation_successful and 'error' not in serializable_global_metrics and 'error' not in serializable_local_metrics else "Failure",
        "global_metrics": serializable_global_metrics,
        "local_metrics": serializable_local_metrics,
    }

    try:
        with open(evaluation_result_path, "w") as f:
            json.dump(evaluation_result, f, indent=4)
        log.info(f"Evaluation results successfully saved to {evaluation_result_path}")
    except TypeError as e:
         log.error(f"TypeError saving results to JSON (Non-serializable data might still be present): {e}", exc_info=True)
    except Exception as e:
        log.error(f"Failed to save evaluation results to JSON: {e}", exc_info=True)


    # --- TSV 요약 출력 (Global metrics 기반) ---
    checkpoints_path: str = os.path.join(args.output, "checkpoints.tsv")
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

            header: List[str] = ["description", "sequence"] + [f"metric:{k}" for k in gm_keys]
            writer.writerow(header)
            row: List[Any] = [input_seq.description, cds_sequence] + gm_values
            writer.writerow(row)
        log.info(f"✅ Evaluation summary successfully saved to TSV: {checkpoints_path}")
    except Exception as e:
        log.error(f"Failed to save evaluation summary to TSV: {e}", exc_info=True)

    log.info("evaluate_only script finished.")


if __name__ == "__main__":
    main()