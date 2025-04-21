# sequence_evaluator.py

import sys
import re
import pylru
from tqdm import tqdm
from concurrent import futures
from collections import Counter, defaultdict
import numpy as np # Keep numpy import as it's used in collect_scores
from .log import hbar_stars, log
from typing import Dict, List, Any, Optional, Tuple

class FoldEvaluator:
    """Handles RNA secondary structure folding using different engines."""
    def __init__(self, engine: str):
        self.engine = engine
        self.pat_find_loops = re.compile(r'\.{2,}')
        self.initialize()

    def initialize(self):
        """Initializes the folding engine."""
        log.debug(f"Initializing folding engine: {self.engine}")
        if self.engine == 'vienna':
            try:
                import RNA
                self._fold = RNA.fold
                log.info("ViennaRNA engine initialized.")
            except ImportError:
                log.error('ViennaRNA module is not available. Try "pip install ViennaRNA" to install.')
                raise
        elif self.engine == 'linearfold':
            try:
                import linearfold as lf
                def linearfold_wrapper(seq):
                     structure = lf.fold(seq)[0]
                     return structure, float('nan')
                self._fold = linearfold_wrapper
                log.info("LinearFold engine initialized (Note: MFE not calculated).")
            except ImportError:
                log.error('LinearFold module is not available. Try "pip install linearfold-unofficial" to install.')
                raise
        else:
            raise ValueError(f'Unsupported RNA folding engine: {self.engine}')

    def __call__(self, seq: str) -> Dict[str, Any]:
        """Performs folding and basic analysis."""
        try:
            folding, mfe = self._fold(seq)
        except Exception as e:
            log.error(f"RNA folding failed for sequence snippet '{seq[:30]}...': {e}")
            return {'folding': '.' * len(seq), 'mfe': float('inf'), 'stems': [], 'loops': {}, 'error': str(e)}

        stems = self.find_stems(folding)
        cleaned_folding = self.unfold_unstable_structure_simple(folding, stems)
        cleaned_stems = self.find_stems(cleaned_folding)
        cleaned_loops = dict(Counter(map(len, self.pat_find_loops.findall(cleaned_folding))))

        return {
            'folding': cleaned_folding,
            'mfe': mfe,
            'stems': cleaned_stems,
            'loops': cleaned_loops,
            'original_folding': folding,
        }

    @staticmethod
    def find_stems(structure: str) -> List[List[List[int]]]:
        """Finds stem structures (paired regions) in dot-bracket notation."""
        stack = []
        stemgroups = []
        for i, s in enumerate(structure):
            if s == '(':
                stack.append(i)
            elif s == ')':
                if not stack: continue
                peer = stack.pop()
                is_contiguous = False
                if stemgroups:
                    # Check for simple adjacency
                    if stemgroups and peer == stemgroups[-1][0][-1] - 1 and i == stemgroups[-1][1][-1] + 1:
                        is_contiguous = True

                if is_contiguous:
                    stemgroups[-1][0].append(peer)
                    stemgroups[-1][1].append(i)
                else:
                    stemgroups.append([[peer], [i]])
        stemgroups.sort(key=lambda x: x[0][0])
        return stemgroups

    @staticmethod
    def unfold_unstable_structure_simple(folding: str, stems: List[List[List[int]]]) -> str:
        """Unfolds structures based on a simple lone pair rule (length 1 stem)."""
        lonepairs = [p for p in stems if len(p[0]) == 1]
        if not lonepairs:
            return folding

        folding_list = list(folding)
        unfolded_indices = set()
        for p5_indices, p3_indices in lonepairs:
            if 0 <= p5_indices[0] < len(folding_list):
                 unfolded_indices.add(p5_indices[0])
            if 0 <= p3_indices[0] < len(folding_list):
                 unfolded_indices.add(p3_indices[0])

        for index in unfolded_indices:
            folding_list[index] = '.'

        return ''.join(folding_list)

class SequenceEvaluator:
    """Evaluates sequences using various scoring functions and folding."""
    folding_cache_size = 8192

    def __init__(self, scoring_funcs: Dict[str, type], scoreopts: Dict[str, Dict],
                 execopts: Any, mutantgen: Any, species: str,
                 length_cds: int, quiet: bool):
        self.scoring_funcs = scoring_funcs
        self.scoreopts = scoreopts
        self.execopts = execopts
        self.length_cds = length_cds
        self.mutantgen = mutantgen
        self.species = species
        self.quiet = quiet
        self.initialize()

    def initialize(self):
        """Initializes FoldEvaluator, caches, and scoring function instances."""
        self.foldeval = FoldEvaluator(self.execopts.folding_engine)
        self.folding_cache = pylru.lrucache(self.folding_cache_size)
        self.scorefuncs_nofolding = []
        self.scorefuncs_folding = []
        self.annotationfuncs = []
        self.penalty_metric_flags = {}
        additional_opts = {'_length_cds': self.length_cds}

        log.info("Initializing scoring functions...")
        if not self.scoring_funcs:
             log.warning("No scoring functions provided or discovered.")
             return

        for funcname, cls in self.scoring_funcs.items():
            opts = self.scoreopts.get(funcname, {})
            is_off = ('weight' in opts and opts['weight'] == 0) or ('off' in opts and opts['off'])
            use_annot_when_off = getattr(cls, 'use_annotation_on_zero_weight', False)

            if is_off and not use_annot_when_off:
                log.debug(f"  - Skipping {funcname} (off and no annotation needed)")
                continue

            current_opts = opts.copy()
            current_opts.update(additional_opts)
            for reqattr in getattr(cls, 'requires', []):
                if hasattr(self, reqattr):
                     current_opts['_' + reqattr] = getattr(self, reqattr)
                else:
                     log.warning(f"Dependency '{reqattr}' not found for scoring function '{funcname}'.")

            try:
                scorefunc_inst = cls(**current_opts)
                log.debug(f"  + Initialized {funcname}{' (annotation only)' if is_off else ''}")
                self.annotationfuncs.append(scorefunc_inst)

                if not is_off:
                    if getattr(cls, 'uses_folding', False):
                        self.scorefuncs_folding.append(scorefunc_inst)
                    else:
                        self.scorefuncs_nofolding.append(scorefunc_inst)
                    self.penalty_metric_flags.update(getattr(cls, 'penalty_metric_flags', {}))

            except Exception as e:
                log.error(f"  ! Failed to initialize scoring function '{funcname}': {e}", exc_info=True)
        log.info(f"Scoring function initialization complete. {len(self.annotationfuncs)} functions loaded for annotation/local metrics.")


    def evaluate(self, seqs: List[str], executor: futures.Executor) -> Tuple[Optional[List[float]], Optional[List[Dict]], Optional[List[Dict]], Optional[List[Dict]], Optional[List[Dict]]]:
        """
        Evaluates a list of sequences. Returns (total_scores, scores, metrics, foldings, local_metrics) or Nones on error.
        """
        with SequenceEvaluationSession(self, seqs, executor) as sess:
            try:
                 sess.evaluate()
            except Exception as e:
                 log.critical(f"Critical error during SequenceEvaluationSession.evaluate(): {e}", exc_info=True)
                 sess.errors.append(f"Session level error: {e}")

            if not sess.errors:
                total_scores = []
                for i, s_dict in enumerate(sess.scores):
                    score_sum = sum(v for k, v in s_dict.items() if isinstance(v, (int, float)) and '_error' not in k)
                    total_scores.append(score_sum)
                return total_scores, sess.scores, sess.metrics, sess.foldings, sess.local_metrics
            else:
                log.error(f"Evaluation session finished with errors: {sess.errors}")
                return None, None, None, None, None

    def get_folding(self, seq: str) -> Dict[str, Any]:
        """Gets folding results from cache or calculates using FoldEvaluator."""
        if seq in self.folding_cache:
            return self.folding_cache[seq]
        try:
            folding_result = self.foldeval(seq)
            if 'error' not in folding_result:
                 self.folding_cache[seq] = folding_result
            return folding_result
        except Exception as e:
            log.error(f"Error during folding calculation for sequence snippet '{seq[:30]}...': {e}", exc_info=True)
            return {'folding': '.' * len(seq), 'mfe': float('inf'), 'stems': [], 'loops': {}, 'error': str(e)}

    def prepare_evaluation_data(self, seq: str) -> Dict[str, Any]:
        """Calculates local metrics and annotations for a single sequence. (May be deprecated)"""
        log.warning("prepare_evaluation_data is called. The main evaluate method now returns local metrics. Consider using evaluate.")
        folding = self.get_folding(seq)
        if 'error' in folding:
             log.error(f"Cannot prepare evaluation data due to folding error for seq '{seq[:30]}...'")
             return {'local-metrics': {'error': 'Folding failed'}, 'annotations': {'error': 'Folding failed'}}

        localmet = {}
        annotations = {}

        for fun in self.annotationfuncs:
            func_name = getattr(fun, 'name', 'unknown_function')
            try:
                 if hasattr(fun, 'evaluate_local'):
                    uses_folding = getattr(fun, 'uses_folding', False)
                    if uses_folding:
                        result = fun.evaluate_local(seq, folding)
                    else:
                        result = fun.evaluate_local(seq)
                    if isinstance(result, dict):
                         localmet.update(result)
                    else:
                         log.warning(f"evaluate_local for {func_name} returned non-dict: {type(result)}")

                 if hasattr(fun, 'annotate_sequence'):
                    uses_folding = getattr(fun, 'uses_folding', False)
                    if uses_folding:
                        result = fun.annotate_sequence(seq, folding)
                    else:
                        result = fun.annotate_sequence(seq)
                    if isinstance(result, dict):
                         annotations.update(result)
                    else:
                         log.warning(f"annotate_sequence for {func_name} returned non-dict: {type(result)}")

            except Exception as e:
                 log.error(f"Error during prepare_evaluation_data for func {func_name}: {e}", exc_info=True)
                 localmet[f"{func_name}_error"] = f"prepare_evaluation_data: {e}"
                 annotations[f"{func_name}_error"] = f"prepare_evaluation_data: {e}"

        return {'local-metrics': localmet, 'annotations': annotations}


class SequenceEvaluationSession:
    """Manages the evaluation process for a batch of sequences."""
    def __init__(self, evaluator: SequenceEvaluator, seqs: list[str],
                 executor: futures.Executor):
        self.evaluator = evaluator
        self.seqs = seqs
        self.executor = executor

        self.scores: List[Dict[str, float]] = [defaultdict(float) for _ in range(len(seqs))]
        self.metrics: List[Dict[str, Any]] = [{} for _ in range(len(seqs))]
        self.local_metrics: List[Dict[str, Any]] = [{} for _ in range(len(seqs))] # Local metrics
        self.foldings: List[Optional[Dict[str, Any]]] = [None] * len(seqs)
        self.errors: List[str] = []

        self.folding_cache = evaluator.folding_cache
        self.foldings_submitted_indices = set()
        self.foldings_completed_count = 0

        self.pbar = None
        self.quiet = evaluator.quiet

        self.num_folding_tasks = len(seqs)
        self.num_scoring_tasks = len(evaluator.scorefuncs_nofolding) + len(evaluator.scorefuncs_folding)
        self.num_local_metric_tasks_approx = len(seqs) # Simplified count for pbar
        self.total_tasks_approx = self.num_folding_tasks + self.num_scoring_tasks + self.num_local_metric_tasks_approx

    def __enter__(self):
        log.info(f'Starting evaluation session for {len(self.seqs)} sequences...')
        self.pbar = tqdm(total=self.total_tasks_approx, disable=self.quiet,
                         file=sys.stderr, unit='task', desc='Initializing...', smoothing=0.1)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if self.pbar is not None:
            final_n = self.pbar.n
            if not self.errors and self.pbar.n < self.pbar.total:
                 final_n = self.pbar.total
            self.pbar.n = final_n

            final_desc = "Finished"
            if self.errors:
                 final_desc = f"Finished with errors ({len(self.errors)})"
                 self.pbar.colour = 'red'
            self.pbar.set_description(final_desc, refresh=True)
            self.pbar.close()

        if self.errors:
             log.error(f"Evaluation session finished with {len(self.errors)} errors.")
             for i, err in enumerate(self.errors[:5]):
                  log.error(f"  Error {i+1}: {err}")
             if len(self.errors) > 5:
                  log.error(f"  ... and {len(self.errors) - 5} more errors (see logs for details).")
        else:
             log.info('Evaluation session finished successfully.')

    def evaluate(self) -> None:
        """Orchestrates folding, scoring, and local metric calculation using futures."""
        if self.pbar: self.pbar.set_description('Submitting tasks')
        jobs = set()

        # --- 1. Submit Folding Tasks ---
        for i, seq in enumerate(self.seqs):
            if seq in self.folding_cache:
                self.foldings[i] = self.folding_cache[seq]
                self.foldings_completed_count += 1
                if self.pbar: self.pbar.update()
            else:
                try:
                     future = self.executor.submit(self.evaluator.get_folding, seq)
                     future._seqidx = i
                     future._type = 'folding'
                     jobs.add(future)
                     self.foldings_submitted_indices.add(i)
                except Exception as e:
                     self.handle_exception(e, context=f"submitting folding for seq {i}")

        # --- 2. Submit Non-Folding Scoring Tasks ---
        for scorefunc in self.evaluator.scorefuncs_nofolding:
            try:
                 future = self.executor.submit(scorefunc.score, self.seqs)
                 future._scorefunc_name = scorefunc.name
                 future._type = 'scoring'
                 jobs.add(future)
            except Exception as e:
                 self.handle_exception(e, context=f"submitting {scorefunc.name}")

        # --- 3. Process Completed Tasks Until All Foldings Are Done ---
        if self.pbar: self.pbar.set_description('Folding/Scoring...')
        processed_in_loop3 = 0
        while self.foldings_completed_count < self.num_folding_tasks and jobs and not self.errors:
            try:
                done, jobs = futures.wait(jobs, timeout=1.0, return_when=futures.FIRST_COMPLETED)
                if not done:
                     pending_folding = self.foldings_submitted_indices - set(i for i, f in enumerate(self.foldings) if f is not None)
                     if not pending_folding and self.foldings_completed_count < self.num_folding_tasks:
                          log.warning("Timeout waiting for folding, but no pending folding jobs detected. Check executor state.")
                     continue
            except Exception as e:
                 self.handle_exception(e, context="futures.wait (folding phase)")
                 break

            for future in done:
                processed_in_loop3 += 1
                if future._type == 'folding':
                    self.collect_folding(future)
                elif future._type == 'scoring':
                    self.collect_scores(future)
                else:
                     log.warning(f"Unknown future type encountered: {future._type}")
        log.debug(f"Folding phase loop processed {processed_in_loop3} futures.")

        # --- 4. Submit Folding-Dependent Scoring Tasks ---
        if not self.errors and self.foldings_completed_count == self.num_folding_tasks:
            log.debug("All folding tasks completed. Submitting folding-dependent scoring...")
            if any(isinstance(f, dict) and 'error' in f for f in self.foldings if f is not None):
                  log.warning("Cannot proceed with folding-dependent scoring due to errors in folding results.")
                  if "Folding errors detected" not in self.errors: self.errors.append("Folding errors detected")
            else:
                  valid_foldings = [f for f in self.foldings if isinstance(f, dict) and 'error' not in f]
                  if len(valid_foldings) == len(self.seqs):
                       for scorefunc in self.evaluator.scorefuncs_folding:
                            try:
                                 future = self.executor.submit(scorefunc.score, self.seqs, self.foldings)
                                 future._scorefunc_name = scorefunc.name
                                 future._type = 'scoring'
                                 jobs.add(future)
                            except Exception as e:
                                 self.handle_exception(e, context=f"submitting {scorefunc.name}")
                  else:
                       log.error(f"Mismatch in valid folding results ({len(valid_foldings)}) vs sequences ({len(self.seqs)}). Skipping folding-dependent scoring.")
                       if "Missing valid folding results" not in self.errors: self.errors.append("Missing valid folding results")
        elif not self.errors:
             log.error(f"Folding not fully completed ({self.foldings_completed_count}/{self.num_folding_tasks}). Skipping folding-dependent scoring.")
             if "Incomplete folding" not in self.errors: self.errors.append("Incomplete folding")


        # --- 5. Process Remaining Scoring Results ---
        if self.pbar: self.pbar.set_description('Scoring...')
        processed_in_loop5 = 0
        while jobs and not self.errors:
            try:
                done, jobs = futures.wait(jobs, timeout=1.0, return_when=futures.FIRST_COMPLETED)
                if not done: continue
            except Exception as e:
                 self.handle_exception(e, context="futures.wait (scoring phase)")
                 break

            for future in done:
                processed_in_loop5 += 1
                if future._type == 'folding':
                    self.collect_folding(future)
                    log.warning("Collected folding result in final processing loop.")
                elif future._type == 'scoring':
                    self.collect_scores(future)
                else:
                     log.warning(f"Unknown future type encountered: {future._type}")
        log.debug(f"Final processing loop processed {processed_in_loop5} futures.")


        # --- 6. Calculate Local Metrics (Serially for simplicity) ---
        if not self.errors:
            if self.pbar: self.pbar.set_description('Local Metrics...')
            log.info('Calculating local metrics...')
            if self.pbar: self.pbar.n = self.num_folding_tasks + self.num_scoring_tasks # Reset progress
            if self.pbar: self.pbar.total = self.num_folding_tasks + self.num_scoring_tasks + self.num_local_metric_tasks_approx # Adjust total

            for i, seq in enumerate(self.seqs):
                folding_result = self.foldings[i]

                if isinstance(folding_result, dict) and 'error' in folding_result:
                    log.warning(f"Skipping local metrics for sequence {i} due to previous folding error: {folding_result['error']}")
                    self.local_metrics[i]['error'] = f"Folding error: {folding_result['error']}"
                    if self.pbar: self.pbar.update(1)
                    continue
                if folding_result is None:
                     log.warning(f"Skipping local metrics for sequence {i} due to missing folding result (None).")
                     self.local_metrics[i]['error'] = 'Missing folding result'
                     if self.pbar: self.pbar.update(1)
                     continue

                current_local_metrics = {}
                functions_with_evaluate_local = [f for f in self.evaluator.annotationfuncs if hasattr(f, 'evaluate_local')]

                for fun in functions_with_evaluate_local:
                    func_name = getattr(fun, 'name', 'unknown_function')
                    try:
                        uses_folding = getattr(fun, 'uses_folding', False)
                        if uses_folding:
                            if isinstance(folding_result, dict) and 'folding' in folding_result:
                                 result = fun.evaluate_local(seq, folding_result)
                            else:
                                 raise ValueError(f"Valid folding dictionary required for {func_name} but not available")
                        else:
                            result = fun.evaluate_local(seq)

                        if isinstance(result, dict):
                             current_local_metrics.update(result)
                        else:
                             log.warning(f"evaluate_local for {func_name} returned non-dict: {type(result)}")

                    except Exception as exc:
                        log.error(f"Error during evaluate_local for {func_name} on sequence {i}: {exc}", exc_info=False)
                        current_local_metrics[f"{func_name}_error"] = str(exc)

                self.local_metrics[i] = current_local_metrics
                if self.pbar: self.pbar.update(1) # Update once per sequence

            log.info('Finished calculating local metrics.')
        else:
             log.warning("Skipping local metrics calculation due to previous errors in the session.")


    def collect_scores(self, future):
        """Collects results from score() method calls and updates scores/metrics."""
        func_name = getattr(future, '_scorefunc_name', 'unknown_function')
        try:
            ret = future.result()
            if ret is None:
                 log.warning(f"Scoring function {func_name} returned None.")
                 self.errors.append(f"{func_name} returned None")
                 return
            if isinstance(ret, tuple) and len(ret) == 2 and isinstance(ret[0], dict) and isinstance(ret[1], dict):
                 scoreupdates, metricupdates = ret
            else:
                 log.warning(f"Scoring function {func_name} returned unexpected data type: {type(ret)}. Expected (dict, dict).")
                 self.errors.append(f"{func_name} returned wrong type")
                 return
        except Exception as exc:
            self.handle_exception(exc, context=f"collecting results from {func_name}")
            return

        num_seqs = len(self.scores)

        # Update scores
        for k, updates in scoreupdates.items():
            if isinstance(updates, list) and len(updates) == num_seqs:
                for idx, s in enumerate(self.scores):
                    if isinstance(updates[idx], (int, float)):
                        s[k] = updates[idx]
                    else:
                         log.warning(f"Non-numeric score value '{updates[idx]}' for '{k}' from {func_name} at index {idx}. Skipping.")
            elif num_seqs == 1 and not isinstance(updates, list):
                value_to_update = None
                if isinstance(updates, (int, float)):
                    value_to_update = updates
                elif isinstance(updates, np.ndarray) and updates.size == 1:
                    value_to_update = updates.item()
                    if not isinstance(value_to_update, (int, float)):
                         log.warning(f"Extracted non-numeric score value '{value_to_update}' from NumPy array for '{k}' from {func_name}. Skipping.")
                         value_to_update = None
                else:
                    log.warning(f"Score update data for '{k}' from {func_name} for single sequence has unexpected type/size: {type(updates)}. Skipping.")

                if value_to_update is not None:
                     self.scores[0][k] = value_to_update
            else:
                 log.warning(f"Score update data for '{k}' from {func_name} has incorrect format or length (Expected list of len {num_seqs}, got {type(updates)}).")

        # Update metrics (Global metrics)
        for k, updates in metricupdates.items():
            if isinstance(updates, list) and len(updates) == num_seqs:
                for idx, m in enumerate(self.metrics):
                    m[k] = updates[idx]
            elif num_seqs == 1 and not isinstance(updates, list):
                 value_to_update = None
                 if isinstance(updates, np.ndarray) and updates.size == 1:
                      value_to_update = updates.item()
                 elif not isinstance(updates, list):
                      value_to_update = updates
                 else:
                      log.warning(f"Metric update data for '{k}' from {func_name} for single sequence has unexpected type/size: {type(updates)}. Skipping.")

                 if value_to_update is not None:
                      self.metrics[0][k] = value_to_update

            else:
                 log.warning(f"Metric update data for '{k}' from {func_name} has incorrect format or length (Expected list of len {num_seqs}, got {type(updates)}).")


    def collect_folding(self, future):
        """Collects results from folding tasks, updates cache and state."""
        seq_idx = getattr(future, '_seqidx', -1)
        if seq_idx < 0 or seq_idx >= len(self.seqs):
             log.error(f"Invalid sequence index ({seq_idx}) associated with folding future.")
             if "Invalid sequence index in folding task" not in self.errors: self.errors.append("Invalid sequence index in folding task")
             self.foldings_completed_count += 1
             if self.pbar: self.pbar.update()
             return

        try:
            folding_result = future.result()

            if folding_result is None:
                log.warning(f"Folding task for sequence {seq_idx} returned None unexpectedly.")
                if f'Folding task returned None (Seq {seq_idx})' not in self.errors: self.errors.append(f'Folding task returned None (Seq {seq_idx})')
                self.foldings[seq_idx] = {'error': 'Folding returned None'}
            elif isinstance(folding_result, dict) and 'error' in folding_result:
                 log.error(f"Folding failed for sequence {seq_idx}: {folding_result['error']}")
                 self.foldings[seq_idx] = folding_result
                 if f"Folding error seq {seq_idx}" not in self.errors: self.errors.append(f"Folding error seq {seq_idx}")
            elif isinstance(folding_result, dict) and 'folding' in folding_result:
                 self.foldings[seq_idx] = folding_result
                 self.folding_cache[self.seqs[seq_idx]] = folding_result
            else:
                 log.error(f"Unexpected folding result type for sequence {seq_idx}: {type(folding_result)}")
                 self.foldings[seq_idx] = {'error': f'Unexpected folding result type: {type(folding_result)}'}
                 if f"Unexpected folding result type seq {seq_idx}" not in self.errors: self.errors.append(f"Unexpected folding result type seq {seq_idx}")

        except Exception as exc:
            self.handle_exception(exc, context=f"folding sequence {seq_idx}")
            self.foldings[seq_idx] = {'error': f'Task execution failed: {exc}'}

        self.foldings_completed_count += 1


    def handle_exception(self, exc: Exception, context: str = "general evaluation"):
        """Handles exceptions during session execution, logs details, and adds to error list."""
        import traceback
        import io
        errormsg = io.StringIO()
        traceback.print_exc(file=errormsg)
        exc_repr = repr(exc)
        msg_list = [hbar_stars, f'ERROR occurred during {context}: {exc_repr}', f'Traceback:\n{errormsg.getvalue()}', hbar_stars]
        log.error('\n'.join(msg_list))
        error_summary = f"{context}: {exc_repr}"
        if error_summary not in self.errors: self.errors.append(error_summary)