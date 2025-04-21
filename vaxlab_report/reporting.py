# reporting.py

import plotly.offline as pyo
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import jinja2
import os
import re
import json
import time
import csv
import argparse
from typing import Dict, Any, List, Optional, Mapping
from vaxlab_report.log import log # Import log

# ExecutionOptions 클래스 import 시도 및 예외 처리
try:
    from vaxlab_report.evolution_chamber import ExecutionOptions
except ImportError:
    log.warning("Could not import ExecutionOptions from vaxlab_report.evolution_chamber. Using placeholder.")
    class ExecutionOptions:
        """Placeholder for ExecutionOptions if import fails."""
        def __init__(self, **kwargs): self._options = kwargs
        def to_dict(self) -> Dict[str, Any]: return self._options

class TemplateFiltersMixin:
    """ Mixin class providing Jinja2 template filters. """
    @staticmethod
    def filter_localtime(timestamp: float) -> str:
        try: return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(float(timestamp)))
        except (ValueError, TypeError): return "Invalid timestamp"
    @staticmethod
    def filter_format_bool(value: bool, truevalue: str, falsevalue: str) -> str:
        return [falsevalue, truevalue][int(bool(value))]
    @staticmethod
    def filter_format_number(value: Any) -> str:
        if isinstance(value, (int, float)):
            if abs(value) < 1e-4 and value != 0: return f"{value:.2e}"
            if isinstance(value, float):
                 formatted = f"{value:.3f}".rstrip('0').rstrip('.')
                 return formatted if '.' in formatted else f"{value:.0f}"
            return str(value)
        return str(value)

    @classmethod
    def set_filters(cls, env: jinja2.Environment):
        for name in dir(cls):
            if name.startswith('filter_'):
                func = getattr(cls, name)
                if callable(func):
                    env.filters[name.replace('filter_', '', 1)] = func


class ReportGenerator:
    def __init__(
        self,
        status: Dict[str, Any],
        args: argparse.Namespace,
        metainfo: Dict[str, Any],
        scoring_options: Dict[str, Any],
        execution_options: ExecutionOptions,
        inputseq: Dict[str, Any],
        outputseq: Dict[str, Any],
        scoring_functions: Dict[str, Any],
        metric_labels: Optional[Mapping[str, str]] = None,
        metric_descriptions: Optional[Mapping[str, str]] = None,
        metric_ranges: Optional[Mapping[str, str]] = None
    ):
        self.status = status
        self.args = args
        self.metainfo = metainfo
        self.scoring_options = scoring_options
        self.execution_options = execution_options
        self.inputseq = inputseq
        self.outputseq = outputseq
        self.scoring_functions = scoring_functions
        self.metric_labels = metric_labels if metric_labels is not None else {}
        self.metric_descriptions = metric_descriptions if metric_descriptions is not None else {}
        self.metric_ranges = metric_ranges if metric_ranges is not None else {}

        self.templates: Dict[str, jinja2.Template] = {}
        self.load_templates()

    def load_templates(self):
        """Loads Jinja2 templates."""
        try:
            template_dir = os.path.join(os.path.dirname(__file__), 'report_template')
            if not os.path.isdir(template_dir):
                 import vaxlab_report
                 template_dir = os.path.join(os.path.dirname(vaxlab_report.__file__), 'report_template')

            if not os.path.isdir(template_dir):
                 raise FileNotFoundError(f"Template directory not found: {template_dir}")

            log.debug(f"Using template directory: {template_dir}")
            template_loader = jinja2.FileSystemLoader(searchpath=template_dir)
            env = jinja2.Environment(loader=template_loader, autoescape=True)
            TemplateFiltersMixin.set_filters(env)
            self.templates['report.html'] = env.get_template('report.html')
            log.debug("Jinja2 templates loaded successfully.")
        except Exception as e:
            log.error(f"Error loading Jinja2 templates: {e}"); raise e

    def prepare_template_data(self) -> Dict[str, Any]:
        """Prepares the dictionary to be passed to the template."""
        evaluations_data = self.status.get('evaluations', {})
        optimized_eval = evaluations_data.get('optimized', {})

        params: Dict[str, Any] = {
            'args': vars(self.args),
            'metainfo': self.metainfo,
            'scoring_options': self.scoring_options,
            'execution_options': self.execution_options.to_dict() if hasattr(self.execution_options, 'to_dict') else {},
            'inputseq': self.inputseq,
            'outputseq': self.outputseq,
            'evaluations': evaluations_data if isinstance(evaluations_data, dict) else {},
            'positional_plot_div': self.status.get('positional_plot_div', None),
            'metric_labels': self.metric_labels,
            'metric_descriptions': self.metric_descriptions,
            'metric_ranges': self.metric_ranges
        }

        if 'optimized' not in params['evaluations'] or not isinstance(params['evaluations']['optimized'], dict):
             params['evaluations']['optimized'] = {}
        if 'global_metrics' not in params['evaluations']['optimized'] or not isinstance(params['evaluations']['optimized'].get('global_metrics'), dict):
             params['evaluations']['optimized']['global_metrics'] = optimized_eval.get('global_metrics', {})
        if 'local_metrics' not in params['evaluations']['optimized']:
             params['evaluations']['optimized']['local_metrics'] = optimized_eval.get('local_metrics', {})

        return params

    def generate(self) -> str:
        """Generates the HTML report."""
        template = self.templates.get('report.html')
        if not template: raise ValueError("Report template 'report.html' not loaded.")
        try:
            params: Dict[str, Any] = self.prepare_template_data()
            log.debug("Prepared template parameters including specific metric metadata.")
            output: str = template.render(**params)
            log.debug("Template rendered successfully.")
            return output
        except Exception as e:
            log.error(f"Error during template rendering: {e}"); raise e

def save_checkpoints(results: List[Dict[str, Any]], output_file: str):
    """Saves checkpoint data to a TSV file."""
    if not results:
        log.warning("No results provided to save_checkpoints.")
        return

    try:
        with open(output_file, 'w', newline='', encoding='utf-8') as f:
            fieldnames = ["description", "sequence"]
            first_result = results[0]
            global_metrics = first_result.get("global_metrics")
            metric_keys = [] # Initialize metric_keys
            if isinstance(global_metrics, dict):
                 metric_keys = sorted(global_metrics.keys())
                 fieldnames += [f"metric:{k}" for k in metric_keys]
            else:
                 log.warning(f"Could not infer metric keys from first result: {first_result}")

            writer = csv.writer(f, delimiter='\t')
            writer.writerow(fieldnames)

            for res in results:
                 if not isinstance(res, dict): continue

                 row_data = { "description": res.get("description", ""), "sequence": res.get("sequence", "") }
                 global_metrics_data = res.get("global_metrics", {})
                 if isinstance(global_metrics_data, dict):
                      for k in metric_keys:
                           row_data[f"metric:{k}"] = global_metrics_data.get(k, "")

                 writer.writerow([row_data.get(field, "") for field in fieldnames])

        log.info(f"✅ Checkpoints saved to: {output_file}") # Use log.info
    except Exception as e:
        log.error(f"Error saving checkpoints to {output_file}: {e}") # Use log.error