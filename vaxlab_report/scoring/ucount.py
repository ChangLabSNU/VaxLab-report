#   
# VaxPress
#
# Copyright 2023 Seoul National University
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# “Software”), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
# NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
# THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#

from . import ScoringFunction
from collections import Counter

class UCountFitness(ScoringFunction):

    name = 'ucount'
    description = 'Uracil Content'
    priority = 10

    use_annotation_on_zero_weight = True

    arguments = [
        ('weight', dict(
            type=float, default=3.0, metavar='WEIGHT',
            help='scoring weight for uracil content (default: 3.0)')),
        ('target', dict(
            type=float, default=0.25, metavar='TARGET',
            help='target ratio of uracil content (default: 0.25)')),
    ]

    def __init__(self, weight, target, _length_cds):
        self.weight = weight
        self.target = target

    def _calculate_percentages(self, seq):
        counts = Counter(seq.upper())
        total_count = len(seq)
        a_percent = (counts.get('A', 0) / total_count) 
        c_percent = (counts.get('C', 0) / total_count) 
        g_percent = (counts.get('G', 0) / total_count) 
        u_percent = (counts.get('U', 0) / total_count) 
        return f"{u_percent:.2f}"

    def score(self, seqs):
        scores = []
        percentages = []
        for seq in seqs:
            percentage_str = self._calculate_percentages(seq)
            percentages.append(percentage_str)
            u_percent = float(percentage_str)  # Extract U percentage
            score = -abs(u_percent - self.target * 100)
            scores.append(score * self.weight)
        return {'ucount': scores}, {'ucount': percentages}

    def evaluate_local(self, seq):
        percentage_str = self._calculate_percentages(seq)
        return {'ucount': percentage_str}