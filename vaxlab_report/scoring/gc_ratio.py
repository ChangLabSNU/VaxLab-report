# gc_ratio.py

from . import ScoringFunction
import numpy as np

def gc_content_sliding_window(seq, winsize, stride):
    chars = np.frombuffer(seq.encode(), dtype=np.uint8)
    isgc = ((chars == ord('G')) + (chars == ord('C')))
    gc = []
    for i in range(0, len(chars) - winsize + 1, stride):
        gc.append(np.mean(isgc[i:i+winsize]))
    return np.array(gc)

def compute_gc_ratio(seq, winsize, stride):
    """Computes the GC ratio instead of a penalty."""
    gc_values = gc_content_sliding_window(seq, winsize, stride)
    return np.mean(gc_values)  # Average GC ratio over the windows

class GCRatioFitness(ScoringFunction):

    name = 'gc'
    description = 'GC Ratio'
    priority = 50

    use_annotation_on_zero_weight = True

    arguments = [
        ('weight', dict(
            type=float, default=3.0, metavar='WEIGHT',
            help='scoring weight for GC ratio (default: 3.0)')),
        ('window-size', dict(
            type=int, default=50, metavar='SIZE',
            help='size of window for GC content calculation (default: 50)')),
        ('stride', dict(
            type=int, default=5, metavar='STRIDE',
            help='size of stride for GC content calculation (default: 5)')),
    ]

    def __init__(self, weight, window_size, stride, _length_cds):
        num_windows = (_length_cds - window_size) // stride + 1
        num_windows = max(num_windows, 1)

        self.weight = weight / num_windows
        self.window_size = window_size
        self.stride = stride

    def score(self, seqs):
        gc_ratios = [compute_gc_ratio(seq, self.window_size, self.stride)
                        for seq in seqs]
        # The score is now directly based on the ratio, no penalty calculation
        scores = [s * self.weight for s in gc_ratios]  
        return {'gc': scores}, {'gc': gc_ratios}  # Changed 'gc_penalty' to 'gc'

    def evaluate_local(self, seq):
        gc = gc_content_sliding_window(seq, self.window_size, self.stride)
        centers = (
            np.arange(0, len(seq) - self.window_size + 1, self.stride) +
            self.window_size / 2)
        return {'gc': (centers, gc)} # Keep 'gc' here