"""Process-local thread limits; never modify the notebook's environment."""
from threadpoolctl import threadpool_limits
import os

_thread_limits = None


def initialize_worker():
    """Avoid nested BLAS/OpenMP parallelism inside a spawned worker."""
    global _thread_limits
    # Only worker processes run this initializer. Cover native libraries that
    # are imported later while unpickling tasks, without touching the parent.
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[name] = '1'
    _thread_limits = threadpool_limits(limits=1)


def initialize_esp_worker():
    initialize_worker()
    import psi4

    psi4.set_num_threads(1)
    # Include libraries loaded by Psi4 as well as NumPy/SciPy.
    initialize_worker()
