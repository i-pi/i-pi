"""Unit tests for batched MACE preprocessing."""

import threading
from concurrent.futures import ThreadPoolExecutor

import pytest

pytest.importorskip("mace")
torch = pytest.importorskip("torch")

from ipi.pes._mace import BatchedMACE


def test_parallel_graph_construction_preserves_structure_order():
    """Graph workers run concurrently without reordering their results."""
    calculator = object.__new__(BatchedMACE)
    calculator._graph_executor = ThreadPoolExecutor(max_workers=2)
    barrier = threading.Barrier(2)
    thread_ids = set()
    thread_ids_lock = threading.Lock()

    def build_graph(config):
        with thread_ids_lock:
            thread_ids.add(threading.get_ident())
        barrier.wait(timeout=5)
        return config * 2

    calculator._config_to_atomic_data = build_graph
    try:
        result = calculator._configs_to_dataset([0, 1, 2, 3])
    finally:
        calculator._graph_executor.shutdown()

    assert result == [0, 2, 4, 6]
    assert len(thread_ids) == 2


def test_batch2natoms_uses_batch_pointers():
    """Atom counts are obtained without copying the positions tensor."""
    batch = {"ptr": torch.tensor([0, 2, 7, 10])}

    assert BatchedMACE.batch2natoms(batch) == [2, 5, 3]
