"""Random-number buffers and checkpoint state must preserve stochastic dynamics."""

import numpy as np
import pytest

from ipi.inputs.prng import InputRandom
from ipi.utils.prng import Random


@pytest.fixture
def random_factory():
    """Close the worker pools created by each test."""
    generators = []

    def make(**kwargs):
        generator = Random(**kwargs)
        generators.append(generator)
        return generator

    yield make
    for generator in generators:
        if hasattr(generator, "executor"):
            generator.executor.shutdown(wait=True)


@pytest.mark.parametrize(
    "n_threads,shape,contiguous",
    [
        (1, (3,), True),
        (1, (67, 3), True),
        (2, (33, 3), True),
        (4, (33, 3), True),
        (2, (67, 3), True),
        (4, (401, 3), True),
        (4, (401, 3), False),
    ],
)
def test_gaussian_fill_updates_the_supplied_array(
    random_factory, n_threads, shape, contiguous
):
    """Both serial and parallel fills must overwrite every output element."""
    generator = random_factory(seed=109, n_threads=n_threads)
    output = np.full(shape, np.nan) if contiguous else np.full(shape[::-1], np.nan).T
    generator.gfill(output)
    assert np.isfinite(output).all()
    assert np.any(output != output.flat[0])


@pytest.mark.parametrize("n_threads", [1, 2, 4])
def test_checkpoint_restores_all_random_streams(random_factory, n_threads):
    """Restarting preserves the next draws, including each worker's state."""
    original = random_factory(seed=312, n_threads=n_threads)
    original.gvec((401, 3))
    checkpoint = InputRandom()
    checkpoint.store(original)
    restored = checkpoint.fetch()
    try:
        # Check every stream explicitly: a small draw would use only stream zero.
        for old, new in zip(original.rng, restored.rng):
            np.testing.assert_array_equal(old.normal(size=11), new.normal(size=11))
        assert restored.n_threads == n_threads
        for draw in (
            lambda rng: rng.u,
            lambda rng: rng.g,
            lambda rng: rng.gamma(2.5, size=7),
            lambda rng: rng.poisson(3.0, size=7),
            lambda rng: rng.uniform(-2, 3, size=7),
            lambda rng: rng.integers(2, 9, size=7, endpoint=True),
        ):
            np.testing.assert_array_equal(draw(original), draw(restored))
        left = np.arange(12).reshape(3, 4)
        right = left.copy()
        original.shuffle(left, axis=1)
        restored.shuffle(right, axis=1)
        np.testing.assert_array_equal(left, right)
        np.testing.assert_array_equal(
            np.sort(left, axis=1), np.arange(12).reshape(3, 4)
        )
    finally:
        if hasattr(restored, "executor"):
            restored.executor.shutdown(wait=True)
