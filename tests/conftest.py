"""pytest setup: BLAS pinned before numpy loads (E31), and the regression tolerance of this platform."""
from ttmt.runtime import pin_blas

pin_blas()

import os  # noqa: E402

import pytest  # noqa: E402

from ttmt import config  # noqa: E402


def pytest_addoption(parser):
    parser.addoption("--match", choices=("auto", "exact", "tol"), default=os.environ.get("TTMT_MATCH", "auto"),
                     help="regression comparison: exact (Linux CI), tol (declared tolerance), auto (exact on Linux)")


@pytest.fixture(scope="session")
def match_mode(request):
    from ttmt.runtime import is_exact_platform
    m = request.config.getoption("--match")
    return ("exact" if is_exact_platform() else "tol") if m == "auto" else m



def pytest_runtest_setup(item):
    for mk in item.iter_markers(name="rawdata"):
        if "TTMT_DATA_ROOT" not in os.environ:
            pytest.skip("raw-data regression: TTMT_DATA_ROOT is not set")
        for ds in mk.args:
            if not config.has_dataset(ds):
                pytest.skip(f"raw dataset {ds} not found under TTMT_DATA_ROOT={config.data_root()}")
