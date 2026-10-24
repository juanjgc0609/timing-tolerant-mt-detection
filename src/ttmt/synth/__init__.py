"""Synthetic generator v2: full-wave layered-medium synthetics independent of the detector's forward operator.

Modules: ``media`` (layered models, own ray tracer for the true arrival times), ``fk`` (wavenumber-integration
Green's functions), ``analytic`` (whole-space reference), ``source`` (mechanisms and source-time functions),
``coda``, ``noise``, ``scenario`` (windows and truth), ``sets`` (dataset build and partitioned reading).
None imports the detector modules (checked by tests/test_synth.py)."""
