# Task 2I records (regenerated)

`task2i_parts.tar.xz` is a **regeneration** of the task 2I primary-arm (`prim`) parts with the frozen code (`OPENBLAS_CORETYPE=Prescott`),
not the original registered records, which were never stored (E49). Equivalence to the registered run is shown by
reproducing every registered number except run times; see `results/rescore_duplicates/README.md`.

Check and extract from `code/` (the frozen scripts read `code/task2i_parts`):

```bash
cd code
(cd ../results/task2i/records && sha256sum -c SHA256SUMS)
tar -xJf ../results/task2i/records/task2i_parts.tar.xz
```
