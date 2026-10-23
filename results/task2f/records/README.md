# Task 2F records (regenerated)

`task2f_parts.tar.xz` is a **regeneration** of the task 2F parts with the frozen code (`OPENBLAS_CORETYPE=Prescott`),
not the original registered records, which were never stored (E49). Equivalence to the registered run is shown by
reproducing every registered number except run times; see `results/rescore_duplicates/README.md`.

Check and extract from `code/` (the frozen scripts read `code/task2f_parts`):

```bash
cd code
(cd ../results/task2f/records && sha256sum -c SHA256SUMS)
tar -xJf ../results/task2f/records/task2f_parts.tar.xz
```
