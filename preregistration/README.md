# Preregistration

Every confirmatory task (2A to 2L, R0, 2M) was written into an append-only log before its run, the code was frozen with
SHA-256 checksums before evaluation, and the result and its independent verification were appended afterwards. After
each append, the SHA-256 of the whole log and the UTC time were added to the hash chain, and from Task 2E on the
registration hash was sent to the author in the conversation before any run.

## Files

| File | Content |
|:--|:--|
| `PREREGISTRATION_LOG_EN.txt` | English translation of the log, entry by entry, for reading. |
| `original_sealed_es/PREREGISTRO_1j_1k.txt` | The log as sealed, in the original Spanish. This is the file the hashes certify; it cannot be translated without breaking the chain. |
| `original_sealed_es/PREREGISTRO_hashes.txt` | The hash chain: UTC time and SHA-256 of the sealed log after each append. |
| `drafts/task2M_registration_DRAFT.txt` | Draft of Task 2M (G1-R with real noise on independent days) as prepared before the author's three protocol decisions; the registered text is the entry "Tarea 2M" of the sealed log (16 October 2026), which resolves them. |

## Verifying the chain

```bash
cd original_sealed_es
tail -1 PREREGISTRO_hashes.txt                 # last recorded hash
sha256sum PREREGISTRO_1j_1k.txt                # must match it
```

Earlier hashes certify earlier states of the same append-only file: each earlier state is a byte prefix of the
final log. The copies of the log that were stored next to the results of each task (sizes 47540, 56925, 57567 and
63478 bytes) were checked to be exact prefixes of the final file and are therefore not duplicated here; their hashes
are in the chain (for example `ca3f62fc...2d38` for the state after the verification of Task 2L).

## Section names

Documents in this repository cite sections of the log by their English names ("Task 2I", "Result of 2I",
"Independent verification of 2I"). In the sealed Spanish file the same sections read "Tarea 2I", "Resultado de 2I"
and "Verificacion independiente de 2I"; the English log keeps the same line numbering as the sealed file, so a line
cited in one is the same line in the other.

## Code freezes

The checksums of the frozen code of each task are in `code/` (`task2d_code_hashes.txt` to `task2l_code_hashes.txt`),
`code/explore_contrib/` (`C11_PLAN.sha256`, `C14_FREEZE.txt`, `C16_FREEZE.txt`), `code/r0/R0_FREEZE.txt` and
`results/task2M/T2M_FREEZE.sha256`. Paths
inside them are those of the original working directory (`/home/claude/work/...`); check them from a tree with the
same layout or by file name.
