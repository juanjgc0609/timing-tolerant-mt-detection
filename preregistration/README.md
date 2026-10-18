# Preregistration

Every confirmatory task (2A to 2L, R0, 2M, R1) was written into an append-only log before its run, and its result and,
where one was made, its verification were appended afterwards. After each append, the SHA-256 of the whole log and the
UTC time were added to the hash chain, and from Task 2E on the registration hash was sent to the author in the
conversation before any run. The chain begins inside the registration of Task 2D, so for 2A to 2C it certifies only
that the text has not changed since 26 September 2026, not that the registration preceded the run; code hash lists
exist from 2D on. Table 3 of `docs/thesis/Thesis.md` gives, task by task, which of these steps each test had.

## Files

| File | Content |
|:--|:--|
| `PREREGISTRATION_LOG_EN.txt` | English translation of the log, entry by entry, for reading. |
| `original_sealed_es/PREREGISTRO_1j_1k.txt` | The log as sealed, in the original Spanish. This is the file the hashes certify; it cannot be translated without breaking the chain. |
| `original_sealed_es/PREREGISTRO_hashes.txt` | The hash chain: UTC time and SHA-256 of the sealed log after each append. |
| `drafts/taskR1_registration_DRAFT.txt` | Draft of Task R1, the registered real-data test on the unread Hi-net swarm that follows the gate met in exploration C18j; not registered, with four decisions left to the author. |
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

## Format of the hash file

`PREREGISTRO_hashes.txt` has, for each append, a line with the UTC time followed by a line with the SHA-256 of the
whole log at that moment; the first hash has no time line, and one time line (5 October 2026) is not followed by a hash; every hash in the
file is nevertheless the digest of a byte prefix of the log, and the last one is the digest of the whole file. The English translation keeps the line
breaks between paragraphs and the number of lines of every paragraph, but within a long paragraph a word can fall
one line earlier or later than in the Spanish file, so a cited line should be read with its neighbours.

## Section names

Documents in this repository cite sections of the log by their English names ("Task 2I", "Result of 2I",
"Independent verification of 2I"). In the sealed Spanish file the same sections read "Tarea 2I", "Resultado de 2I"
and "Verificacion independiente de 2I"; the English log keeps the same line numbering as the sealed file, so a line
cited in one is the same line in the other.

## Code freezes

The checksums of the frozen code of each task are in `code/` (`task2d_code_hashes.txt` to `task2l_code_hashes.txt`),
`code/explore_contrib/` (`C11_PLAN.sha256`, `C14_FREEZE.txt`, `C16_FREEZE.txt`), `code/r0/R0_FREEZE.txt` and
`results/task2M/T2M_FREEZE.sha256` and `code/r1/FREEZE_R1.sha256`. Paths
inside them are those of the original working directory (`/home/claude/work/...`); check them from a tree with the
same layout or by file name.
