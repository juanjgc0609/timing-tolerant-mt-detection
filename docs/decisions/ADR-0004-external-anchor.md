# ADR-0004: External timestamp of registrations

Date: 18 October 2026. Status: accepted; revised on 18 October 2026 (the first version is kept at the end).

Context. The hash chain of the sealed log shows that stored entries and hashes agree, but its timing rested on a
private chat record; no registration was in the public repository before its run. The first version of this record
took the GitHub commit time as the external timestamp. That is not a third-party time: the author and committer dates
of a git commit are written by the machine that makes the commit and can be set to any value, GitHub does not show
the time of a push in the history, and its event API keeps push events for about 90 days only. A registration and its
result can also arrive in the same push.

Decision.
1. Every registration, amendment and freeze list is committed and pushed, and then anchored by an **annotated tag
   and a GitHub Release** named `reg/<task>/<kind>-<n>` (kind is `registration`, `amendment` or `freeze`) on that
   commit. The `created_at` time of a Release is set by the GitHub server and is readable by anyone through
   `GET /repos/{owner}/{repo}/releases/tags/{tag}`; that time is the external timestamp.
2. The SHA-256 of the registration file is also stamped with OpenTimestamps (`ots stamp`), and the `.ots` proof is
   committed next to it once it is complete. This anchor does not depend on GitHub.
3. Confirmatory tests on real data (phase 6) are, in addition, registered on OSF Registries, which issues a DOI and a
   frozen, dated copy. This needs the author's OSF account.
4. Synthetic evaluation seeds are derived after the anchor from a value nobody could know before it: the first 8
   bytes of the SHA-256 of the commit hash of the registration tag concatenated with the Release `created_at` string.
   The derivation is written in the registration.
5. The run starts only after the Release exists. The result entry of the log cites the tag, the Release time and, when
   available, the OpenTimestamps proof.

Consequences. The order "registered, then run" becomes checkable by a third party for every future test: the
Release time precedes the first archived output, and for synthetic tests the seeds could not have been chosen before
the registration was public. For real data the anchor proves when the analysis was fixed, not that the data were never
looked at; the data partitions of the protocol (Section 2) cover that part. Earlier registrations (2E to R1) keep their
original status and are not re-anchored after the fact.

---

First version (18 October 2026, superseded): "Every registration, amendment and freeze list is committed and pushed to
the public GitHub repository before the run starts; the result entry cites those commits."
