# ADR-0004: External timestamp of registrations

Date: 18 October 2026. Status: accepted.

Context. The hash chain of the sealed log shows that stored entries and hashes agree, but its timing rested on a
private chat record; no registration was in the public repository before its run.

Decision. Every registration, amendment and freeze list is committed and pushed to the public GitHub repository
before the run starts; the result entry cites those commits. The sealed log and the hash chain continue unchanged.

Consequences. The order "registered, then run" becomes verifiable by a third party from the commit history.
