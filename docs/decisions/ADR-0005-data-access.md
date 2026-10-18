# ADR-0005: Data access

Date: 18 October 2026. Status: proposed; requires a decision of the author.

Context. The working environment reaches only GitHub (including release assets and LFS) and PyPI. The large public
archives (SCEDC and NCEDC on AWS, IRIS FDSN, Zenodo, Hugging Face, the SeisBench store) are refused by the egress
proxy. No unused real dataset remains in the repository. The allowed domains of the environment can be extended by the
account owner in the claude.ai settings (Capabilities).

Decision options. (a) The author allows the domains listed in `docs/data_access.md`; (b) the author downloads the data
on a personal computer and attaches it or links that computer; (c) only GitHub-hosted data are used. Option (a) is
recommended: the confirmatory phase needs several unseen sequences from different networks.

Consequences. Until a decision is taken, phases 1 to 3 proceed with the data already held and GitHub-hosted examples.
