# Data card: Alto Tiberina 2010, noise days

- Source: repository `avuan/PyMPA37`, commit `4a14b134afec3c1f527b738ca3bae4cdc7cbbbf6`, folder `main.pympa.dir/24h/`
  (GPL-3.0). Original data: INGV, network IV (Alto Tiberina Near Fault Observatory); cite the network and the PyMPA
  authors.
- Record: continuous three-component data decimated to 25 Hz; 2010-07-23 (14 stations) and 2010-02-10, 2010-02-11,
  2010-02-12, 2010-12-31; station inventory `inv.IV`, `stations.list`; TABOO catalogue `main.pympa.dir/templates.zmap`
  (not the copies in the `input.*` folders).
- Checksums: `data/checksums/altotiberina_20100723_C16.sha256`, `altotiberina_other_days_task2M.sha256`.
- Partitions: the day 2010-07-23 was the development day (C16); the other four days carried the injected events of
  Task 2M (registered), with the primary day containing transients and two near-Gaussian days as the secondary
  replicate. The calibration seed of 2M was used by exploration C18f.
- Uses: C16 (exploratory), 2M (registered, run once; E48 for one resumed lane).
- Status: development only from now on.
