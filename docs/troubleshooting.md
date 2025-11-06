## Troubleshooting and FAQ

### Common issues

- Manifest paths incorrect: update absolute paths in your manifest files.
- Docker runs out of memory: increase Docker memory to 8 GB+.
- Deblur and sample names: avoid underscores; use alphanumerics, dashes, periods.
- Long QZV loads (Empress): allow >10 minutes in browser.

### Tips

- Use Snakemake `--printshellcmds` and `--dryrun` to inspect planned execution.
- Delete individual outputs to force regeneration, or target a file/rule directly.
- To archive/compare runs, change `run_name` or copy output directories before re-running.

### Performance

- Set `--cores` according to machine capacity; certain steps parallelize well (e.g., taxonomy).
- For large datasets, test truncation/trim parameters on a subset first.

### Where to ask questions

- Tourmaline GitHub issues: `https://github.com/aomlomics/tourmaline/issues`
- Project wiki and README for additional examples.


