# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

himut detects single-molecule somatic single-base substitutions (SBS) from PacBio CCS (HiFi) reads aligned with minimap2 (`--cs` tag and `@RG SM:` tag are required). Python `>=3.8,<3.10`; core libraries are pysam, cyvcf2, pytabix, pyfastx, numpy, scipy, plotnine/matplotlib.

## Commands

```sh
bash install.sh                 # python setup.py install --record files.txt
pip install -e .                # editable install alternative
himut <subcommand> -h           # CLI entry point: himut.__main__:main
python -m himut <subcommand>    # same, without installing the console script
```

Release: `pyproject.toml` (Poetry) is the source of truth; `setup.py` is generated from it (`poetry2setup > setup.py`, see `upload.sh`) and must be regenerated when version or dependencies change. The version string is also hardcoded in `src/himut/__main__.py` (`__version__`) and must be bumped there too. Upload with `twine upload dist/*`.

Tests/lint: `noxfile.py` and `tests/test_main.py` are leftovers from a project template (the test uses `click.testing.CliRunner`, but the CLI is argparse-based), so there is no meaningful test suite. Verify changes by running the relevant `himut` subcommand on a small BAM/VCF region (`--region chr`).

Note: `scipy` is imported by `phaselib.py` but is not declared in `pyproject.toml`/`setup.py`.

## Architecture

- `parse_args.py` defines all argparse subcommands (`call`, `phase`, `sbs52`, `sbs96`, `sbs1536`, `burden`, `tricount`, `normcounts`); `__main__.py` dispatches each to a single top-level function, passing every option positionally. Adding a CLI option means touching `parse_args.py`, the call in `__main__.py`, and the long positional signature of the target function (and usually the worker-arg tuple it builds).
- `caller.py` (`call`): `call_somatic_substitutions` loads contigs/chunks (`util.load_loci`), optional phased het SNPs (`vcflib.load_phased_hetsnps`), read-length/mismatch thresholds (`bamlib.get_thresholds`), builds the VCF header, then runs `get_somatic_substitutions` once per chromosome via `multiprocessing.Pool.starmap`. Workers write results into `Manager().dict()`s keyed by chromosome, which are written out in chromosome order by `vcflib.dump_sbs` / `dump_phased_sbs`.
- `normcounts.py` (`normcounts`) re-walks the BAM with the **same read/base filters as `caller.py`** to count callable trinucleotide contexts and normalise SBS96 counts. Filter logic changes in the caller must be mirrored here, and users are expected to run both with identical parameters.
- `phaselib.py` (`phase`): builds a graph of het SNP pairs linked by CCS reads (binomial test on allele co-occurrence) to produce haplotype blocks; same per-chromosome Pool/Manager pattern.
- Helper libraries: `bamlib` (BAM metadata, thresholds, trimming/mismatch-window checks), `cslib` (parses minimap2 `cs` tag into substitutions/indels/positions), `haplib` (assigns reads to haplotypes via het-SNP hamming distance), `gtlib` (germline genotype likelihoods/GQ), `somlib` (somatic vs. non-somatic likelihoods), `vcflib` (VCF I/O, headers, germline priors for `--non_human_sample`), `mutlib` (SBS52/96/1536 counting, plotting, burden), `reflib` (reference trinucleotide counts), `util` (input validation, loci loading, PoN default params via `load_pon_params`).
- `--create_panel_of_normal` overrides user thresholds with relaxed values from `util.load_pon_params()`.

`scripts/` holds standalone analysis utilities (spectrum plotting, HDP/SigProfiler matrix prep, FASTQ/BAM stats, etc.). They are not part of the installed package and are run directly with `python scripts/<name>.py`.
