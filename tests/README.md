# Tests

```
tests/
├── run_tests.sh              the integration suite (needs tools, databases, SLURM)
├── submit_tests.sh           run it as a SLURM job
└── data/
    ├── rawfastq/             two paired-end samples, sub-sampled ECHO reads
    ├── single_end_rawfastq/  one single-end file, plus the bioBakery tutorial demo sample
    └── tutorial_output/      reference output from the upstream tutorial run (KneadData,
                              MetaPhlAn, HUMAnN, BAQLaVa) for HD32R1_subsample
```

`tests/results/` and `tests/work/` are created by a run and gitignored.

## What is checked where

| Check | Needs | Where |
|---|---|---|
| Registry, component files, orchestrator and docs agree | Python | `bin/make_diagrams.py --check`, CI |
| No site profile drops a resource request | Nextflow | `bin/check_profile_resources.py`, CI |
| Every `--components` / `--only` selection builds a valid graph | Nextflow | `nextflow run -preview`, CI |
| Components produce the right output on real data | tools, databases, SLURM | `tests/run_tests.sh` |

## The integration suite

```bash
cd /path/to/HVP-nextflow
bash tests/run_tests.sh                                  # on a login node, submits SLURM jobs
mkdir -p tests/results && sbatch tests/submit_tests.sh   # or as one batch job
```

The cases are listed at the top of `run_tests.sh`: BAQLaVa with and without
KneadData, the placeholders (they run, warn, and produce nothing), the assembly
component, `--only` dependency resolution, and an unknown component name
being refused at launch. The assembly cases need the simulated MAG fixture
(`$ASSEMBLY_FIXTURE`, default `~/biobakery_assembly_test`) and are skipped
without it.
