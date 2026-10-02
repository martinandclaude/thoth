# thoth

<https://martinandclaude.github.io/thoth/>

A pipeline-neutral GIAB benchmarking harness. Give it a run file listing each
sample's VCFs, per caller; it scores every caller against every applicable GIAB
truth set — small variants with hap.py/vcfeval, stratified; structural variants
with truvari — and emits one long-format metrics table plus a provenance record.

Named for the Egyptian god of measurement and record-keeping, who weighed and
wrote down the result. It treats every caller equally and is tied to none.

## Requirements

- Linux x86_64 with Apptainer on the submit host and every compute node (the
  images are linux/amd64; CI uses Apptainer 1.4.5).
- SLURM, through the profile below, or a single host (see the end of Quick start).
- Python 3.12 for the driver environment (Snakemake and the script rules).
- On the host: `curl`, `sha256sum` and `tar`, which the fetch and unpack rules
  call outside any container, and `git` to record the code's commit.
- Outbound HTTPS from the submit host, for `fetch_containers.sh` and the fetch
  rules; without it, set `fetch_on_login_node: false` in `config/site.yaml`.
- About 14 GB of shared disk: the GRCh38 truth sets and stratifications
  (`check_manifest` logs the download total), the unpacked stratifications, the
  reference SDF, hap.py's copy of the reference FASTA and the SIFs.
- 128 GB of RAM per hap.py job (`small_variants.happy_mem_mb`).
- The reference FASTA your VCFs were called against.

## Quick start

```bash
git clone https://github.com/martinandclaude/thoth && cd thoth
python3.12 -m venv .venv && . .venv/bin/activate
pip install --require-hashes -r requirements-driver.lock

cp config/site.example.yaml config/site.yaml   # then set resources_dir and references
cp -r profiles/slurm.example profiles/slurm    # then set slurm_account, slurm_partition, apptainer-args

workflow/scripts/fetch_containers.sh --test         # download SIFs by digest, verify, test (PENDING: see Containers)
snakemake --profile profiles/slurm check_manifest   # every URL, before the first download
snakemake --profile profiles/slurm resources_all    # truth sets, stratifications

cp runs/run.example.yaml runs/my-run.yaml      # then set run_id and queries
snakemake --profile profiles/slurm -n --configfile runs/my-run.yaml   # dry run: read the plan and notes
snakemake --profile profiles/slurm --configfile runs/my-run.yaml
```

Run everything from the repository root, in the activated environment. The
Snakefile loads `config/config.yaml` and `config/site.yaml` itself; pass only
the run file on the command line. `config/config.yaml` deliberately carries no
site keys: Snakemake applies command-line config files *after* the Snakefile's,
so an empty `resources_dir:` there would silently overwrite yours.

Beware one CLI trap: `--configfile` takes several values, so
`snakemake --configfile run.yaml check_manifest` reads `check_manifest` as a
second config file. Put targets first.

Without SLURM, on one host, replace `--profile profiles/slurm` with
`--cores N --software-deployment-method apptainer --apptainer-args "--bind /your/shared/fs"`.

## Layout

```
config/config.yaml          global settings and the benchmark matrix
config/site.example.yaml    template for config/site.yaml, this cluster's paths (gitignored)
profiles/slurm.example/     template for profiles/slurm/, executor settings (gitignored)
runs/run.example.yaml       template for runs/*.yaml, one file per run (gitignored)
resources/manifest.tsv      every truth set and stratification: URL and sha256
containers/<image>/         the recipe for every image thoth runs
workflow/
  Snakefile
  containers.tsv            tool -> GHCR image, SIF digest, smoke test
  rules/*.smk
  scripts/                  script rules and fetch_containers.sh
tests/                      planted-truth fixtures and the end-to-end test
docs/                       the Pages site and the scripts that write its data
.github/workflows/          containers (build, test, attest, publish), tests (end to end)
requirements-driver.txt     driver environment, direct pins; .lock is the hashed closure
```

## The run file

One entry per sample on a reference, with any number of callers:

```yaml
run_id: my-run
queries:
  - name: sample-HG002
    sample: HG002
    reference: GRCh38
    vcf_small:            # small variants -> hap.py
      deepvariant: /path/deepvariant.vcf.gz
      clair3: /path/clair3.vcf.gz
    vcf_sv:               # structural variants -> truvari
      sniffles: /path/sniffles.vcf.gz
      cutesv: /path/cutesv.vcf
```

Each (query, caller) VCF is scored against every truth set the manifest has for
that sample; the metrics table carries a `caller` column. `sample` is the GIAB
ID (HG001-HG007). Plain `.vcf` is fine for either class, but the header must
declare its contigs: a contig the reference lacks fails the run. Unknown keys
fail at parse time, in the run file and in the config, so a misspelt `vcf_svs:`
is an error, not a run that silently has no SV rows. See
`runs/run.example.yaml` for every key.

Per-arm outputs live in `results/<query>/<caller>/<truth>/` and are shared by
every run file that uses that query name, so a query name must always mean the
same VCFs. Give different inputs a new name.

## Outputs

```
results/report/<run_id>.metrics.tsv      the result, long format
results/report/<run_id>.runinfo.json     provenance record (below)
results/share/<run_id>/                  what may leave this machine
results/<query>/<caller>/prep/           contig check, prepped VCF, variant counts
results/<query>/<caller>/<truth>/happy/  hap.py output, incl. the annotated VCF
results/<query>/<caller>/<truth>/truvari/  truvari bench (and refine) output
results/report/manifest_check.tsv        from check_manifest
```

`metrics.tsv` has one value per row, in the columns run_id, query, caller,
sample, reference, origin, share, score_field, truth, strat, strat_subset, tool,
type, subtype, subset, filter, genotype, qq_field, metric, value.

`results/share/<run_id>/` holds the metrics table and `MANIFEST.txt` (code
version, config_sha256, container digests) for every query, plus the annotated
hap.py VCF, with every small-variant TP, FP and FN, for queries marked
`share: public`. SV per-variant files and `runinfo.json`, which records user,
host and input paths, are never staged.

## Rule graph

```
resources/manifest.tsv
        ├─ check_manifest ─────────► report/manifest_check.tsv
        └─ fetch_resource ─┬─ truth VCF (+ fetch_index: .tbi) + confident BED
                           │    └─ prep_sv_truth (drop ALT=*) ─► SV truth for truvari
                           └─ stratification tarball
                                └─ extract_stratifications
                                     └─ absolutize_strat_tsv ─► <subset>.<sample>.abs.tsv

reference FASTA ─┬─ rtg_format (hap.py image) ─► <ref>.sdf
                 ├─ faidx_reference ──────────► <fasta>.fai
                 └─ happy_reference ──────────► <ref>.iupacN.fa (ambiguity codes as N, for hap.py)

run file (queries: sample × callers)
        ├─ per caller in vcf_small: check_contigs ─► prep_query ─► prepped.vcf.gz
        └─ per caller in vcf_sv:    check_contigs ─► prep_sv_query (sort + index)
                                                          │                │
                               per (query × caller × truth)                │
                                                          ▼                ▼
                                          truvari / truvari_refine       happy
                                          (HG002 only: v5.0q with      (vcfeval, own-sample
                                           refine, CMRG)                stratifications)
                                                          └► tidy_metrics ◄┘
                                                                 └► runinfo
                                                                      └► share_bundle
```

`resources_all` fetches every truth row in the manifest, and the
stratifications of each reference that has truth rows (so far GRCh38 only).

## What thoth changes from the defaults

Refine for GIAB v5.0q, `dup-to-ins` and `max-resolve` in truvari, `ALT=*`
removed from SV truth, stratification by the query's own sample only,
half-calls such as `./1` kept, and hap.py given the reference with its IUPAC
ambiguity codes as N, as the truth sets write them. Each is covered by a
planted-truth test in `tests/`, and each is explained, with real results, on the
site: <https://martinandclaude.github.io/thoth/#settings>.

Small variants are scored on chr1-22 by default (`small_variants.regions`); SVs
on every chromosome in the truth BED.

Two things to know before running: **SV truth sets exist for GIAB HG002 only**
(SV calls for other samples are skipped with a note), and **insertion precision
is a size-and-position precision** under truvari's default sequence-similarity
threshold.

## Containers

Every variant tool runs from a SIF published by this repository to
`ghcr.io/martinandclaude/thoth/<image>` and pinned by digest in
`workflow/containers.tsv`: `hts` (bcftools, samtools) and `truvari` are built
here from `containers/`; `happy` is the community image `jmcdani20/hap.py`,
pinned by digest in `containers/happy/SOURCE` and converted unchanged. Its
bundled RTG runs vcfeval and also builds the reference SDF. Downloads, unpacking
and file copies use the host's tools; Snakemake and the script rules run in the
driver environment.

The containers workflow converts each image to a SIF on GitHub and publishes it
as an ORAS artifact: a single registry blob whose digest is the file's own
sha256. `workflow/scripts/fetch_containers.sh` downloads it over plain HTTPS
and refuses any file that does not match, so the digest is both pin and
integrity check, and no node ever pulls or converts an image, a step many
clusters cannot do.

**Before the first publish** every digest in `workflow/containers.tsv` reads
`PENDING`: `fetch_containers.sh` reports the image as not yet published and
exits 1, the tests workflow fails, and any job needing that image fails on its
missing path. The containers workflow builds and tests recipe changes on pull
requests and publishes only from `main`. After it has run on `main`, copy each
image's `digest` and `built_from` (the commit the recipe was built from) from
its job summary into `workflow/containers.tsv` in a reviewed commit. Each SIF
also carries that commit as the label `org.opencontainers.image.revision`
(`apptainer inspect <file>.sif`), and the provenance record keeps both.

`fetch_containers.sh` downloads anonymously, so each GHCR package (`hts`,
`truvari`, `happy`) must be public; GHCR creates new packages private, so change
the visibility in each package's settings after the first publish. Each SIF's
build provenance is attested; to check a downloaded file:

```bash
gh attestation verify <file>.sif --repo martinandclaude/thoth
```

The smoke tests are the `test_cmd` column of `workflow/containers.tsv`, so CI
and `fetch_containers.sh --test` run the same checks. Recipes pin base images by
digest, sources by sha256 and truvari's Python packages by hash; Debian
packages come from the live archive, so a rebuild matches in function
(`tests/run.sh`), not byte for byte, and the SIF digest is the exact pin. The
Apptainer version CI uses is pinned in `.github/actions/apptainer/action.yml`.

## Provenance

`resources/manifest.tsv` pins every truth and stratification file by sha256.
`fetch_resource` writes the sha256 of what it downloaded beside the file, and a
mismatch with the pin fails the run at parse time and again when the provenance
record is written; an empty pin is reported as unpinned.

Each run writes `results/report/<run_id>.runinfo.json`: the git commit and dirty
flag and a code_sha256 that identifies the code without git, the merged config
and its config_sha256, every image digest, the sha256 of each truth and
stratification file this run used (and whether the manifest pins it), the .fai
checksum of each reference, per-caller input sizes, and hap.py's own command
line per arm.

## Tests

`tests/run.sh` runs the real pipeline, in the real containers, on planted data
whose correct answers are known exactly (`tests/make_fixtures.py`), and asserts
every count (`tests/check_results.py`). It needs Apptainer and the driver
environment:

```bash
THOTH_CONTAINERS_DIR=tests/work/sifs workflow/scripts/fetch_containers.sh --test
tests/run.sh
```

CI runs it on every push to `main` and every pull request. To see the plan
without Apptainer: `python3 tests/make_fixtures.py`, then
`snakemake -n --configfile tests/work/test.yaml`.

## The site's data

`docs/index.html` reads `docs/data/evidence.json`: Oxford Nanopore's public
GIAB 2025.01 calls scored by thoth, beside ONT's published numbers and a rerun
of ONT's SV recipe. With the truvari SIF from `workflow/containers.tsv`:

```bash
apptainer exec <containers_dir>/truvari-<digest>.sif python docs/seqsim_baseline.py
python docs/extract_evidence.py --runinfo <run>.runinfo.json --repro <repro> <run>.metrics.tsv
```

`docs/data/README.md` describes the comparison run and `<repro>`. Only
aggregate rates and counts for public GIAB samples leave the metrics table;
review the JSON before committing.

## Not finished, deliberately

- **hap.py from source.** Converted unchanged from a pinned community image for now.
- **GRCh38 truth sets only** in the manifest. GIAB publishes GRCh37 and CHM13v2.0
  v5.0q as well — in the same `v5.0q/` directory, which, unlike every other GIAB
  release, has no `NIST` prefix and no per-reference subdirectory.
- **No plots.** `metrics.tsv` is long-format and plots in a few lines of ggplot or
  seaborn; chart code would make this a reporting tool rather than a scoring one.
- **hap.py's memory is headroom.** `happy_mem_mb` was set above a measured
  failure, not from a measured peak; read MaxRSS with `sacct` after a run and
  tighten it. The refine rule's resources come from a measured genome-wide run
  (see `workflow/rules/bench_sv.smk`).

## Citing

Cite thoth from `CITATION.cff` (GitHub: "Cite this repository"), and the truth
sets and tools whose numbers it reports:

- GIAB v4.2.1: Wagner J et al. Benchmarking challenging small variants with
  linked and long reads. *Cell Genomics* 2, 100128 (2022).
  <https://doi.org/10.1016/j.xgen.2022.100128>
- GIAB CMRG v1.00: Wagner J et al. Curated variation benchmarks for challenging
  medically relevant autosomal genes. *Nat Biotechnol* 40, 672–680 (2022).
  <https://doi.org/10.1038/s41587-021-01158-1>
- GIAB v5.0q: unpublished; its release README asks users to cite the T2T HG002
  Q100 project: Hansen NF et al. A complete diploid human genome benchmark for
  personalized genomics. bioRxiv (2025).
  <https://doi.org/10.1101/2025.09.21.677443>
- GIAB stratifications v3.6: Dwarshuis N et al. The GIAB genomic
  stratifications resource for human reference genomes. *Nat Commun* 15, 9029
  (2024). <https://doi.org/10.1038/s41467-024-53260-y>
- hap.py: Krusche P et al. Best practices for benchmarking germline
  small-variant calls in human genomes. *Nat Biotechnol* 37, 555–560 (2019).
  <https://doi.org/10.1038/s41587-019-0054-x>
- RTG vcfeval: Cleary JG et al. Comparing variant call files for performance
  benchmarking of next-generation sequencing variant calling pipelines.
  bioRxiv 023754 (2015). <https://doi.org/10.1101/023754>
- truvari: English AC et al. Truvari: refined structural variant comparison
  preserves allelic diversity. *Genome Biol* 23, 271 (2022).
  <https://doi.org/10.1186/s13059-022-02840-6>

## Licence

thoth is released under the MIT licence (`LICENSE`). Logo generated with ChatGPT
(OpenAI).

The published images keep their own licences:

- `hts`: bcftools, samtools and htslib, MIT/Expat (htslib's CRAM code: modified
  BSD); the texts are in the image under `/usr/share/doc/<pkg>/LICENSE`.
- `truvari`: truvari, MIT, and the Python packages in
  `containers/truvari/requirements.lock`, each with its licence in its dist-info.
- `happy`: Illumina's hap.py, simplified BSD, in the community image
  `jmcdani20/hap.py`, republished unchanged; its bundled RTG Tools is BSD
  2-clause, and other bundled components keep their own licences.

Base-OS package copyrights are under `/usr/share/doc/` in each image.
