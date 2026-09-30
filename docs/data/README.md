# Data behind the Pages site

| File | What it holds | Written by |
|---|---|---|
| `evidence.json` | Every number `index.html` shows | `../extract_evidence.py` |
| `seqsim_baseline.json` | truvari's sequence similarity for unrelated random sequence pairs | `../seqsim_baseline.py` |
| `ont_giab_2025.01_published.tsv` | ONT's own published GIAB 2025.01 benchmark numbers, with the URL of each | copied by hand from ONT's public outputs |

Only aggregate counts and rates for public GIAB samples leave a metrics table.
Query and caller names, users, hosts and paths do not; `extract_evidence.py`
takes provenance from a runinfo file through a whitelist (metrics and config
checksums, the thoth commit, container digests).

Most of this file is about one section of the page: ONT's published benchmark,
reproduced and then rescored with thoth's settings. The other sections are
summarised briefly at the end.

---

## ONT's GIAB 2025.01 benchmark: reproduced, then scored thoth's way

### What was compared, and why

Oxford Nanopore released GIAB sequencing (2025.01, sup basecalling, two flow
cells per sample for HG001-HG007), its wf-human-variation calls (v2.4.1: Clair3
for small variants, Sniffles2 for SVs), and its own benchmark of one flow cell
per sample: hap.py against GIAB v4.2.1 for all seven samples, truvari for the
HG002 SVs.

A benchmarking pipeline that claims its settings matter should first be able to
give the published numbers from the published calls. So the comparison has three
steps:

1. Score ONT's public call files with thoth and check that the small-variant
   results match ONT's.
2. Recover ONT's SV recipe (ONT publishes no truvari parameters) and rerun it,
   to show that the scoring, not the calls, explains the SV differences.
3. Score the same SV calls with thoth's settings, on the truth ONT used and on
   GIAB v5.0q.

### Sources

| What | Where |
|---|---|
| Release page | https://epi2me.nanoporetech.com/giab-2025.01/ |
| ONT's calls | `https://ont-open-data.s3.amazonaws.com/giab_2025.01/analysis/wf-human-variation/sup/<SAMPLE>/<FLOWCELL>/output/SAMPLE.wf_snp.vcf.gz` and `SAMPLE.wf_sv.vcf.gz` |
| ONT's hap.py results | `https://ont-open-data.s3.amazonaws.com/giab_2025.01/analysis/happy-benchmark/sup/<SAMPLE>_<FLOWCELL>/happy_GIABv4.2.1.summary.csv` (and the annotated VCF beside it) |
| ONT's truvari results | `https://ont-open-data.s3.amazonaws.com/giab_2025.01/analysis/truvari-benchmark/hg002/SUP_PAW70337/` (`summary.json`, `refine.variant_summary.json`, the bench VCFs, `phab_bench/`) |
| SV truth ONT used | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/ |
| GIAB v5.0q | the URLs in `resources/manifest.tsv` |

`ont_giab_2025.01_published.tsv` holds ONT's PASS rows, read on 2026-09-29, with
the exact source file of each row. The flow cells ONT benchmarked are HG001
PAW81754, HG002 PAW70337, HG003 PAY87794, HG004 PAY88428, HG005 PAW88001, HG006
PBA16846 and HG007 PAY12990.

### Step 1: small variants

ONT's `SAMPLE.wf_snp.vcf.gz` for each benchmarked flow cell went into a thoth
run as a `vcf_small` query (hap.py 0.3.12 with vcfeval, PASS only, chr1-22,
GIAB v4.2.1 and, for HG002, v5.0q and CMRG). Against ONT's 14 published rows
(7 samples, SNP and indel):

- every truth total is identical;
- recall and precision differ by at most 0.012 percentage points (HG006 indel
  precision); every other gap is under 0.01 point;
- TP, FN and FP counts differ by at most 68 per sample and type, except HG004
  SNPs, where thoth finds 231 fewer false positives (and 20 fewer true
  positives).

The page's table shows every row.

**Where the residue comes from.** On HG004 PAY88428 chr20-22, where it was
tested, the difference is in the input, not in the scoring: ONT's hap.py ran
(VCF `fileDate` 2024-12-03) on an earlier Clair3 call set than the
`SAMPLE.wf_snp.vcf.gz` it published (concatenated and annotated 2025-01-10).
The other samples were not tested this way; their hap.py outputs are also dated
early December 2024. What was established on HG004 chr20-22:

- *ONT's benchmarked query, rebuilt.* hap.py's output VCF keeps every query
  record in its QUERY column. Taking the records with an ALT allele in their
  genotype, and hap.py's QQ as the score, rebuilds the call set ONT actually
  scored (318,023 records on chr20-22). Scored with thoth's exact hap.py path,
  it reproduces ONT's decisions record for record: 0 of 318,023 differ, and
  every SNP and indel TP, FN, FP, UNK and total equals ONT's. The same rebuild
  from thoth's own hap.py output gives thoth's numbers exactly (the control).
- *A hotspot.* In chr20:26,072,136-26,102,512 (30 kb), ONT's hap.py VCF has 778
  SNP false positives. 393 of them are at positions with no record at all in
  the published VCF. All 389 of thoth's false positives there are in it.
- *Chunk pattern.* The two call sets differ in whole Clair3 5-Mb chunks.
  Seven 5-Mb chunks on chr20-22 differ (6-39 % of shared records differ in QUAL
  or genotype; for example chr20:25-30 Mb 19 %, chr21:40-45 Mb 39 %), against
  0.3-1.0 % in all the others. This is consistent with a partial re-run; why
  ONT's run differs is not known.
- *What does not matter.* On the same query, none of these changed a single
  decision: hap.py 0.3.12 with RTG 3.10.1 against RTG 3.11 (and hap.py 0.3.14
  and 0.3.15 with RTG 3.11); the no-alt reference against the hs38d1 analysis
  set with ambiguity codes as N; stratification; `--roc QUAL` against `--roc
  GQ`.

Genome-wide, the net differences are smaller than on chr20-22 (HG004 SNP FP
231 fewer genome-wide against 388 fewer on chr20-22), so other chunks go the
other way. That is inferred from the totals; it was not checked contig by
contig. Scoring ONT's published VCFs gives thoth's numbers; on HG004 the
published summary describes an earlier call set.

**Ambiguity codes.** hap.py rewrites a query's REF from the reference. GRCh38
keeps 94 IUPAC ambiguity bases (B, K, M, R, S, W, Y) where truth sets and
callers write N, and vcfeval stops at the first mismatch: ONT's HG002 calls
against v5.0q fail at chr3:16902879 (`TGTGB` against `TGTGN`). thoth gives
hap.py a copy of the reference with those bases as N (rule `happy_reference`);
every other step reads the original FASTA. Without it, those calls cannot be
scored against v5.0q at all.

**ROC field.** ONT's hap.py QQ equals the call's QUAL, so ONT ran `--roc QUAL`;
thoth uses GQ (`small_variants.default_score_field`). This changes ROC curves
only, not counts.

### Step 2: ONT's SV recipe, recovered and rerun

ONT's truvari outputs record no parameters (`log.txt` is empty, there is no
`params.json`). The recipe was recovered from the outputs themselves and
confirmed by reproducing the counts:

- **Truth: the GIAB draft V0.018** (`NIST_HG002_DraftBenchmark_defrabbV0.018-20240716`),
  unmodified, not v5.0q. Every bcftools provenance line in the header of ONT's
  `tp-base.vcf.gz` matches this draft; none matches the November 2024 draft.
- **truvari 4.3.1.** The set of output files, the 100-bp padding of every
  refined region and the absence of `ALT=*` records fit 4.3.1 and no other
  release; truvari 5.4.0 cannot reproduce the counts with any flags tried.
- **bench `--pick ac --passonly`, everything else default** (refdist 500,
  chunksize 1000, pctseq 0.7, pctsize 0.7, pctovl 0, sizemin 50, sizefilt 30,
  sizemax 50000). Changing any one of pick mode, `--passonly`, refdist or
  sizefilt, or dropping the reference, changes the counts; chunksize 1000 is
  confirmed by the byte-identical `candidate.refine.bed`.
- **refine `--recount --use-region-coords --use-original-vcfs --align mafft`**
  over bench's `candidate.refine.bed`: GIAB's refine recipe for this draft.

The commands, run with the biocontainers truvari 4.3.1 image
(`https://depot.galaxyproject.org/singularity/truvari:4.3.1--pyhdfd78af_0`,
sha256 `e1ae92b51307f7c454816f43e0cdfd19dc8b09e9c9574241599b203ef4c08ba2`;
pysam 0.22.1 and an x86 mafft 7.526, which refine's counts depend on):

```sh
truvari bench -b GRCh38_HG2-T2TQ100-V1.1.vcf.gz -c SAMPLE.wf_sv.vcf.gz -o bench \
  --includebed GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed -f GRCh38.fa \
  --pick ac --passonly
truvari refine --recount --use-region-coords --use-original-vcfs --align mafft \
  --reference GRCh38.fa --regions bench/candidate.refine.bed -t 8 bench
```

Run it as `apptainer exec --env TMPDIR=$PWD/tmp --bind $PWD --bind <reference dir>
truvari_4.3.1--pyhdfd78af_0.sif truvari ...`, with `tmp/` created first:
truvari 4.3.1 passes its temporary directory to `bcftools sort --temp-dir`,
which treats it as a name prefix, so the default `/tmp` would make bcftools
write into the container's read-only root.

Inputs and their sha256:

| File | sha256 |
|---|---|
| ONT `HG002/PAW70337/output/SAMPLE.wf_sv.vcf.gz` | `ab4108ff572d378c5bbd01b11f1c74b3fce96cb09369385cabc9b36cace7cb02` |
| its `.tbi` | `90e389050efe8b05e980b738c1d1d71d0d110fbcbd06946ef0324511119a2af9` |
| GIAB V0.018 `GRCh38_HG2-T2TQ100-V1.1.vcf.gz` | `25e2a5f9dadbd69ab7321fe0068938d18ce45cdd70b6987f66622dd216c90f89` |
| its `.tbi` | `c2e9340274ff4e13f9a508865da5c49190fa9621901fd66aedd7b096effeafa1` |
| GIAB V0.018 `GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed` | `1c5b600d448953d4350c23528c4d4a4ede833006827c3f348bd966382375c5cb` |

ONT used the GRCh38 no-alt analysis set; the rerun used the full_plus_hs38d1
analysis set, which gives identical bench counts and an identical
`candidate.refine.bed`.

**Result.**

| | truth SVs | found | missed | false calls | recall | precision |
|---|---|---|---|---|---|---|
| bench, ONT published | 28,188 | 21,062 | 7,126 | 1,957 | 74.72 % | 90.55 % |
| bench, rerun | 28,188 | 21,062 | 7,126 | 1,957 | 74.72 % | 90.55 % |
| refine, ONT published | 22,966 | 22,072 | 894 | 234 | 96.11 % | 98.82 % |
| refine, rerun | 22,960 | 22,066 | 894 | 233 | 96.11 % | 98.83 % |

Bench matches every count and rate in ONT's `summary.json`; its `tp-base` and
`fn` VCFs match ONT's record for record. Only the phased cells of `gt_matrix`
differ, for the reason below. Refine lands 6 truth SVs and 1 false call away.

**Why refine is not exact: a later version of the call file.** ONT benchmarked a
`SAMPLE.wf_sv.vcf.gz` dated `fileDate=2024/12/18` in its outputs; the public file
says `2025/01/09`. The 20,708 calls ONT scored are all in the public file with
the same position, alleles, filter and length, but 2,242 of them are phased the
other way round (0|1 against 1|0). Bench ignores phase; refine builds
haplotypes from it. With ONT's phasing restored (from ONT's own output VCFs,
and inferred from phase-set neighbours for calls ONT never wrote out),
`refine.variant_summary.json`, `refine.region_summary.json` and
`refine.regions.txt` became byte-identical to ONT's. With the public file the
region summary already matches exactly and the variant summary is 6 truth SVs
off. The mafft build also matters: the same mafft 7.526 built for arm64 moves
refine's truth total by about 20 (22-23 truth SVs, 36-39 calls).

### What `--recount` does to the refine truth total

Bench scores ONT's calls against 28,188 truth SVs; ONT's refine reports 22,966.
With `--recount`, truvari 4.3.1's refine reports two parts added together:

| | truth SVs |
|---|---|
| bench truth SVs outside the 2,955 refined regions, counted as bench counted them | 19,192 |
| inside those regions: the variants phab re-derives from the realigned haplotypes (ONT's `phab_bench`; 3,768 in the rerun) | 3,774 |
| **refine truth total** | **22,966** |

Inside the refined regions bench had 8,996 truth SVs (28,188 − 19,192); after
realignment they are represented by 3,774 harmonised variants. ONT's 96.11 %
refine recall is therefore a share of 22,966, not of the 28,188 truth SVs bench
scores. The inside part also depends on the call set: the rerun changed only the
calls' phasing, and the truth total moved from 22,966 to 22,960.

truvari 5.4.0, which thoth runs, has no `--recount`; its `refine.variant_summary.json`
keeps every original truth record, so thoth's refine truth total equals bench's
(28,188 on V0.018, 28,123 on v5.0q).

`extract_evidence.py` derives the split from the rerun directory: the rerun's
refine total minus its `phab_bench` total is the part outside (19,192), bench's
total minus that is the part inside (8,996), and `refine.regions.txt` gives the
number of refined regions. The 19,192 is the same in ONT's run: in every rerun,
with any phasing or mafft build, refine minus `phab_bench` was exactly 19,192,
and 22,966 − 3,774 (ONT's own `phab_bench`) is 19,192.

### Step 3: the same calls on thoth's settings

thoth scored ONT's HG002 PAW70337 Sniffles calls with its v5.0q settings, on
V0.018 and on v5.0q: truvari 5.4.0, `--refdist 2000 --chunksize 5000 --pick ac
--passonly --pctseq 0.7 --pctsize 0.7 --pctovl 0 --sizemin 50 --sizemax 50000
--dup-to-ins --max-resolve 50000`, `ALT=*` records removed from the truth, then
`truvari refine` with its defaults.

| truth | comparison | truth SVs | found | missed | false calls | recall | precision |
|---|---|---|---|---|---|---|---|
| V0.018 | bench | 28,188 | 21,404 | 6,784 | 1,712 | 75.93 % | 91.75 % |
| V0.018 | refine | 28,188 | 24,425 | 3,763 | 884 | 86.65 % | 95.74 % |
| v5.0q | bench | 28,123 | 21,354 | 6,769 | 1,693 | 75.93 % | 91.82 % |
| v5.0q | refine | 28,123 | 24,375 | 3,748 | 865 | 86.67 % | 95.82 % |

**Where the bench gain comes from.** Measured on HG002 PAW70337 against
V0.018, in truth SVs found (TP-base), from ONT's 21,062 to thoth's 21,404
(+342):

| change | TP-base | step |
|---|---|---|
| ONT's recipe (truvari 4.3.1) | 21,062 | |
| truvari 5.4.0, ONT's settings | 21,214 | +152 |
| ... of which: sequence similarity compared on REF/ALT directly (with rolling), not through the reference | | +144 |
| ... of which: symbolic and breakend calls 5.x keeps and 4.3.1 drops (35 calls: 16 `<INS>`, 9 `<DUP>`, 2 `<INV>`, 8 BND; 8 of the `<INS>` match) | | +8 |
| `-r 2000 -C 5000` (GIAB's recommended window) | 21,401 | +187 |
| `--dup-to-ins`, `--max-resolve 50000`, `ALT=*` removed from the truth | 21,404 | +3 |

The +144 is measured independently: truvari 4.3.1 run without a reference gives
21,206, and 5.4.0 with the symbolic and breakend calls filtered out gives
exactly the same counts. No 5.x flag restores 4.3.1's reference-based
similarity.

### Why thoth scores it this way

- **The window GIAB recommends.** GIAB's v5.0q README recommends `truvari bench
  --refine --pick ac -r 2000 -C 5000`; the V0.018 draft's README gives the same
  window (`--pick ac --passonly -r 2000 -C 5000`, then the refine recipe above).
  ONT left refdist and chunksize at 500 and 1000. On these calls the
  wider window finds about 190 more truth SVs; 286 of ONT's misses had a call
  with at least 70 % sequence and size similarity 500-1,999 bp away.
- **The current truvari.** truvari 5 is the maintained release; 4.3.1's
  `--recount` and `--use-region-coords` are gone from it. Both versions' refine
  need a pinned pysam (below 0.24). Its sequence comparison finds 144 more
  matches on these calls than 4.3.1's; that is a version difference, not
  evidence of better accuracy.
- **A refine recall over every truth SV.** thoth's refine recall is the share of
  the whole truth set found, with the same denominator as bench, so the gain
  from haplotype-aware comparison reads directly (75.93 % to 86.67 % on v5.0q).
  A `--recount` recall is a legitimate measure of agreement after realignment,
  but its denominator is smaller than the truth set and moves with the calls
  being scored and with the aligner build, so it cannot be set beside a bench
  recall or compared across call sets on the same footing. This is a statement
  about what the number means, not about which calls are better: ONT's calls
  are the same in every row.

What this comparison does not show: that thoth's settings are optimal for every
caller or truth set, or that a refine recall computed either way is the "true"
recall. The planted-truth tests behind each thoth setting are described in the
settings table on the page and in the repository's tests.

### Versions

| | ONT | thoth |
|---|---|---|
| small variants | hap.py, RTG 3.11 (per ONT's output header), `--roc QUAL` | hap.py 0.3.12 (image `jmcdani20/hap.py`, bundled RTG 3.10.1), vcfeval, `--roc GQ`, reference with ambiguity codes as N |
| SVs | truvari 4.3.1, mafft aligner | truvari 5.4.0 (this run: the biocontainers image `5.4.0--pyhdfd78af_0`; thoth's own image pins the same dependency versions), refine with its default aligner |
| truth | v4.2.1; SVs V0.018 | v4.2.1; SVs V0.018 and v5.0q |

The container digests of the thoth run are in `evidence.json` under
`ont.provenance`, with the sha256 of its metrics table.

### Reproducing this section

1. Download ONT's calls for the benchmarked flow cells (URLs above).
2. Add the V0.018 draft to `resources/manifest.tsv` and a benchmark entry with
   the same `truvari:` block as v5.0q in `config/config.yaml`:

   ```
   HG002	GRCh38	V0.018	sv_vcf	truth/HG002/GRCh38/V0.018/GRCh38_HG2-T2TQ100-V1.1.vcf.gz	yes	https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/GRCh38_HG2-T2TQ100-V1.1.vcf.gz	25e2a5f9dadbd69ab7321fe0068938d18ce45cdd70b6987f66622dd216c90f89
   HG002	GRCh38	V0.018	sv_bed	truth/HG002/GRCh38/V0.018/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed	no	https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/analysis/NIST_HG002_DraftBenchmark_defrabbV0.018-20240716/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed	1c5b600d448953d4350c23528c4d4a4ede833006827c3f348bd966382375c5cb
   ```

   ```yaml
     - truth: V0.018
       strat: v3.6
       strat_subset: all
       truvari:
         refine: true
         pick: ac
         refdist: 2000
         chunksize: 5000
   ```

3. Write a run file with one query per sample (caller `ont-<flow cell, lower
   case>` under `vcf_small`; for HG002 also `ont-sniffles-paw70337` under
   `vcf_sv`) and run thoth.
4. Rerun ONT's SV recipe with the commands in step 2 into a directory `<repro>`,
   so that `<repro>/bench/` holds bench's and refine's outputs.
5. Regenerate `evidence.json` (next section).

---

## Regenerating `evidence.json`

```sh
python docs/extract_evidence.py \
  --query <query> --caller <small-variant caller> \
  --runinfo <run>.runinfo.json \
  --ont-metrics <comparison run>.metrics.tsv \
  --ont-repro <repro> \
  --ont-runinfo <comparison run>.runinfo.json \
  <run>.metrics.tsv
```

`--ont-metrics` and `--ont-repro` go together and add the `ont` key; without
them the other keys are written exactly as before. `--ont-runinfo` is optional.
Review the JSON before committing.

`ont` holds:

- `small`: per sample and variant type, ONT's `published` hap.py summary and
  `thoth`'s (recall, precision, truth total, TP, FN, FP);
- `sv`: HG002 PAW70337 SV rows (`scored_by` ONT as published, ONT's recipe
  rerun, or thoth; `truth`; `method` bench or refine; recall, precision, truth
  total `base`, `tp`, `fn`, `fp`);
- `recount_split`: the refine truth total split described above;
- `provenance`: whitelisted checksums of the comparison run.

Rates in `ont` carry six digits, so a two-decimal percentage is rounded once,
from the source value (ONT's 98.8247 % refine precision would otherwise show as
98.83 %); the other keys carry five.

## The other evidence on the page

- `truth_sets`, `strata`: one DeepVariant HG002 call set scored by thoth against
  GIAB v4.2.1, v5.0q and CMRG v1.00, and broken down by GIAB v3.6
  stratifications (own-sample strata only).
- `sv_bench`, `sv_refine_bench`, `sv_refine`: Sniffles and cuteSV HG002 calls,
  truvari 5.4.0 bench and refine against v5.0q (and bench against CMRG).
- `seqsim_baseline`: truvari's sequence similarity for pairs of unrelated
  uniform-random sequences, by length (`../seqsim_baseline.py` describes how it
  is computed).
- `provenance`: checksums and container digests of that run.
