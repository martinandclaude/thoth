"""Build the planted-truth fixture set for the end-to-end test.

Deterministic (fixed seeds), so the expected counts in check_results.py are
facts about these files, not estimates. Writes everything under tests/work/:

  ref/chr1.fa                    60 kb random chr1, with one IUPAC code (B at
                                 7004), as GRCh38 has 94
  giab/                          truth sets + a stratification tarball, served
                                 to the pipeline through file:// URLs
  manifest.tsv                   points at giab/, same schema as the real one
  calls/                         one small-variant caller, one SV caller
  test.yaml                      run config: site values, test matrix, queries

Planted small variants (HG002, truth TESTSMALL):
  10 SNPs + a 3 bp deletion + a 3 bp insertion in the truth. The query calls 9
  SNPs right, one with the wrong ALT (FN, and FP at the allele level), adds one
  novel SNP (FP), calls the deletion and misses the insertion.
  Plus, in truth and query alike, a 4 bp deletion at 7000 whose REF ends on the
  reference's B and writes it N, as GIAB and callers do; hap.py must read a
  reference that says N there too, or vcfeval stops on the REF mismatch.
  -> SNP TP 9 / FN 1 / FP 2; INDEL TP 2 / FN 1.
  Plus a ./1 half-call at 59000, outside the confident region (0-58000), which
  prep must keep.

Queries: e2e-HG002 (public, both classes); e2e-internal (share: internal, the
same small-variant calls, which the share bundle must withhold); e2e-HG003 (SV
calls only, for a sample the manifest has small-variant truth for but no SV
truth, so its SV caller is skipped with a note).

Planted SVs (HG002):
  truth TESTSV (compared with refine, like v5.0q):
    a 300 bp deletion the caller writes as two adjacent 150 bp deletions on
    the same haplotype -> bench misses it, refine matches it
  truth TESTSV2 (plain bench, like CMRG):
    an 800 bp tandem duplication the truth writes as a sequence-resolved
    insertion and the caller as symbolic <DUP> -> matches only with
    dup-to-ins; plus a 30 kb deletion called as symbolic <DEL> -> sequence
    compared only with max-resolve >= 30000
"""
import hashlib
import io
import random
import tarfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
W = HERE / "work"
N = 60000
IUPAC = 7004  # the reference's one ambiguity base; inside the small-variant BED, clear of the SVs


def main():
    rng = random.Random(7)
    seq = "".join(rng.choice("ACGT") for _ in range(N))
    for d in ("ref", "giab/truth", "giab/strat", "calls"):
        (W / d).mkdir(parents=True, exist_ok=True)

    ref = seq[:IUPAC - 1] + "B" + seq[IUPAC:]
    (W / "ref/chr1.fa").write_text(">chr1\n" + "\n".join(ref[i:i + 60] for i in range(0, N, 60)) + "\n")
    (W / "ref/chr1.fa.fai").write_text(f"chr1\t{N}\t6\t60\t61\n")

    def b(pos):
        return seq[pos - 1]

    def alt(base):
        return {"A": "G", "C": "T", "G": "A", "T": "C"}[base]

    # ---------------------------------------------------------------- small variants
    hdr = ["##fileformat=VCFv4.2", f"##contig=<ID=chr1,length={N}>",
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
           '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002"]

    def rec(pos, ref, a, gq=50):
        return f"chr1\t{pos}\t.\t{ref}\t{a}\t50\tPASS\t.\tGT:GQ\t0/1:{gq}"

    snps = [500, 1500, 2500, 4000, 6000, 8000, 11000, 13000, 15000, 17000]
    truth = [(p, b(p), alt(b(p))) for p in snps]
    truth.append((9000, seq[8999:9003], b(9000)))            # 3 bp deletion
    truth.append((12000, b(12000), b(12000) + "ACG"))         # 3 bp insertion
    truth.append((IUPAC - 4, seq[IUPAC - 5:IUPAC - 1] + "N", b(IUPAC - 4)))  # deletion ending on the B
    truth.sort()
    query = []
    for p, r, a in truth:
        if p == 13000:
            query.append((p, r, next(x for x in "ACGT" if x not in (r, a))))   # wrong ALT
        elif p == 12000:
            continue                                                         # missed insertion
        else:
            query.append((p, r, a))
    query.append((5000, b(5000), alt(b(5000))))                               # novel FP
    query.sort()
    _write_vcf(W / "giab/truth/small.vcf", hdr, [rec(*t) for t in truth])
    (W / "giab/truth/small.bed").write_text("chr1\t0\t58000\n")
    halfcall = f"chr1\t59000\t.\t{b(59000)}\t{alt(b(59000))}\t50\tPASS\t.\tGT:GQ\t./1:30"
    (W / "calls/snvcaller.vcf").write_text(
        "\n".join(hdr + [rec(*t, gq=60 - i) for i, t in enumerate(query)] + [halfcall]) + "\n")

    # ---------------------------------------------------------------- SVs
    svh = ["##fileformat=VCFv4.2", f"##contig=<ID=chr1,length={N}>",
           '##ALT=<ID=DEL,Description="Deletion">', '##ALT=<ID=DUP,Description="Duplication">',
           '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">',
           '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">',
           '##INFO=<ID=END,Number=1,Type=Integer,Description="End">',
           '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
           "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002"]

    def dele(pos, n, label, gt="0|1"):
        return f"chr1\t{pos}\t{label}\t{seq[pos-1:pos+n]}\t{b(pos)}\t60\tPASS\tSVTYPE=DEL;SVLEN=-{n};END={pos+n}\tGT\t{gt}"

    def ins(pos, s, label):
        return f"chr1\t{pos}\t{label}\t{b(pos)}\t{b(pos)}{s}\t60\tPASS\tSVTYPE=INS;SVLEN={len(s)};END={pos}\tGT\t0|1"

    def sym(pos, kind, n, label):
        svlen = -n if kind == "DEL" else n
        return f"chr1\t{pos}\t{label}\t{b(pos)}\t<{kind}>\t60\tPASS\tSVTYPE={kind};SVLEN={svlen};END={pos+n}\tGT\t0|1"

    # TESTSV (refine): one 300 bp deletion, plus an ALT=* placeholder that prep_sv_truth must drop
    _write_vcf(W / "giab/truth/sv_refine.vcf", svh,
               [f"chr1\t3000\tt_altstar\t{b(3000)}\t*\t60\tPASS\t.\tGT\t0|1", dele(10000, 300, "t_del300")])
    (W / "giab/truth/sv_refine.bed").write_text("chr1\t8000\t12000\n")
    # TESTSV2 (plain bench): tandem dup as resolved INS, and a 30 kb deletion
    P1, n1, P2, n2 = 20000, 800, 25000, 30000
    _write_vcf(W / "giab/truth/sv_plain.vcf", svh,
               [ins(P1, seq[P1:P1 + n1], "t_dup_as_ins"), dele(P2, n2, "t_del30kb")])
    (W / "giab/truth/sv_plain.bed").write_text(f"chr1\t19000\t{N}\n")
    # the SV caller: written UNSORTED and uncompressed on purpose -- prep must cope
    calls = [sym(P2, "DEL", n2, "q_sym_del30kb"), dele(10150, 150, "q_del150b"),
             sym(P1, "DUP", n1, "q_sym_dup"), dele(10000, 150, "q_del150a")]
    (W / "calls/svcaller.vcf").write_text("\n".join(svh + calls) + "\n")

    # ---------------------------------------------------------------- stratifications
    # Shaped like a GIAB v3.6 <ref>@all tarball: one top-level directory, a
    # stratification TSV with paths relative to it, and GenomeSpecific strata
    # for two samples so the own-sample filter has something to drop.
    strata = {
        "LowComplexity/GRCh38_left.bed": "chr1\t0\t30000\n",
        "LowComplexity/GRCh38_right.bed": "chr1\t30000\t60000\n",
        "GenomeSpecific/GRCh38_HG002_v4.2.1_testregion.bed": "chr1\t0\t10000\n",
        "GenomeSpecific/GRCh38_HG007_v4.2.1_testregion.bed": "chr1\t0\t10000\n",
    }
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        def add(name, data):
            info = tarfile.TarInfo(f"GRCh38@all/{name}")
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))
        for name, bed in strata.items():
            add(name, bed.encode())
        tsv = "".join(f"{Path(n).stem.replace('.bed', '').replace('GRCh38_', '')}\t{n}\n" for n in strata)
        add("GRCh38-all-stratifications.tsv", tsv.encode())
    (W / "giab/strat/strat.tar.gz").write_bytes(buf.getvalue())

    # ---------------------------------------------------------------- manifest + config
    g = f"file://{W / 'giab'}"
    rows = [
        ("HG002", "GRCh38", "TESTSMALL", "smvar_vcf", "truth/small.vcf.gz", "yes", f"{g}/truth/small.vcf.gz"),
        ("HG002", "GRCh38", "TESTSMALL", "smvar_bed", "truth/small.bed", "no", f"{g}/truth/small.bed"),
        ("HG003", "GRCh38", "TESTSMALL", "smvar_vcf", "truth/HG003_small.vcf.gz", "yes", f"{g}/truth/small.vcf.gz"),
        ("HG003", "GRCh38", "TESTSMALL", "smvar_bed", "truth/HG003_small.bed", "no", f"{g}/truth/small.bed"),
        ("HG002", "GRCh38", "TESTSV", "sv_vcf", "truth/sv_refine.vcf.gz", "yes", f"{g}/truth/sv_refine.vcf.gz"),
        ("HG002", "GRCh38", "TESTSV", "sv_bed", "truth/sv_refine.bed", "no", f"{g}/truth/sv_refine.bed"),
        ("HG002", "GRCh38", "TESTSV2", "sv_vcf", "truth/sv_plain.vcf.gz", "yes", f"{g}/truth/sv_plain.vcf.gz"),
        ("HG002", "GRCh38", "TESTSV2", "sv_bed", "truth/sv_plain.bed", "no", f"{g}/truth/sv_plain.bed"),
        ("-", "GRCh38", "vT", "strat_tar", "strat/GRCh38/vT/genome-stratifications-GRCh38-vT.tar.gz", "no", f"{g}/strat/strat.tar.gz"),
    ]
    # Pinned where the file exists now. The truth VCFs are bgzipped later, by
    # run.sh, so they stay unpinned: the run must report that, not fail on it.
    def sha256(url):
        path = Path(url.removeprefix("file://"))
        return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else ""

    (W / "manifest.tsv").write_text(
        "sample\treference\trelease\tkind\tdest\tindex\turl\tsha256\n"
        + "".join("\t".join((*r, sha256(r[-1]))) + "\n" for r in rows))

    (W / "test.yaml").write_text(f"""# Written by tests/make_fixtures.py -- do not edit.
run_id: e2e
results_dir: {W}/results
resources_dir: {W}/resources
containers_dir: {W}/sifs
manifest: {W}/manifest.tsv
references:
  GRCh38:
    fasta: {W}/ref/chr1.fa
small_variants:
  regions: chr1
  happy_mem_mb: 4000
benchmarks:
  - truth: TESTSMALL
    strat: vT
    strat_subset: all
  - truth: TESTSV
    strat: vT
    strat_subset: all
    truvari:
      refine: true
      pick: ac
      refdist: 2000
      chunksize: 5000
  - truth: TESTSV2
    strat: vT
    strat_subset: all
queries:
  - name: e2e-HG002
    sample: HG002
    reference: GRCh38
    share: public
    vcf_small:
      snvcaller: {W}/calls/snvcaller.vcf
    vcf_sv:
      svcaller: {W}/calls/svcaller.vcf
  - name: e2e-internal
    sample: HG002
    reference: GRCh38
    share: internal
    vcf_small:
      withheld: {W}/calls/snvcaller.vcf
  - name: e2e-HG003
    sample: HG003
    reference: GRCh38
    vcf_sv:
      svcaller: {W}/calls/svcaller.vcf
""")
    print(f"fixtures written to {W}")


def _write_vcf(path, header, records):
    # Left uncompressed: tests/run.sh bgzips and tabix-indexes it inside the hts
    # image (tabix needs BGZF), so the host needs no compression tools.
    path.write_text("\n".join(header + records) + "\n")


if __name__ == "__main__":
    main()
