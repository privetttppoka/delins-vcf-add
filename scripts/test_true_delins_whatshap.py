#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path

import pysam

from true_delins_common import (
    Candidate,
    GT,
    SnpVariant,
    add_info_fields,
    build_contig_name_map,
    candidate_record_id,
    gt_class,
    gt_has_alt,
    gt_label,
    index_vcf,
    is_biallelic_snp,
    is_pass_record,
    join_values,
    normalize_gt,
    select_sample,
    write_vcf_with_renamed_contigs,
)


mapq = 20
min_run_length = 2
min_vcf_alt_depth = 5

tsv_columns = [
    "candidate_id",
    "sample",
    "chrom",
    "start",
    "end",
    "n_snps",
    "positions",
    "ref_hap",
    "alt_hap",
    "genotypes",
    "phased",
    "phase_sets",
    "alt_haplotypes",
    "dp",
    "alt_depth",
    "gq",
    "decision",
    "reasons",
    "evidence",
]


@dataclass(frozen=True)
class PhaseSite:
    gt: GT
    phased: bool
    phase_set: int | str | None
    gt_class: str
    alt_haplotype: int | None
    dp: int | None
    alt_depth: int | None
    gq: int | None


@dataclass(frozen=True)
class PhaseResult:
    candidate: Candidate
    sites: tuple[PhaseSite, ...]
    decision: bool
    reasons: tuple[str, ...]
    evidence: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Use WhatsHap to find adjacent SNPs that form true delins/MNV events."
    )
    parser.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ with split SNP calls.")
    parser.add_argument("--bam", required=True, help="BAM/CRAM used for WhatsHap phasing.")
    parser.add_argument("--out-vcf", required=True, help="Output unphased VCF with confirmed SNP groups replaced by MNV.")
    parser.add_argument("--out-tsv", required=True, help="Output TSV log with all candidates and decisions.")
    parser.add_argument("--summary-json", help="Optional JSON file with runtime and summary counters.")
    parser.add_argument("--sample", help="Sample name. Default: first sample in VCF.")
    parser.add_argument("--reference", help="Optional reference FASTA for WhatsHap. If omitted, --no-reference is used.")
    return parser.parse_args()


def timed_command(command: list[str]) -> float:
    started = time.perf_counter()
    subprocess.run(command, check=True)
    return time.perf_counter() - started


def run_whatshap(args: argparse.Namespace, input_vcf: str, phased_vcf: str) -> float:
    command = [
        "whatshap",
        "phase",
        "--output",
        phased_vcf,
        "--mapping-quality",
        str(mapq),
        "--only-snvs",
        "--ignore-read-groups",
    ]
    command.extend(["--reference", args.reference] if args.reference else ["--no-reference"])
    if args.sample:
        command.extend(["--sample", args.sample])
    command.extend([input_vcf, args.bam])
    return timed_command(command)


def alt_depth_from_ad(ad: object) -> int | None:
    if ad is None:
        return None
    try:
        values = tuple(ad)
    except TypeError:
        return None
    return values[1] if len(values) > 1 else None


def alt_haplotype_from_gt(gt: GT) -> int | None:
    if gt is None:
        return None

    alt_indexes = [index for index, allele in enumerate(gt) if allele is not None and allele > 0]
    if len(alt_indexes) != 1:
        return None
    return alt_indexes[0] + 1


def parse_phase_record(
    record: pysam.VariantRecord,
    sample: str | None,
    contig_map: dict[str, str],
) -> tuple[SnpVariant, PhaseSite] | None:
    if not is_pass_record(record) or not is_biallelic_snp(record):
        return None

    gt = gq = dp = alt_depth = None
    phased = False
    phase_set = None
    if sample is not None:
        sample_data = record.samples[sample]
        gt = normalize_gt(sample_data.get("GT"))
        if not gt_has_alt(gt):
            return None
        gq = sample_data.get("GQ")
        dp = sample_data.get("DP")
        alt_depth = alt_depth_from_ad(sample_data.get("AD"))
        phased = bool(sample_data.phased)
        phase_set = sample_data.get("PS")

    variant = SnpVariant(
        chrom=contig_map.get(record.chrom, record.chrom),
        pos=record.pos,
        ref=record.ref.upper(),
        alt=record.alts[0].upper(),
        var_id=record.id,
        qual=record.qual,
        gt=gt,
        gq=gq,
    )
    site = PhaseSite(
        gt=gt,
        phased=phased,
        phase_set=phase_set,
        gt_class=gt_class(gt),
        alt_haplotype=alt_haplotype_from_gt(gt),
        dp=dp,
        alt_depth=alt_depth,
        gq=gq,
    )
    return variant, site


def decide_by_phase(sites: tuple[PhaseSite, ...]) -> tuple[bool, tuple[str, ...], str]:
    reasons: list[str] = []
    if any(site.alt_depth is None or site.alt_depth < min_vcf_alt_depth for site in sites):
        reasons.append("low_vcf_alt_depth")

    classes = {site.gt_class for site in sites}
    if classes == {"hom_alt"}:
        return not reasons, tuple(reasons), "all_hom_alt" if not reasons else "phase_failed"

    if classes != {"het"}:
        reasons.append("inconsistent_gt_class")

    het_sites = [site for site in sites if site.gt_class == "het"]
    if not het_sites:
        reasons.append("no_het_phase_information")
        return False, tuple(reasons), "phase_failed"
    if any(not site.phased for site in het_sites):
        reasons.append("unphased_het")

    phase_sets = {site.phase_set for site in het_sites}
    alt_haps = {site.alt_haplotype for site in het_sites}
    if len(phase_sets) != 1 or None in phase_sets:
        reasons.append("different_or_missing_phase_set")
    if len(alt_haps) != 1 or None in alt_haps:
        reasons.append("alt_on_different_haplotypes")

    if reasons:
        return False, tuple(reasons), "phase_failed"
    return True, tuple(), f"PS_{next(iter(phase_sets))}|ALT_HP_{next(iter(alt_haps))}"


def append_candidate_result(
    results: list[PhaseResult],
    run_variants: list[SnpVariant],
    run_sites: list[PhaseSite],
) -> None:
    if len(run_variants) < min_run_length:
        return

    candidate = Candidate(tuple(run_variants))
    sites = tuple(run_sites)
    decision, reasons, evidence = decide_by_phase(sites)
    results.append(PhaseResult(candidate, sites, decision, reasons, evidence))


def starts_new_run(variant: SnpVariant, run_variants: list[SnpVariant]) -> bool:
    return bool(run_variants) and (
        variant.chrom != run_variants[-1].chrom or variant.pos != run_variants[-1].pos + 1
    )


def analyze_phased_vcf(
    phased_vcf: str,
    sample_name: str | None,
    contig_map: dict[str, str],
) -> tuple[list[PhaseResult], str | None]:
    results: list[PhaseResult] = []
    run_variants: list[SnpVariant] = []
    run_sites: list[PhaseSite] = []

    with pysam.VariantFile(phased_vcf) as vcf:
        sample = select_sample(vcf, sample_name)
        for record in vcf:
            parsed = parse_phase_record(record, sample, contig_map)
            if parsed is None:
                continue

            variant, site = parsed
            if starts_new_run(variant, run_variants):
                append_candidate_result(results, run_variants, run_sites)
                run_variants = []
                run_sites = []

            run_variants.append(variant)
            run_sites.append(site)

    append_candidate_result(results, run_variants, run_sites)
    return results, sample


def result_row(result: PhaseResult, sample: str | None) -> dict[str, str | int]:
    candidate = result.candidate
    return {
        "candidate_id": candidate.candidate_id,
        "sample": sample or ".",
        "chrom": candidate.chrom,
        "start": candidate.start,
        "end": candidate.end,
        "n_snps": len(candidate.variants),
        "positions": join_values(candidate.positions),
        "ref_hap": candidate.ref_hap,
        "alt_hap": candidate.alt_hap,
        "genotypes": join_values(gt_label(site.gt) for site in result.sites),
        "phased": join_values("1" if site.phased else "0" for site in result.sites),
        "phase_sets": join_values(site.phase_set for site in result.sites),
        "alt_haplotypes": join_values(site.alt_haplotype for site in result.sites),
        "dp": join_values(site.dp for site in result.sites),
        "alt_depth": join_values(site.alt_depth for site in result.sites),
        "gq": join_values(site.gq for site in result.sites),
        "decision": "TRUE_DELINS" if result.decision else "NO",
        "reasons": "." if result.decision else ";".join(result.reasons),
        "evidence": result.evidence,
    }


def write_tsv(results: list[PhaseResult], sample: str | None, path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=tsv_columns, delimiter="\t")
        writer.writeheader()
        for result in results:
            writer.writerow(result_row(result, sample))


def add_mnv_header(header: pysam.VariantHeader) -> None:
    add_info_fields(
        header,
        {
            "WHATSHAP_TRUE_DELINS": ("0", "Flag", "MNV record added by WhatsHap true delins detection."),
        },
    )


def make_mnv_record(header: pysam.VariantHeader, result: PhaseResult, sample: str | None) -> pysam.VariantRecord:
    candidate = result.candidate
    quals = [variant.qual for variant in candidate.variants if variant.qual is not None]
    gqs = [variant.gq for variant in candidate.variants if variant.gq is not None]
    gt = result.sites[0].gt
    if gt is not None and None not in gt:
        gt = tuple(sorted(gt))

    record = header.new_record(
        contig=candidate.chrom,
        start=candidate.start - 1,
        stop=candidate.end,
        id=candidate_record_id(candidate),
        qual=min(quals) if quals else None,
        alleles=(candidate.ref_hap, candidate.alt_hap),
    )
    record.filter.add("PASS")
    record.info["WHATSHAP_TRUE_DELINS"] = True

    if sample is not None and sample in header.samples:
        if "GT" in header.formats and gt is not None:
            record.samples[sample]["GT"] = gt
            record.samples[sample].phased = False
        if "GQ" in header.formats and gqs:
            record.samples[sample]["GQ"] = min(gqs)
    return record


def write_corrected_vcf(args: argparse.Namespace, results: list[PhaseResult], sample: str | None) -> None:
    confirmed = [result for result in results if result.decision]
    by_start = {(result.candidate.chrom, result.candidate.start): result for result in confirmed}
    source_keys = {
        (result.candidate.chrom, variant.pos)
        for result in confirmed
        for variant in result.candidate.variants
    }

    with pysam.VariantFile(args.vcf) as in_vcf:
        header = in_vcf.header.copy()
        add_mnv_header(header)
        mode = "wz" if args.out_vcf.endswith(".gz") else "w"
        with pysam.VariantFile(args.out_vcf, mode, header=header) as out_vcf:
            for record in in_vcf:
                key = (record.chrom, record.pos)
                if key in by_start:
                    out_vcf.write(make_mnv_record(out_vcf.header, by_start[key], sample))
                if key in source_keys:
                    continue
                out_vcf.write(record)
    index_vcf(args.out_vcf)


def write_summary(
    args: argparse.Namespace,
    sample: str | None,
    results: list[PhaseResult],
    phase_sec: float,
    analyze_sec: float,
    total_sec: float,
) -> None:
    if not args.summary_json:
        return
    payload = {
        "vcf": args.vcf,
        "bam": args.bam,
        "sample": sample,
        "out_vcf": args.out_vcf,
        "n_candidates": len(results),
        "n_true_delins": sum(result.decision for result in results),
        "elapsed_sec": {
            "whatshap_phase": phase_sec,
            "analyze_and_write": analyze_sec,
            "total": total_sec,
        },
        "parameters": {
            "mapq": mapq,
            "min_run_length": min_run_length,
            "min_vcf_alt_depth": min_vcf_alt_depth,
            "only_snvs": True,
            "ignore_read_groups": True,
            "strict_adjacent_snps": True,
        },
    }
    with open(args.summary_json, "w", encoding="utf-8") as out:
        json.dump(payload, out, indent=2, sort_keys=True)
        out.write("\n")


def main() -> int:
    args = parse_args()
    started = time.perf_counter()

    out_dir = Path(args.out_vcf).resolve().parent
    with tempfile.TemporaryDirectory(prefix="whatshap_delins.", dir=out_dir) as tmpdir:
        with pysam.VariantFile(args.vcf) as vcf, pysam.AlignmentFile(args.bam) as bam:
            vcf_to_bam_contigs = build_contig_name_map(
                vcf.header.contigs,
                bam.references,
                source_label="VCF",
                target_label="BAM",
            )
        bam_to_vcf_contigs = {bam_contig: vcf_contig for vcf_contig, bam_contig in vcf_to_bam_contigs.items()}

        whatshap_vcf = str(Path(tmpdir, "input.bam_contigs.vcf.gz"))
        write_vcf_with_renamed_contigs(args.vcf, whatshap_vcf, vcf_to_bam_contigs)
        index_vcf(whatshap_vcf)

        phased_vcf = str(Path(tmpdir, "phased.vcf.gz"))
        phase_sec = run_whatshap(args, whatshap_vcf, phased_vcf)
        index_vcf(phased_vcf)

        analyze_started = time.perf_counter()
        results, sample = analyze_phased_vcf(phased_vcf, args.sample, bam_to_vcf_contigs)
        write_tsv(results, sample, args.out_tsv)
        write_corrected_vcf(args, results, sample)
        analyze_sec = time.perf_counter() - analyze_started

    total_sec = time.perf_counter() - started
    write_summary(args, sample, results, phase_sec, analyze_sec, total_sec)

    print(
        f"WhatsHap candidates {len(results)}; confirmed {sum(result.decision for result in results)} true delins; "
        f"phase {phase_sec:.2f}s; total {total_sec:.2f}s.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
