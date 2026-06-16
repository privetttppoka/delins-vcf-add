#!/usr/bin/env python3
"""Find true adjacent-SNP delins/MNV events directly from BAM read haplotypes."""

from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from collections import Counter
from dataclasses import dataclass

import pysam

from true_delins_common import (
    Candidate,
    add_info_fields,
    build_candidates,
    candidate_record_id,
    fmt_float,
    gt_class,
    gt_label,
    index_vcf,
    join_values,
    load_snp_variants,
)


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
    "decision",
    "reasons",
    "n_reads_seen",
    "n_reads_filtered",
    "n_partial_reads",
    "n_full_reads",
    "n_ref_hap",
    "n_alt_hap",
    "n_mixed_hap",
    "alt_fraction",
    "mixed_fraction",
    "observed_haplotypes",
    "min_qual",
    "min_gq",
    "elapsed_sec",
]


min_run_length = 2  # Minimum number of adjacent SNPs required to form a candidate.
min_mapq = 20  # Minimum read mapping quality used for BAM evidence.
min_baseq = 20  # Minimum base quality at every candidate SNP position.
min_full_reads = 5  # Minimum reads spanning all SNP positions in the candidate.
min_alt_hap_reads = 5  # Minimum full-spanning reads carrying the complete ALT haplotype.
het_min_af = 0.25  # Minimum ALT haplotype fraction expected for heterozygous candidates.
het_max_af = 0.75  # Maximum ALT haplotype fraction expected for heterozygous candidates.
hom_alt_min_af = 0.80  # Minimum ALT haplotype fraction expected for homozygous ALT candidates.
max_mixed_fraction = 0.20  # Maximum allowed fraction of reads carrying neither REF nor ALT haplotype.


@dataclass
class ReadHaplotypeStats:
    candidate: Candidate
    sample: str | None
    hap_counts: Counter[str]
    n_reads_seen: int = 0
    n_reads_filtered: int = 0
    n_partial_reads: int = 0
    elapsed_sec: float = 0.0
    decision: bool = False
    reasons: list[str] | None = None

    @property
    def n_ref_hap(self) -> int:
        return self.hap_counts[self.candidate.ref_hap]

    @property
    def n_alt_hap(self) -> int:
        return self.hap_counts[self.candidate.alt_hap]

    @property
    def n_full_reads(self) -> int:
        return sum(self.hap_counts.values())

    @property
    def n_mixed_hap(self) -> int:
        return self.n_full_reads - self.n_ref_hap - self.n_alt_hap

    @property
    def alt_fraction(self) -> float:
        return self.n_alt_hap / self.n_full_reads if self.n_full_reads else 0.0

    @property
    def mixed_fraction(self) -> float:
        return self.n_mixed_hap / self.n_full_reads if self.n_full_reads else 0.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Find read-backed true delins/MNV events from adjacent SNPs.")
    parser.add_argument("--vcf", required=True, help="Input VCF/VCF.GZ with SNP calls.")
    parser.add_argument("--bam", required=True, help="Coordinate-sorted BAM/CRAM with index.")
    parser.add_argument("--out-tsv", required=True, help="Output TSV with candidate haplotype statistics.")
    parser.add_argument("--out-vcf", required=True, help="Output unphased VCF with confirmed SNP groups replaced by one MNV.")
    parser.add_argument("--summary-json", help="Optional JSON file with runtime and summary counters.")
    parser.add_argument("--sample", help="Sample name. Default: first sample in VCF.")
    return parser.parse_args()


def read_is_usable(read: pysam.AlignedSegment) -> bool:
    return not (
        read.is_unmapped
        or read.is_secondary
        or read.is_supplementary
        or read.is_duplicate
        or read.is_qcfail
        or read.mapping_quality < min_mapq
    )


def read_haplotype(read: pysam.AlignedSegment, positions0: tuple[int, ...]) -> str | None:
    if read.query_sequence is None:
        return None

    wanted = set(positions0)
    query_by_ref = {
        ref_pos: query_pos
        for query_pos, ref_pos in read.get_aligned_pairs(matches_only=False)
        if query_pos is not None and ref_pos in wanted
    }

    bases: list[str] = []
    for ref_pos in positions0:
        query_pos = query_by_ref.get(ref_pos)
        if query_pos is None:
            return None
        if read.query_qualities is not None and read.query_qualities[query_pos] < min_baseq:
            return None
        bases.append(read.query_sequence[query_pos].upper())
    return "".join(bases)


def count_candidate_haplotypes(
    bam: pysam.AlignmentFile,
    candidate: Candidate,
) -> tuple[Counter[str], int, int, int]:
    positions0 = tuple(pos - 1 for pos in candidate.positions)
    hap_counts: Counter[str] = Counter()
    seen = filtered = partial = 0

    for read in bam.fetch(candidate.chrom, candidate.start - 1, candidate.end):
        seen += 1
        if not read_is_usable(read):
            filtered += 1
            continue

        hap = read_haplotype(read, positions0)
        if hap is None:
            partial += 1
            continue

        hap_counts[hap] += 1

    return hap_counts, seen, filtered, partial


def decide(stats: ReadHaplotypeStats) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    gt_labels = {gt_label(variant.gt) for variant in stats.candidate.variants}

    if len(gt_labels) > 1:
        reasons.append("inconsistent_gt")
    if stats.n_full_reads < min_full_reads:
        reasons.append("low_full_read_depth")
    if stats.n_alt_hap < min_alt_hap_reads:
        reasons.append("low_alt_haplotype_support")
    if stats.mixed_fraction > max_mixed_fraction:
        reasons.append("high_mixed_haplotype_fraction")

    if len(gt_labels) == 1:
        klass = gt_class(stats.candidate.variants[0].gt)
        if klass == "het" and not (het_min_af <= stats.alt_fraction <= het_max_af):
            reasons.append("het_alt_fraction_out_of_range")
        elif klass == "hom_alt" and stats.alt_fraction < hom_alt_min_af:
            reasons.append("hom_alt_fraction_too_low")
        elif klass in {"hom_ref", "unknown"}:
            reasons.append(f"unsupported_gt_{klass}")

    return not reasons, reasons


def analyze_candidate(
    bam: pysam.AlignmentFile,
    candidate: Candidate,
    sample: str | None,
) -> ReadHaplotypeStats:
    started = time.perf_counter()
    hap_counts, seen, filtered, partial = count_candidate_haplotypes(bam, candidate)
    stats = ReadHaplotypeStats(
        candidate=candidate,
        sample=sample,
        hap_counts=hap_counts,
        n_reads_seen=seen,
        n_reads_filtered=filtered,
        n_partial_reads=partial,
        elapsed_sec=time.perf_counter() - started,
    )
    stats.decision, stats.reasons = decide(stats)
    return stats


def row(stats: ReadHaplotypeStats) -> dict[str, str | int | float]:
    candidate = stats.candidate
    observed = ";".join(f"{hap}:{count}" for hap, count in stats.hap_counts.most_common())
    quals = [variant.qual for variant in candidate.variants if variant.qual is not None]
    gqs = [variant.gq for variant in candidate.variants if variant.gq is not None]
    return {
        "candidate_id": candidate.candidate_id,
        "sample": stats.sample or ".",
        "chrom": candidate.chrom,
        "start": candidate.start,
        "end": candidate.end,
        "n_snps": len(candidate.variants),
        "positions": join_values(candidate.positions),
        "ref_hap": candidate.ref_hap,
        "alt_hap": candidate.alt_hap,
        "genotypes": join_values(gt_label(variant.gt) for variant in candidate.variants),
        "decision": "TRUE_DELINS" if stats.decision else "NO",
        "reasons": "." if stats.decision else ";".join(stats.reasons or []),
        "n_reads_seen": stats.n_reads_seen,
        "n_reads_filtered": stats.n_reads_filtered,
        "n_partial_reads": stats.n_partial_reads,
        "n_full_reads": stats.n_full_reads,
        "n_ref_hap": stats.n_ref_hap,
        "n_alt_hap": stats.n_alt_hap,
        "n_mixed_hap": stats.n_mixed_hap,
        "alt_fraction": fmt_float(stats.alt_fraction),
        "mixed_fraction": fmt_float(stats.mixed_fraction),
        "observed_haplotypes": observed or ".",
        "min_qual": fmt_float(min(quals)) if quals else ".",
        "min_gq": min(gqs) if gqs else ".",
        "elapsed_sec": fmt_float(stats.elapsed_sec),
    }


def write_tsv(stats: list[ReadHaplotypeStats], path: str) -> None:
    with open(path, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=tsv_columns, delimiter="\t")
        writer.writeheader()
        for item in stats:
            writer.writerow(row(item))


def add_mnv_header(header: pysam.VariantHeader) -> None:
    add_info_fields(
        header,
        {
            "TRUE_DELINS": ("0", "Flag", "MNV record added by true delins detection."),
        },
    )


def make_mnv_record(header: pysam.VariantHeader, stats: ReadHaplotypeStats, sample: str | None) -> pysam.VariantRecord:
    candidate = stats.candidate
    quals = [variant.qual for variant in candidate.variants if variant.qual is not None]
    gqs = [variant.gq for variant in candidate.variants if variant.gq is not None]
    record = header.new_record(
        contig=candidate.chrom,
        start=candidate.start - 1,
        stop=candidate.end,
        id=candidate_record_id(candidate),
        qual=min(quals) if quals else None,
        alleles=(candidate.ref_hap, candidate.alt_hap),
    )
    record.filter.add("PASS")
    record.info["TRUE_DELINS"] = True

    if sample is not None and sample in header.samples:
        if "GT" in header.formats and candidate.variants[0].gt is not None:
            record.samples[sample]["GT"] = candidate.variants[0].gt
        if "DP" in header.formats:
            record.samples[sample]["DP"] = stats.n_full_reads
        if "AD" in header.formats:
            record.samples[sample]["AD"] = (stats.n_ref_hap, stats.n_alt_hap)
        if "VAF" in header.formats:
            record.samples[sample]["VAF"] = (stats.alt_fraction,)
        if "GQ" in header.formats and gqs:
            record.samples[sample]["GQ"] = min(gqs)
    return record


def write_mnv_vcf(args: argparse.Namespace, stats: list[ReadHaplotypeStats], sample: str | None) -> None:
    confirmed = [item for item in stats if item.decision]
    by_start = {(item.candidate.chrom, item.candidate.start): item for item in confirmed}
    source_keys = {
        (item.candidate.chrom, variant.pos)
        for item in confirmed
        for variant in item.candidate.variants
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


def write_summary(args: argparse.Namespace, stats: list[ReadHaplotypeStats], sample: str | None, elapsed_sec: float) -> None:
    if not args.summary_json:
        return
    payload = {
        "vcf": args.vcf,
        "bam": args.bam,
        "sample": sample,
        "n_candidates": len(stats),
        "n_true_delins": sum(item.decision for item in stats),
        "elapsed_sec": elapsed_sec,
        "parameters": {
            "min_run_length": min_run_length,
            "strict_adjacent_snps": True,
            "min_mapq": min_mapq,
            "min_baseq": min_baseq,
            "min_full_reads": min_full_reads,
            "min_alt_hap_reads": min_alt_hap_reads,
            "het_min_af": het_min_af,
            "het_max_af": het_max_af,
            "hom_alt_min_af": hom_alt_min_af,
            "max_mixed_fraction": max_mixed_fraction,
        },
    }
    with open(args.summary_json, "w", encoding="utf-8") as out:
        json.dump(payload, out, indent=2, sort_keys=True)
        out.write("\n")


def main() -> int:
    args = parse_args()

    started = time.perf_counter()
    variants, sample = load_snp_variants(args.vcf, args.sample)
    candidates = build_candidates(variants, min_run_length)

    with pysam.AlignmentFile(args.bam) as bam:
        stats = [analyze_candidate(bam, candidate, sample) for candidate in candidates]

    write_tsv(stats, args.out_tsv)
    write_mnv_vcf(args, stats, sample)

    elapsed = time.perf_counter() - started
    write_summary(args, stats, sample, elapsed)
    print(
        f"Processed {len(candidates)} candidates; confirmed {sum(item.decision for item in stats)} true delins; "
        f"elapsed {elapsed:.2f}s.",
        file=sys.stderr,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
