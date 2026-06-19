"""Shared helpers for true-delins/MNV detection scripts."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from typing import Iterable

import pysam


GT = tuple[int | None, ...] | None


@dataclass(frozen=True)
class SnpVariant:
    chrom: str
    pos: int
    ref: str
    alt: str
    var_id: str | None
    qual: float | None
    gt: GT
    gq: int | None


@dataclass(frozen=True)
class Candidate:
    variants: tuple[SnpVariant, ...]

    @property
    def chrom(self) -> str:
        return self.variants[0].chrom

    @property
    def start(self) -> int:
        return self.variants[0].pos

    @property
    def end(self) -> int:
        return self.variants[-1].pos

    @property
    def positions(self) -> tuple[int, ...]:
        return tuple(variant.pos for variant in self.variants)

    @property
    def ref_hap(self) -> str:
        return "".join(variant.ref for variant in self.variants)

    @property
    def alt_hap(self) -> str:
        return "".join(variant.alt for variant in self.variants)

    @property
    def candidate_id(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}:{self.ref_hap}>{self.alt_hap}"


def normalize_gt(gt: GT) -> GT:
    return None if gt is None else tuple(gt)


def gt_has_alt(gt: GT) -> bool:
    return bool(gt) and any(allele is not None and allele > 0 for allele in gt)


def gt_label(gt: GT) -> str:
    if gt is None:
        return "."
    return "/".join("." if allele is None else str(allele) for allele in gt)


def gt_class(gt: GT) -> str:
    if gt is None or any(allele is None for allele in gt):
        return "unknown"
    alt_count = sum(1 for allele in gt if allele and allele > 0)
    if alt_count == 0:
        return "hom_ref"
    if alt_count == len(gt):
        return "hom_alt"
    return "het"


def is_pass_record(record: pysam.VariantRecord) -> bool:
    filters = tuple(record.filter.keys())
    return not filters or filters == ("PASS",)


def is_biallelic_snp(record: pysam.VariantRecord) -> bool:
    if len(record.ref) != 1 or not record.alts or len(record.alts) != 1:
        return False
    return len(record.alts[0]) == 1 and record.ref.upper() in "ACGTN" and record.alts[0].upper() in "ACGTN"


def select_sample(vcf: pysam.VariantFile, requested: str | None) -> str | None:
    samples = list(vcf.header.samples)
    if not samples:
        return None
    if requested is None:
        return samples[0]
    if requested not in samples:
        raise SystemExit(f"Sample {requested!r} is not present in VCF. Available: {', '.join(samples)}")
    return requested


def snp_from_record(record: pysam.VariantRecord, sample: str | None) -> SnpVariant:
    gt = None
    gq = None
    if sample is not None:
        sample_data = record.samples[sample]
        gt = normalize_gt(sample_data.get("GT"))
        gq = sample_data.get("GQ")
    return SnpVariant(
        chrom=record.chrom,
        pos=record.pos,
        ref=record.ref.upper(),
        alt=record.alts[0].upper(),
        var_id=record.id,
        qual=record.qual,
        gt=gt,
        gq=gq,
    )


def load_snp_variants(vcf_path: str, sample_name: str | None) -> tuple[list[SnpVariant], str | None]:
    variants: list[SnpVariant] = []
    with pysam.VariantFile(vcf_path) as vcf:
        sample = select_sample(vcf, sample_name)
        for record in vcf:
            if not is_pass_record(record):
                continue
            if not is_biallelic_snp(record):
                continue

            variant = snp_from_record(record, sample)
            if sample is not None and not gt_has_alt(variant.gt):
                continue
            variants.append(variant)
    return variants, sample


def build_candidates(variants: Iterable[SnpVariant], min_run_length: int) -> list[Candidate]:
    candidates: list[Candidate] = []
    run: list[SnpVariant] = []

    for variant in variants:
        if run and variant.chrom == run[-1].chrom and variant.pos == run[-1].pos + 1:
            run.append(variant)
            continue
        if len(run) >= min_run_length:
            candidates.append(Candidate(tuple(run)))
        run = [variant]

    if len(run) >= min_run_length:
        candidates.append(Candidate(tuple(run)))
    return candidates


def fmt_float(value: float) -> str:
    return f"{value:.6g}"


def join_values(values: Iterable[object | None]) -> str:
    return ",".join("." if value is None else str(value) for value in values)


def candidate_record_id(candidate: Candidate) -> str | None:
    ids = [variant.var_id for variant in candidate.variants if variant.var_id not in {None, "."}]
    return ";".join(ids) if ids else None


def add_info_fields(header: pysam.VariantHeader, fields: dict[str, tuple[str, str, str]]) -> None:
    for key, (number, field_type, description) in fields.items():
        if key not in header.info:
            header.info.add(key, number, field_type, description)


def contig_aliases(contig: str) -> tuple[str, ...]:
    aliases = [contig]
    if contig.startswith("chr"):
        aliases.append(contig[3:])
        if contig == "chrM":
            aliases.extend(["MT", "M"])
    else:
        aliases.append(f"chr{contig}")
        if contig in {"MT", "M"}:
            aliases.append("chrM")
    return tuple(dict.fromkeys(aliases))


def build_contig_name_map(
    source_contigs: Iterable[str],
    target_contigs: Iterable[str],
    source_label: str = "source",
    target_label: str = "target",
) -> dict[str, str]:
    target_set = set(target_contigs)
    contig_map: dict[str, str] = {}
    missing: list[str] = []

    for contig in dict.fromkeys(source_contigs):
        for alias in contig_aliases(contig):
            if alias in target_set:
                contig_map[contig] = alias
                break
        else:
            missing.append(contig)

    if missing:
        shown = ", ".join(missing[:10])
        raise SystemExit(
            f"{source_label} contigs are absent from {target_label}: {shown}. "
            "This usually means the files use different chromosome naming styles."
        )
    return contig_map


def rename_contig_header_line(line: str, contig_map: dict[str, str]) -> str:
    """Rename contig ID inside a VCF ##contig header line."""
    prefix = "##contig=<ID="
    if not line.startswith(prefix):
        return line

    end = line.find(",", len(prefix))
    if end == -1:
        end = line.find(">", len(prefix))
    if end == -1:
        return line

    contig = line[len(prefix):end]
    renamed = contig_map.get(contig)
    if renamed is None:
        return line
    return f"{line[:len(prefix)]}{renamed}{line[end:]}"


def write_vcf_with_renamed_contigs(input_vcf: str, output_vcf: str, contig_map: dict[str, str]) -> None:
    """Write a temporary VCF with contig names changed to match another file."""
    opener = gzip.open if input_vcf.endswith(".gz") else open
    writer = pysam.BGZFile if output_vcf.endswith(".gz") else open

    with opener(input_vcf, "rt") as inp, writer(output_vcf, "w") as out:
        for line in inp:
            if line.startswith("##contig=<ID="):
                line = rename_contig_header_line(line, contig_map)
            elif not line.startswith("#"):
                fields = line.split("\t", 1)
                fields[0] = contig_map.get(fields[0], fields[0])
                line = "\t".join(fields)
            out.write(line.encode() if output_vcf.endswith(".gz") else line)


def index_vcf(path: str, no_index: bool = False) -> None:
    if path.endswith(".gz") and not no_index:
        pysam.tabix_index(path, preset="vcf", force=True)
