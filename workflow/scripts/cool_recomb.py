import pysam
import pandas as pd
from pathlib import Path

MAP_BED = pysam.TabixFile(
    "/private/groups/russelllab/cade/wolb-mix/p3p6/results/mappability/min_1_map.bed.gz",
    "r",
)


def get_annotation(chrom, start, end):
    if chrom == "NC_012416.1":
        file = "/private/groups/russelllab/cade/wolb-mix/recomb/sorted_wri.gtf.gz"
    else:
        file = "/private/groups/russelllab/cade/wolb-mix/recomb/sorted_wmel.gtf.gz"

    gtf = pysam.TabixFile(file, parser=pysam.asGTF())

    entries = gtf.fetch(chrom, start, end)
    prods = set()

    for e in entries:
        attrs = e.attributes.split(";")
        prod = [p for p in attrs if "product" in p]
        for p in prod:
            prods.add(p.split("=")[1])
    return prods


def percentage_overlap_with_bed(
    chrom,
    read_start,
    read_end,
    bedfile=MAP_BED,
):
    """
    Determine the percentage of a given read (specified by chrom, read_start, and read_end)
    that falls into a BED interval.

    Args:
    - bedfile: Path to the indexed BED file (with tabix index).
    - chrom: Chromosome name.
    - read_start: Read start position.
    - read_end: Read end position.

    Returns:
    - Percentage of the read that overlaps with the BED intervals.
    """

    total_overlap = 0

    for region in bedfile.fetch(chrom, read_start, read_end):
        bed_chrom, bed_start, bed_end = region.split("\t")[0:3]
        bed_start, bed_end = int(bed_start), int(bed_end)

        # Calculate overlap between read and this BED interval
        overlap_start = max(read_start, bed_start)
        overlap_end = min(read_end, bed_end)
        overlap = max(0, overlap_end - overlap_start)

        total_overlap += overlap

    read_length = read_end - read_start
    percentage_overlap = total_overlap / read_length

    return round(percentage_overlap, 3)


def get_timepoint(file: Path) -> int:
    t = int(file.name.split("-")[4].replace("T", ""))
    return t


def has_reads(files):
    has_reads = []
    for f in files:
        s = pysam.AlignmentFile(f, "rb")
        reads = list(s.fetch())
        s.close()
        # samp_name = f.name.split(".bam")[0]
        # wmel_cov = coverage.loc[coverage["sample_id"] == samp_name]["wmel"].values[0]
        # wri_cov = coverage.loc[coverage["sample_id"] == samp_name]["wri"].values[0]
        # print(f"{f.name=}, {len(reads)=}, {wmel_cov=}, {wri_cov=}")
        if len(reads) > 0:
            has_reads.append(f)
    return has_reads


def does_overlap(read1: pysam.AlignedSegment, read2: pysam.AlignedSegment) -> bool:
    """
    Check if two reads overlap.

    Args:
        read1 (pysam.AlignedSegment): First read.
        read2 (pysam.AlignedSegment): Second read.

    Returns:
        bool: True if the reads overlap, otherwise False.
    """
    if read1.reference_name == read2.reference_name:
        if (
            read1.reference_end > read2.reference_start
            and read1.reference_start < read2.reference_end
        ):
            return True
    return False


def find_overlaps(files: list[Path]):
    overlaps = []

    for i in range(len(files)):
        # first search file_i for overlaps.
        seen = set()  # keep track of read_names so we dont consider the pair.
        bam_i = pysam.AlignmentFile(files[i], "rb")
        for read_i in bam_i.fetch():
            if read_i.query_name in seen:
                continue
            seen.add(read_i.query_name)
            read_i_mate = bam_i.mate(read_i)
            o = [(files[i].name, read_i, read_i_mate)]
            # fetch reads in current (i) file that overlap with read_i
            for other_read_i in bam_i.fetch(
                read_i.reference_name, read_i.reference_start, read_i.reference_end
            ):
                if other_read_i.query_name in seen:
                    continue
                seen.add(other_read_i.query_name)
                other_read_i_mate = bam_i.mate(other_read_i)
                if does_overlap(read_i_mate, other_read_i_mate):
                    o.append((files[i].name, other_read_i, other_read_i_mate))
            # finished checking file_i, now we check other files.
            for j in range(i + 1, len(files)):
                bam_j = pysam.AlignmentFile(files[j], "rb")
                for read_j in bam_j.fetch(
                    read_i.reference_name, read_i.reference_start, read_i.reference_end
                ):
                    if read_j.query_name in seen:
                        continue
                    seen.add(read_j.query_name)
                    read_j_mate = bam_j.mate(read_j)
                    if does_overlap(read_i_mate, read_j_mate):
                        o.append((files[j].name, read_j, read_j_mate))
            if len(o) > 1:
                overlaps.append(o)
    return overlaps


def print_overlaps(overlaps):
    all_annotations = set()
    out = []
    for t in overlaps:
        print(f"found {len(t)} read pairs ")
        out.append(len(t))
        for p in t:
            print(p[0])
            r1_ref, r1_start, r1_end, r1_mapq, r1_qname, r1_strand = (
                p[1].reference_name,
                p[1].reference_start,
                p[1].reference_end,
                p[1].mapq,
                p[1].query_name,
                p[1].is_reverse,
            )
            r1_ratio = percentage_overlap_with_bed(r1_ref, r1_start, r1_end)
            r1_annotation = get_annotation(r1_ref, r1_start, r1_end)
            all_annotations.update(r1_annotation)
            r2_ref, r2_start, r2_end, r2_mapq, r2_strand = (
                p[2].reference_name,
                p[2].reference_start,
                p[2].reference_end,
                p[2].mapq,
                p[2].is_reverse,
            )
            r2_ratio = percentage_overlap_with_bed(r2_ref, r2_start, r2_end)
            r2_annotation = get_annotation(r2_ref, r2_start, r2_end)
            all_annotations.update(r2_annotation)
            print(
                r1_ref,
                r1_start,
                r1_end,
                r1_ratio,
                r1_mapq,
                r1_annotation,
                r1_strand,
                sep=",",
            )
            print(
                r2_ref,
                r2_start,
                r2_end,
                r2_ratio,
                r2_mapq,
                r2_annotation,
                r2_strand,
                sep=",",
            )
            # print("+++++++++++++++++++++++++++++++++")
        print("======================================")
    return all_annotations, out


# get Mix-S2_wMel-wRi-1_100
s2_files_100 = list(
    Path(
        "/private/groups/russelllab/cade/wolb-mix/p3p6/recomb_results/best_chimeras"
    ).glob("Mix-S2_wMel-wRi-1_100*.bam")
)
s2_a_files = sorted([f for f in s2_files_100 if f.name[-5] == "A"], key=get_timepoint)
s2_b_files = sorted([f for f in s2_files_100 if f.name[-5] == "B"], key=get_timepoint)
s2_c_files = sorted([f for f in s2_files_100 if f.name[-5] == "C"], key=get_timepoint)

s2_files_1k = list(
    Path(
        "/private/groups/russelllab/cade/wolb-mix/p3p6/recomb_results/best_chimeras"
    ).glob("Mix-S2_wMel-wRi-1_1K*.bam")
)
for f in s2_files_1k:
    # fix swapped business, kinda hacky but w/e
    f = str(f)
    if "Mix-S2_wMel-wRi-1_1K-T7-1805" in f:
        f = f.replace("Mix-S2_wMel-wRi-1_1K-T7-1805", "Mix-S2wMel-DOX-1_1-T7-1805")
    f = Path(f)

s2_a_files_1k = sorted([f for f in s2_files_1k if f.name[-5] == "A"], key=get_timepoint)
s2_b_files_1k = sorted([f for f in s2_files_1k if f.name[-5] == "B"], key=get_timepoint)
s2_c_files_1k = sorted([f for f in s2_files_1k if f.name[-5] == "C"], key=get_timepoint)


JW18_files_100 = list(
    Path(
        "/private/groups/russelllab/cade/wolb-mix/p3p6/recomb_results/best_chimeras"
    ).glob("Mix-JW18_wMel-wRi-1_100*.bam")
)
JW18_a_files = sorted(
    [f for f in JW18_files_100 if f.name[-5] == "A"], key=get_timepoint
)
JW18_b_files = sorted(
    [f for f in JW18_files_100 if f.name[-5] == "B"], key=get_timepoint
)
JW18_c_files = sorted(
    [f for f in JW18_files_100 if f.name[-5] == "C"], key=get_timepoint
)

JW18_files_1k = list(
    Path(
        "/private/groups/russelllab/cade/wolb-mix/p3p6/recomb_results/best_chimeras"
    ).glob("Mix-JW18_wMel-wRi-1_1K*.bam")
)
JW18_a_files_1k = sorted(
    [f for f in JW18_files_1k if f.name[-5] == "A"], key=get_timepoint
)
JW18_b_files_1k = sorted(
    [f for f in JW18_files_1k if f.name[-5] == "B"], key=get_timepoint
)
JW18_c_files_1k = sorted(
    [f for f in JW18_files_1k if f.name[-5] == "C"], key=get_timepoint
)


allgroups = {
    "s2_100a": s2_a_files,
    "s2_100b": s2_b_files,
    "s2_100c": s2_c_files,
    "s2_1ka": s2_a_files_1k,
    "s2_1kb": s2_b_files_1k,
    "s2_1kc": s2_c_files_1k,
    "jw_100a": JW18_a_files,
    "jw_100b": JW18_b_files,
    "jw_100c": JW18_c_files,
    "jw_1ka": JW18_a_files_1k,
    "jw_1kb": JW18_b_files_1k,
    "jw_1kc": JW18_c_files_1k,
}

d = {"JW18":[], "S2":[]}
for name, fg in allgroups.items():
    has = has_reads(fg)
    overlaps = find_overlaps(has)
    an, num = print_overlaps(overlaps)
    if "jw" in name:
        d["JW18"].extend(num)
    if "s2" in name:
        d["S2"].extend(num)
    print(name, an)

# with open("recombs.txt", "w") as f:
    
print("Infection", "Number of alignments", sep=",")

for i,j in zip(d["JW18"], d["S2"]):
    print(f"JW18,{i}\nS2,{j}")
        
    
