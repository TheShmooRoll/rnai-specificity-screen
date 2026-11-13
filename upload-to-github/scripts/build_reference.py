#!/usr/bin/env python3
"""Prepare the versioned Bristol N2 reference used by the browser."""

import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import shutil
from urllib.parse import unquote
from urllib.request import urlretrieve
from zipfile import ZipFile


ACCESSION = "GCF_000002985.6"
ANNOTATION = "WormBase WS298"
DOWNLOAD = (
    f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{ACCESSION}/download"
    "?include_annotation_type=GENOME_GFF,GENOME_FASTA,RNA_FASTA"
)
REPORT_URL = f"https://www.ncbi.nlm.nih.gov/datasets/genome/{ACCESSION}/"
# Pin the uncompressed input files because the ZIP wrapper can be regenerated.
SOURCE_HASHES = {
    "GCF_000002985.6_WBcel235_genomic.fna":
        "739d1fb774ef9f5a2de12fd16de3e559201d00d69502b7013d820584d9369220",
    "genomic.gff":
        "b6296b398ce09b944696d7519d5013886e3ee2172760a6a6849774afc116334d",
    "rna.fna":
        "d1b4cbd014f5efcb3e09cccd1e1c1f10735e500f880912ecd81872ccedc23a9c",
}
BASE_CODES = bytes(b"ACGT".find(bytes([value])) % 5 for value in range(256))
COMPLEMENT = bytes.maketrans(b"ACGT", b"TGCA")
FEATURES_WITHOUT_RNA = {"region", "sequence_feature", "exon", "CDS"}


def fasta(source):
    name = None
    chunks = []
    for raw in source:
        if raw.startswith(b">"):
            if name is not None:
                yield name, b"".join(chunks).upper()
            name = raw[1:].split()[0].decode()
            chunks = []
        else:
            chunks.append(raw.strip())
    if name is not None:
        yield name, b"".join(chunks).upper()


def attributes(field):
    return dict(item.split("=", 1) for item in field.split(";") if "=" in item)


def stable_gene_id(attrs):
    refs = attrs.get("Dbxref", "").split(",")
    for prefix in ("WormBase:WBGene", "GeneID:"):
        for ref in refs:
            if ref.startswith(prefix):
                return ref.removeprefix("WormBase:")
    return attrs["ID"]


def parse_annotation(source):
    genes = {}
    features = {}
    exons = defaultdict(set)
    coding = defaultdict(set)
    chromosomes = {}
    for raw in source:
        if raw.startswith(b"#"):
            continue
        fields = raw.decode().rstrip("\n").split("\t")
        if len(fields) != 9:
            continue
        accession, _, kind, start, end, _, strand, _, attr_text = fields
        attrs = attributes(attr_text)
        start, end = int(start) - 1, int(end)
        if kind == "region":
            chromosomes[accession] = attrs.get("chromosome", attrs.get("Name", accession))
        elif kind in ("gene", "pseudogene"):
            symbol = unquote(attrs.get("gene", attrs.get("Name", attrs["ID"])))
            locus = attrs.get("locus_tag", "").removeprefix("CELE_")
            genes[attrs["ID"]] = (
                stable_gene_id(attrs), symbol, locus, attrs.get("gene_biotype", kind),
                accession, start, end, strand
            )
        elif kind in ("exon", "CDS"):
            records = exons if kind == "exon" else coding
            for parent in attrs.get("Parent", "").split(","):
                records[parent].add((start, end))
        elif kind not in FEATURES_WITHOUT_RNA and attrs.get("Parent") and "ID" in attrs:
            accession_id = attrs.get("transcript_id", "")
            name = unquote(attrs.get("standard_name", accession_id or attrs.get("Name", attrs["ID"])))
            features[attrs["ID"]] = (
                attrs["Parent"], name, accession_id, kind, accession, start, end, strand
            )
    return genes, features, exons, coding, chromosomes


def owner_of(feature, features, genes):
    while feature not in genes:
        if feature not in features:
            return None
        feature = features[feature][0]
    return feature


def coding_coordinates(exon_spans, coding_spans, strand):
    output = []
    offset = 0
    for left, right in exon_spans:
        for cds_left, cds_right in coding_spans:
            start, end = max(left, cds_left), min(right, cds_right)
            if end > start:
                if strand == "+":
                    output.append((offset + start - left, offset + end - left))
                else:
                    output.append((offset + right - end, offset + right - start))
        offset += right - left
    return [value for pair in sorted(output) for value in pair]


def gzip_json(path, value):
    with path.open("wb") as out:
        with gzip.GzipFile(filename="", mode="wb", fileobj=out, mtime=0) as gz:
            gz.write(json.dumps(value, ensure_ascii=True, separators=(",", ":")).encode())


def convert(archive_path, destination):
    assert len(BASE_CODES) == 256
    with ZipFile(archive_path) as archive:
        package = f"ncbi_dataset/data/{ACCESSION}/"
        report = json.loads(archive.read("ncbi_dataset/data/assembly_data_report.jsonl").splitlines()[0])
        if report["accession"] != ACCESSION or report["annotationInfo"]["name"] != ANNOTATION:
            raise ValueError("The download does not contain the pinned assembly and annotation.")
        for filename, expected in SOURCE_HASHES.items():
            with archive.open(package + filename) as source:
                actual = hashlib.file_digest(source, "sha256").hexdigest()
            if actual != expected:
                raise ValueError(f"{filename} differs from the pinned WS298 input.")
        with archive.open(package + ACCESSION + "_WBcel235_genomic.fna") as source:
            genome = dict(fasta(source))
        with archive.open(package + "rna.fna") as source:
            rna = dict(fasta(source))
        with archive.open(package + "genomic.gff") as source:
            genes, features, exons, coding, chromosome_labels = parse_annotation(source)

    genome_order = sorted(genome)
    genome_indices = {key: i for i, key in enumerate(genome_order)}
    transcripts_by_gene = defaultdict(list)
    used_accessions = set()
    counts = Counter()
    mature_mirna_parents = {row[0] for row in features.values() if row[3] == "miRNA"}
    for feature_id, row in features.items():
        parent, name, rna_id, kind, chromosome, start, end, strand = row
        gene = owner_of(parent, features, genes)
        if gene is None:
            continue
        if kind == "primary_transcript" and feature_id in mature_mirna_parents:
            continue
        ordered = sorted(exons.get(feature_id) or [(start, end)], reverse=strand == "-")
        pieces = [genome[chromosome][left:right] for left, right in ordered]
        if strand == "-":
            pieces = [piece.translate(COMPLEMENT)[::-1] for piece in pieces]
        genomic_sequence = b"".join(pieces)
        sequence = rna.get(rna_id, genomic_sequence)
        if rna_id in rna:
            used_accessions.add(rna_id)
        if genomic_sequence != sequence:
            counts["rnaSequenceDiffersFromExons"] += 1
        cds = coding_coordinates(ordered, coding.get(feature_id, ()), strand)
        if len(genomic_sequence) != len(sequence):
            counts["rnaLengthDiffersFromExons"] += 1
            cds = []
        record = [
            rna_id or feature_id, name, kind, len(sequence),
            [point for pair in ordered for point in pair], cds
        ]
        transcripts_by_gene[gene].append((record, sequence))
        counts[kind] += 1

    # The mitochondrial genes and a few pseudogenic RNAs can have exons directly
    # under the gene feature in RefSeq. Use that mature model where present.
    for gene_key, gene in genes.items():
        if gene_key in transcripts_by_gene or gene_key not in exons:
            continue
        stable_id, symbol, _, kind, chromosome, _, _, strand = gene
        ordered = sorted(exons[gene_key], reverse=strand == "-")
        pieces = [genome[chromosome][left:right] for left, right in ordered]
        if strand == "-":
            pieces = [piece.translate(COMPLEMENT)[::-1] for piece in pieces]
        sequence = b"".join(pieces)
        cds = coding_coordinates(ordered, coding.get(gene_key, ()), strand)
        record = [stable_id, symbol, kind, len(sequence), [x for pair in ordered for x in pair], cds]
        transcripts_by_gene[gene_key].append((record, sequence))
        counts["geneLevelMatureModels"] += 1

    destination.mkdir(parents=True, exist_ok=True)
    chromosome_rows = []
    for index, accession in enumerate(genome_order):
        asset = f"chromosome-{index}.bin.gz"
        sequence = genome[accession]
        (destination / asset).write_bytes(gzip.compress(sequence.translate(BASE_CODES), compresslevel=6, mtime=0))
        chromosome_rows.append([chromosome_labels.get(accession, accession), accession, asset, len(sequence)])

    gene_rows = []
    transcript_rows = []
    total_bases = 0
    with (destination / "transcripts.bin.gz").open("wb") as out:
        with gzip.GzipFile(filename="", mode="wb", fileobj=out, compresslevel=6, mtime=0) as gz:
            for key, gene in sorted(genes.items(), key=lambda row: (row[1][1].casefold(), row[1][0])):
                stable_id, symbol, locus, kind, chromosome, start, end, strand = gene
                mature = sorted(transcripts_by_gene.get(key, []), key=lambda row: row[0][1])
                first = len(transcript_rows)
                gene_index = len(gene_rows)
                for record, sequence in mature:
                    transcript_rows.append([*record, gene_index, total_bases])
                    gz.write(sequence.translate(BASE_CODES))
                    total_bases += len(sequence)
                gene_rows.append([
                    stable_id, symbol, locus, kind, genome_indices[chromosome],
                    start, end, strand, first, len(mature)
                ])

    metadata = {
        "schema": 1,
        "assembly": report["assemblyInfo"]["assemblyName"],
        "accession": ACCESSION,
        "annotation": ANNOTATION,
        "annotationDate": report["annotationInfo"]["releaseDate"],
        "species": "Caenorhabditis elegans",
        "strain": "Bristol N2",
        "source": REPORT_URL,
        "download": DOWNLOAD,
        "sourceSha256": SOURCE_HASHES,
        "geneCount": len(gene_rows),
        "transcriptCount": len(transcript_rows),
        "transcriptBases": total_bases,
        "featureCounts": dict(sorted(counts.items())),
        "refseqRnaRecords": len(rna),
        "unusedRefseqRnaRecords": sorted(rna.keys() - used_accessions),
    }
    gzip_json(destination / "catalog.json.gz", {
        "metadata": metadata, "chromosomes": chromosome_rows, "genes": gene_rows,
        "transcripts": transcript_rows
    })
    metadata["assets"] = {
        path.name: {"bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
        for path in sorted(destination.glob("*.gz"))
    }
    (destination / "reference.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"{metadata['annotation']}: {len(gene_rows):,} genes, {len(transcript_rows):,} mature RNA models, "
          f"{total_bases:,} bases. {len(metadata['unusedRefseqRnaRecords'])} RefSeq RNA records not used.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    parser.add_argument("--output", type=Path, default=Path("src/data"))
    parser.add_argument("--download", action="store_true", help="Fetch the pinned assembly package")
    options = parser.parse_args()
    if options.download:
        options.archive.parent.mkdir(parents=True, exist_ok=True)
        temporary = options.archive.with_suffix(".download")
        urlretrieve(DOWNLOAD, temporary)
        shutil.move(temporary, options.archive)
    convert(options.archive, options.output)
