import argparse
import csv
import os
import re
import subprocess

from collections import defaultdict
from pyfaidx import Fasta


# ============================================================
# 1. READ UNANNOTATED GENE IDS
# ============================================================

def read_targets(ids_file=None, results_file=None,
                 id_column="gene_id",
                 name_column="gene_name"):

    gene_ids = []

    if ids_file:

        with open(ids_file) as f:
            gene_ids = [
                line.strip()
                for line in f
                if line.strip()
            ]

    elif results_file:

        with open(results_file, newline="") as f:

            sample = f.read(4096)
            f.seek(0)

            delimiter = "\t" if "\t" in sample else ","

            reader = csv.DictReader(f, delimiter=delimiter)

            if id_column not in reader.fieldnames:
                raise ValueError(
                    f"ID column '{id_column}' not found. "
                    f"Available columns: {reader.fieldnames}"
                )

            if name_column not in reader.fieldnames:
                raise ValueError(
                    f"Annotation column '{name_column}' not found. "
                    f"Available columns: {reader.fieldnames}"
                )

            for row in reader:

                gene_id = row[id_column].strip()
                gene_name = (row[name_column] or "").strip()

                missing_names = {
                    "", "NA", "N/A", "NaN", "nan", "None"
                }

                if gene_name in missing_names or gene_name == gene_id:
                    if gene_id:
                        gene_ids.append(gene_id)

    else:
        raise ValueError("Provide either an ID file or a DESeq2 results file.")

    # Remove duplicates while preserving order
    gene_ids = list(dict.fromkeys(gene_ids))

    print(f"Found {len(gene_ids)} unannotated gene IDs.")

    return gene_ids


# ============================================================
# 2. PARSE GTF ATTRIBUTES
# ============================================================

def parse_attributes(attribute_string):

    attributes = {}

    matches = re.finditer(
        r'(\S+)\s+"([^"]*)"\s*;?',
        attribute_string
    )

    for match in matches:
        attributes[match.group(1)] = match.group(2)

    return attributes


# ============================================================
# 3. EXTRACT EXON COORDINATES FROM GTF
# ============================================================

def get_exons(gtf_file, target_ids):

    target_lookup = {
        gene_id.split(".")[0]: gene_id
        for gene_id in target_ids
    }

    transcripts = defaultdict(lambda: {
        "gene_id": None,
        "transcript_id": None,
        "gene_name": "",
        "chromosome": None,
        "strand": None,
        "exons": set()
    })

    with open(gtf_file) as f:

        for line in f:

            if line.startswith("#"):
                continue

            columns = line.rstrip("\n").split("\t")

            if len(columns) != 9:
                continue

            chromosome = columns[0]
            feature = columns[2]

            if feature != "exon":
                continue

            start = int(columns[3])
            end = int(columns[4])
            strand = columns[6]

            attributes = parse_attributes(columns[8])

            gene_id = attributes.get("gene_id")
            transcript_id = attributes.get("transcript_id")

            if not gene_id or not transcript_id:
                continue

            clean_id = gene_id.split(".")[0]

            if clean_id not in target_lookup:
                continue

            original_id = target_lookup[clean_id]

            key = (original_id, transcript_id)

            transcripts[key]["gene_id"] = original_id
            transcripts[key]["transcript_id"] = transcript_id

            transcripts[key]["gene_name"] = attributes.get(
                "gene_name", ""
            )

            transcripts[key]["chromosome"] = chromosome
            transcripts[key]["strand"] = strand

            transcripts[key]["exons"].add((start, end))

    print(f"Found {len(transcripts)} transcripts.")

    return transcripts


# ============================================================
# 4. REVERSE COMPLEMENT
# ============================================================

def reverse_complement(sequence):

    complement = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn"
    )

    return sequence.translate(complement)[::-1]


# ============================================================
# 5. EXTRACT SPLICED TRANSCRIPT SEQUENCES
# ============================================================

def extract_sequences(transcripts, genome_file,
                      output_fasta, output_csv):

    genome = Fasta(genome_file)

    genome_contigs = set(genome.keys())

    extracted = 0

    with open(output_fasta, "w") as fasta_out, \
         open(output_csv, "w", newline="") as csv_out:

        writer = csv.writer(csv_out)

        writer.writerow([
            "Gene_ID",
            "Transcript_ID",
            "Gene_Name",
            "Chromosome",
            "Start",
            "End",
            "Strand",
            "Exon_Count",
            "Sequence_Length"
        ])

        for key, info in transcripts.items():

            gene_id = info["gene_id"]
            transcript_id = info["transcript_id"]

            chromosome = info["chromosome"]
            strand = info["strand"]

            if chromosome not in genome_contigs:

                print(
                    f"WARNING: Chromosome {chromosome} "
                    f"not found for {gene_id}"
                )

                continue

            # Sort exons in genomic order
            exons = sorted(info["exons"])

            sequence_parts = []

            for start, end in exons:

                # GTF: 1-based inclusive
                # pyfaidx: 0-based, end-exclusive

                exon_sequence = str(
                    genome[chromosome][start-1:end]
                )

                sequence_parts.append(exon_sequence)

            # Concatenate exons
            sequence = "".join(sequence_parts).upper()

            # Reverse complement entire transcript
            # for negative-strand genes

            if strand == "-":
                sequence = reverse_complement(sequence)

            if not sequence:
                continue

            # Unique transcript identifier in FASTA
            fasta_id = f"{gene_id}|{transcript_id}"

            fasta_out.write(f">{fasta_id}\n")

            for i in range(0, len(sequence), 60):
                fasta_out.write(sequence[i:i+60] + "\n")

            starts = [x[0] for x in exons]
            ends = [x[1] for x in exons]

            writer.writerow([
                gene_id,
                transcript_id,
                info["gene_name"],
                chromosome,
                min(starts),
                max(ends),
                strand,
                len(exons),
                len(sequence)
            ])

            extracted += 1

    print(f"Extracted {extracted} transcript sequences.")

    return extracted


# ============================================================
# 6. IDENTIFY GENES WITHOUT EXON RECORDS
# ============================================================

def save_missing_ids(target_ids, transcripts, output_file):

    found_ids = {
        info["gene_id"]
        for info in transcripts.values()
    }

    missing = set(target_ids) - found_ids

    with open(output_file, "w") as f:

        for gene_id in sorted(missing):
            f.write(gene_id + "\n")

    print(f"{len(missing)} genes had no exon records in the GTF.")


# ============================================================
# 7. RUN REMOTE NCBI BLAST
# ============================================================

def run_blast(input_fasta, output_file):

    command = [

        "blastn",

        "-query", input_fasta,

        "-db", "nt",

        "-remote",

        "-out", output_file,

        "-outfmt",
        "6 qseqid saccver pident length qlen qcovs evalue bitscore stitle",

        "-evalue", "1e-5",

        "-max_target_seqs", "10"

    ]

    print("\nSubmitting sequences to NCBI BLAST...")

    subprocess.run(command, check=True)

    print("BLAST completed.")


# ============================================================
# 8. MAIN
# ============================================================

def main():

    parser = argparse.ArgumentParser(
        description="Extract spliced gene sequences and run remote BLAST."
    )

    inputs = parser.add_mutually_exclusive_group(required=True)

    inputs.add_argument("--ids", help="File containing Ensembl IDs")

    inputs.add_argument("--results", help="DESeq2 CSV/TSV results file")

    parser.add_argument(
        "--id-column",
        default="gene_id",
        help="Gene ID column in DESeq2 results"
    )

    parser.add_argument(
        "--name-column",
        default="gene_name",
        help="Gene annotation column in DESeq2 results"
    )

    parser.add_argument(
        "--gtf",
        required=True,
        help="GTF annotation file"
    )

    parser.add_argument(
        "--genome",
        required=True,
        help="Reference genome FASTA"
    )

    parser.add_argument(
        "--output",
        default="output",
        help="Output directory"
    )

    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Output paths

    fasta_file = os.path.join(
        args.output,
        "spliced_transcripts.fasta"
    )

    coordinates_file = os.path.join(
        args.output,
        "transcript_coordinates.csv"
    )

    missing_file = os.path.join(
        args.output,
        "missing_gene_ids.txt"
    )

    blast_file = os.path.join(
        args.output,
        "blast_results.tsv"
    )

    # Step 1: Read IDs

    print("\nSTEP 1: Reading gene IDs")

    gene_ids = read_targets(
        ids_file=args.ids,
        results_file=args.results,
        id_column=args.id_column,
        name_column=args.name_column
    )

    # Step 2: Find exons

    print("\nSTEP 2: Finding exons in GTF")

    transcripts = get_exons(
        args.gtf,
        gene_ids
    )

    # Step 3: Report missing IDs

    save_missing_ids(
        gene_ids,
        transcripts,
        missing_file
    )

    if not transcripts:
        raise RuntimeError("No matching exon records found.")

    # Step 4: Extract sequences

    print("\nSTEP 3: Extracting spliced sequences")

    extracted = extract_sequences(
        transcripts,
        args.genome,
        fasta_file,
        coordinates_file
    )

    if extracted == 0:
        raise RuntimeError("No sequences were extracted.")

    # Step 5: BLAST

    print("\nSTEP 4: Running remote BLAST")

    run_blast(
        fasta_file,
        blast_file
    )

    print("\nPIPELINE COMPLETE")

    print("\nOutput files:")
    print(f"FASTA: {fasta_file}")
    print(f"Coordinates: {coordinates_file}")
    print(f"Missing IDs: {missing_file}")
    print(f"BLAST results: {blast_file}")


if __name__ == "__main__":
    main()
