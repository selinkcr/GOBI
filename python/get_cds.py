import os
import csv


genome_gff_paths = {
    "Bombyx Mori": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Bombyx Mori/GCF_030269925.1/GCF_030269925.1_ASM3026992v2_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Bombyx Mori/GCF_030269925.1/genomic.gff"),
    "Bombyx Mandarina": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/BombyxMandarina/ncbi_dataset/data/GCF_003987935.1/GCF_003987935.1_ASM398793v1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/BombyxMandarina/ncbi_dataset/data/GCF_003987935.1/genomic.gff"),
    "Plutella Xylostella": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Plutella xylostella/ncbi_dataset/data/GCF_932276165.1/GCF_932276165.1_ilPluXylo3.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Plutella xylostella/ncbi_dataset/data/GCF_932276165.1/genomic.gff"),
    "Galleria Mellonella": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Galleria mellonella/ncbi_dataset/data/GCF_026898425.1/GCF_026898425.1_CSIRO_AGI_GalMel_v1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Galleria mellonella/ncbi_dataset/data/GCF_026898425.1/genomic.gff"),
    "Achroia Grisella": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Achroia grisella/ncbi_dataset/data/GCF_030625045.1/GCF_030625045.1_CSIRO_AGI_Agris_v1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Achroia grisella/ncbi_dataset/data/GCF_030625045.1/genomic.gff"),
    "Manduca Sexta": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Manduca sexta/ncbi_dataset/data/GCF_014839805.1/GCF_014839805.1_JHU_Msex_v1.0_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Manduca sexta/ncbi_dataset/data/GCF_014839805.1/genomic.gff"),
    "Cydia Amplana": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Cydia amplana/ncbi_dataset/data/GCF_948474715.1/GCF_948474715.1_ilCydAmpl1.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Cydia amplana/ncbi_dataset/data/GCF_948474715.1/genomic.gff"),
    "Vanessa Cardui": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Vanessa Cardui/ncbi_dataset/data/GCF_905220365.1/GCF_905220365.1_ilVanCard2.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Vanessa Cardui/ncbi_dataset/data/GCF_905220365.1/genomic.gff"),
    "Maniola Jurtina": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Maniola Jurtina/ncbi_dataset/data/GCF_905333055.1/GCF_905333055.1_ilManJurt1.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Maniola Jurtina/ncbi_dataset/data/GCF_905333055.1/genomic.gff"),
    "Bombus Impatiens": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Bombus Impatiens/ncbi_dataset/data/GCF_000188095.3/GCF_000188095.3_BIMP_2.2_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Bombus Impatiens/ncbi_dataset/data/GCF_000188095.3/genomic.gff"),
    "Apis Mellifera": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Appis Mellifera/ncbi_dataset/data/GCF_003254395.2/GCF_003254395.2_Amel_HAv3.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Appis Mellifera/ncbi_dataset/data/GCF_003254395.2/genomic.gff"),
    "Monomorium Pharaonis": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Monomorium pharaonis/ncbi_dataset/data/GCF_013373865.1/GCF_013373865.1_ASM1337386v2_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Monomorium pharaonis/ncbi_dataset/data/GCF_013373865.1/genomic.gff"),
    "Leptopilina Boulardi": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Leptopilina boulardi/ncbi_dataset/data/GCF_019393585.1/GCF_019393585.1_ZJU_Lbou_2.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Leptopilina boulardi/ncbi_dataset/data/GCF_019393585.1/genomic.gff"),
    "Vespula Vulgaris": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Vespula vulgaris/ncbi_dataset/data/GCF_905475345.1/GCF_905475345.1_iyVesVulg1.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Vespula vulgaris/ncbi_dataset/data/GCF_905475345.1/genomic.gff"),
    "Cotesia Glomerata": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Cotesia glomerata/ncbi_dataset/data/GCF_020080835.1/GCF_020080835.1_MPM_Cglom_v2.3_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Cotesia glomerata/ncbi_dataset/data/GCF_020080835.1/genomic.gff"),
    "Nilaparvata Lugens": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Nilaparvata lugens/ncbi_dataset/ncbi_dataset/data/GCF_014356525.2/GCF_014356525.2_ASM1435652v1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Nilaparvata lugens/ncbi_dataset/ncbi_dataset/data/GCF_014356525.2/genomic.gff"),
    "Homalodisca Vitripennis": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Homalodisca vitripennis/ncbi_dataset/ncbi_dataset/data/GCF_021130785.1/GCF_021130785.1_UT_GWSS_2.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Homalodisca vitripennis/ncbi_dataset/ncbi_dataset/data/GCF_021130785.1/genomic.gff"),
    "Planococcus Citri": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Planococcus citri/ncbi_dataset/ncbi_dataset/data/GCF_950023065.1/GCF_950023065.1_ihPlaCitr1.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Planococcus citri/ncbi_dataset/ncbi_dataset/data/GCF_950023065.1/genomic.gff"),
    "Cimex Lectularius": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Cimex lectularius/ncbi_dataset/data/GCF_000648675.2/GCF_000648675.2_Clec_2.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Cimex lectularius/ncbi_dataset/data/GCF_000648675.2/genomic.gff"),
    "Halyomorpha Halys": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Halyomorpha halys/ncbi_dataset/data/GCF_000696795.3/GCF_000696795.3_Hhal_1.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Halyomorpha halys/ncbi_dataset/data/GCF_000696795.3/genomic.gff"),
    "Ischnura Elegans": ("/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Ischnura elegans /ncbi_dataset/data/GCF_921293095.1/GCF_921293095.1_ioIscEleg1.1_genomic.fna", "/Users/peanut/Flanking_Genes/Flanking_Genes_New_2/Ischnura elegans /ncbi_dataset/data/GCF_921293095.1/genomic.gff"),
}

# 🔹 Define gene IDs for each species
gene_dict = {
    "Bombyx Mori": [693030, 101738200, 134200129],
    "Bombyx Mandarina": [114246424, 114247220, 114239327],
    "Plutella Xylostella": [105393633, 105392941, 119694615],
    "Galleria Mellonella": [113519866, 113513958],
    "Achroia Grisella": [131851741, 131848063],
    "Manduca Sexta": [115443893, 115449923, 115447212, 115442891],
    "Cydia Amplana": [134655470, 134653699, 134647545],
    "Vanessa Cardui": [124538299, 124537616, 124535929],
    "Maniola Jurtina": [123875302, 123871537],
    "Bombus Impatiens": [100741543, 105681168],
    "Apis Mellifera": [726750],
    "Monomorium Pharaonis": [105830070],
    "Leptopilina Boulardi": [127279011],
    "Vespula Vulgaris": [127062096, 127067120],
    "Cotesia Glomerata": [123262670, 123274246],
    "Nilaparvata Lugens": [111057060, 111050149, 111046429, 111045123],
    "Homalodisca Vitripennis": [124361079, 124363296, 124363764, 124363763],
    "Planococcus Citri": [135844163, 135836629],
    "Cimex Lectularius": [106670165, 106670162, 106670148, 106661923],
    "Halyomorpha Halys": [106679084],
    "Ischnura Elegans": [124173615, 124161124, 124159079, 124174215, 124153170],
}
output_file = "cds_extraction.csv"

# CSV Header
csv_columns = [
    "species", "gff_seqid", "gff_start", "gff_end", "gff_strand",
    "gff_gene_id", "gff_attribute_ID", "sequence"
]

debug = True


def reverse_complement(seq):
    """Return the reverse complement of a DNA sequence."""
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]


def load_fasta(fasta_file):
    """Load genome sequences from a FASTA file into a dictionary."""
    fasta_sequences = {}
    with open(fasta_file, "r") as file:
        seq_id = ""
        sequence = []
        for line in file:
            if line.startswith(">"):
                if seq_id:
                    fasta_sequences[seq_id] = "".join(sequence)
                seq_id = line.strip().split()[0][1:]  # Get first part after ">"
                sequence = []
            else:
                sequence.append(line.strip())

        if seq_id:
            fasta_sequences[seq_id] = "".join(sequence)

    return fasta_sequences


def extract_sequence(fasta_dict, seqid, start, end, strand):
    """Extract sequence from FASTA file based on coordinates."""
    if seqid not in fasta_dict:
        return "Sequence not found"

    sequence = fasta_dict[seqid][start-1:end]  # Adjust for 1-based indexing

    if strand == "-":
        sequence = reverse_complement(sequence)

    return sequence


def parse_gff_cds(gff_file, fasta_dict, gene_id):
    """Extract CDS and their sequences for a specific gene ID."""
    cds_info = []
    mRNA_ids = {}  # Store all mRNA IDs linked to the gene

    with open(gff_file, "r") as file:
        for line in file:
            if line.startswith("#"):
                continue  # Skip comments

            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue  # Skip malformed lines

            seqid, source, feature_type, start, end, score, strand, phase, attributes = fields
            start, end = int(start), int(end)

            # Convert attributes to dictionary
            attr_dict = {kv.split("=")[0]: kv.split("=")[1] for kv in attributes.split(";") if "=" in kv}
            dbxref = attr_dict.get("Dbxref", "")

            # Find the correct mRNA linked to this gene
            if feature_type == "mRNA" and f"GeneID:{gene_id}" in dbxref:
                mRNA_ids[attr_dict.get("ID", "")] = True  # Save mRNA ID

            # Extract CDS associated with the correct mRNA
            if feature_type == "CDS" and "Parent" in attr_dict:
                parent_id = attr_dict["Parent"]
                if parent_id in mRNA_ids:
                    sequence = extract_sequence(fasta_dict, seqid, start, end, strand)
                    cds_info.append([seqid, start, end, strand, gene_id, attr_dict.get("ID", ""), sequence])

    return cds_info


if __name__ == '__main__':
    with open(output_file, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(csv_columns)  # Write CSV header

        for species, gene_ids in gene_dict.items():
            if species not in genome_gff_paths:
                print(f"Skipping {species} (No genome/GFF paths provided)")
                continue

            fasta_file, gff_file = genome_gff_paths[species]

            if not os.path.exists(fasta_file) or not os.path.exists(gff_file):
                print(f"ERROR: Missing files for {species}. Check paths:\nFASTA: {fasta_file}\nGFF: {gff_file}")
                continue

            # Load genome sequences (without Biopython)
            fasta_dict = load_fasta(fasta_file)

            for gene_id in gene_ids:
                if debug:
                    print(f"Processing {species}, Gene ID: {gene_id}")

                cds_info = parse_gff_cds(gff_file, fasta_dict, gene_id)

                # Write CDS data
                for row in cds_info:
                    writer.writerow([species] + row)

    print(f"CDS extraction complete. Results saved in {output_file}")
