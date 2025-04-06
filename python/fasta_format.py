def parse_fasta_and_filter(input_path, output_path):
    species_seen = set()
    selected_records = []

    with open(input_path, "r") as infile:
        header = ""
        sequence = []

        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if header and sequence:
                    species = get_species(header)
                    if species and species not in species_seen:
                        species_seen.add(species)
                        new_header = f">{species.replace(' ', '_')}|{get_transcript(header)}|{get_geneid(header)}"
                        selected_records.append((new_header, ''.join(sequence)))
                header = line
                sequence = []
            else:
                sequence.append(line)

        # Handle the last record
        if header and sequence:
            species = get_species(header)
            if species and species not in species_seen:
                species_seen.add(species)
                new_header = f">{species.replace(' ', '_')}|{get_transcript(header)}|{get_geneid(header)}"
                selected_records.append((new_header, ''.join(sequence)))

    with open(output_path, "w") as outfile:
        for header, seq in selected_records:
            outfile.write(header + "\n")
            for i in range(0, len(seq), 60):
                outfile.write(seq[i:i + 60] + "\n")


def get_species(header):
    tag = "[organism="
    if tag in header:
        return header.split(tag)[1].split("]")[0].strip()
    return None


def get_transcript(header):
    tag = "[transcript="
    if tag in header:
        return header.split(tag)[1].split("]")[0].strip()
    return "no_transcript"


def get_geneid(header):
    tag = "[GeneID="
    if tag in header:
        return header.split(tag)[1].split("]")[0].strip()
    return "no_geneid"


# Example usage
input_path = "/Users/peanut/Library/CloudStorage/OneDrive-Personal/gobi/Flanking_Genes_New_2/single copy/ncbi_dataset/data/cds.fna"
output_path = "/Users/peanut/Library/CloudStorage/OneDrive-Personal/gobi/Flanking_Genes_New_2/single copy/ncbi_dataset/data/cleaned_cds.fna"
parse_fasta_and_filter(input_path, output_path)
