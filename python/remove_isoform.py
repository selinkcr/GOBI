import re
from collections import defaultdict

def filter_isoforms_and_reformat_headers(input_fasta, output_fasta):
    species_count = defaultdict(int)

    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        keep = False
        for line in infile:
            if line.startswith('>'):
                # Check isoform status
                isoform_match = re.search(r'\[isoform=(X\d+)\]', line)
                if isoform_match:
                    keep = isoform_match.group(1) == 'X1'
                else:
                    keep = True

                if keep:
                    # Extract organism name
                    organism_match = re.search(r'\[organism=([^\]]+)\]', line)
                    if organism_match:
                        organism = organism_match.group(1).replace(" ", "_")
                        species_count[organism] += 1
                        count = species_count[organism]

                        # Remove original [organism=...] part
                        line = re.sub(r'\[organism=[^\]]+\]', '', line)

                        # Build new header with numbered organism
                        numbered_organism = f"{count}_{organism}" if count > 1 else organism
                        rest_of_header = line[1:].strip()  # Remove original '>' and strip newline
                        line = f">{numbered_organism}|{rest_of_header}\n"
            if keep:
                outfile.write(line)


# Example usage
filter_isoforms_and_reformat_headers("/Users/peanut/Library/CloudStorage/OneDrive-Personal/gobi/Flanking_Genes_New_2/proteinNew.faa", "/Users/peanut/Library/CloudStorage/OneDrive-Personal/gobi/Flanking_Genes_New_2/ProteinNew.fna")
