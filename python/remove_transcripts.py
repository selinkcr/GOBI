import re

def filter_and_reformat_headers(input_fasta, output_fasta):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        keep = False
        for line in infile:
            if line.startswith('>'):
                # Check transcript status
                transcript_match = re.search(r'\[transcript=(X\d+)\]', line)
                if transcript_match:
                    keep = transcript_match.group(1) == 'X1'
                else:
                    keep = True

                if keep:
                    # Extract organism
                    organism_match = re.search(r'\[organism=([^\]]+)\]', line)
                    if organism_match:
                        organism = organism_match.group(1).replace(" ", "_")
                        # Remove the [organism=...] part
                        line = re.sub(r'\[organism=[^\]]+\]', '', line)
                        # Prepend organism to header
                        line = f">{organism}|{line[1:]}"
            if keep:
                outfile.write(line)

filter_and_reformat_headers("/Users/peanut/Library/CloudStorage/OneDrive-Personal/gobi/Flanking_Genes_New_2/single copy/ncbi_dataset/data/cds.fna", "/Users/peanut/Library/CloudStorage/OneDrive-Personal/gobi/Flanking_Genes_New_2/single copy/ncbi_dataset/data/cds.fna")
