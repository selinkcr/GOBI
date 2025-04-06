gene_positions = {}

with open("/Users/peanut/Library/CloudStorage/OneDrive-Personal/vscode/regions_new_result/halysFlanking.tsv") as file:
    for line in file:
        cols = line.strip().split("\t")
        gene = cols[1]
        start = int(cols[6])
        end = int(cols[7])
        pos = min(start, end)  # Use the lower coordinate (handles reverse strand)

        if gene not in gene_positions:
            gene_positions[gene] = pos
        else:
            gene_positions[gene] = min(gene_positions[gene], pos)

# Sort genes by their genomic position
sorted_genes = [gene for gene, _ in sorted(gene_positions.items(), key=lambda x: x[1])]

# Print each gene on its own line
for gene in sorted_genes:
    print(gene)
