import pandas as pd

def filter_top_blast_hits_preserve_format(blast_file, output_file):
    # Read BLAST TSV (assumes no header, tab-separated)
    df = pd.read_csv(blast_file, sep="\t", header=None, comment="#")

    # Extract species from column 1 (sseqid), which is at index 1
    df["species"] = df[1].str.extract(r"(^[^|]+)")

    # Convert e-value (column 10, index 10) to numeric
    df[10] = pd.to_numeric(df[10], errors="coerce")

    # Sort by e-value and keep first (best) entry per species
    best_hits = df.sort_values(10).groupby("species", as_index=False).first()

    # Drop the extra 'species' column
    best_hits = best_hits.drop(columns=["species"])

    # Save result with original format
    best_hits.to_csv(output_file, sep="\t", header=False, index=False)
    print(f"Saved top hits per species to: {output_file}")

# Example usage
if __name__ == "__main__":
    filter_top_blast_hits_preserve_format("regions_new_result/achroia.tsv", "regions_new_gff/achroia.tsv")
