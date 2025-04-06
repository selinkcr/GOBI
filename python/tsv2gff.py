import csv

with open("/Users/peanut/Library/CloudStorage/OneDrive-Personal/vscode/regions_result/CDScotesia.tsv") as infile, open("/Users/peanut/Library/CloudStorage/OneDrive-Personal/vscode/regions_gff/CDScotesia.gff", "w") as outfile:
    reader = csv.reader(infile, delimiter="\t")
    outfile.write("##gff-version 3\n")
    count = 1
    for row in reader:
        query_id = row[0]
        target_id = row[1]
        target_start = int(row[8])
        
        target_end = int(row[9])
        score = row[11]

        # Ensure start < end
        start = min(target_start, target_end)
        end = max(target_start, target_end)

        # Determine strand
        strand = "+" if target_start <= target_end else "-"

        # GFF fields
        seqid = target_id
        source = "blast"
        feature_type = "match"
        phase = "."
        attributes = f"ID=match{count};Name={query_id}"

        gff_line = f"{seqid}\t{source}\t{feature_type}\t{start}\t{end}\t{score}\t{strand}\t{phase}\t{attributes}\n"
        outfile.write(gff_line)
        count += 1
