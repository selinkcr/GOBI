import csv

input_file = "/Users/peanut/Library/CloudStorage/OneDrive-Personal/vscode/regions_new_result/Cotesia.tsv"
output_file = "//Users/peanut/Library/CloudStorage/OneDrive-Personal/vscode/regions_new_gff/Cotesia.gff"

gff_lines = []
with open(input_file, "r") as infile:
    tsv_reader = csv.reader(infile, delimiter="\t")
    match_id = 1
    for row in tsv_reader:
        query = row[0]
        subject = row[1]
        identity = row[2]
        qstart = int(row[6])
        qend = int(row[7])
        
        start = min(qstart, qend)
        end = max(qstart, qend)
        strand = "+" if qend >= qstart else "-"
        score = identity  # or use int(float(identity)) for rounded

        gff_line = [
            query,            # seqid
            "blast",          # source
            "match",          # type
            str(start),       # start
            str(end),         # end
            score,            # score
            strand,           # strand
            ".",              # phase
            f"ID=match{match_id};Name={subject}"  # attributes
        ]
        gff_lines.append("\t".join(gff_line))
        match_id += 1

with open(output_file, "w") as outfile:
    for line in gff_lines:
        outfile.write(line + "\n")

print(f"Converted {match_id - 1} entries to GFF format.")
