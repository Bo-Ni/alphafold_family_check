file_h = open("../temp_files/nodata_analysis/nodata_alignment_old.fasta")
file = file_h.read()
file_h.close()

seq = file.split(">")
seq = [">" + i for i in seq]

file_h = open("../info/summary2_filtered.txt")
file = file_h.read()
file_h.close()

hits = file.split("Group")[1].split("\n")[1:-1]

hits = [i.split("---")[1].split("_")[0] for i in hits]

data = []
for i in hits:
    for i1 in seq:
        if i in i1:
            if i1 not in data:
                data.append(i1)


newfile = open("../temp_files/nodata_analysis/nodata_alignment_filtered.fasta", "w")
for i1 in data:
    newfile.write(i1)
newfile.close()
