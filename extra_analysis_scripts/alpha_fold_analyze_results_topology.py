from alpha_fold_support_scripts import DownloadAlphaFoldSequences

data_h = open("../info/summary2_filtered.txt")
data = data_h.readlines()
data_h.close()
uni_ids = [i.split("---")[1].split("_")[0] for i in data if "---" in i]
unique_ids = []
for i in uni_ids:
    if uni_ids not in unique_ids:
        unique_ids.append(i)

print(unique_ids)
DownloadAlphaFoldSequences(unique_ids, "../temp_files/full_alignment_filtered.fasta")
