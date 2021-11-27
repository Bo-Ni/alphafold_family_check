from alpha_fold_support_scripts import DownloadAlphaFoldSequences
from alpha_fold_analyze_results_families import get_non_human_knot_info, get_human_knot_info

knot_list = get_human_knot_info()
knot_list.extend(get_non_human_knot_info())

uni_ids = knot_list = [i[0] for i in knot_list]
unique_ids = []
for i in uni_ids:
    if i not in unique_ids:
        unique_ids.append(i)

DownloadAlphaFoldSequences(unique_ids, "../temp_files/full_alignment.fasta")

print(unique_ids)