import os
from alpha_fold_support_scripts import DecodeSavefile, DownloadAlphaFoldSequences

to_download = []
organism_files = os.listdir("../temp_files/organism_data")
organism_files = [i for i in organism_files if "alpha" in i]
for organism in organism_files:
    organism_data = DecodeSavefile(f"../temp_files/organism_data/{organism}")
    for single_data in organism_data:
        if single_data[1] == "nodata":
            to_download.append(single_data[0][0])

DownloadAlphaFoldSequences(to_download, "../temp_files/nodata_alignment1.fasta")
