from tqdm import tqdm
from urllib import request
import os
from alpha_fold_analyze_results import DecodeSavefile


to_download = []
organism_files = os.listdir("../temp_files/organism_data")
organism_files = [i for i in organism_files if "alpha" in i]
for organism in organism_files:
    organism_data = DecodeSavefile(f"../temp_files/organism_data/{organism}")
    for single_data in organism_data:
        if single_data[1] == "nodata":
            to_download.append(single_data[0][0])

batch_size = 50
data = [to_download[i:i + batch_size] for i in range(0, len(to_download), batch_size)]

full_ali = ""
for i in tqdm(data):
    url = f"https://www.uniprot.org/uniprot/?query=id:{'%20OR%20id:'.join(i)}&format=fasta&sort=score"
    site = request.urlopen(url).read().decode("utf-8")
    newfile = open("../temp_files/nodata_alignment.fasta", "a")
    newfile.write(site)
    newfile.close()



