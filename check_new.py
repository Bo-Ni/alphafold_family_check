from tqdm import tqdm
from urllib import request


def check_id(sim_info):
    query = "id:" + "%20OR%20id:".join(sim_info[1][:100])
    new_url = f"https://www.uniprot.org/uniprot/?query={query}&format=tab&columns=id,database(Pfam)&sort=score"
    try:
        temp_h = request.urlopen(new_url)
    except:
        print(len(sim_info[1]))
        print(new_url)
    temp = temp_h.readlines()
    temp = [i.decode("utf-8").rstrip().split("\t") for i in temp[1:]]
    temp = [i for i in temp if len(i) == 2]  # No family

    try:
        temp = [[i[0], "_".join([i1 for i1 in i[1].split(";") if i1])] for i in temp]
    except:
        print(temp)
    temp_d = {}

    for i1 in temp:
        if i1[1] in temp_d:
            temp_d[i1[1]] += 1
        else:
            temp_d[i1[1]] = 1

    fin_result = [sim_info[0], len(sim_info[1]), []]
    for i in temp_d:
        fin_result[-1].append([i, temp_d[i]])

    return fin_result


def decision_func(data):
    all = data[1]
    if data[2][0][1] / all > 0.6:
        data[2] = data[2][0]

    return data


data_h = open("savefilec.txt")
data = data_h.readlines()
data_h.close()
data_s = [i.rstrip().split(";") for i in data]
data_s = [i for i in data_s if not i[1]]
data_s = [i[0].split("_")[0] for i in data_s]

urls = [f"https://www.uniprot.org/uniref/UniRef90_{i}.list" for i in data_s]

result_list = []
for i in tqdm(range(len(urls))):
    temp_h = request.urlopen(urls[i])
    temp = temp_h.readlines()
    temp = [i.decode("utf-8").rstrip() for i in temp]
    search_object = [data_s[i], temp]
    result = check_id(search_object)
    result = decision_func(result)
    print(result)
    result_list.append(result)







