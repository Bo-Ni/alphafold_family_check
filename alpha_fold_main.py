import os
from tqdm import tqdm
import asyncio
from alpha_fold_domain_check import AlphaFoldKnotDomainCheck
from alpha_fold_structure_check import AlphaFoldStructureCheck


def AlphaFoldWorkflow(data, savefile=None):
    def continue_search(data, savefile_name):
        temp_h = open(savefile_name)
        temp = temp_h.readlines()
        temp_h.close()

        newdata = []
        for i in data:
            t_i = "_".join([str(i0) for i0 in i])
            state = False
            for i1 in temp:
                if t_i in i1:
                    state = True
            if not state:
                newdata.append(i)

        return newdata

    async def MAIN_WORKER_FUNC(batch):
        batch_str = [i0[0] for i0 in batch]

        path1 = asyncio.create_task(AlphaFoldKnotDomainCheck(batch))
        path2 = asyncio.create_task(AlphaFoldStructureCheck(batch_str))
        MAIN_WORKER = await asyncio.gather(path1, path2)

        return MAIN_WORKER


    if not savefile:
        savefile_name = "savefile.txt"
    elif not os.path.exists(savefile):
        savefile_name = savefile
    elif os.path.exists(savefile):
        savefile_name = savefile
        data = continue_search(data, savefile)

    result = []
    batch_size = 1
    data = [data[i:i + batch_size] for i in range(0, len(data), batch_size)]
    for i in tqdm(range(len(data))):
        temp_batch_string_domain = []
        # ASYNCIO PART ---------------------------
        domain_result, struct_result = asyncio.run(MAIN_WORKER_FUNC(data[i]))
        # GATHER ---------------------------------
        for i11 in domain_result:
            for i12 in struct_result:
                if i12[0] in i11[0]:
                    temp_result = i11[:]
                    temp_result.extend(i12[1:])
                    result.append(temp_result)
                    temp_string = ";".join(["_".join([str(i0) for i0 in i11[0]]),
                                            "_".join([str(i0) for i0 in i11[1]]),
                                            ",".join(["_".join([str(i0) for i0 in i]) for i in i11[2]]),
                                            i12[1],
                                            str(i12[2]),
                                            ",".join(str(i12[3]))]) + "\n"
                    temp_batch_string_domain.append(temp_string)

        temp_file = open(savefile_name, "a")
        temp_file.write("".join(temp_batch_string_domain))
        temp_file.close()

    return result


if __name__ == "__main__":
    res = AlphaFoldWorkflow([["Q9FLD5", 1, 100], ["Q6PL18", 1, 100], ["Q5T9A4", 1, 100]])
    # print(res)
