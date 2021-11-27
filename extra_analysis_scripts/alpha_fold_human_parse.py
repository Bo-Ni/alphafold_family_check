from alpha_fold_main import AlphaFoldWorkflow


if __name__ == "__main__":
    data_h = open("../input_files/Knots on Alphafold proteins - Sulkowska LabTable.csv")
    data = data_h.readlines()
    data_h.close()
    UNI_IDs = [[i.split(",")[1].strip('"').rstrip("F1").rstrip("-"),
                int(i.split(",")[5].strip('"').split("-")[0]),
                int(i.split(",")[5].strip('"').split("-")[1])] for i in data[1:]]

    result = AlphaFoldWorkflow(UNI_IDs, savefile="temp_files/HUMAN_alpha_fold.txt")

    for i in result:
        print(i)