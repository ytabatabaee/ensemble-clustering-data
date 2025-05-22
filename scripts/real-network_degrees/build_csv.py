import graph_tool.all as gt

dataset_arr = []
dataset_sizes = {}
with open("./dataset_tuple.list", "r") as f:
    for line in f:
        current_dataset = line.strip().split(",")[0]
        current_dataset_size = line.strip().split(",")[1]
        if current_dataset == "dataset":
            continue
        dataset_arr.append(current_dataset)
        dataset_sizes[current_dataset] = int(current_dataset_size)

with open("./output/average_degree", "w") as f:
    f.write(f"dataset,size,average_degree\n")
    for current_dataset in dataset_arr:
        current_dataset_size = dataset_sizes[current_dataset]
        current_dataset_average_degree = gt.collection.ns_info[current_dataset]["analyses"]["average_degree"]
        f.write(f"{current_dataset},{current_dataset_size},{current_dataset_average_degree}\n")
