import os
import sys

res = []
cwd = os.getcwd()
filenames = sorted(os.listdir(cwd + "/gen"))[-int(sys.argv[-1]):]
print(filenames)
for filename in filenames:
    with open(cwd + "/gen/" + filename, "r") as file:
        lines = file.readlines()
        for i in range(len(lines)):
            if "++ ls" in lines[i]:
                res.append({})
            elif "mpi_time:" in lines[i]:
                res[-1][int(lines[i - 1].rstrip().split(" ")[-1].replace(".in", ""))] = lines[i].rstrip().split(":")[1][:-1]

for d in res:
    keys = sorted(d.keys())
    for k in keys:
        print(d[k])
    print("\r")