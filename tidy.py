import os
new_dir = "C:\\Users\\DuYih\\Desktop\\metadata"
os.chdir(new_dir)
# 定义列表名称
sample_names = [
    "CK0001", "CK0002", "CK0003", "CK0601", "CK0602", "CK0603", "CK1201", "CK1202", "CK1203",
    "CK1801", "CK1802", "CK1803", "LDPE0601", "LDPE0602", "LDPE0603", "LDPE1201", "LDPE1202", "LDPE1203",
    "LDPE1801", "LDPE1802", "LDPE1803", "PET0601", "PET0602", "PET0603", "PET1201", "PET1202", "PET1203",
    "PET1801", "PET1802", "PET1803", "PP0601", "PP0602", "PP0603", "PP1201", "PP1202", "PP1203", "PP1801",
    "PP1802", "PP1803", "SLDPE0601", "SLDPE0602", "SLDPE0603", "SLDPE1201", "SLDPE1202", "SLDPE1203",
    "SLDPE1801", "SLDPE1802", "SLDPE1803", "SPET0601", "SPET0602", "SPET0603", "SPET1201", "SPET1202",
    "SPET1203", "SPET1801", "SPET1802", "SPET1803", "SPP0601", "SPP0602", "SPP0603", "SPP1201", "SPP1202",
    "SPP1203", "SPP1801", "SPP1802", "SPP1803"
]

# 读取DNA文本文件
with open("C:\\Users\\DuYih\\Desktop\\metadata\\meta.txt", "r") as f:
    lines = f.readlines()
# 创建文件夹
if not os.path.exists("output"):
    os.makedirs("output")

# 根据名称行将序列保存为单独的文件
for i in range(len(lines)):
    if lines[i].startswith(">"):
        name = lines[i].strip()[1:].split("_")[0]
        
        # if any(sample_name in name for sample_name in sample_names):
        #     seq = ""
        #     j = i + 1
        #     while j < len(lines) and not lines[j].startswith(">"):
        #         seq += lines[j].strip()
        #         j += 1
        #     file_name = os.path.join("output", name + ".txt")
        #     with open(file_name, "a") as f:
        #         f.write(seq + "\n")
        # name = lines[i].strip()[1:]
        if name in sample_names:
            seq = ""
            j = i + 1
            while j < len(lines) and not lines[j].startswith(">"):
                seq += lines[j].strip()
                j += 1
                # print(seq)
            file_name = os.path.join("output", name + ".txt")
            
            with open(file_name, "a") as f:
                f.write(name + "\n" + seq + "\n")