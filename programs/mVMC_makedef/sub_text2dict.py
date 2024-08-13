import json
import os

# https://stackoverflow.com/questions/5608702/how-can-i-convert-a-string-to-either-int-or-float-with-priority-on-int
def int_or_float(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s

def make_text2dict():
    file_txt = "dat_input_makedef.txt"
    file_json = "dat_output_makedef.json"
    data_dict = {}
    for line in open(file_txt):
        li = line.strip()
        if not li.startswith("#"): # skip "^#" lines
            if line.rstrip(): # skip empty line
                # print(line.rstrip())
                key, value = line.strip().split()
                data_dict[key] = int_or_float(value)
    # print(data_dict)
    dir_def = data_dict["dir_def"]
    if not os.path.exists(dir_def):
        os.makedirs(dir_def)
    #with open(dir_def+"/"+file_json,"w") as json_file:
    with open(file_json,"w") as json_file:
        json.dump(data_dict,json_file)
    #with open(dir_def+"/"+file_json,"r") as json_file:
    with open(file_json,"r") as json_file:
        loaded_dict = json.load(json_file)
    # print(loaded_dict)
    # print(loaded_dict.keys())
    # print(list(loaded_dict))
    # print(loaded_dict.values())
    # print()
    # print(list(loaded_dict)[0])
    # print(loaded_dict[list(loaded_dict)[0]])
    # print(loaded_dict["Lx"])

if __name__ == "__main__":
    make_text2dict()
