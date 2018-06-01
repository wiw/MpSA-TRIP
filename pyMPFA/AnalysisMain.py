#C:\Python27\python.exe
#!/usr/bin/env python
# encoding: utf-8

import os, pickle, re, copy, json

CONFIG = {
    "experiment_dir": "/home/anton/backup/input/trip/RUN_2018-05-10/results/bcRead_1__bcmutProb_80/Lib_33-40",
    "content": ("control_e", "control_n", "expression", "normalization"),
    "exception": {
      "experiment_dir": "/home/anton/backup/MCBDataCloud/LCD TRIP/Experiment_2018-05-10/Results/bcRead_3__bcmutProb_80/Lib_33-40",
      "content": ("control_m", "mapping")
    },
    "control": {
        "wt-bc1": "TTCCAAGTGCAGGTTAGGCG",
        "wt-bc2": "TGTGTACGGCTTGCTCTCAA",
        "deltaC-bc3": "GAGCCCGGATCCACTCCAAG",
        "deltaC-bc4": "TGTCACGTCAGCTAACCCAC"
    },
    "output_for_R": "/home/anton/backup/input/trip/RUN_2018-05-10/results/statistics"
}

def dump_to_json(obj, output_file):
    with open(output_file, "w") as handle:
        json.dump(obj, handle)

def SaveDictToPy(dictVar, filename):
    with open(filename + ".py", "wb") as handle:
        for k, v in dictVar.items():
            if type(v) == str:
                handle.write(str(k) + " = '" + str(v) + "'\n")
            else:
                handle.write(str(k) + " = " + str(v) + "\n")

def load_pickle(CONFIG):
    all_data = {}
    control_data = {}
    def pickle_opener(i):
        with open(i, "rb") as handle:
            unserialized_data = pickle.load(handle)
        return unserialized_data
    for item in CONFIG["content"]:
        regex0 = re.compile(".*" + item)
        path_to = filter(regex0.match, os.listdir(CONFIG["experiment_dir"]))[0]
        path_to = os.path.join(CONFIG["experiment_dir"], path_to, "Dump")
        if item in ("control_m", "mapping"):
            what_is_pickle = "resultDict"
        else:
            what_is_pickle = "bcDict"
        regex1 = re.compile(".*" + what_is_pickle + ".*")
        pickles = filter(regex1.match, os.listdir(path_to))
        for f in pickles:
            filename, ext = os.path.splitext(f)
            fst, snd, repl, dct = filename.split("_")
            del(fst, snd, ext, dct)
            filename = item + "-" + repl[1:]
            f = os.path.join(path_to, f)
            if re.search("control", item) is not None:
                control_data[filename] = pickle_opener(f)
            else:
                all_data[filename] = pickle_opener(f)
    if CONFIG.get("exception") is not None and bool(CONFIG["exception"]):
        exception_all_data, exception_control_data = load_pickle(CONFIG["exception"])
        all_data.update(exception_all_data)
        control_data.update(exception_control_data)
    return all_data, control_data

def get_control_count(CONFIG, control):
    compiled_control = {}
    for item in control:
        compiled_control.setdefault(item, copy.deepcopy(CONFIG["control"]))
        for alias, bc in CONFIG["control"].items():
            counter = sum([count for mutated_bc, count in control[item][bc]])
            compiled_control[item][alias] = counter
    return compiled_control

# >>> data.keys()
# ['expression-1', 'normalization-2', 'normalization-1', 'expression-2', 'mapping-1', 'mapping-2']



def align_experiments(CONFIG, data):
    pass
    
def main(CONFIG):
    data, control = load_pickle(CONFIG)
    output_control = get_control_count(CONFIG, control)
    if not os.path.exists(CONFIG["output_for_R"]): os.makedirs(CONFIG["output_for_R"])
    dump_to_json(output_control, os.path.join(CONFIG["output_for_R"], "control.json"))


if __name__ == '__main__':
    main()