#!/usr/bin/env python3

import os
from glob import glob

from multiprocessing import Pool
from subprocess import check_call, CalledProcessError, Popen, PIPE, STDOUT

from termcolor import cprint


BASE = "LungCancer/10x/11_CNV/LeMoNe_v2.5"
JAVA = "lemone.jar"
MATLAB = "matlab"
CLUSTER = "lemone_select_clusters"

__dir__ = os.path.abspath(os.path.dirname(__file__))

QUITE = True


def call(cmd, cwd=__dir__):
    if QUITE:
        w = open(os.devnull, "w+")
        check_call(cmd, shell=True, stderr=w, stdout=w, timeout=600, cwd=cwd)
    else:
        cprint(cmd, "white", "on_blue")
        check_call(cmd, shell=True, timeout=600, cwd=cwd)


def ganesh(input_file, output_dir, n_jobs=10):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    input_file = os.path.abspath(input_file)
    output_dir = os.path.abspath(output_dir)
    
    num = 0
    with open(input_file) as r:
        for line in r:
           num += 1 
           
    if not num // 2:
        raise ValueError("%s only %d rows" % (input_file, num))
    
    cmds = []
    for i in range(20):
        out_name = os.path.join(output_dir, f"net{i}")
        cmds.append(f"java -Xmx1000m -jar {BASE}/{JAVA}  -task ganesh -data_dir {os.path.dirname(input_file)} -data_file {os.path.basename(input_file)} -output_file {out_name} -init_num_clust 2500")
        
    with Pool(n_jobs) as p:
        p.map(call, cmds)
        
        
def fuzzy(input_dir, output_file):
    
    call(f"java -Xmx1000m -jar {BASE}/{JAVA} -task fuzzy -data_dir {input_dir} -output_file {output_file}")
    

def matlab(fuzzy, output_dir):
    cmd = f"matlab -r \"addpath('{BASE}/{MATLAB}');[s,data,idx] = makepairwiseprob('{fuzzy}');[s,n1,n2,nc] = pcutoff(s,0.2);p = matrixClustSym(s);[xopt, qopt,n1,n2,nc] = optCutoff(p,s);writeclusters(qopt,idx, '{output_dir}/tight_clusters.txt');disp('job done');quit;\""
    
    p = Popen(cmd, stdout=PIPE, stderr=STDOUT, shell=True)
    
    while True:
        try:
            line = p.stdout.readline().decode("utf-8").strip()

            # print(line)
            if "Error in optCutof" in line:
                print("%s - %s" % (cmd, line))
                p.terminate()
                return
            elif "job done" in line:
                print("%s - %s" % (cmd, line))
                break
        except Exception:
            p.kill()
        
    call(f"perl {BASE}/{CLUSTER} 10 < {output_dir}/tight_clusters.txt > {output_dir}/data/tc4.txt")
    
def regulator(input_dir, tf_mtx, tf_list, out_file, cluster):
    call(f"java -Xmx1000m -jar {BASE}/{JAVA} -task regulators -data_dir {input_dir} -data_file {tf_mtx} -reg_file {tf_list} -cluster_file {cluster} -output_file {out_file} -num_clust 200")
    

def print_reg(input_file, output_file):
    call(f"java -Xmx1000m -jar {BASE}/{JAVA} -task print_regulators -data_dir {input_file} -output_file {output_file} -all_regulators")
    
def figures(input_file):
    
    with open(os.path.join(input_file, "all_genes_map"), "w+") as w:
        with open(os.path.join(input_file, "data/tf.txt")) as r:
            for line in r:
                line = "\t".join(line.split()[:2])
                w.write(line + "\n")
        
    data = []        
    with open(os.path.join(input_file, "reg_all.txt")) as r:
        for line in r:
            line = line.split()
            line[-1] = float(line[-1])
            data.append(line)
    
    data = sorted(data, key=lambda x: x[-1], reverse=True)
    
    data = data[:round(len(data) * .1)]
    
    with open(os.path.join(input_file, "reg_top.txt"), "w+") as w:
        for line in data:
            w.write("\t".join([str(x) for x in line]) + "\n")
            
    for x in glob(os.path.join(input_file, "*.eps")):
        os.remove(x)
    
    call(f"java -Xmx1000m -jar {BASE}/{JAVA} -task figures -data_dir {input_file}/reg_xml/ -top_regulators reg_top.txt -map_file all_genes_map", cwd=input_file)


def pipeline(args):
    
    if os.path.exists(args["xml"]):
        fuzzy(args["xml"], args["fuzzy_file"])
        cprint("%s fuzzy finished" % args["key"], "white", "on_green")
        
    if os.path.exists(args["fuzzy_file"] + ".txt"):
        matlab(args["fuzzy_file"] + ".txt", args["indir"])
        cprint("%s matlab finished" % args["key"], "white", "on_green")

    if os.path.exists(args["tf"]) and os.path.exists(args["tf_list"]) and os.path.exists(args["indir"] + "/data/tc4.txt"):
        if not os.path.exists(os.path.join(args["indir"], "reg_xml")):
            os.makedirs(os.path.join(args["indir"], "reg_xml"))

        regulator(
            args["indir"], 
            "data/" + os.path.basename(args["tf"]), 
            "data/" + os.path.basename(args["tf_list"]), 
            os.path.join(args["indir"], "reg_xml/reg"),
            cluster="data/tc4.txt"
        )

        print_reg(os.path.join(args["indir"], "reg_xml/"), os.path.join(args["indir"], "reg_all.txt"))
        figures(args["indir"])
        
        cprint("%s regulator finished" % args["key"], "white", "on_green")


def main(input_dir, n_jobs=20):
    
    cells = [os.path.join(input_dir, x) for x in os.listdir(input_dir)]
    cells = [x for x in cells if os.path.isdir(x)]
    
    cmds = []
    
    global QUITE 
    QUIET = False
    
    indir = os.path.join(input_dir, "LeMoNe")

    rnas = {}
    tf_names = {}
    tfs = {}

    if os.path.exists(indir):
        for i in os.listdir(indir):
            if os.path.exists(os.path.join(indir, i, "data/rna.txt")):
                rnas[i] = os.path.join(indir, i, "data/rna.txt")
            
            if os.path.exists(os.path.join(indir, i, "data/tf.txt")):
                tfs[i] = os.path.join(indir, i, "data/tf.txt")
                
            if os.path.exists(os.path.join(indir, i, "data/tf_list.txt")):
                tf_names[i] = os.path.join(indir, i, "data/tf_list.txt")
    
    for i in glob(os.path.join(input_dir, "rna_*.txt")):
        key = os.path.basename(i).replace(".txt", "").replace("rna_", "")
        rnas[key] = i
        
    
    for i in glob(os.path.join(input_dir, "tf_name_*.txt")):
        key = os.path.basename(i).replace(".txt", "").replace("tf_name_", "")
        tf_names[key] = i
        
    
    for i in glob(os.path.join(input_dir, "tf_*.txt")):
        if i not in tf_names.values():
            key = os.path.basename(i).replace(".txt", "").replace("tf_", "")
            tfs[key] = i
            
    for key in set(rnas.keys()) & set(tf_names.keys()) & set(tfs.keys()):
        print(key)
        # 整理文件
        rna = rnas[key]
        tf = tfs[key]
        tf_list = tf_names[key]
        indir = os.path.join(input_dir, "LeMoNe", key)
        xml = os.path.join(indir, "xml")
        
        data = os.path.join(indir, "data")
        if not os.path.exists(data):
            os.makedirs(data)
        
        if os.path.exists(rna) and rna != os.path.join(data, "rna.txt"):
            os.rename(rna, os.path.join(data, "rna.txt"))
        rna = os.path.join(data, "rna.txt") 
        if os.path.exists(tf) and tf != os.path.join(data, "tf.txt"):
            os.rename(tf, os.path.join(data, "tf.txt"))
        tf = os.path.join(data, "tf.txt")
        if os.path.exists(tf_list) and tf_list != os.path.join(data, "tf_list.txt"):
            os.rename(tf_list, os.path.join(data, "tf_list.txt"))
        tf_list = os.path.join(data, "tf_list.txt")
            
        fuzzy_file = os.path.join(indir, "fuzzy")
        
        if os.path.exists(rna):
            # print("ganesh")
            try:
                ganesh(os.path.join(data, "rna.txt"), xml, n_jobs)
                if not os.path.exists(xml):
                    os.makedirs(xml)
                    
                cmds.append({
                    "indir": indir, "data": data, 
                    "tf": tf, "tf_list": tf_list, 
                    "rna": rna, "xml": xml,
                    "fuzzy_file": fuzzy_file,
                    "key": key
                })
                
                pipeline(cmds[-1])
            except (CalledProcessError, ValueError) as err:
                print(err)
    
    # for cell in sorted(cells):
    #     for disease in ["LUAD", "LUSC"]:
    #         for key in ["Malignant", "Non_Malignant", "Mix"]:
    #             # 指定文件
    #             rna = os.path.join(cell, disease, "rna_" + key + ".txt")
                
    #             tf = os.path.join(cell, disease, "tf_" + key + ".txt")
    #             tf_list = os.path.join(cell, disease, "tf_name_" + key + ".txt")
                
    #             if key == "Non_Malignent":
    #                 key = "Non_Malignant"
                    
    #             # print(cell, disease, key)
                
    #             # 整理文件
    #             indir = os.path.dirname(rna)
    #             indir = os.path.join(indir, "LeMoNe", key)
    #             xml = os.path.join(indir, "xml")
   
    #             data = os.path.join(indir, "data")
    #             if not os.path.exists(data):
    #                 os.makedirs(data)
                
    #             if os.path.exists(rna):
    #                 os.rename(rna, os.path.join(data, "rna.txt"))
    #             rna = os.path.join(data, "rna.txt") 
    #             if os.path.exists(tf):
    #                 os.rename(tf, os.path.join(data, "tf.txt"))
    #             tf = os.path.join(data, "tf.txt")
    #             if os.path.exists(tf_list):
    #                 os.rename(tf_list, os.path.join(data, "tf_list.txt"))
    #             tf_list = os.path.join(data, "tf_list.txt")
                    
    #             fuzzy_file = os.path.join(indir, "fuzzy")
                
    #             if os.path.exists(rna):
    #                 # print("ganesh")
    #                 try:
    #                     ganesh(os.path.join(data, "rna.txt"), xml, n_jobs)
    #                     if not os.path.exists(xml):
    #                         os.makedirs(xml)
                            
    #                     cmds.append({
    #                         "indir": indir, "data": data, 
    #                         "tf": tf, "tf_list": tf_list, 
    #                         "rna": rna, "xml": xml,
    #                         "fuzzy_file": fuzzy_file,
    #                         "key": "%s-%s-%s" % (os.path.basename(cell), disease, key)
    #                     })
                        
    #                     # pipeline(cmds[-1])
    #                 except (CalledProcessError, ValueError) as err:
    #                     print(err)
                                 
    # with Pool(n_jobs) as p:
    #     p.map(pipeline, cmds)
        
        
if __name__ == '__main__':
    from fire import Fire
    Fire(main)