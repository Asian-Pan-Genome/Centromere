import os
import glob


def get_chain_lst(hicat_dir):
    files = glob.glob("*.chain.xls")
    sorted_files = sorted(files, key=lambda x: int(x.split('.')[0]))
    return sorted_files


def get_block_info(chainfile):
    with open(chainfile, 'r') as inf:
        first_line = inf.readline().strip()
        block_info_lst = first_line[1:].split(',')
    return block_info_lst

def file_exists_and_not_empty(file_path):
    return os.path.isfile(file_path) and os.path.getsize(file_path) > 0

def haploid_chains_merge(haploid_hicat_dir, sample, hap):

    script_dir = os.getcwd()

    os.chdir(haploid_hicat_dir)
    sorted_files = get_chain_lst(haploid_hicat_dir)

    outf = open(f"HiCAT.hor.summary.xls","w")
    for chainfile in sorted_files:
        block_info_lst = get_block_info(chainfile)
        blockid = block_info_lst[0]
        bstart = block_info_lst[2]
        bend = block_info_lst[3]
        blen = block_info_lst[4]
        bnum = block_info_lst[5]
        hor_exist_flag = file_exists_and_not_empty(f"{blockid}.all_layer.xls.reorder.xls") 
        if hor_exist_flag: 
            with open(f"{blockid}.all_layer.xls.reorder.xls", "r") as inf:
                for line in inf:
                    line = line.strip()
                    outf.write(f"{sample}\t{hap}\t{blockid}\t{bstart}\t{bend}\t{blen}\t{bnum}\t{line}\n")
        else:
            with open(f"{blockid}.chain.xls") as inf:
                for line in inf:
                    if line.startswith('#'):
                        ichr = line.strip().split(',')[1]
            tmp = "\t".join(['Na']*8)
            outf.write(f"{sample}\t{hap}\t{blockid}\t{bstart}\t{bend}\t{blen}\t{bnum}\t{ichr}\t{tmp}\n")

    os.chdir(script_dir)

def all_sample_update():
    with open("rerun.xls", "r") as inf:
        for line in inf:
            tokens = line.strip().split('\t')
            sample=tokens[0]
            hap = tokens[1]
            haploid_hicat_dir = tokens[2]
            haploid_chains_merge(haploid_hicat_dir, sample, hap)
            print(f"{sample} {hap} done!")

if __name__ == "__main__":
    all_sample_update()
