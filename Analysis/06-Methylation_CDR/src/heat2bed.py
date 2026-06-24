import argparse
import sys
import textwrap
import gzip
import math
import seaborn as sns

parser = argparse.ArgumentParser(prog='',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
heat to bed
'''))

parser.add_argument('-bed', metavar='str', help='get bed')
parser.add_argument('-bedlist', metavar='str', help='get bed list')
parser.add_argument('-data', metavar='str', help='get data')
parser.add_argument('-pal', metavar='str', default= "Spectral", help='color pallele')
parser.add_argument('-color', default=False, action='store_true', help='Add color')
parser.add_argument('-sub', default=False, action='store_true', help='get sub heat')
parser.add_argument('-eq', default=False, action='store_true', help='get average heat')
parser.add_argument('-num', metavar='<INT>', type=int, default=100, help='seperate number')
parser.add_argument('-bd', metavar='<INT>', type=int, default=4, help='local number')
parser.add_argument('-intval', metavar='<INT>', type=int, default=5000, help='interval number')

if len(sys.argv) == 1:
    parser.print_help()
    parser.exit()

args = parser.parse_args()

def read_bed(fh, bound):
    out = {}
    basecount = 0
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        lining = line.split("\t")
        if line.startswith("#"):
            continue
        org_chrom = lining[0]
        qst = int(lining[1])
        qed = int(lining[2])
        tst = int(lining[4])
        ted = int(lining[5])
        score = float(lining[6])
        if qst == tst:
            basecount = basecount + 1
        if (tst - qst) / args.intval != 0 and (tst - qst) / args.intval <= bound and (tst - qst) / args.intval >= -bound:
            if (org_chrom, qst, qed) in out:
                out[(org_chrom, qst, qed)].append((score, tst, ted))
            else:
                out[(org_chrom, qst, qed)] = [(score, tst, ted)]
    return out, basecount

def outbed(out):
    for info in out.keys():
        org_chrom = info[0]
        scores = [i[0] for i in out[info]]
        qst = info[1]
        qed = info[2]
        chrom = org_chrom.split(":")[0]
        org_s = int(org_chrom.split(":")[1].split("-")[0])
        org_e = int(org_chrom.split(":")[1].split("-")[1])
        rst = org_s + qst
        red = org_s + qed
        out_score = sum(scores) / len(scores)
        print("%s\t%i\t%i\theat\t%.10f" % (chrom, rst, red, out_score))

def bounds_getter(numbers, intv):
    cutnum = len(numbers) / intv
    bounds = []
    i = 0
    while i * cutnum < len(numbers):
        k = int(i * cutnum)
        if k == 0:
            i = i + 1
            continue
        bounds.append(numbers[k])
        i = i + 1
    return bounds

def equal_getter(numbers, intv):
    i = 0
    nl = []
    minm = min(numbers)
    maxm = max(numbers)
    iv = (maxm - minm) / intv
    while i < intv:
        nl.append(minm + iv * i)
        i = i + 1
    nl.append(maxm)
    return nl
        

def split_data(fh, intv):
    numbers = []
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        numbers.append(float(line))
    numbers.sort()
    if args.eq:
        bounds = equal_getter(numbers, intv)
    else:
        bounds = bounds_getter(numbers, intv)
   
    return bounds

def getnum(id, bounds):
    i = 0
    while i < len(bounds):
        if id <= bounds[i]:
            break
        i = i + 1
    return i
    

def addcolor(fh, bounds):
    num = len(bounds) + 1
    color_l = sns.color_palette(args.pal, num).as_hex()
    color_l.reverse()
    for line in fh:
        line = line.rstrip()
        line = line.replace("\n", "")
        lining = line.split("\t")
        id = float(lining[4])
        index = getnum(id, bounds)
        color = color_l[index]
        print("%s\t+\t%s\t%s\t%s" % ("\t".join(lining), lining[1], lining[2], color))

def addcolor_v2(fh, bounds, outname):
    num = len(bounds) + 1
    color_l = sns.color_palette(args.pal, num).as_hex()
    color_l.reverse()
    out, basecount = read_bed(fh, bd)
    outseq = []
    for info in out.keys():
        qst = info[1]
        qed = info[2]
        o = out[info]
        chrom, orst = get_st(info[0])
        for qu in o:
            id = qu[0]
            tst = qu[1]
            ted = qu[2]
            if (tst - qst) / args.intval >= 0 and (tst - qst) / args.intval <= args.bd:
                index = getnum(id, bounds)
                color = color_l[index]
                outseq.append("%s\t%i\t%i\t%s\t%i\t%i\t%.10f\t%s" % (chrom, qst + orst, qed + orst, chrom, tst + orst, ted + orst, id, color))
                #print("%s\t%i\t%i\t%s\t%i\t%i\t%.10f\t%s" % (chrom, qst, qed, chrom, tst, ted, id, color))
    return outseq

def get_st(chrom):
    if ":" in chrom:
        ch = chrom.split(":")[0]
        st = int(chrom.split(":")[1].split("-")[0])
        return ch, st
    else:
        return chrom, 0

if args.color:
    with open(args.data, "r") as fh:
        bounds = split_data(fh, args.num)

if args.color:
    with open(args.bed, "rt") as fh:
        addcolor(fh, bounds)
else:
    if args.sub:
        bd = math.inf
    else:
        bd = args.bd
        
    if args.sub:
        numbers = []
        with open(args.data, "rt") as fh:
            bounds = split_data(fh, args.num)
        filelist = []
        with open(args.bedlist, "r") as fh:
            for line in fh:
                line = line.rstrip()
                line = line.replace("\n", "")
                filelist.append(line)
        for filename in filelist:
            naming = filename.split("/")[-1]
            outname = naming.split(".")[0] 
            with gzip.open(filename, "rt") as fh:
                outseq = addcolor_v2(fh, bounds, "%s.heat" % (outname))
            with open(outname, "w") as out:
                out.write("\n".join(outseq))
    else:
        with gzip.open(args.bed, "rt") as fh:
            out, basecount = read_bed(fh, bd)
            outbed(out)

