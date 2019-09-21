import sys
import argparse
ap=argparse.ArgumentParser()
ap.add_argument("--input", help="input file.", default="/dev/stdin")




args=ap.parse_args()
inFile = open(args.input)

counts={}
#count=[]
c=0

for line in inFile:
    line=line.rstrip()+"\n"
    if c==0:
        state_prev=line
        c=1
        continue
    state_cur=line

    if state_prev not in counts:
        counts[state_prev]={}
    if state_cur not in counts[state_prev]:
        counts[state_prev][state_cur]=0

    counts[state_prev][state_cur]=counts[state_prev][state_cur]+1
    state_prev=state_cur


for i in counts.keys():
    j=counts[i].keys()
    for jx in j:
        i_s=i.rstrip()
        jx_s=jx.rstrip()
        sys.stdout.write("%s\t%s\t%s\n" % (i_s,jx_s,counts[i][jx]) )
