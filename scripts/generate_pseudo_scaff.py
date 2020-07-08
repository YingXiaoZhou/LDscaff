import vcf
import sys
import os
import random

cf = sys.argv[1] # input chromosome vcf file
od = sys.argv[2] # output prefix for pseudo scaffolds
n  = int(sys.argv[3]) # number of scaffolds
cl = int(sys.argv[4]) # chromosome length
minor = float(sys.argv[5]) # loci minor ratio cutoff

# reset chromosome length by loci pos
vcf_reader = vcf.Reader(open(cf, 'r'))
for record in vcf_reader:
    cl = record.POS+1

vcf_reader = vcf.Reader(open(cf, 'r'))
step   = cl / n
temps = range(step, cl, step)
breaks = []
for t in temps:
    breaks.append(random.randint(-(step/10), +(step/10)) + t)
# breaks = sorted(random.sample(range(cl), n))
# breaks = sorted([random.randint(0, cl) for i in range(n)])
pos    = 0
epos   = breaks[0]
cur    = 0
idx    = 0

vcf_writer = vcf.Writer(open("{}_scaffold_{}.vcf".format(od, idx), "w"), vcf_reader)
for record in vcf_reader:
    if record.POS < epos:# breaks[cur]:
        record.POS = record.POS - pos
        vcf_writer.write_record(record)
    else:
        cur += 1
        idx += 1
        # assert(cur <= len(breaks))
        if cur < len(breaks):
            pos = breaks[cur-1]
            epos = breaks[cur]
        elif cur == len(breaks):
            pos = breaks[cur-1]
            epos = cl
        else:
            print "break index error, pls check"
            exit(1)
        vcf_writer.close()
        vcf_writer = vcf.Writer(open("{}_scaffold_{}.vcf".format(od, idx), "w"), vcf_reader)
        record.POS = record.POS - pos
        vcf_writer.write_record(record)
