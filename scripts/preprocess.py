import vcf
import sys
import os
import random

cf = sys.argv[1] # input chromosome vcf file
od = sys.argv[2] # output prefix for pseudo scaffolds
ns = int(sys.argv[3]) # number of scaffolds
nl = int(sys.argv[4]) # number of markers
cl = int(sys.argv[5]) # chromosome length
minor = float(sys.argv[6]) # loci minor allele ratio cutoff
gap = int(sys.argv[7]) # gap between scaffolds

# temp
minor = 0

def get_gt(sample):
    return sample.gt_type

def check_minor(record, cutoff):
    A = 0
    a = 0
    cnt = 0
    for sample in record.samples:
        if get_gt(sample) == 0:
            A += 2
            cnt += 1
        elif get_gt(sample) == 1:
            A += 1
            a += 1
            cnt += 1
        elif get_gt(sample) == 2:
            a += 2
            cnt += 1
        else:
            return False
    if float(A) / float(2 * cnt) < cutoff or float(a) / float(2 * cnt) < cutoff:
        return False
    return True

# vcf_reader = vcf.Reader(open(cf, 'r'))
vcf_reader = vcf.Reader(filename=cf)
step  = cl / ns
temps = range(step, cl, step)[:ns-1]
breaks = []
for t in temps:
    breaks.append(random.randint(-(step/10), +(step/10)) + t)
breaks.append(cl)
spos = 0
epos = breaks[0]
cur  = 0

if not os.path.isdir(od):
    os.mkdir(od)

chrom = vcf_reader.next().CHROM

logf1 = vcf.Writer(open(os.path.join(od, "snps.info.1"), 'w'), vcf_reader) # snp info for */*.snps.1
logf2 = vcf.Writer(open(os.path.join(od, "snps.info.2"), 'w'), vcf_reader) # snp info for */*.snps.2

while 1:
    if not os.path.isdir(os.path.join(od, "scaffold_{}".format(cur))):
        os.mkdir(os.path.join(od, "scaffold_{}".format(cur)))
    # snps.1
    front_snps = {}
    for record in [r for r in vcf_reader.fetch(chrom, spos, epos) if check_minor(r, minor)][:nl]:
        logf1.write_record(record)
        for sample in record.samples:
            if sample.sample not in front_snps:
                front_snps[sample.sample] = [get_gt(sample)]
            else:
                front_snps[sample.sample].append(get_gt(sample))
    for ind, ss in front_snps.items():
        with open(os.path.join(od, "scaffold_{}".format(cur), "{}.snps.1".format(ind)), "w") as ofs:
            for i, s in enumerate(ss):
                if i < nl - 1:
                    ofs.write("{} ".format(s))
                else:
                    ofs.write("{}\n".format(s))
    # snps.2
    back_snps = {}
    sp = epos - 5000
    while len([r for r in vcf_reader.fetch(chrom, sp, epos) if check_minor(r, minor)]) < nl:
        sp -= 5000
        if sp < 0:
            sp = 0
    for record in [r for r in vcf_reader.fetch(chrom, sp, epos)if check_minor(r, minor)][-nl:]:
        logf2.write_record(record)
        for sample in record.samples:
            if sample.sample not in back_snps:
                back_snps[sample.sample] = [get_gt(sample)]
            else:
                back_snps[sample.sample].append(get_gt(sample))
    for ind, ss in back_snps.items():
        with open(os.path.join(od, "scaffold_{}".format(cur), "{}.snps.2".format(ind)), "w") as ofs:
            for i, s in enumerate(ss):
                if i < nl - 1:
                    ofs.write("{} ".format(s))
                else:
                    ofs.write("{}\n".format(s))
    cur += 1
    if cur == len(breaks):
        break
    spos = breaks[cur-1] + gap
    epos = breaks[cur]

logf1.close()
logf2.close()
