import vcf
import sys
import os

sf = sys.argv[1]
# sc = sys.argv[2]
n  = int(sys.argv[2])
od = sys.argv[3]

# cands = []
# with open(sc) as ifs:
#     for line in ifs.readlines():
#         cands.append(line.replace("\n", ""))

vcf_reader = vcf.Reader(open(sf, 'r'))

cur = ""
snps = {}
pos = 0
parsed_chrom = []
# valid = 0
merge = {}
for record in vcf_reader:
    if record.CHROM != cur:
        if cur != "": # and valid == 1:
            # if not os.path.isdir(od):
            #     os.mkdir(od)
            # for ind, ss in snps.items():
            #     if len(ss) < 2 * n:
            #         continue
            #     # with open(os.path.join(od, "{}.snps.1".format(ind)), "w") as ofs:
            #     #     for idx, s in enumerate(ss[:n]):
            #     #         if idx < n - 1:
            #     #             ofs.write("{} ".format(s))
            #     #         else:
            #     #             ofs.write("{}".format(s))
            #     # with open(os.path.join(od, "{}.snps.2".format(ind)), "w") as ofs:
            #     #     for idx, s in enumerate(ss[-n:]):
            #     #         if idx < n - 1:
            #     #             ofs.write("{} ".format(s))
            #     #         else:
            #     #             ofs.write("{}".format(s))
            merge[cur] = snps
        assert(record.CHROM not in parsed_chrom)
        parsed_chrom.append(record.CHROM)
        cur = record.CHROM
        pos = 0
        # get genotype label
        # snps = [genotype]
        snps = {}
        for sample in record.samples:
            if sample.gt_type == None:
                snps[sample.sample] = [-1]
            else:
                snps[sample.sample] = [sample.gt_type]
        # if cur in cands:
        #     valid = 1
        # else:
        #     valid = 0
    else:
        # if valid == 0:
        #     continue
        # snps.append(genotype)
        for sample in record.samples:
            if sample.gt_type == None:
                snps[sample.sample].append(-1)
            else:
                snps[sample.sample].append(sample.gt_type)
        assert(record.POS > pos)
        pos = record.POS
# # if cur in cands:
# if not os.path.isdir(od):
#     os.mkdir(od)
# for ind, ss in snps.items():
#     # if len(ss) < 2 * n:
#     #     continue
#     with open(os.path.join(od, "{}.snps.1".format(ind)), "w") as ofs:
#         for idx, s in enumerate(ss[:n]):
#             if idx < n - 1:
#                 ofs.write("{} ".format(s))
#             else:
#                 ofs.write("{}".format(s))
#     with open(os.path.join(od, "{}.snps.2".format(ind)), "w") as ofs:
#         for idx, s in enumerate(ss[-n:]):
#             if idx < n - 1:
#                 ofs.write("{} ".format(s))
#             else:
#                 ofs.write("{}".format(s))
if cur not in merge:
    merge[cur] = snps

# output
scaf_out = open(os.path.join(od, 'scaffold.list'), 'w')
samp_out = open(os.path.join(od, 'sample.list'), 'w')
samp_list = list(merge.values()[0].keys())
scaf_list = []
handlers = {}
for sample in samp_list:
    handlers[sample] = open(os.path.join(od, "{}.snps".format(sample)), 'w')
for scaf, sss in merge.items():
    no = len(sss.values()[0])
    mask = [1] * no
    for i in range(no):
        tot = len(sss.keys())
        mis = 0
        for ind, ss in sss.items():
            if ss[i] == -1:
                mis += 1
        if mis > tot * 0.3:
            mask[i] = 0
    if sum(mask) < 2 * n:
        continue
    scaf_list.append(scaf)
    frnt_idxs = []
    back_idxs = []
    for i, m in enumerate(mask):
        if m == 1:
            frnt_idxs.append(i)
            if len(frnt_idxs) == n:
                break
    for i, m in enumerate(mask[::-1]):
        if m == 1:
            back_idxs.append(i)
            if len(back_idxs) == n:
                break
    for ind, ss in sss.items():
        for idx, ii in enumerate(frnt_idxs):
            if idx < n - 1:
                handlers[ind].write("{} ".format(ss[ii]))
            else:
                handlers[ind].write("{}\n".format(ss[ii]))
        for idx, ii in enumerate(back_idxs):
            if idx < n - 1:
                handlers[ind].write("{} ".format(ss[ii]))
            else:
                handlers[ind].write("{}\n".format(ss[ii]))
for sample in samp_list:
    handlers[sample].close()
for samp in samp_list:
    samp_out.write("{}\n".format(samp))
for scaf in scaf_list:
    scaf_out.write("{}\n".format(scaf))
scaf_out.close()
samp_out.close()
