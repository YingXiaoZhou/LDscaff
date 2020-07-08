import vcf
import sys

vf = sys.argv[1]
of = sys.argv[2]

vcf_reader = vcf.Reader(open(vf, 'r'))
with open(of, 'w') as ofs:
    for sample in vcf_reader.samples:
        ofs.write("{}\n".format(sample))
