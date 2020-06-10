import argparse

parser.add_argument('-paths', help='each line contains a path')
parser.add_argument('-contigs', help='names of contigs used')
parser.add_argument('-i', help='input fasta filename')
parser.add_argument('-o', help='outout fasta filename')
args = parser.parse_args()

print(args.i)
