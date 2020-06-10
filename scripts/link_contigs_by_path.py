from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import itertools
import os,sys

def read_path(input_file):
	with open(input_file) as f:
		paths = [li.rstrip('\n') for li in f.readlines()]
	return paths

def build_contig_names_dic(contigs_names):
	l = len(contigs_names)*2
	odd_idx = list(range(0,l,2))
	even_idx = [i+1 for i in odd_idx]
	idx = odd_idx + even_idx
	contigs_names = contigs_names + contigs_names
	contig_names_dic = dict(zip(idx, contigs_names))
	return contig_names_dic

#contig_names_dic : idx:name
def get_contigname_and_order(paths, contig_names_dic):
	contig_names_paths = {}
	count = 0
	for path in paths:
		path = list(map(int, path.split()))
		contig_end_num = len(path)
		odd_idxs = [i for i in range(contig_end_num) if i%2==0]
		path_sim = [path[odd_idx] for odd_idx in odd_idxs]
		contig_fullname_data = []
		contig_name = [contig_names_dic[path[i]] for i in odd_idxs]
		contig_orientation = list(map(int, [path[i]<path[i+1] for i in odd_idxs]))
		contig_names_paths[count] = (contig_name, contig_orientation)
		count += 1
	return contig_names_paths

def record_filtered_iterator(output_name, output_iterator):
	with open(output_name, "w") as output_handle:
		SeqIO.write(output_iterator, output_handle, "fasta")

def record_filtered_iterator_v1(output_iterator, output_handle):
	SeqIO.write(output_iterator, output_handle, "fasta")


def assign_seq_to_path(contig_names_paths, seq_id):
	for key in contig_names_paths:
		contig_name = contig_names_paths[key][0]
		if seq_id in contig_name: 
			contig_orientation = contig_names_paths[key][1]
			idx = contig_name.index(seq_id)
			orientation = contig_orientation[idx]
			return key, idx, orientation
	return -1,-1,-1 #sequences that fail

def build_empty_iterator(contig_names_paths):
	initial_cycle_iterator = {}
	for key in contig_names_paths:
		contig_names = contig_names_paths[key][0]
		initial_cycle_iterator[key] = [0 for i in range(len(contig_names))]
	return initial_cycle_iterator

def build_Seq_record(seq, seq_id):
	Seq_record = SeqRecord(Seq(seq), id=seq_id, name='', description='')
	return Seq_record

def build_iterator(contig_names_paths, record_iterator):
	initial_cycle_iterator = build_empty_iterator(contig_names_paths)
	for record in record_iterator:
		cycle_num, seq_idx, seq_orientation = assign_seq_to_path(contig_names_paths, record.id)
		if cycle_num == -1: 
			continue
		if seq_orientation==1:
			record.seq = record.seq[::-1]
		initial_cycle_iterator[cycle_num][seq_idx] = record
	return initial_cycle_iterator

def construct_sequence(initial_cycle_iterator):
	cycle_finished_iterator = []
	l = len(initial_cycle_iterator.keys())
	for key in range(l):
		gap = ''.join(['N' for i in range(100)])
		cycle_seq_raw = initial_cycle_iterator[key]
		cycle_seq = gap.join([str(contig.seq) for contig in initial_cycle_iterator[key]])
		cycle_id = 'cycle_'+str(key)
		Cycle_record = build_Seq_record(cycle_seq, cycle_id)
		cycle_finished_iterator.append(Cycle_record)
	return cycle_finished_iterator

def generate_iterator_for_a_path(contig_dic, record_iterator):
	contig_names = contig_dic.keys()
	output_iterator = []
	for record in record_iterator:
		if record.id in contig_names:  
			orientation = contig_dic[record.id]
			if orientation==1:
				record.seq = record.seq[::-1]
				output_iterator.append(record)
	return output_iterator

def generate_sequences_for_a_path_v1(record_iterator, contig_names, contig_orientation ):
	seq = ['' for i in range(len(contig_names))]
	output_iterator = []
	for record in record_iterator:
		if record.id in contig_names:
			idx = contig_names.index(record.id)
			orientation = contig_orientation[idx]
			if orientation==1:
				# record = record.seq[::-1]
				seq[idx] = SeqRecord(record.seq[::-1], id=record.id, name='', description='')
			else:
				seq[idx] = record
	gap = ''.join(['N' for i in range(100)])
	return seq

def write_per_path(contig_names_paths, input_fasta, contig_names_dic):
	l = len(contig_names_paths)
	xx = []
	with open(input_fasta, "r") as input_handle:
		x = SeqIO.parse(input_handle, "fasta")
		xx.append(itertools.tee(x, l))
		for i in range(l):
			record_iterator = xx[0][i]
			contig_name, contig_orientation = contig_names_paths[i]
			output_name = 'cycle_'+str(i)+'.fasta'
			output_iterator = generate_iterator_for_a_path(contig_dic, record_iterator)
			seq = generate_sequences_for_a_path(record_iterator, contig_name, contig_orientation)
			record_filtered_iterator(output_name, output_iterator)

def write_per_path_v1(contig_names_paths, input_fasta, contig_names_dic, output_name):
	l = len(contig_names_paths) # cycle num
	with open(input_fasta, "r") as input_handle, open(output_name, "w") as output_handle:
		record_iterator = SeqIO.parse(input_handle, "fasta")
		initial_cycle_iterator = build_iterator(contig_names_paths, record_iterator)
		cycle_finished_iterator = construct_sequence(initial_cycle_iterator)
		SeqIO.write(cycle_finished_iterator, output_handle, "fasta")

def main():
	# path_inputfile = sys.argv[1]
	# contigs_names_file = sys.argv[2]
	# input_fasta = sys.argv[3]
	# output_fasta = sys.argv[4]
	parser.add_argument('-paths', help='each line contains a path')
	parser.add_argument('-contigs', help='names of contigs used')
	parser.add_argument('-i', help='input fasta filename')
	parser.add_argument('-o', help='outout fasta filename')
	args = parser.parse_args()

	input_fasta = args.i
	output_fasta = args.o
	paths = read_path(args.paths)
	contigs_names = read_path(args.contigs)
	contig_names_dic = build_contig_names_dic(contigs_names)
	contig_names_paths = get_contigname_and_order(paths, contig_names_dic)
	write_per_path_v1(contig_names_paths, input_fasta, contig_names_dic, output_fasta)

main()


