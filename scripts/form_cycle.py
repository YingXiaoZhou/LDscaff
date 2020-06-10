import os,sys
import numpy as np

def read_scaffold_name(scaffold_list='scaf_rerun.list'):
	with open(scaffold_list) as f:
		lines = f.readlines()
	scaf_names = [li.rstrip('\n') for li in lines]
	return scaf_names

def read_contig_name_from_node(node_num, scaf_names):
	idx = node_num/2
	return scaf_names[idx]

def read_scaffold_length(scaffold_list='scaf_rerun.list', scafoold_length='scaf_rerun.length'):
	scaf_names = read_scaffold_name(scaffold_list)
	length_dic = {}
	with open(scafoold_length) as f:
		for li in f:
			items = li.rstrip('\n').split()
			length_dic[items[0]] = int(items[1])
	scaf_lengths = []
	for si in scaf_names:
		scaf_lengths.append(length_dic[si])
	return scaf_lengths

def get_length_of_a_scaffold(node_num):
	if node_num%2==0: #front side
		another_node = node_num+1
	else:
		another_node = node_num-1
	return another_node

def find_N50(lengths):
	length_sum = sum(lengths)
	length_sum_half = float(length_sum/2)
	lengths.sort(reverse=True)
	s = 0
	for i in lengths:
		s+=i
		if s>length_sum_half:
			return i
	print('check your code again')
	return None

def cycle_get_length(path, scaf_lengths):
	valid_nodes = path[:-1]
	power_idx = list(range(0,len(valid_nodes),2))
	power_idx = [i+1 for i in power_idx]
	nodes_sets = [valid_nodes[i] for i in power_idx]
	contig_lengths = []
	contig_lengths = [scaf_lengths[int(i/2)] for i in nodes_sets]
	length_super_scaffold = 100*(int(len(path)/2)-1) + sum(contig_lengths)
	return length_super_scaffold


def linear_get_length(path, scaf_lengths):
	valid_nodes = path
	power_idx = list(range(0,len(valid_nodes),2))
	power_idx = [i+1 for i in power_idx]
	nodes_sets = [valid_nodes[i] for i in power_idx]
	contig_lengths = [scaf_lengths[int(i/2)] for i in nodes_sets]
	length_super_scaffold = 100*(int(len(path)/2)-1) + sum(contig_lengths)
	return length_super_scaffold

def print_length_scaffold(length_super_scaffold_lengths, tag):
	for li in length_super_scaffold_lengths:
		print('Length: ', tag, li)

def get_super_scaffold_length(paths, scaf_lengths):
	super_scaffold_lengths = []
	for path in paths:
		l = len(path)
		if l%2==0:
			length_super_scaffold = linear_get_length(path, scaf_lengths)
		else:
			length_super_scaffold = cycle_get_length(path, scaf_lengths)
		super_scaffold_lengths.append(length_super_scaffold)
	return super_scaffold_lengths

def get_statistic(best_links_file):
	strength_dic = {}
	couple_dic = {}
	from_strength_to_links = {}
	with open(best_links_file) as f:
		for li in f:
			items = li.rstrip('\n').split()
			head = int(items[0])
			tail = int(items[1])
			if abs(head-tail)==1 and tail%2==1: 
				print('Head:{} Tail:{}'.format(head, tail))
				continue
			LD = float(items[2])
			strength_dic[head] = LD
			strength_dic[tail] = LD
			if not LD in from_strength_to_links:
				from_strength_to_links[LD] = [head, tail]
			else:
				from_strength_to_links[LD].append(head)
				from_strength_to_links[LD].append(tail)
			couple_dic[head] = tail
			couple_dic[tail] = head
	return strength_dic, couple_dic, from_strength_to_links

def get_another_node_of_a_scaffold(node_num):
	if node_num%2==0: #front side
		another_node = node_num+1
	else:
		another_node = node_num-1
	return another_node


def Get_strength_from_a_path(path, strength_dic):
	if path[0]==path[-1]:
		valid_nodes = path[:-1]
	else:
		valid_nodes = path
	valid_nodes = path[:-1]
	power_idx = list(range(0,len(valid_nodes),2))
	power_idx = [i+1 for i in power_idx]
	if path[0]!=path[-1]:
		power_idx = power_idx[:-1]
	nodes_sets = [valid_nodes[i] for i in power_idx]
	LD = [strength_dic[i] for i in nodes_sets]
	return LD

def break_a_cycle(path, strength_dic, LD):
	LD_idx = np.array(range(len(LD)))
	threshold = sorted(LD)[0]
	LD_array = np.array(LD)
	break_point = list(LD_idx[LD_array==threshold])[0]
	new_paths = [] 
	break_point = 2*(break_point+1)
	path1 = path[:break_point]
	path2 = path[break_point:]
	newpath = path2 + path1[1:]
	return newpath

def break_a_path(threshold, path, strength_dic, LD):
	LD_idx = np.array(range(len(LD)))
	LD_array = np.array(LD)
	break_points = LD_idx[LD_array<threshold]
	if len(break_points) < 1:
		return path, 0
	else:
		new_paths = [] 
		break_points = [2*(i+1) for i in break_points]
		break_points.insert(0,0)
		break_points.append(len(path))
		for i in range(len(break_points)-1):
			pos_start = break_points[i]
			pos_end = break_points[i+1]
			new_paths.append(path[pos_start:pos_end])
		return new_paths, 1

def break_each_cycle(paths, strength_dic):
	paths_after_break_cycles = []
	for pi in paths:
		LD = Get_strength_from_a_path(pi, strength_dic)
		path = break_a_cycle(pi, strength_dic, LD)	
		paths_after_break_cycles.append(path)
	return paths_after_break_cycles

def break_path_main(threshold, paths, strength_dic):
	new_paths = []	
	for pi in paths:
		LD = Get_strength_from_a_path(pi, strength_dic)
		new_path, tag = break_a_path(threshold, pi, strength_dic, LD)
		if tag == 0:
			new_paths.append(new_path)
		else:
			for ni in new_path:
				new_paths.append(ni)
	return new_paths

def get_test_threshold(LDs):
	all_links_strength = sorted(list(set(LDs)))
	test_thresholds = all_links_strength[:100]
	return test_thresholds

def form_cycle(couple_dic):
	#start from a node
	#add next node into path until it become a cycle
	#stops when all nodes are used
	nodes_used_record = []
	paths = []
	nodes_not_used_record = list(couple_dic.keys())

	while(len(nodes_not_used_record)>0):
		path = []
		start_node = nodes_not_used_record[0]
		path.append(start_node)
		start_node_sister_node = get_another_node_of_a_scaffold(start_node)
		path.insert(0, start_node_sister_node)
		nodes_not_used_record.remove(start_node)
		nodes_not_used_record.remove(start_node_sister_node)
		while (path[0]!=path[-1]):
			couple_node = couple_dic[start_node]
			path.append(couple_node)
			if not couple_node in nodes_not_used_record: 
				break
			sister_node = get_another_node_of_a_scaffold(couple_node)
			path.append(sister_node)
			nodes_not_used_record.remove(couple_node)
			nodes_not_used_record.remove(sister_node)			
			start_node = sister_node
		paths.append(path)
	return paths

def collect_LD_for_paths(paths, strength_dic):
	LDs = []
	count = 0
	for pi in paths:
		LD = Get_strength_from_a_path(pi, strength_dic)
		LDs += LD
	return LDs

def calculate_length(paths):
	path_length = 0
	for pi in paths:
		path_length += len(pi)
	return path_length

def check_marked_node(marked_node, paths):
	# print(marked_node)
	for pi in paths:
		if marked_node in pi:
			return pi
	print("check node num")
	return None

def Get_N50(paths):
	scaf_lengths = read_scaffold_length()
	super_scaffold_lengths = get_super_scaffold_length(paths, scaf_lengths)
	N50 = find_N50(super_scaffold_lengths)
	return N50

def calc_raw_N50():
	scaf_lengths = read_scaffold_length(scaffold_list='scaf_rerun.list', scafoold_length='scaf_rerun.length')
	N50 = find_N50(scaf_lengths)
	return N50

def get_superScaf_sum_length(paths):
	scaf_lengths = read_scaffold_length()
	super_scaffold_lengths = get_super_scaffold_length(paths, scaf_lengths)
	return sum(super_scaffold_lengths)

def record_path(paths, output_file):
	with open(output_file, 'w') as f:
		for pi in paths:
			l = ' '.join(list(map(str, pi)))
			f.write('%s\n'%l)

def test(strength_dic, paths, from_strength_to_links):
	np.set_printoptions(precision=2)
	LDs_before = collect_LD_for_paths(paths, strength_dic)
	paths_after_break_cycles = break_each_cycle(paths, strength_dic)
	LDs_after = collect_LD_for_paths(paths_after_break_cycles, strength_dic)
	test_thresholds = get_test_threshold(LDs_after)
	test_thresholds = list(np.linspace(0,1,21)) # test_thresholds = np.linspace(0,1,21)
	for ti in test_thresholds:
		if ti > 0.45:
			continue
		new_paths = break_path_main(ti, paths_after_break_cycles, strength_dic)
		output_file = str(ti) + '.cutoff.path'
		record_path(new_paths, output_file)
		Old_N50 = Get_N50(paths)
		New_N50 = Get_N50(new_paths)
		raw_N50 = 23779253
		Len_sum = get_superScaf_sum_length(new_paths)
		line = map(str, [ti, len(new_paths), len(paths), raw_N50, New_N50])

length_record_file = 'scaf_rerun.length'
scaffold_list = 'scaf_rerun.list'
best_links_file = 'best.links'
strength_dic, couple_dic, from_strength_to_links = get_statistic(best_links_file)
scaffold_lengths = read_scaffold_length()
raw_N50= find_N50(scaffold_lengths)
paths = form_cycle(couple_dic)
test(strength_dic, paths, from_strength_to_links)
