# Author: Yuan Meng
# Time: Aug 2020

import math
import numpy as np
import sys

# ----------------------- customizable perf models ----------------------- #
def im2col_OS_latency(a,b,c,P1,P2): #a: H^2, b:K^2*C_in, c:C_out
	return math.ceil(a/P1)*math.ceil(c/P2)*b/freq + (a*b+b*c)/bandwidth

def im2col_WS_latency(a,b,c,P1,P2):
	return math.ceil(b/P1)*math.ceil(c/P2)*a/freq + (a*b+b*c)/bandwidth

def kn2row_OS_latency(a,b,c,k1,k2,P1,P2): #a: H^2, b:C_in, c:C_out
	return math.ceil(a/P1)*math.ceil(c/P2)*b*k1*k2/freq 

def kn2row_WS_latency(a,b,c,k1,k2,P1,P2):
	return math.ceil(b/P1)*math.ceil(c/P2)*a*k1*k2/freq

def wino_OS_latency(a,b,c,k,m,r,P1,P2): #a: H^2, b:C_in, c:C_out
	return (math.ceil((a/m**2)/P1)*math.ceil(c/P2)*b*(m+r-1)**2 *(k/r)**2 + (k/r)**2*(m+r-1)**2*max(m,r))/freq

def wino_WS_latency(a,b,c,k,m,r,P1,P2):
	return (math.ceil(b/P1)*math.ceil(c/P2)*(a/m**2)*(m+r-1)**2 *(k/r)**2 + (k/r)**2*(m+r-1)**2*max(m,r))/freq

# o_last=H_next, cout_last=cin_nextt
def transition_LD_cost(bool_db,P1,Cout_last,o1_last,o2_last,k1_next,k2_next,algo_last, algo_next, m,r):
	load_store_latency=0
	if (algo_last=="im2col"):
		if (algo_next=="im2col"):
			load_store_latency= o1_last*o2_last * k1_next*k2_next * Cout_last / bandwidth
		elif (algo_next=="kn2row"):
			load_store_latency= o1_last*o2_last*Cout_last / bandwidth
		elif (algo_next=="wino"):
			load_store_latency= o1_last*o2_last*((m+r-1)/m)**2 * Cout_last / (Cout_last*bandwidth/(Cout_last+m**2/o1_last*o2_last))
	elif (algo_last=="kn2row"):
		if (algo_next=="im2col"):
			load_store_latency= o1_last*o2_last * k1_next*k2_next * Cout_last / bandwidth
		elif (algo_next=="kn2row"):
			load_store_latency= o1_last*o2_last*Cout_last / bandwidth     
		elif (algo_next=="wino"):
			load_store_latency= o1_last*o2_last*((m+r-1)/m)**2 * Cout_last / (Cout_last*bandwidth/(Cout_last+m**2/o1_last*o2_last))    
	elif (algo_last=="wino"):
		if (algo_next=="im2col"):
			load_store_latency= o1_last*o2_last*k1_next*k2_next* Cout_last / bandwidth
		elif (algo_next=="kn2row"):
			load_store_latency= o1_last*o2_last*Cout_last / bandwidth
		elif (algo_next=="wino"):
			load_store_latency= o1_last*o2_last*((m+r-1)/m)**2 * Cout_last / (m**2*bandwidth)
	if(bool_db==1): # double buffering, only store-side latency+transformation paid
		return load_store_latency
	else:
		return 2*load_store_latency
# ------------------------------------------------------------------------ #

# num_layer: used to label node
def single_conv_node_cost_write(file_x,cost_vectors_per_conv,num_layer):
	file_x.write("%d \n" %len(cost_vectors_per_conv))
	for i in range(len(cost_vectors_per_conv)):
		file_x.write('%f ' %cost_vectors_per_conv[i])
	file_x.write("\n")

# single_conv_node_cost_write(f,cost_vectors_np[0],0)
# single_conv_node_cost_write(f,cost_vectors_np[1],1)
# single_conv_node_cost_write(f,cost_vectors_np[2],2)

if __name__ == "__main__":
	if (len(sys.argv)!=3):
		print ("Command Arg Error! required format: python cnn-const.py <your CNN spec input>.in <DAG generated>.in")
		sys.exit()

	input_name=sys.argv[1]
	output_name=sys.argv[2]
	print(f"Constructing graph from {sys.argv[1]}.in ... ")

	print ("=======================================")
	print ("==========Device information===========")
	print ("=======================================")
	bandwidth = input("DDR bandwidth (GB/s): ")  #10^9 B/s
	if(bandwidth==''):
		bandwidth=64
	bandwidth=int(bandwidth)*10**6 #B/ms

	print ("=======================================")
	print ("=========Hardware information==========")
	print ("=======================================")
	freq = input("Frequncy (MHz): ") #10^6 cycle/s
	if(freq==''):
		freq=200
	freq=int(freq)
	freq=freq*10**3 #cycle/ms
	P1 = input("Number of systolic array rows (input dimension): ") 
	P2 = input("Number of systolic array columns (kernel dimension):") 
	if(P1==''):
		P1=92
	if(P2==''):
		P2=64
	P1=int(P1)
	P2=int(P2)
	opt_db=input("Double Buffering (Yes-1/No-0): ") 
	if(opt_db==''):
		opt_db='1'
	bool_db=int(opt_db)
	if (int(opt_db)==1):
		if ((P1+P2)*freq>bandwidth):
			bool_db=0
		else:
			bool_db=1
	
	print ("=======================================")
	print ("===========Network Meta-data===========")
	print ("=======================================")
	#We currently support F(2x2,3x3) winograd fast algorithm
	m=2
	r=3
	# file_name = input("Model name: ") 
	N = input("Number of CONV & Concat layers: ") 
	if(N==''):
		N=8
	N=int(N)
	v_node_cost = [[] for _ in range(N)]
	branch_node_count=input("Number of branching layers: ")    # number of forking nodes in the graph
	if(branch_node_count==''):
		branch_node_count=1
	branch_node_count=int(branch_node_count)
	print(f"Generating test case file {sys.argv[2]}.in ... ")

	# Keeping track of branching, keep current (branching) node as vc, append vs later
	# inserting edge from vs->vc and replace others->vc into others->vs
	fork_vs_id=[]
	fork_vc_id=[]

	# Keeping track of algo index at each layer to facilitate transition matrix format
	algo_recorder=[[] for _ in range(N)]
	metadata_recorder=[[] for _ in range(N)]
	type_recorder=[[] for _ in range(N)]
	join_tracking=[]
	fr = open(input_name, "r")
	for i in range(N):
		v_node_cost[i]=[]
		temp_v=fr.readline().split(" ")
		type_recorder[i]=temp_v[1].strip("\n")
		# print(temp_v)
		metadata_v_s=fr.readline().split(" ")
		metadata_v_s=list(map(str.strip, metadata_v_s))
		metadata_v = [int(numeric_string) for numeric_string in metadata_v_s]
		metadata_recorder[i]=metadata_v
		if (temp_v[1]=="c\n"):
			# print(metadata_v)
			# --------------------------------- Customizable node algorithmic configurations --------------------------------- #
			if (metadata_v[2]==1): #K=1, only im2col
				algo_recorder[i].append("im2col")
				v_node_cost[i].append(min(im2col_OS_latency(metadata_v[0]*metadata_v[1],metadata_v[2]*metadata_v[3]*metadata_v[4],metadata_v[5],P1,P2),
					im2col_WS_latency(metadata_v[0]*metadata_v[1],metadata_v[2]*metadata_v[3]*metadata_v[4],metadata_v[5],P1,P2)))
			elif (metadata_v[2]==2 or metadata_v[2]!=metadata_v[3]): #K=2 or irregular shaped kernel, im2col or kn2row
				algo_recorder[i].append("im2col")
				algo_recorder[i].append("kn2row")
				v_node_cost[i].append(min(im2col_OS_latency(metadata_v[0]*metadata_v[1],metadata_v[2]*metadata_v[3]*metadata_v[4],metadata_v[5],P1,P2),
					im2col_WS_latency(metadata_v[0]*metadata_v[1],metadata_v[2]*metadata_v[3]*metadata_v[4],metadata_v[5],P1,P2)))
				v_node_cost[i].append(min(kn2row_OS_latency(metadata_v[0]*metadata_v[1],metadata_v[4],metadata_v[5],metadata_v[2],metadata_v[3],P1,P2),
					kn2row_WS_latency(metadata_v[0]*metadata_v[1],metadata_v[4],metadata_v[5],metadata_v[2],metadata_v[3],P1,P2)))
			else: #K>=3, all three
				v_node_cost[i].append(min(im2col_OS_latency(metadata_v[0]*metadata_v[1],metadata_v[2]*metadata_v[3]*metadata_v[4],metadata_v[5],P1,P2),
					im2col_WS_latency(metadata_v[0]*metadata_v[1],metadata_v[2]*metadata_v[3]*metadata_v[4],metadata_v[5],P1,P2)))
				v_node_cost[i].append(min(kn2row_OS_latency(metadata_v[0]*metadata_v[1],metadata_v[4],metadata_v[5],metadata_v[2],metadata_v[3],P1,P2),
					kn2row_WS_latency(metadata_v[0]*metadata_v[1],metadata_v[4],metadata_v[5],metadata_v[2],metadata_v[3],P1,P2)))
				v_node_cost[i].append(min(wino_OS_latency(metadata_v[0]*metadata_v[1],metadata_v[4],metadata_v[5],metadata_v[2],m,r,P1,P2),
					wino_WS_latency(metadata_v[0]*metadata_v[1],metadata_v[4],metadata_v[5],metadata_v[2],m,r,P1,P2)))
				algo_recorder[i].append("im2col")
				algo_recorder[i].append("kn2row")
				algo_recorder[i].append("wino")
			# --------------------------------------------------------------------------------------------------------------- #
		elif(temp_v[1]=="d\n"):
			v_node_cost[i].append(0) #depthwise concat has a computation cost of 0
			algo_recorder[i]=["im2col","kn2row","wino"]
		elif(temp_v[1]=="d" and temp_v[2]=="fork"):
			# v_node_cost[i]=[0]*3*int(temp_v[3].strip("\n")) #forking-depthwise concat has a computation cost of 0 but switch to all following layer algo combinations
			fork_vc_id.append(i) #track vc
			# branch_node_count+=1
		elif(temp_v[1]=="d" and temp_v[2]=="join"):
			join_tracking.append(i)
			v_node_cost[i].append(0) #joining-depthwise concat has a computation cost of 0
			algo_recorder[i]=["im2col","kn2row","wino"]
	
	for i in range(branch_node_count):
		v_node_cost.append([0]) #attach branched vs
		fork_vs_id.append(len(v_node_cost)-1) 
	if (len(fork_vs_id)!=len(fork_vc_id)): #1-to-1 match on vs->vc edge
		print("ERROR: branch_node_count",branch_node_count," does not match input specification",len(fork_vc_id))
		sys.exit()
	
	fr.readline()
	fr.readline()
	edges=[] #COO Format: [[node_x, node_y],[node_y, node_z], ...]
	while (True):
		try:
			line=fr.readline()
			if not line.strip():
				break
			else:
				edges.append(line.strip().split(" "))
		except StopIteration:
			break
	
	# Specifying all vc node vector length (total alg. choice for all downstream layers) and costs (0) for node splitting
	# specify vc algorithmic layout formats
	for branch_index in fork_vc_id:
		v_node_cost[branch_index]=[0]*sum(len(v_node_cost[int(ed[1])])*(int(ed[0])==branch_index) for ed in edges)
		concat_vc_algos=[]
		for ed in edges:
			concat_vc_algos+=algo_recorder[int(ed[1])]*(int(ed[0])==branch_index) 
		algo_recorder[branch_index]=concat_vc_algos #concat of all downstream layers
		# print (branch_index)
	
	# print (v_node_cost)
	print("Populating node costs: Done! :D")

	edge_list = [[int(x) for x in ed] for ed in edges]
	# adding new edges form node splitting
	for i in range(len(fork_vc_id)):
		edge_list.append([fork_vs_id[i],fork_vc_id[i]]) 
		type_recorder.append('d')
		algo_recorder.append('kn2row') # only 1 layout possible for fork-vs: tensor layout

	# print(edge_list)
	# print(algo_recorder)
	# print(type_recorder)
	# Populating transition matrices
	edge_costs=[]
	# print(metadata_recorder)
	for ed in edge_list:
		if (ed[0] in fork_vs_id): #(8,0) in gn_incep case
			# vs->vc costs
			if (type_recorder[ed[0]]=='d'):
			#depth concat: loop over vc then vs vs - only 1 possible layout - tensor(kn2row in)
				edge_cost=[]
				for algo_item in algo_recorder[ed[1]]:
					edge_cost.append(transition_LD_cost(bool_db,P1,
						metadata_recorder[ed[1]][2],metadata_recorder[ed[1]][0],metadata_recorder[ed[1]][1],1,1,
						"kn2row", algo_item, m,r))
				edge_costs.append(edge_cost)
			else:
			#fork conv: vs has multiple possible layouts
				edge_cost=[]
				for algo_item_end in algo_recorder[ed[0]]:
					for algo_item_begin in ["im2col","kn2row","wino"]:
						edge_cost.append(transition_LD_cost(bool_db,P1,
							metadata_recorder[ed[1]][2],metadata_recorder[ed[1]][0],metadata_recorder[ed[1]][1],1,1,
							algo_item_begin, algo_item_end,m,r))	
				edge_costs.append(edge_cost)

		elif (ed[1] in join_tracking and type_recorder[ed[1]]=='d'): #end is a join concat node, only 1 possible layout
			edge_cost=[]
			for algo_item in algo_recorder[ed[0]]:
				edge_cost.append(transition_LD_cost(bool_db,P1,
					metadata_recorder[ed[1]][2],metadata_recorder[ed[1]][0],metadata_recorder[ed[1]][1],1,1,
					"kn2row", algo_item, m,r))	
			edge_costs.append(edge_cost)

		else:
			edge_cost=[]
			for algo_item_end in algo_recorder[ed[1]]:
				for algo_item_begin in algo_recorder[ed[0]]:
					edge_cost.append(transition_LD_cost(bool_db,P1,
						metadata_recorder[ed[1]][4],metadata_recorder[ed[1]][0],metadata_recorder[ed[1]][1],
						metadata_recorder[ed[1]][2],metadata_recorder[ed[1]][3],
						algo_item_begin, algo_item_end, m,r))			
			edge_costs.append(edge_cost)		

	# print(edge_costs)
	# print(len(edge_costs))

	print("Populating edge costs: Done! :D")

	fr.close()
	fw = open("dump/"+output_name, "w")

	fw.write("%d \n" %len(v_node_cost)) #num nodes
	fw.write("%d \n" %len(edge_list))   #num edges

	# write vertices
	for i in range(len(v_node_cost)):
		fw.write("%d \n" %(len(v_node_cost[i])))
		for j in range(len(v_node_cost[i])):
			fw.write('%f ' %v_node_cost[i][j])
		fw.write("\n")

	# write edges
	fw.write("\n")
	for i in range(len(edge_list)):
		fw.write("%d %d\n" %(edge_list[i][0],edge_list[i][1]))
		fw.write("%d %d\n" %(len(algo_recorder[edge_list[i][0]]),len(algo_recorder[edge_list[i][1]])))
		for j in range(len(edge_costs[i])):
			fw.write('%f ' %edge_costs[i][j])
		fw.write("\n")
	fw.close()

	print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print ("~~~~~~~CNN GRAPH CONSTRUCTION COMPLETE~~~~~~~~")
	print ("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

