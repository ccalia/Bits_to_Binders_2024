
import os, sys

usage = '\npython Filter_CD_AF_output_fasta_RFD_pdbs.py path_to_CD_fasta_file path_to_RFD_pdbs path_to_AF_pdbs where_to_put_pdbs_passing_filters MIN_plddt MIN_iptm MAX_ipae MAX_rmsd add_to_label\nSet add_to_label to None to keep original labels without an extension\n'

if len(sys.argv) != 10:
	print(usage)
	sys.exit()


path_to_CD_fasta_file = sys.argv[1]
path_to_RFD_pdbs = sys.argv[2]
path_to_AF_pdbs = sys.argv[3]
where_to_put_pdbs_passing_filters = sys.argv[4] #output folder that will have passing_RFD_pdbs and passing_AF_pdbs subfolders, and a filtered.txt file
MIN_plddt = float(sys.argv[5])
MIN_iptm = float(sys.argv[6])
MAX_ipae = float(sys.argv[7])
MAX_rmsd = float(sys.argv[8])


if sys.argv[9] == 'None':
	add_to_label = ''
else:
	add_to_label = sys.argv[9]



if os.path.exists(where_to_put_pdbs_passing_filters) != True:
	os.system(f'mkdir {where_to_put_pdbs_passing_filters}')


if os.path.exists(os.path.join(where_to_put_pdbs_passing_filters, 'passing_RFD_pdbs')) != True:
	passrfdpath = os.path.join(where_to_put_pdbs_passing_filters, 'passing_RFD_pdbs')
	os.system(f'mkdir {passrfdpath}')


if os.path.exists(os.path.join(where_to_put_pdbs_passing_filters, 'passing_AF_pdbs')) != True:
	passafpath = os.path.join(where_to_put_pdbs_passing_filters, 'passing_AF_pdbs')
	os.system(f'mkdir {passafpath}')


with open(path_to_CD_fasta_file, 'r') as df:
	datalines = [line for line in df]

datadict = {}

for i in range(len(datalines) - 1):
	if i % 2 == 0:
		datadict[datalines[i]] = datalines[i + 1] #keys are fasta header lines, values are target/binder seqs


def parse_header(header):
	split_header = header.split('|')
	label = split_header[0][1:]
	plddt = float(split_header[2].split(':')[-1])
	iptm = float(split_header[3].split(':')[-1])
	ipae = float(split_header[4].split(':')[-1])
	rmsd = float(split_header[5].split(':')[-1].strip())
	return label, plddt, iptm, ipae, rmsd



with open(os.path.join(where_to_put_pdbs_passing_filters, 'filtered.txt'), 'a') as outfile:
	outfile.write(f'Filtered RFD pdbs in {path_to_RFD_pdbs} and AF pdbs in {path_to_AF_pdbs} based on metrics in {path_to_CD_fasta_file}\n')
	outfile.write(f'Criteria: ipae < {MAX_ipae}, rmsd < {MAX_rmsd}, plddt > {MIN_plddt}, iptm > {MIN_iptm}\n')
	outfile.write(f'label,binder_seq,ipae,rmsd,plddt,iptm\n')


pass_count = 0
total_count = 0

for headerline in datadict:
	label, plddt, iptm, ipae, rmsd = parse_header(headerline)
	if (plddt > MIN_plddt) and (iptm > MIN_iptm) and (ipae < MAX_ipae) and (rmsd < MAX_rmsd):
		rfdpath = os.path.join(path_to_RFD_pdbs, label.split('_seq')[0] + '.pdb')
		os.system(f'cp {rfdpath} {passrfdpath}')
		afpath = os.path.join(path_to_AF_pdbs, label + '_af2.pdb')
		afdest = os.path.join(passafpath, label + add_to_label + '_af2.pdb')
		os.system(f'cp {afpath} {afdest}')
		with open(os.path.join(where_to_put_pdbs_passing_filters, 'filtered.txt'), 'a') as outfile:
			seq = datadict[headerline].split('/')[-1].strip()
			flabel = label + add_to_label
			outfile.write(f'{flabel},{seq},{ipae},{rmsd},{plddt},{iptm}\n')
		pass_count += 1
	total_count += 1


print(f'\n{pass_count} / {total_count} passing\n')


