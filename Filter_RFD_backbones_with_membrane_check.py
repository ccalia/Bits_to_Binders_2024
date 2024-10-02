
import Bio.PDB
import numpy as np
import os, glob, sys, re


usage = '\nUsage: python Filter_RFD_backbones_with_membrane_check.py RFD_pdbs_folder_path filtered_pdbs_folder_path require_contact forbid_contact min_num_target_res_in_contact min_num_helices max_num_helices max_helix_length min_num_sheets path_to_membrane_pdb first_align_res last_align_res max_num_clashes\n\nBinder must be chain A, target must be chain B\nrequire_contact and forbid_contact count from 0\nIf there are at least min_num_sheets, min_num_helices is not considered\nTo not use a criterium set it to -1\nfirst_align_res/last_align_res are residue indices from 0 in target that align to the start/end of IKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQY (the I and Y)\nIf max_num_clashes is negative, path_to_membrane_pdb and first_align_res/last_align_res are not used\n'

if len(sys.argv) != 14:
	print(usage)
	sys.exit()

RFD_pdbs_folder_path = sys.argv[1]
filtered_pdbs_folder_path = sys.argv[2]
require_contact = int(sys.argv[3])
forbid_contact = int(sys.argv[4])
min_num_target_res_in_contact = int(sys.argv[5])
min_num_helices = int(sys.argv[6])
max_num_helices = int(sys.argv[7])
max_helix_length = int(sys.argv[8])
min_num_sheets = int(sys.argv[9])
path_to_membrane_pdb = sys.argv[10]
first_align_res = int(sys.argv[11]) #IKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQY
last_align_res = int(sys.argv[12]) #IKISHFLKMESLNFIRAHTPYINIYNCEPANPSEKNSPSTQY
max_num_clashes = int(sys.argv[13])


contact_cutoff_distance = 10.0


if os.path.exists(filtered_pdbs_folder_path) != True:
	os.system(f'mkdir {filtered_pdbs_folder_path}')


parser = Bio.PDB.PDBParser(QUIET = True)


#Credit for contact map: https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/

def calc_residue_dist(residue_one, residue_two) :
	"""Returns the C-alpha distance between two residues"""
	diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
	return np.sqrt(np.sum(diff_vector * diff_vector))


def calc_dist_matrix(chain_one, chain_two) :
	"""Returns a matrix of C-alpha distances between two chains"""
	answer = np.zeros((len(chain_one), len(chain_two)), float)
	for row, residue_one in enumerate(chain_one) :
		for col, residue_two in enumerate(chain_two) :
			answer[row, col] = calc_residue_dist(residue_one, residue_two)
	return answer


def get_A_B_contact_map(path_to_pdb, cutoff_distance):
	structure = parser.get_structure('scratch', path_to_pdb)
	model = structure[0]
	dist_matrix = calc_dist_matrix(model['A'], model['B'])
	contact_map = dist_matrix < cutoff_distance
	return contact_map


#n from 0, ctmapT = transpose of ctmap
def is_res_n_of_B_in_contact_with_A(ctmapT, n):
	if True in ctmapT[n]:
		return True
	else:
		return False


#ctmapT = transpose of ctmap
def num_B_res_in_contact_with_A(ctmapT):
	contacts = [i for i in range(len(ctmapT)) if True in ctmapT[i]]
	return len(contacts)


def count_helices(path_to_pdb):
	structure = parser.get_structure('scratch', path_to_pdb)
	model = structure[0]
	dssp = Bio.PDB.DSSP(model, path_to_pdb)
	secstruct = ''
	for i in range(len(model['A'])): #from https://github.com/lamm-mit/ProteinDiffusionGenerator/blob/main/ProteinDiffusionGenerator_Model_B.ipynb
		a_key = list(dssp.keys())[i]
		secstruct += dssp[a_key][2]
	helices = re.findall('H' + '+', secstruct)
	sheets = re.findall('E' + '+', secstruct)
	if len(helices) == 0:
		max_hlx_len = 0
	else:
		max_hlx_len = max([len(h) for h in helices])
	return len(helices), max_hlx_len, len(sheets)


#Checking clashes - from https://www.blopig.com/blog/2023/05/checking-your-pdb-file-for-clashing-atoms/

atom_radii = {
#	"H": 1.20,  # Who cares about hydrogen??
	"C": 1.70, 
	"N": 1.55, 
	"O": 1.52,
	"S": 1.80,
#	"F": 1.47, 
	"P": 1.80, 
}

def count_clashes(structure, clash_cutoff=0.63):
	# Set what we count as a clash for each pair of atoms
	clash_cutoffs = {i + "_" + j: (clash_cutoff * (atom_radii[i] + atom_radii[j])) for i in atom_radii for j in atom_radii}
	# Extract atoms for which we have a radii
	atoms = [x for x in structure.get_atoms() if x.element in atom_radii]
	coords = np.array([a.coord for a in atoms], dtype="d")
	# Build a KDTree (speedy!!!)
	kdt = Bio.PDB.kdtrees.KDTree(coords)
	# Initialize a list to hold clashes
	clashes = []
	# Iterate through all atoms
	for atom_1 in atoms:
		# Find atoms that could be clashing
		kdt_search = kdt.search(np.array(atom_1.coord, dtype="d"), max(clash_cutoffs.values()))
		# Get index and distance of potential clashes
		potential_clash = [(a.index, a.radius) for a in kdt_search]
		for ix, atom_distance in potential_clash:
			atom_2 = atoms[ix]
			# Exclude clashes from atoms in the same residue
			if atom_1.parent.id == atom_2.parent.id:
				continue
			# Exclude clashes from peptide bonds
			elif (atom_2.name == "C" and atom_1.name == "N") or (atom_2.name == "N" and atom_1.name == "C"):
				continue
			# Exclude clashes from disulphide bridges
			elif (atom_2.name == "SG" and atom_1.name == "SG") and atom_distance > 1.88:
				continue
			if atom_distance < clash_cutoffs[atom_2.element + "_" + atom_1.element]:
				clashes.append((atom_1, atom_2))
	return len(clashes) // 2


#Based on https://combichem.blogspot.com/2013/08/aligning-pdb-structures-with-biopython.html
#Tested alignment with test15/FILTERED_test15RFD3_outputs/RFD3_test15_228.pdb, super_imposer.rms from this alignment with first/last_align_res = 29/70 matched rmsd from Chimera for: match #1:110-151.B@CA #0:97-138.A@CA
def align_to_struc_in_membrane(ref_structure, sample_structure):
	ref_model = ref_structure[0]
	sample_model = sample_structure[0]
	ref_atoms = []
	sample_atoms = []
	for ref_res in ref_model['A']: #get CA atoms from CD20 in membrane
		ref_atoms.append(ref_res['CA'])
	for sample_res in sample_model['B']: #get CA atoms from target struc
		sample_atoms.append(sample_res['CA'])
	super_imposer = Bio.PDB.Superimposer()
	super_imposer.set_atoms(ref_atoms[96:138], sample_atoms[first_align_res:last_align_res + 1])
	super_imposer.apply(sample_model.get_atoms())
	#return super_imposer.rms
	ref_model.detach_child('A') #get rid of cd20 in membrane
	ref_model.add(sample_model['A']) #add binder to membrane structure
	return ref_structure #now just binder and membrane



def check_clashes_with_membrane_structure(path_to_pdb):
	target_binder_structure = parser.get_structure('scratch', path_to_pdb)
	cd20_in_memb_structure = parser.get_structure('scratchm', path_to_membrane_pdb)
	binder_with_memb_structure = align_to_struc_in_membrane(cd20_in_memb_structure, target_binder_structure)
	return count_clashes(binder_with_memb_structure, clash_cutoff = 0.63)




def check_pdb(path_to_pdb):
	num_helices, max_hlx_len, num_sheets = count_helices(path_to_pdb)
	sheets_criteria_met = False
	if min_num_sheets >= 0:
		if num_sheets >= min_num_sheets:
			sheets_criteria_met = True
	if max_helix_length >= 0:
		if max_hlx_len > max_helix_length:
			return False
	if (min_num_helices >= 0) and (not sheets_criteria_met):
		if num_helices < min_num_helices:
			return False
	if max_num_helices >= 0:
		if num_helices > max_num_helices:
			return False
	ctm = get_A_B_contact_map(path_to_pdb, contact_cutoff_distance)
	ctmT = np.transpose(ctm) #Now rows = res in chain B
	if min_num_target_res_in_contact >= 0:
		if num_B_res_in_contact_with_A(ctmT) < min_num_target_res_in_contact:
			return False
	if require_contact >= 0:
		if not is_res_n_of_B_in_contact_with_A(ctmT, require_contact):
			return False
	if forbid_contact >= 0:
		if is_res_n_of_B_in_contact_with_A(ctmT, forbid_contact):
			return False
	if max_num_clashes >= 0:
		nclashes = check_clashes_with_membrane_structure(path_to_pdb)
		if nclashes > max_num_clashes:
			return False
	return True


count = 0
total = 0

for pdbpath in glob.glob(os.path.join(RFD_pdbs_folder_path, '*.pdb')):
	if check_pdb(pdbpath) == True:
		os.system(f'cp {pdbpath} {filtered_pdbs_folder_path}')
		count += 1
	total += 1


with open(os.path.join(filtered_pdbs_folder_path, 'filter.txt'), 'a') as logfile:
	logfile.write(f'Filtered with: RFD_pdbs_folder_path: {RFD_pdbs_folder_path}, filtered_pdbs_folder_path: {filtered_pdbs_folder_path}, require_contact: {require_contact}, forbid_contact: {forbid_contact}, min_num_target_res_in_contact: {min_num_target_res_in_contact}, min_num_helices: {min_num_helices}, max_num_helices: {max_num_helices}, max_helix_length: {max_helix_length}, min_num_sheets: {min_num_sheets}, path_to_membrane_pdb: {path_to_membrane_pdb}, first_align_res: {first_align_res}, last_align_res: {last_align_res}, max_num_clashes: {max_num_clashes}\n')
	logfile.write(f'{count} / {total} pdbs met the criteria.\n\n')


print(f'\n{count} / {total} pdbs met the criteria.\n')


