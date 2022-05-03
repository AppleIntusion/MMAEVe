'''

DISCLAIMER: All comments are made to appease the one and only Gubbin Eel; I do not approve

'''

import manipulate_pdb as manp
import manipulate_files as manf
import molecule_structures as mols
import geom_shapes as gs

import numpy as np
import copy

class Nanodisc(object):
	'''IN PROGRESS
	'''
	
	def __init__(self, struc = [], lines = []):
		self.struc = struc
		self.lines = lines

	''' Purpose:   Convert the struc attribute from a list of Lipid class
	               objects to a list lists where each list contains the
	               lines of the pdb file as strings. And assign it to the
	               Micelle instance "lines" attribute.
	    Arguments: (1) Micelle class instance specified using the class
	                   method format.
	'''
	def struc_to_lines(self):
		ii = 1 # Index variable to track the atom number.
		jj = 1 # Index variable to track the residue number.
		for lipid in self.struc:
			for atom in lipid.atom_info:
				atom = lipid.atom_info[atom]
				atom = atom.to_pdb_file(atom_serial = str(ii), resn_seq = str(jj))
				self.lines.append(atom)
				ii += 1
			self.lines.append("TER\n")
			jj += 1

	''' Purpose:   To write the file lines in the "lines" attribute of a Micelle
	               instance to a file.
	    Arguments: (1) A Micelle instance specified in the class method format.
	               (2) The name of the file to export the pdb structure to.
	                   Specify the full path in the name.
	'''
	def struc_to_lines_protein(self):	
		'''
		Purpose:  To write the file lines in the "lines" attribute of a nanodisc instance to a file.

		Arguments:
		'''
		
		ii = 0 # loop index
		for pro in self.struc:
		
			end = list(pro.atom_info.items())[-1][0]	
			for atom in pro.atom_info:
				atom_num = pro.atom_info[atom].serial
				residue_num = pro.atom_info[atom].resi
				atom = pro.atom_info[atom]
				atom = atom.to_pdb_file(atom_serial = str(atom_num + pro.atom_info[end].serial * ii), resn_seq = str(residue_num + pro.atom_info[end].resi * ii))
				self.lines.append(atom)
		
			self.lines.append("TER\n")
			ii +=1
	
	def struc_to_lines_composite(self):
			
		ii = 1 # Index variable to track the atom number.
		jj = 1 # Index variable to track the residue number.
		counter = 0
		prot = 0
		lip = 0	
		index= 0
		for unit in self.struc:
			index = list(self.struc).index(unit)
			
			elem = list(unit.atom_info.items())[0][0]
			if unit.atom_info[elem].resn != "POC":
				prot += 1	
				end = list(unit.atom_info.items())[-1][0]	
			
			if unit.atom_info[elem].resn == "POC":
				lip += 1	
			for atom in unit.atom_info:
				if unit.atom_info[atom].resn == "POC":
					lipid_atom = unit.atom_info[atom]
					lipid_atom = lipid_atom.to_pdb_file(atom_serial = str(ii), resn_seq = str(jj))

					self.lines.append(lipid_atom)
				else:
					
					atom_num = unit.atom_info[atom].serial
					residue_num = unit.atom_info[atom].resi
					print(residue_num)
					protein_atom = unit.atom_info[atom]
					#atom = atom.to_pdb_file(atom_serial = str(atom_num + unit.atom_info[end].serial * prot), resn_seq = str(residue_num + unit.atom_info[end].resi * prot))
					protein_atom = protein_atom.to_pdb_file(atom_serial = str(ii), resn_seq = str(residue_num + lip + prot * unit.atom_info[end].resi))
					self.lines.append(protein_atom)
				ii += 1
				
			self.lines.append("TER\n")
			jj += 1
			counter +=1
		print(prot)	
	def save_lines(self, file_name):
		manf.write_file(file_name, self.lines)

	''' Purpose:   Build a micelle of variable size and composition.
	    Arguments: (1) A Micelle class object specified in the class method
	                   format.
	               (2) A list of strings. Each string is the lipid name or
	                   abbreviation preceeding the ".pdb" suffix of structures in
	                   the "lipids" folder.
	               (3) A list of fractions correspoinding the the proportion
	                   of each lipid present in the micelle. Must sum to 1.
	               (4) Radius of the micelle. Float.
	               (5) Number of lipids in the micelle. Integer.
	               (6) The name of the output file with its full path specified.
	'''
	def build_top_nanodisc(self, lipids, frac, radius, points):
		# Create a dictionary of lipid structures.
		struc = dict()
		for ii in range(0, len(lipids)):
			struc[ii] = manp.PdbFile(lipid_name = lipids[ii]).to_lipid()

		# Number of entries to replace in each list.
		lipid_num = [int(prob * points) for prob in frac]

		# A list of points containing the lipid ID to be placed at each
		# point. When initialized it will contain only a single id.
		lipid_id = [0 for ii in range(0, points)]

		# A list from 0 to the number of lipids in the micelle. Represents the
		# indicies of entries in the lipid_id list.
		lipid_id_index = [ii for ii in range(0, points)]

		# num_to_replace is the number of lipids to replace with a given lipid srtucture
		# index, begining with the second index "1".
		lip_struc_id = 1 # ID of the structure to replace.
		for num_to_replace in lipid_num[1:]:
			for ii in range(0, num_to_replace):
				# Pick a random index from the lipid index list.
				index = np.random.choice(lipid_id_index)
				# Assign the structure ID to replace the default ID at a
				# given index.
				lipid_id[index] = lip_struc_id
				# Remove the index from the list so that it cannot be used again.
				lipid_id_index.remove(index)
			# Increase ID so that it refelects the next structure ID.
			lip_struc_id += 1

		#print(lipid_id)
		#print(lipid_num)

		# Coordinates on the surface of the sphere which lipids will be mapped to.
		circ_coord = gs.sunflower(points, 2, radius = radius)
		for i in range(0, len(circ_coord)):
			lipid = copy.deepcopy(struc[lipid_id[i]])
			extrema = lipid.lipid_extrema()
		
			# Translate lipids to origin
			lipid.trans_axis_h_to_point(extrema, [0, 0, 0])
			
			lipid.align_axis_o_to_vec(extrema, [1, 1, 0])

			lipid.rotate(axis = 'z', theta = 3.14* 0.25)
			
			lipid.rotate(axis = 'x', theta = 3.14* 0.5)
			# Move lipids
			lipid.trans_axis_o_to_point(extrema, circ_coord[i])
			
			self.struc.append(lipid)	
	
#		self.struc_to_lines()


	def build_bottom_nanodisc(self, lipids, frac, radius, points):
		# Create a dictionary of lipid structures.
		struc = dict()
		for ii in range(0, len(lipids)):
			struc[ii] = manp.PdbFile(lipid_name = lipids[ii]).to_lipid()

		# Number of entries to replace in each list.
		lipid_num = [int(prob * points) for prob in frac]

		# A list of points containing the lipid ID to be placed at each
		# point. When initialized it will contain only a single id.
		lipid_id = [0 for ii in range(0, points)]

		# A list from 0 to the number of lipids in the micelle. Represents the
		# indicies of entries in the lipid_id list.
		lipid_id_index = [ii for ii in range(0, points)]

		# num_to_replace is the number of lipids to replace with a given lipid srtucture
		# index, begining with the second index "1".
		lip_struc_id = 1 # ID of the structure to replace.
		for num_to_replace in lipid_num[1:]:
			for ii in range(0, num_to_replace):
				# Pick a random index from the lipid index list.
				index = np.random.choice(lipid_id_index)
				# Assign the structure ID to replace the default ID at a
				# given index.
				lipid_id[index] = lip_struc_id
				# Remove the index from the list so that it cannot be used again.
				lipid_id_index.remove(index)
			# Increase ID so that it refelects the next structure ID.
			lip_struc_id += 1

		#print(lipid_id)
		#print(lipid_num)

		# Coordinates on the surface of the sphere which lipids will be mapped to.
		circ_coord = gs.sunflower(points, 2, radius = radius)
		for ii in range(0, len(circ_coord)):
			lipid = copy.deepcopy(struc[lipid_id[ii]])
			extrema = lipid.lipid_extrema()
		
			# Translate lipids to origin
			lipid.trans_axis_h_to_point(extrema, [0, 0, 0])
			
			lipid.align_axis_o_to_vec(extrema, [1, 1, 0])

			lipid.rotate(axis = 'z', theta = 3.14* 0.25)
			
			lipid.rotate(axis = 'x', theta = 3.14* 0.5)
			# Move lipids
			lipid.trans_axis_o_to_point(extrema, circ_coord[ii])
			
			lipid.rotate(axis = 'y', theta = 3.14 )
			lipid.translate(t=[0, 0, -50])
	
			self.struc.append(lipid)	
		
	def build_protein(self, protein, frac, radius, points):
		# Create a dictionary of lipid structures.
		struc = dict()
		for ii in range(0, len(protein)):
			struc[ii] = manp.PdbFile(file_name = protein[ii]).to_protein()
	#	protein = copy.deepcopy(struc[protein_id[0]])
		#protein = copy.deepcopy(struc[0])
		#self.struc.append(protein)		
		#self.struc_to_lines_protein()
		#self.save_lines(outfile)
		# Number of entries to replace in each list.
		pro_num = [int(prob * points) for prob in frac]

		# A list of points containing the lipid ID to be placed at each
		# point. When initialized it will contain only a single id.
		protein_id = [0 for ii in range(0, points)]


		# A list from 0 to the number of lipids in the micelle. Represents the
		# indicies of entries in the lipid_id list.
		protein_id_index = [ii for ii in range(0, points)]

		# num_to_replace is the number of lipids to replace with a given lipid srtucture
		# index, begining with the second index "1".
		pro_struc_id = 1 # ID of the structure to replace.
		for num_to_replace in pro_num[1:]:
			for ii in range(0, num_to_replace):
				# Pick a random index from the lipid index list.
				index = np.random.choice(pro_id_index)
				# Assign the structure ID to replace the default ID at a
				# given index.
				pro_id[index] = pro_struc_id
				# Remove the index from the list so that it cannot be used again.
				pro_id_index.remove(index)
			# Increase ID so that it refelects the next structure ID.
			pro_struc_id += 1
	
		#print(lipid_id)
		#print(lipid_num)
#		protein = copy.deepcopy(struc[protein_id[0]])
#		self.struc.append(protein)		
#		self.struc_to_lines_protein()
#		self.save_lines(outfile)

		# Coordinates on the surface of the sphere which lipids will be mapped to.
		circ_coord = gs.pro_circle(points, radius)
		for i in range(0, len(circ_coord)):
			protein = copy.deepcopy(struc[protein_id[i]])
			print(circ_coord[i])	
			extrema = protein.lipid_extrema()
		
			# Translate lipids to origin
			protein.trans_axis_h_to_point(extrema, [0, 0, 0])
			
			protein.align_axis_o_to_vec(extrema, [1, 1, 0])

			protein.rotate(axis = 'z', theta = 3.14* 0.25)
			
			protein.rotate(axis = 'x', theta = 3.14* 0.5)
			# Move lipids
			protein.trans_axis_o_to_point(extrema, circ_coord[i])
			
			self.struc.append(protein)		


	def build_protein_belt(self, protein, frac, radius, points):
		# Create a dictionary of lipid structures.
		struc = dict()
		for ii in range(0, len(protein)):
			struc[ii] = manp.PdbFile(file_name = protein[ii]).to_protein()
	#	protein = copy.deepcopy(struc[protein_id[0]])
		#protein = copy.deepcopy(struc[0])
		#self.struc.append(protein)		
		#self.struc_to_lines_protein()
		#self.save_lines(outfile)
		# Number of entries to replace in each list.
		pro_num = [int(prob * points) for prob in frac]

		# A list of points containing the lipid ID to be placed at each
		# point. When initialized it will contain only a single id.
		protein_id = [0 for ii in range(0, points)]


		# A list from 0 to the number of lipids in the micelle. Represents the
		# indicies of entries in the lipid_id list.
		protein_id_index = [ii for ii in range(0, points)]

		# num_to_replace is the number of lipids to replace with a given lipid srtucture
		# index, begining with the second index "1".
		pro_struc_id = 1 # ID of the structure to replace.
		for num_to_replace in pro_num[1:]:
			for ii in range(0, num_to_replace):
				# Pick a random index from the lipid index list.
				index = np.random.choice(pro_id_index)
				# Assign the structure ID to replace the default ID at a
				# given index.
				pro_id[index] = pro_struc_id
				# Remove the index from the list so that it cannot be used again.
				pro_id_index.remove(index)
			# Increase ID so that it refelects the next structure ID.
			pro_struc_id += 1
	
		#print(lipid_id)
		#print(lipid_num)
#		protein = copy.deepcopy(struc[protein_id[0]])
#		self.struc.append(protein)		
#		self.struc_to_lines_protein()
#		self.save_lines(outfile)

		# Coordinates on the surface of the sphere which lipids will be mapped to.
		circ_coord = gs.pro_circle(points, radius)
		for i in range(0, len(circ_coord)):
			protein = copy.deepcopy(struc[protein_id[i]])
			print(circ_coord[i])	
			extrema = protein.lipid_extrema()
		
			# Translate lipids to origin
			protein.trans_axis_h_to_point(extrema, [0, 0, 0])
			
			protein.align_axis_o_to_vec(extrema, [1, 1, 0])

			protein.rotate(axis = 'z', theta = 3.14* 0.25 + 3.14*i*2/len(circ_coord))
			
			
			# Move lipids
			protein.trans_axis_o_to_point(extrema, circ_coord[i])
			protein.translate(t=[0,0, -5])
			
			self.struc.append(protein)		
	def composite_nanodisc(self, lipids, frac, radius, lipid_points, protein_points, protein, outfile):
		
		self.build_top_nanodisc(lipids, frac, radius, lipid_points)
		self.build_bottom_nanodisc(lipids, frac, radius, lipid_points)
		self.build_protein(protein, frac, radius + 10, protein_points)
		self.struc_to_lines_composite()
		self.save_lines(outfile)

	def composite_nanodisc_belt(self, lipids, frac, radius, lipid_points, protein_points, protein, outfile):
		
		self.build_top_nanodisc(lipids, frac, radius, lipid_points)
		self.build_bottom_nanodisc(lipids, frac, radius, lipid_points)
		self.build_protein_belt(protein, frac, radius + 10, protein_points)
		self.struc_to_lines_composite()
		self.save_lines(outfile)
if __name__ == "__main__":
	
	nanodisc = Nanodisc()

	file_name = "/home/sanchezw/MMAEVe/builder/nano_test.pdb"
	file_name2 = "/home/sanchezw/MMAEVe/builder/pep_test.pdb"
	file_name3 = "/home/sanchezw/MMAEVe/builder/disc_test.pdb"
	file_name4 = "/home/sanchezw/MMAEVe/builder/belt_test.pdb"
	
	#nanodisc.build_nanodisc(["POC"], [1,0] , 100, 200,["apoa1_peptide"], file_name)
#	nanodisc.build_protein(["../proteins/apoa1_peptide.pdb"], [1] , 100, 3, file_name2)				
#	nanodisc.build_nanodisc(["POC"], [1,0] , 100, 200,["apoa1_peptide"], file_name2)
#	nanodisc.composite_nanodisc(["POC"], [1], 100, 200, 30, ["../proteins/apoa1_peptide.pdb"], file_name3)
	nanodisc.composite_nanodisc_belt(["POC"], [1], 100, 200, 15, ["../proteins/apoa1_peptide.pdb"], file_name4)

