'''
Title:   Build Micelle

Author:  Gubbin Eel (Samuel Lindsay, East Carolina University)

Purpose: Builds a micelle out of lipids.

NOTE:    Intended improvements are denoted by the #%# flag.
         Searching for that flag will identify code that needs to be
         updated or will later be improved.
'''

import manipulate_pdb      as manp
import manipulate_files    as manf
import molecule_structures as mols
import geom_shapes         as gs

import numpy as np
import copy

class Micelle(object):
	''' Has two lists as attributes. One is of lipid objects that make up the 
	    structure. The other is a list of those lipids as strings in the 
	    format of pdb file lines.

		Attributes: struc and lines
	'''

	''' Initialize class instance.
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
				atom = atom.to_pdb_file(atom_serial = str(jj), resn_seq = str(ii))
				self.lines.append(atom)
				jj += 1
			self.lines.append("TER\n")
			ii += 1

	''' Purpose:   To write the file lines in the "lines" attribute of a Micelle
	               instance to a file.
	    Arguments: (1) A Micelle instance specified in the class method format.
	               (2) The name of the file to export the pdb structure to.
	                   Specify the full path in the name.
	'''
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
	def build_micelle(self, lipids, frac, radius, points, outfile):
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
		sph_coord = gs.fib_sphere(points, radius = radius)
		print(sph_coord)
		for ii in range(0, len(sph_coord)):
			# Copy of lipid
			lipid = copy.deepcopy(struc[lipid_id[ii]])

			# Determine which oxygen and hydrogen atom are farthest from each other. 
			extrema = lipid.lipid_extrema()

			# Translate the Lipid such that the position of the hydrogen extreme is 
			# at the origin.
			lipid.trans_axis_h_to_point(extrema, [0, 0, 0])

			# Rotate lipid to align with the position vector on the sphere's surface.
			lipid.align_axis_o_to_vec(extrema, sph_coord[ii])
			

			# Translate the lipid to the position on the sphere.
			lipid.trans_axis_o_to_point(extrema, sph_coord[ii])

			# Add lipid to list
			self.struc.append(lipid)

		# Convert familythe list of micelle lipids to a list of lists containing strings for 
		# each atom.
		self.struc_to_lines()

		# Write the micelle to a file.
		self.save_lines(outfile)

if __name__ == "__main__":
	# Initialize Micelle class instance.
	micelle = Micelle()

	# Output file.
	file_name = "/home/gubbin/Documents/mdProjects/programs/MMAEVe/lipids/test1.pdb"
	micelle.build_micelle(["POC"], [1.0], 100, 800, file_name)
