'''
Title:   Manipulate PDB

Author:  Gubbin Eel (Samuel Lindsay, East Carolina University)

Purpose: Defines functions used to import, format, and export .pdb
         files.

NOTE:    Intended improvements are denoted by the #%# flag.
         Searching for that flag will identify code that needs to be
         updated or will later be improved.
'''

import manipulate_files    as manf
import molecule_structures as mols

class PdbFile():
	'''
	#%# Description of the class.
	
	'''

	#%# Update to include a count of file lines and fields. At the
	#   moment I just need it to import ATOM information.
	def __init__(self, file_name = "None", lipid_name = "None"):
		# Variable to control whether class instance is initialized.
		run = False
		if file_name == "None" and lipid_name == "None":
			print("Either file_path or lipid_name must be specified. Both" + \
			      "cannot be unspecified.")
		elif file_name != "None" and lipid_name != "None":
			print("Either file_path or lipid_name must be specified. Both" + \
			      "cannot be specified.")
		else:
			run = True

		# File name is generated if a lipid name is specified instead of a file
		# name.
		if run == True and file_name != "None":
			pass
		elif run == True and lipid_name  != "None":
			file_name = "../lipids/" + lipid_name + ".pdb"

		# If name is specified properly then initialize the class instance.
		if run == True:
			file_lines = manf.import_file_contents(file_name)
			atom_stuff = []
			for line in file_lines:
				if line[0:4] == "ATOM":
					'''
					Atom fields with their index number in the list.
					0:  "ATOM  "        1:  Atom serial num 2:  Atom name 
					3:  Alternate loc   4:  Residue Name    5:  Chain ID 
					6:  Residue seq num 7:  iCode           8:  x-coord 
					9:  y-coord         10: z-coord         11: Occupancy 
					12: Temp factor     13: Elem sym        14: Atomic charge
					'''
					# All line entries for an ATOM entry in a .pdb file.
					line      = [line[0:6],   line[6:11],  line[12:16],
					             line[16],    line[17:20], line[21],
				                 line[22:26], line[26],    line[30:38],
					             line[38:46], line[46:54], line[54:60],
					             line[60:66], line[76:78], line[78:80]]

					# Stripping whitespace and formatting numeric entries to be
					# either integers or floats, whichever is most appropriate.
					line[0]   = line[0].strip(manf.whitespace)
					line[1]   = int(line[1].strip(manf.whitespace))
					line[2]   = line[2].strip(manf.whitespace)
					line[3]   = line[3].strip(manf.whitespace)
					line[4]   = line[4].strip(manf.whitespace)
					line[6]   = int(line[6].strip(manf.whitespace))
					line[5]   = line[7].strip(manf.whitespace)
					line[8]   = float(line[8].strip(manf.whitespace))
					line[9]   = float(line[9].strip(manf.whitespace))
					line[10]  = float(line[10].strip(manf.whitespace))
					line[11]  = float(line[11].strip(manf.whitespace))
					line[12]  = float(line[12].strip(manf.whitespace))
					line[13]  = line[13].strip(manf.whitespace)
					try:
						line[14]  = float(line[14].strip(manf.whitespace))
					except ValueError:
						line[14] = ''

					atom_stuff.append(line)

			self.line      = file_lines
			self.numline   = len(file_lines)
			self.atom_info = atom_stuff
		else:
			print("Nothing happened. Congratulations.")

	''' Purpose: Take a PdbFile class instance and convert it to a Lipid
	             class instance.
	'''
	def to_lipid(self):
		molecule_atoms = dict() 
		for atom in self.atom_info:
			atom_mass = mols.atomic_masses[atom[13]]
			molecule_atoms[atom[1]] = mols.Atom(serial = atom[1],   x      = atom[8],  
			                                    y      = atom[9],   z      = atom[10], 
			                                    mass   = atom_mass, resi   = 1,        
			                                    resn   = atom[4],   elem   = atom[13], 
			                                    name   = atom[2],   charge = atom[14])

		return mols.Lipid(atom_dict = molecule_atoms)


