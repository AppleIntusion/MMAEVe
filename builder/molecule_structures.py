'''
Title:   Molecule Structures

Author:  Gubbin Eel (Samuel Lindsay, East Carolina University)

Purpose: Defines classes needed to represent molecules along with
         the functions and methods needed to effectively manipulate
         them.

NOTE:    Intended improvements are denoted by the #%# flag.
         Searching for that flag will identify code that needs to be
         updated or will later be improved.
'''

'''-----------------
| Required Modules |
-----------------'''
import geom_shapes as gs

import copy
import numpy as np

'''
Variables
'''
####################################################################
# Atomic masses of the most abundant isotope of each of the most   #
# common elements in biological molecules                          #
#                                                                  #
# Mass data:                                                       #
#     G. Audi, A. H. Wapstra Nucl. Phys A. 1993, 565, 1-65         #
#     G. Audi, A. H. Wapstra Nucl. Phys A. 1995, 595, 409-480      #
#                                                                  #
# Abundance:                                                       #
#     K. J. R. Rosman, P. D. P. Taylor, Pure Appl. Chem. 1999, 71, #
#     1593-1607                                                    #
####################################################################
atomic_masses = {'H': 1.007825, 'He': 4.002603, 'Li': 7.016004, \
                 'Be': 9.012182, 'B': 11.009305, 'C': 12.000000, \
                 'N': 14.003074, 'O': 15.994915, 'F': 18.998403, \
                 'P': 30.973762, 'S': 31.972071}
class Atom(object):
	'''
	Contains important information needed to describe atoms.

	Attributes: serial, x, y, z, mass, elem, name, resn, and resi.
	'''

	''' Initialize class instance. Numerica values are 0 and string values are
	    "NA" unless specified.
	'''
	def __init__(self, serial = 0, x = 0, y = 0, z = 0, mass = 0, resi = 0, resn = "NA", \
                     chain = " ", elem = "NA", name = "NA", charge = 0):
		self.serial = serial
		self.x      = x   
		self.y      = y   
		self.z      = z   
		self.mass   = mass
		self.resi   = resi
		self.resn   = resn
		self.chain  = chain
		self.elem   = elem
		self.name   = name
		self.charge = charge

	''' Purpose:   Translates an atom in the x, y, and z-directions by 
	               amounts specified in a provided list.
	    Arguments: (1) An atom object specified using the method format.
	               (2) A list containing the values each coordinate should
	                   be modified by.
	'''
	def translate(self, t = [0, 0, 0]):
		self.x += t[0]
		self.y += t[1]
		self.z += t[2]

	''' Purpose: To calculate the distance between two atoms.
		Arguments: (1) An atom object apecified using the method format.
	               (2) A second atom object.
	'''
	def distance(self, other):
		dist = (((self.x - other.x) ** 2)  +          \
		        ((self.y - other.y) ** 2)  +          \
		        ((self.z - other.z) ** 2)) ** (1 / 2)
		return dist

	''' Purpose:   To generate a vector in the form of a list that points
	               toward another atom.
	    Arguments: (1) An atom object specified using the method format.
	               (2) A second atom object supplied as an argument.
	'''
	def distance_vec_atom(self, other):
		dist_vec = [other.x - self.x,
					other.y - self.y,
					other.z - self.z]
		return dist_vec

	''' Purpose:   To generate a vector in the form of a list that points
	               toward a point.
	    Arguments: (1) An atom object specified using the method format.
	               (2) A list or other indexable data structure with x, y,
	                   and z coordinates specified. 
	'''
	def distance_vec_point(self, point):
		dist_vec = [point[0] - self.x,
					point[1] - self.y,
					point[2] - self.z]
		return dist_vec

	''' Purpose:   To rotate an Atom object about the x-axis.
		Arguments: (1) An atom object specified using the method format.
	               (2) Angle by which to rotate the molecule, in radians.
	'''
	def rot_x(self, thetax = 0):
		y = self.y
		z = self.z
		self.y = (np.cos(thetax) * y) + (-1 * (np.sin(thetax) * z))
		self.z = (np.sin(thetax) * y) + (np.cos(thetax) * z)

	''' Purpose:   To rotate an Atom object about the y-axis.
		Arguments: (1) An atom object specified using the method format.
	               (2) Angle by which to rotate the molecule, in radians.
	'''
	def rot_y(self, thetay = 0):
		x = self.x
		z = self.z
		self.x = (np.cos(thetay) * x) + (np.sin(thetay) * z)
		self.z = (np.sin(thetay) * x * -1) + (np.cos(thetay) * z)

	''' Purpose:   To rotate an Atom object about the z-axis.
		Arguments: (1) An atom object specified using the method format.
	               (2) Angle by which to rotate the molecule, in radians.
	'''
	def rot_z(self, thetaz = 0):
		x = self.x
		y = self.y
		self.x = (np.cos(thetaz) * x) - (np.sin(thetaz) * y) 
		self.y = (np.sin(thetaz) * x) + (np.cos(thetaz) * y)

	''' Purpose:    To retrieve the coordinates of an Atom object as a list.
		Arguments: (1) An atom object specified using the method format.
	'''
	def to_pos_vec(self):
		return [self.x, self.y, self.z]

	''' Purpose:   To format an atom object as a list of PDB file lines.
	    Arguments: (1) An instance of the lipid class specified using the
	               class method format.
	               (2) The atom serial number as an integer. If it is 0
	                   the the value already assigned will be used.
	               (3) The residue sequence number as an integer. If it is
	                   0 the value already assigned will be used.
	'''
	def to_pdb_file(self, atom_serial = '0', resn_seq = '0', chain=''):
		# All the required entries for an ATOM field entry in a pdb file.
		record_name = "ATOM"
		if atom_serial == '0':
			atom_serial = str(self.serial)
		else:
			pass
		if int(atom_serial) > 99999:
			atom_serial = "99999"
		atom_name = self.name
		alternate_loc = ' '
		resn_name = self.resn 
		chain = self.chain
		if resn_seq == '0':
			resn_seq = str(self.resi)
		else:
			pass
		inser_code = ' '
		x = str(format(round(self.x, 3), ".3f"))
		y = str(format(round(self.y, 3), ".3f"))
		z = str(format(round(self.z, 3), ".3f"))
		occ = '1.00'
		tmp_fac = '0.00'
		seg_iden = ' '
		elem = self.elem
		charge = str(self.charge)
		line = []
		# Record name
		line.append(record_name + ((6 - len(record_name)) * ' '))
		# Atom serial number
		line.append(((5 - len(atom_serial)) * ' ') + atom_serial)
		# Atom name
		line.append(' ' + atom_name + ((4 - len(atom_name)) * ' ')) #DEBUG line.append('  ' + atom_name + ((3 - len(atom_name)) * ' '))
		# Alternate location
		line.append(alternate_loc)
		# Residue name
		line.append(resn_name + ' ')
		# Chain identifier
		line.append(chain)
		# Residue sequence number
		line.append(((4 - len(resn_seq)) * ' ') + resn_seq)
		# Code for the insertion of residues
		line.append(inser_code)
		# x-coordinate
		line.append('   ' + ((8 - len(x)) * ' ') + x)
		# y-coordinate
		line.append(((8 - len(y)) * ' ') + y)
		# z-coordinate
		line.append(((8 - len(z)) * ' ') + z)
		# Occupancy
		line.append(((6 - len(occ)) * ' ') + occ)
		# Temperature factor
		line.append(((6 - len(tmp_fac)) * ' ') + tmp_fac)
		# Segment identifier
		line.append('   ' + seg_iden + ((4 - len(seg_iden)) * ' '))
		# Elemental symbol
		line.append(((2 - len(elem)) * '    ') + elem)
		# Charge
		line.append(((2 - len(charge)) * ' ') + charge)
		line.append('\n')
		line = ''.join(line)
		return line

class Lipid(object):
	'''
	Contains the constituent atoms of single molecule along with their
	important identifiers.	

	Attributes: A dictionary composed of Atom instances that represents the constituent
	            atoms of a molecule.
	'''

	''' Initializes class instance assigning an empty dictionary. 
	'''
	def __init__(self, atom_dict = dict()):
		self.atom_info = atom_dict

	''' Purpose:   Translates a molecule in the x, y, and z directions
	               using a provided list.
	    Arguments: (1) An instance of the lipid class specified using the
	                   class method format.
	               (2) A list of three floats used to translate the lipid 
	                   molecule.
	'''
	def translate(self, t = [0, 0, 0]):
		for atom in self.atom_info:
			self.atom_info[atom].translate(t = t)

	''' Purpose:   To rotate a molecule by a specified amount around either
	               the x, y, or z-axis.
		Arguments: (1) An instance of the lipid class specified using the
	                   class method format.
	               (2) The axis to rotate about: 'x', 'y', or 'z'. 
	               (3) How far to rotate, in radians; default of 0.
	'''
	def rotate(self, axis = 'x', theta = 0):
			for atom in self.atom_info:
				if axis == 'x':
					self.atom_info[atom].rot_x(thetax = theta)
				if axis == 'y':
					self.atom_info[atom].rot_y(thetay = theta)
				if axis == 'z':
					self.atom_info[atom].rot_z(thetaz = theta)

	''' Purpose:   To determine which hydrogen and oxygen atom in the lipid
	               molecule are farthest from each other and provide a
	               dictionary that contains the key to identify them in the
	               molecule dictionary for the Lipid instance.
	    Arguments: (1) An instance of the lipid class specified using the
	                   class method format.
	'''
	def lipid_extrema(self): 
		o_atoms = [] # List to hold Atom class instances for oxygen atoms.
		h_atoms = [] # List to hold Atom class instances for hydrogen atoms.
		# Retrieve class instances from the dictionary.
		for atom in self.atom_info:
			if self.atom_info[atom].elem == "O":
				o_atoms.append(copy.deepcopy(self.atom_info[atom]))
			elif self.atom_info[atom].elem == "H":
				h_atoms.append(copy.deepcopy(self.atom_info[atom]))
			else:
				pass

		# List of lists to contain the atomic serial numbers and their distance.
		distance_comp = []
		o_len, h_len = len(o_atoms), len(h_atoms)
		for oxy in o_atoms:
			for hyd in h_atoms:
				dist    = oxy.distance(hyd) # Distance between h and o atoms.
				o_index = oxy.serial # Molecule dictionary key.
				h_index = hyd.serial # Molecule dictionary key.
				distance_comp.append([o_index, h_index, dist])

		# Sort the list from largest to shortest distance.
		distance_comp = sorted(distance_comp, key = lambda x: x[2])
		distance_comp = distance_comp[-1]
		distance_comp = {"O":    distance_comp[0], "H": distance_comp[1], 
		                 "dist": distance_comp[2]}
		return distance_comp

	''' Purpose:   Translate the axis hydrogen of a lipid to a point. The origin
	               when building micelles and vesicles.
	    Arguments: (1) Lipid class instance specified in the class method 
	                   format.
	               (2) Extrema dictionary provided by the lipid_extrema lipid
	                   class method.
	               (3) Point to translate to.
	'''
	def trans_axis_h_to_point(self, extrema, point):
		# Get axis hydrogen index.
		h_ext = extrema["H"]
		# Distance vector required to translate a given atom to a point.
		distance_vec = self.atom_info[h_ext].distance_vec_point(point) 
		# Translate the lipid along a provided vector.
		self.translate(t = distance_vec)

	''' Purpose:   Align a lipid with a given position vector. The axis
                   hydrogen remains in its current position. Usually the
	               origin when building micells or vesicles.
                   surface of the sphere. 
	    Arguments: (1) Lipid class instance specified in the class method 
	                   format.
	               (2) Extrema dictionary provided by the lipid_extrema lipid
	                   class method.
	               (3) Point to translate to.
	'''
	def align_axis_h_to_vec(self, extrema, point):
		# Get axis oxygen index.
		h_ext = extrema["H"]
		# Rotate the vector about the z-axis so that it is in the same plane as the
		# sphere vector.
		# Find the angle between the lipid vector and the sphere vector when
		# they are both projected onto the xy-plane. As well as the angles made
		# with the x_axis.
		sphere_vec = [point[0], point[1], 0]
		lipid_vec = [self.atom_info[h_ext].x, self.atom_info[h_ext].y, 0]
		angle = gs.directional_angle(lipid_vec, sphere_vec, [0, 0, 1])
		# Rotate by the specified amount about the y-axis.
		self.rotate(axis = 'z', theta = angle)

		# Rotate the vector about the y-axis so that it is in the same plane as the
		# sphere vector.
		# Find the angle between the lipid vector and the sphere vector when
		# they are both projected onto the xz-plane. As well as the angles made with
		# the x-axis.
		sphere_vec = [point[0], 0, point[2]]
		lipid_vec = [self.atom_info[h_ext].x, 0, self.atom_info[h_ext].z]
		angle = gs.directional_angle(lipid_vec, sphere_vec, [0, 1, 0])
		# Rotate by the specified amount about the y-axis.
		self.rotate(axis = 'y', theta = angle)

		# Rotate the vector about the x-axis so that it is in the same plane as the
		# sphere vector.
		# Find the angle between the lipid vector and the sphere vector when
		# they are both projected onto the yz-plane. As well as the angles made with
		# the y-axis.
		sphere_vec = [0, point[1], point[2]]
		lipid_vec = [0, self.atom_info[h_ext].y, self.atom_info[h_ext].z]
		angle = gs.directional_angle(lipid_vec, sphere_vec, [1, 0, 0])
		# Rotate by the specified amount about the y-axis.
		self.rotate(axis = 'x', theta = angle)

	''' Purpose:   Align a lipid with a given position vector. The axis
                   hydrogen remains in its current position. Usually the
	               origin when building micells or vesicles.
                   surface of the sphere. 
	    Arguments: (1) Lipid class instance specified in the class method 
	                   format.
	               (2) Extrema dictionary provided by the lipid_extrema lipid
	                   class method.
	               (3) Point to translate to.
	'''
	def align_axis_o_to_vec(self, extrema, point):
		# Get axis oxygen index.
		o_ext = extrema["O"]
		# Rotate the vector about the z-axis so that it is in the same plane as the
		# sphere vector.
		# Find the angle between the lipid vector and the sphere vector when
		# they are both projected onto the xy-plane. As well as the angles made
		# with the x_axis.
		sphere_vec = [point[0], point[1], 0]
		lipid_vec = [self.atom_info[o_ext].x, self.atom_info[o_ext].y, 0]
		angle = gs.directional_angle(lipid_vec, sphere_vec, [0, 0, 1])
		# Rotate by the specified amount about the y-axis.
		self.rotate(axis = 'z', theta = angle)

		# Rotate the vector about the y-axis so that it is in the same plane as the
		# sphere vector.
		# Find the angle between the lipid vector and the sphere vector when
		# they are both projected onto the xz-plane. As well as the angles made with
		# the x-axis.
		sphere_vec = [point[0], 0, point[2]]
		lipid_vec = [self.atom_info[o_ext].x, 0, self.atom_info[o_ext].z]
		angle = gs.directional_angle(lipid_vec, sphere_vec, [0, 1, 0])
		# Rotate by the specified amount about the y-axis.
		self.rotate(axis = 'y', theta = angle)

		# Rotate the vector about the x-axis so that it is in the same plane as the
		# sphere vector.
		# Find the angle between the lipid vector and the sphere vector when
		# they are both projected onto the yz-plane. As well as the angles made with
		# the y-axis.
		sphere_vec = [0, point[1], point[2]]
		lipid_vec = [0, self.atom_info[o_ext].y, self.atom_info[o_ext].z]
		angle = gs.directional_angle(lipid_vec, sphere_vec, [1, 0, 0])
		# Rotate by the specified amount about the y-axis.
		self.rotate(axis = 'x', theta = angle)

	''' Purpose:   Translate the axis oxygen of a lipid to a point.
	    Arguments: (1) Lipid class instance specified in the class method 
	                   format.
	               (2) Extrema dictionary provided by the lipid_extrema lipid
	                   class method.
	               (3) Point to translate to.
	'''
	def trans_axis_o_to_point(self, extrema, point):
		# Get axis oxygen index.
		o_ext = extrema["O"]
		# Distance vector required to translate oxygen to a point.
		distance_vec = self.atom_info[o_ext].distance_vec_point(point)
		# Translate the lipid to the position on the sphere.
		self.translate(t = distance_vec)

class Protein(Lipid):
	'''
	Inherited Lipid Class
	'''

	
