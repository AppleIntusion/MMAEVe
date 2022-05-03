import manipulate_pdb as manp
import geom_shapes         as gs
import build_micelle as bmi

import copy

class Vesicle(bmi.Micelle):

	def build_outer_surface(self, outer_lipids, outer_frac, outer_radius, outer_points):
		# Create a dictionary of lipid structures.
		struc = dict()
		for ii in range(0, len(outer_lipids)):
			struc[ii] = manp.PdbFile(lipid_name = outer_lipids[ii]).to_lipid()

		# Number of entries to replace in each list.
		lipid_num = [int(prob * outer_points) for prob in outer_frac]

		# A list of points containing the lipid ID to be placed at each
		# point. When initialized it will contain only a single id.
		lipid_id = [0 for ii in range(0, outer_points)]

		# A list from 0 to the number of lipids in the micelle. Represents the
		# indicies of entries in the lipid_id list.
		lipid_id_index = [ii for ii in range(0, outer_points)]

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
		sph_coord = gs.fib_sphere(outer_points, radius = outer_radius)

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

	def build_inner_surface(self, inner_lipids, inner_frac, inner_radius, inner_points):
		# Create a dictionary of lipid structures.
		struc = dict()
		for ii in range(0, len(inner_lipids)):
			struc[ii] = manp.PdbFile(lipid_name = inner_lipids[ii]).to_lipid()

		# Number of entries to replace in each list.
		lipid_num = [int(prob * inner_points) for prob in inner_frac]

		# A list of points containing the lipid ID to be placed at each
		# point. When initialized it will contain only a single id.
		lipid_id = [0 for ii in range(0, inner_points)]

		# A list from 0 to the number of lipids in the micelle. Represents the
		# indicies of entries in the lipid_id list.
		lipid_id_index = [ii for ii in range(0, inner_points)]

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
		sph_coord = gs.fib_sphere(inner_points, radius = inner_radius)

		for ii in range(0, len(sph_coord)):
			# Copy of lipid
			lipid = copy.deepcopy(struc[lipid_id[ii]])

			# Determine which oxygen and hydrogen atom are farthest from each other. 
			extrema = lipid.lipid_extrema()

			# Translate the Lipid such that the position of the hydrogen extreme is 
			# at the origin.
			lipid.trans_axis_o_to_point(extrema, [0, 0, 0])

			# Rotate lipid to align with the position vector on the sphere's surface.
			lipid.align_axis_h_to_vec(extrema, sph_coord[ii])

			# Translate the lipid to the position on the sphere.
			lipid.trans_axis_h_to_point(extrema, sph_coord[ii])

			# Add lipid to list
			self.struc.append(lipid)

		# Convert familythe list of micelle lipids to a list of lists containing strings for 
		# each atom.
#		self.struc_to_lines()

	def build_vesicle(self, \
                      outer_lipids, outer_frac, outer_radius, outer_points, \
                      inner_lipids, inner_frac, inner_radius, inner_points, \
                      outfile):
		self.build_outer_surface(outer_lipids, outer_frac, outer_radius, outer_points)
		self.build_inner_surface(inner_lipids, inner_frac, inner_radius, inner_points)
		self.struc_to_lines()
		self.save_lines(outfile)
		





if __name__ == "__main__":
	vesicle = Vesicle()
	#vesicle_file = "/home/gubbin/Documents/mdProjects/programs/MMAEVe/lipids/test2.pdb"
	test_file = "/home/sanchezw/MMAEVe/lipids/test2.pdb"
	vesicle.build_vesicle(["POC"], [1.0], 100, 100, ["POC"], [1.0], 75, 50, test_file)

	
