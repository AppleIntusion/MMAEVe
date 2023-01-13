'''---------------------------------------------------*
| Title:   Lipid Structure                            |
|                                                     |
| Author:  Gubbin Eel (Satanic Overlord of the Swamp) |
|                                                     |
| Purpose: Defines methods common to all structure    |
|          building procedures.                       |
*---------------------------------------------------'''

# Import program modules
import manipulate_files as manf
import geom_shapes      as gs

# Import outside modules
import numpy as np
import copy
import scipy


class lipidStructure(object):
    ''' 
    Description:
        Used to contain the Lipids and Proteins that compose the 
        structure. Intend to expand it incrementally to include 
        other categories of biomolecules such as Sugars, DNA, RNA, 
        and Small Ligands. None are supported at this time. Has two 
        lists as attributes. One is of lipid objects that make up the 
        structure. The other is a list of those lipids as strings in 
        the format of pdb file lines.

    Attributes: 
        lipids) 1 x N NP Array of Lipid instances. Used to store all
        of the Lipid instances that compose the lipid portion of the 
        system.
        proteins) 1 x N NP Array of Protein instances. Used to store 
        all of the Protein instances that compose the protein portion
        of the system.
        leaf_1_points) Points used to place the Lipid instances of the 
        upper or outer membrane leaflets.
        leaf_2_points) Points used to place the Lipid instances of the
        lower of inner membrane leaflets.
        protein_points) Points used to place Protein instances.
        leaf_1_structures) Lipid instances that are used to build the
        upper or outer membrane leaflets. 
        leaf_2_structures) Lipid instances that are used to build the
        lower or inner membrane leaflets.
        protein_structures) Protein instances that are used during
        system construction.
        leaf_1_ratios) Fraction of the leaf_1_points that correspond
        to each of the leaf_1_structures.
        leaf_2_ratios) Fraction of the leaf_2_points that correspond
        to each of the leaf_2_structures.
        protein_ratios) Fraction of the protein_points that correspond
        to each of the protein_structures.
        leaf_1_ids) Index corresponding to the member of 
        leaf_1_structures that will be used to place a point.
        leaf_2_ids) Index corresponding to the member of 
        leaf_2_structures that will be used to place a point.
        protein_ids) Index corresponding to the member of 
        protein_structures that will be used to place a point.
        lines)
    '''

    def __init__(self, 
                 lipids             = [], 
                 proteins           = [], 
                 leaf_1_points      = np.array([]), 
                 leaf_2_points      = np.array([]), 
                 protein_points     = np.array([]), 
                 leaf_1_structures  = [], 
                 leaf_2_structures  = [], 
                 protein_structures = [], 
                 leaf_1_ratios      = np.array([]), 
                 leaf_2_ratios      = np.array([]), 
                 protein_ratios     = np.array([]), 
                 leaf_1_ids         = [],
                 leaf_2_ids         = [],
                 protein_ids        = [],
                 lines              = []):
        ''' 
        Initialize class instance.
        '''
        self.lipids             = lipids            
        self.proteins           = proteins          
        self.leaf_1_points      = leaf_1_points     
        self.leaf_2_points      = leaf_2_points     
        self.protein_points     = protein_points    
        self.leaf_1_structures  = leaf_1_structures 
        self.leaf_2_structures  = leaf_2_structures 
        self.protein_structures = protein_structures
        self.leaf_1_ratios      = leaf_1_ratios     
        self.leaf_2_ratios      = leaf_2_ratios     
        self.protein_ratios     = protein_ratios    
        self.leaf_1_ids         = leaf_1_ids        
        self.leaf_2_ids         = leaf_2_ids        
        self.protein_ids        = protein_ids       
        self.lines              = lines             

    def collect_lipid_blocks(self, leaf_1 = dict(), leaf_2 = dict()):
        '''
        Purpose:   To import and store the lipid structures and
                   ratios.
        Arguments: self) lipidStructure instance.
                   leaf_1) Dictionary. Key is string and Value is a 
                           Float <=1.0. The key is the name of the 
                           lipid while the value is the proportion of 
                           the Bilayer that the lipid should make up. 
                           All values should add to 1.0. Represents the 
                           entire Micelle, the outer leaflet of a 
                           Vesicle, and the upper leaflet of Nanodiscs 
                           and Bilayers.
                   leaf_2) Dictionary. Key is string and Value is a 
                           Float <=1.0. The key is the name of the 
                           lipid while the value is the proportion of 
                           the Bilayer that the lipid should make up. 
                           All values should add to 1.0. Not used for 
                           Micelles, the inner leaflet of a Vesicle, 
                           and the lower leaflet of Nanodiscs and 
                           Bilayers.
        Requires:  Nothing
        Updates:   leaf_1_structures, leaf_2_structures
        Returns:   Nothing
        '''
        self.leaf_1_ratios = np.array(list(leaf_1.values()))
        self.leaf_2_ratios = np.array(list(leaf_2.values()))
        self.leaf_1_structures = []
        self.leaf_2_structures = []
    
        for lipid in leaf_1:
            pdb =   manf.PdbFile(lipid_name = lipid)
            lipid = pdb.to_lipid()
            lipid.lipid_extrema()
            lipid.add_radius()
            self.leaf_1_structures.append(lipid)
    
        for lipid in leaf_2:
            pdb =   manf.PdbFile(lipid_name = lipid)
            lipid = pdb.to_lipid()
            lipid.lipid_extrema()
            lipid.add_radius()
            self.leaf_2_structures.append(lipid)

    def collect_protein_blocks(self, proteins = dict()):
        '''
        Purpose:   To import and store the protein structures and
                   ratios.
        Arguments: self) lipidStructure instance.
                   proteins) Dictionary. Key is String and Value is 
                   another dictionary. The nested dictionary should 
                   have three Keys "Fraction", "Head", "Tail". The 
                   Keys of the main ditionary are the names of the 
                   proteins that you want to use. For the secondary 
                   Dictionary, "Fraction" should be a Float <= 1.0 
                   with all "Fractions" summing to 1.0. "Head" and 
                   "Tail" are strings. They should be formatted as 
                   described in the protein_extrema method 
                   associated with the the Protein class. 
                   An example is provided below.
                   proteins = {"APO" : {"Fraction" : 0.5 
                                        "Head"     : "123-12-ARG"
                                        "Tail"     : "444-87-GLY"}
                               "ANX" : {"Fraction" : 0.5
                                        "Head"     : "336-38-TYR"
                                        "Tail"     : "901-99-ASP"}}
        Requires:  Nothing
        Modifies:  protein_ratios, protein_structures
        Returns:   Nothing
        '''
        for name in proteins:
            self.protein_ratios = np.append(self.protein_ratios,
                                            proteins[name]["Fraction"])
            pdb = manf.PdbFile(protein_name = name)
            pro = pdb.to_protein()
            pro.generate_lookup()
            pro.protein_extrema(proteins[name]["Head"], 
                                proteins[name]["Tail"])
            pro.add_radius()
            self.protein_structures.append(pro)

    def generate_protein_points(self,                protein_radius = 10.0, 
                                protein_number = 10, pro_width = 10.0,   
                                pro_length = 10.0,   pro_height = 10.0,    
                                structure = "Micelle"):
        '''
        Purpose:   Assign protein points for the various structures.
        Arguments: self) lipidStructure instance.
                   protein_radius) Float. Radius of the sphere, 
                   circular grid, or circle that proteins will be 
                   placed on.
                   protein_number) Integer. Number of points for 
                   proteins to be placed at.
                   pro_width) Float. Width of the rectangular gird
                   for proteins to be placed on.
                   pro_length) Float. Length of the rectangular grid
                   that proteins will be placed on.
                   pro_height) Float. Distance to shift protein 
                   after placement on a circular or rectangular grid.
                   structure) String. The shape of the point 
                   distribution. Options include: "Sphere", "Disc", 
                   "Grid", "Circle"
        Requires:  None
        Modifies:  protein_points
        Returns:   None
        '''
        if structure == "Sphere":
            protein_points = gs.fib_sphere(protein_number, 
                                           protein_radius)
        elif structure == "Disc":
            protein_points = gs.sunflower(protein_number, 
                                          protein_radius)
        elif structure == "Circle":
            protein_points = gs.circle(protein_number, 
                                       protein_radius,
                                       pro_height)
        elif structure == "Grid":
            protein_points = gs.grid(protein_number, pro_length,
                                     pro_width, height = height)
        self.protein_points = protein_points

    def generate_lipid_points(self,                 leaf_1_radius = 10.0, 
                              leaf_2_radius = 10.0, leaf_1_number = 10,   
                              leaf_2_number = 10,   width = 10.0, 
                              length = 10.0,        height = 10.0,         
                              structure = "Micelle"):
        '''
        Purpose:   Assign protein points for the various possible 
                   structures.
        Arguments: self) lipidStructure instance.
                   leaf_1_radius) Float. The radius of the micelle,
                   vesicle outer-leaf, or nanodisc upper-leaf.
                   leaf_2_radius) Float. The radius of the vesicle 
                   inner-leaf or nanodisc lower-leaf.
                   leaf_1_number) Integer. The number of lipids in
                   the vesicle outer-leaf or nanodisc upper-leaf.
                   leaf_2_number) Integer. The number of lipids in
                   the vesicle inner-leaf or nanodisc lower-leaf.
                   width) Float. Width of a bilayer's upper and lower
                   leaflets.
                   length) Float. Length of a bilayer's upper and 
                   lower leaflets.
                   height) Float. How far to elevate the upper level
                   of bilayer points above the lower bilayer points.
                   structure) String. The structure that the 
                   generated points should correspond to: 
        Requires:  None
        Modifies:  leaf_1_points, leaf_2_points 
        Returns:   None
        '''
        if structure == "Micelle":
            leaf_1_points  = gs.fib_sphere(leaf_1_number, 
                                           leaf_1_radius)
            leaf_2_points  = gs.fib_sphere(leaf_2_number, 
                                           leaf_2_radius)
        elif structure == "Vesicle":
            leaf_1_points  = gs.fib_sphere(leaf_1_number, 
                                           leaf_1_radius)
            leaf_2_points  = gs.fib_sphere(leaf_2_number, 
                                           leaf_2_radius)
        elif structure == "Nanodisc":
            leaf_1_points  = gs.sunflower(leaf_1_number, 
                                          leaf_1_radius,
                                          height = height)
            leaf_2_points  = gs.sunflower(leaf_2_number, 
                                          leaf_2_radius)
        elif structure == "Bilayer":
            leaf_1_points  = gs.grid(leaf_1_number, length,
                                     width, height = height)
            leaf_2_points  = gs.grid(leaf_2_number, length,
                                     width)
        self.leaf_1_points =  leaf_1_points
        self.leaf_2_points =  leaf_2_points

    def remove_cylinder(self, axis_points):
        '''
        Purpose:   Assign lipid points for the various structures.
        Arguments: 1) lipidStructure instance.
                   2) N x 4 NP Array of floats. Contains axis 
                      identifier as well as 2D planar coordinates and 
                      the radius of the cylinder. Where 0.0, 1.0, 2.0 
                      are the axis identifiers for x, y, and z, 
                      respectively. For example, [[1.0, 1.0, 2.0, 4.0]] 
                      specifies that the cylinder will be parallel to 
                      the y-axis and pass through the point 
                      (1.0, 0.0, 2.0). All lipids within the specified 
                      distance, 4.0 A of the cylinder, in this case, 
                      will be removed. The two coordinates should be
                      specified in the order x, y, z. For example, If 
                      the point is in the xz-plane then the 
                      x-coordinate should preceed the z-coordinate.
        Requires:     leaf_1_points, leaf_2_points, protein_points
        Modifies:     leaf_1_points, leaf_2_points, protein_points
        Returns:      None.
        '''
        # Used too access the appropriate frames.
        axes = {0.0 : [1, 2], 1.0 : [0, 2], 2.0 : [0, 1]}
        point_info = [[self.leaf_1_points,  self.leaf_1_ids], 
                      [self.leaf_2_points,  self.leaf_2_ids], 
                      [self.protein_points, self.protein_ids]]
        attribute_to_update = [["leaf_1_points",  "leaf_1_ids"], 
                               ["leaf_2_points",  "leaf_2_ids"], 
                               ["protein_points", "protein_ids"]]
        for ii in range(len(point_info)):
            for jj in range(len(axis_points)):
                # Retrieve indicies for coordinate selection
                projection_index  = axes[axis_points[jj, 0]]
                # Point in the plane of projection to use for distance
                # calculation
                plane_point       = axis_points[jj, [1, 2]]
                # Points should be farther than this from the plain point
                radius            = axis_points[jj, 3]
                # Projection of points onto the appropriate plane
                projection_points = point_info[ii][0][:, projection_index]

                # Calculate distance and remove points that are within
                # the defined cylinder.
                distances = projection_points - plane_point
                distances = np.square(distances)
                distances = np.sum(distances, axis = 1)
                distances = np.sqrt(distances)

                point_info[ii][0] = point_info[ii][0][distances > radius]
                point_info[ii][1] = point_info[ii][1][distances > radius]

            setattr(self, attribute_to_update[ii][0], point_info[ii][0])
            setattr(self, attribute_to_update[ii][1], point_info[ii][1])

    def gen_ids(self):
        ''' 
        Purpose:   Assign an ID for each point generated by the 
                   assign_lipid_points or assign_protein_points 
                   methods. These IDs will correspond to which 
                   Protein or Lipid instance will be placed at the 
                   corresponding point.
        Arguments: self) lipidStructure instance.
        Requires:  leaf_1_points, leaf_1_ratios, leaf_1_ids, 
                   leaf_2_points, leaf_2_ratios, leaf_2_ids
        Modifies:  leaf_1_ids, leaf_2_ids, protein_ids
        Returns:   Nothing
        '''
        point_info = [[self.leaf_1_points, self.leaf_1_ratios, 
                       self.leaf_1_ids], 
                      [self.leaf_2_points, self.leaf_2_ratios,
                       self.leaf_2_ids], 
                      [self.protein_points, self.protein_ratios,
                       self.protein_ids]]
        attribute_to_update = ["leaf_1_ids", "leaf_2_ids", 
                               "protein_ids"]
        ii = 0
        for point in point_info:
            number_to_add = len(point[0]) * point[1]
            number_to_add = number_to_add.astype(int)
            point[2] = np.zeros(len(point[0]))
            for jj in range(len(number_to_add)):
                length = number_to_add[jj]
                to_add = np.zeros(length) + jj
                point[2] = np.append(point[2], to_add)
                point[2] = point[2][length:]
            np.random.shuffle(point[2])
            point[2] = point[2].astype(int)
            setattr(self, attribute_to_update[ii], point[2])
            ii += 1

    def struc_to_lines(self):
        ''' 
        Purpose:   Convert the lipids and proteins to pdb file lines
                   and add the lines to the "lines" attribute of
                   a lipidStructure instance.
        Arguments: self) lipidStructure instance.
        Requires:  lipids, proteins
        Modifies:  lines
        Returns:   Nothing
        '''
        atom_number    = 1
        residue_number = 1
        for lipid in self.lipids:
            self.lines += lipid.to_pdb_lines(atom_start = atom_number, 
                                             res_num = residue_number)
            self.lines += ["TER\n"]
            atom_number += len(lipid.serial)
            residue_number += 1

        for protein in self.proteins:
            self.lines += protein.to_pdb_file(
                                     atom_start = atom_number,
                                     resi_start = residue_number)
            self.lines += ["TER\n"]
            protein_atoms = 0
            for res in protein.residues: 
                protein_atoms += len(res.serial)
            atom_number += protein_atoms
            residue_number += len(protein.residues)

    def save_lines(self, file_name):
        ''' 
        Purpose:   To write the file lines in the "lines" attribute of 
                   a lipidStructure instance to a file.
        Arguments: self) lipidStructure instance.
                   file_name) The name of the file to export the .pdb 
                   structure to. Specify the full path in the name.
        Requires:  lines
        Modifies:  Nothing
        Returns:   Nothing. Writes the current pdb file lines for the 
                   system to the specified file.
        '''
        manf.write_file(file_name, self.lines)

    def get_protein_spatial_array(self):
        '''
        Purpose:   To get the xyz coordinates and atomic radius of 
                   every Protein instance in a single N x 3 NP array. 
                   Its intended purpose is to use when checking for 
                   Lipid-Protein overlap. Since only Lipids will be 
                   removed, protein coordinates do not need to be 
                   updated.
        Arguments: self) lipidStructure instance.
        Requires:  proteins
        Modifies:  Nothing
        Returns:   N x 4 NP Array of Floats. The first three floats
                   represent the x, y, and z coordinates, 
                   respectively, of a given atom. The final float
                   represents the atomic radius.
        '''
        xyz_combined = np.array([[0., 0., 0.]])
        radii_combined = np.array([0.])
        for protein in self.proteins:
            for res in protein.residues:
                xyz_combined = np.append(xyz_combined, res.xyz, 
                                         axis = 0)
                radii_combined = np.append(radii_combined, 
                                           res.radius)
        xyz_combined = np.delete(xyz_combined, 0, axis = 0)
        radii_combined = np.delete(radii_combined, 0)

        return np.insert(xyz_combined,   len(xyz_combined[0]), 
                         radii_combined, axis = 1)

    def get_lipid_spatial_array(self):
        '''
        Purpose:   To get the xyz coordinates, atomic radius, and 
                   current index of every Lipid instance in a single 
                   N x 5 NP array. Its intended purpose is to use 
                   when checking for Lipid-Protein overlap.
        Arguments: self) lipidStructure instance.
        Requires:  lipids
        Modifies:  Nothing
        Returns:   N x 5 NP Array of Floats. The first three floats
                   represent the x, y, and z coordinates, 
                   respectively, of a given atom. The fourth float
                   represents the atomic radius. The fifth float
                   is the index of the lipid in the current state of.
                   the 
        '''
        lipid_xyz    = []
        lipid_radius = []
        lipid_num    = []
        ii = 0
        for lipid in self.lipids:
            lipid_xyz.append(lipid.xyz)
            lipid_radius.append(lipid.radius)
            lipid_num.append(lipid.radius - lipid.radius + ii)
            ii += 1

        lipid_xyz    = np.concatenate(lipid_xyz)
        lipid_radius = np.concatenate(lipid_radius)
        lipid_num    = np.concatenate(lipid_num)

        lipid_info = np.insert(lipid_xyz, len(lipid_xyz[0]), 
                               lipid_radius, axis = 1)
        lipid_info = np.insert(lipid_info, len(lipid_info[0]), 
                               lipid_num, axis = 1)

        return(lipid_info)

    def remove_overlap(self, cutoff = 3.0):
        '''
        Purpose:   Remove any Lipid instances that overlap with
                   Protein instances by a specified amount. Overlap
                   is calculated using the atomic center of mass
                   rather than the van der Waals radius. While 
                   slightly suboptimal, it was easier to implement.
                   The minimization module of any Molecular Dynamics 
                   software worth its salt should be able to handle
                   any residual overlap. Lipid-lipid overlap is not 
                   considered as the user has control over the lipid
                   density. If a reasonable density was specified then
                   minimization should take care of any Lipid-Lipid
                   overlap.
        Arguments: self) lipidStructure instance.
                   cutoff) Float. Lipids with atoms within this cutoff
                           threshold from a protein atom will be 
                           removed.
        Requires:  proteins, lipids
        Modifies:  lipids
        Returns:   Nothing
        '''
        # Generate protein coordinate KDTree lookup
        protein_spatial = self.get_protein_spatial_array()
        protein_xyz     = protein_spatial[:, 0:3]
        protein_lookup  = scipy.spatial.KDTree(protein_xyz)

        # Generate lipid coordinate KDTree lookup
        lipid_spatial = self.get_lipid_spatial_array()
        lipid_xyz     = lipid_spatial[:, 0:3]
        lipid_index   = lipid_spatial[:, 4]
        lipid_lookup  = scipy.spatial.KDTree(lipid_xyz)

        # Generate sparse distance matrix
        sparse = protein_lookup.sparse_distance_matrix(lipid_lookup, 3.0, 
                                                       output_type = "ndarray")

        # Remove lipids that overlap with proteins
        hit_index     = list(set([dist[1] for dist in sparse]))
        removal_index = list(set(lipid_index[hit_index]))
        removal_index = [int(ii) for ii in removal_index]
        self.lipids   = np.delete(self.lipids, removal_index, axis = 0)
