# This script implements the 4-Point Congruent Sets (4PlCS) algorithm for coarse registration of a point cloud (PCD) to an IFC model.
# The process is divided into four main steps:
# 1. Load a PCD, extract planar patches, and identify valid 4PlCSs.
# 2. Load an IFC model, extract planar patches, and identify valid 4PlCSs.
# 3. Match the 4PlCSs from the PCD and IFC to determine the most likely transformation.
# 4. Refine the transformation using the Iterative Closest Point (ICP) algorithm.
#
# Key assumptions:
# - The geometry of both the PCD and IFC models is in meters.
# - The Z-axis of the point cloud is oriented vertically upwards, similar to the IFC model.

# Import required libraries
import ifcopenshell
import ifcopenshell.geom
import ifcopenshell.util
import ifcopenshell.util.shape
import open3d as o3d
import numpy as np
import math
import os
import copy
import time
import re
import itertools
from typing import List, Tuple
from collections import defaultdict
from scipy.spatial.transform import Rotation as R
import argparse
np.set_printoptions(precision=2, suppress=True)

def run_4plcs(ifc_path, pcd_path):
    # Set the type of PCD. This allows for different parameter sets.
    PCDTYPE = "TLS"
    #PCDTYPE = "Other" 
    print(f"PCDTYPE is set to '{PCDTYPE}'")


    ###############################################################################
    # # 1. Extract PCD 4PlCSs
    ###############################################################################

    # Read Input pointcloud
    pcd = o3d.io.read_point_cloud(pcd_path)
    print(f"Input pcd {pcd}")

    # Get the base name of the pcd file for the output file name
    pcd_filename = os.path.splitext(os.path.basename(pcd_path))[0]

    # Downsample pcd
    if PCDTYPE == "TLS":
        pcd_voxel_sampling_size = 0.05
    else:
        pcd_voxel_sampling_size = 0.05

    downpcd_tmp = pcd.voxel_down_sample(pcd_voxel_sampling_size)
    print(f"Point pcd voxel sampling size: {pcd_voxel_sampling_size}.")
    print(f"Input downsampled PCD: {downpcd_tmp}.")

    sampling_ratio = float(len(downpcd_tmp.points)) / float(len(pcd.points))
    sampling_count = int(round(1/sampling_ratio))
    downpcd_before_transform = pcd.uniform_down_sample(sampling_count)
    print(f"Point pcd sampling ratio: {sampling_ratio} (sampling cout = {sampling_count}).")
    print(f"Input uniformly downsampled PCD: {downpcd_before_transform}.")


    # Calculate pcd normals
    if PCDTYPE == "TLS":
        my_radius = 0.20
    else:
        my_radius = 0.20

    downpcd_before_transform.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=my_radius, max_nn=30))

    assert (downpcd_before_transform.has_normals()), "Point cloud should have normals, but doesn't"

    # o3d.visualization.draw_geometries([downpcd_before_transform], window_name="Downsampled Point Cloud")


    ##########################################################################
    # ## Known Expected Location:
    # Do we know the expected location? 
    # If no, then set `KNOWN_TARGET_LOCATION` to `False`.
    # If yes, then set `KNOWN_TARGET_LOCATION` to `True` and give the location coordinate to `pcd_origin_target`.
    
    KNOWN_TARGET_LOCATION = False

    print(f"KNOWN_TARGET_LOCATION is set to '{KNOWN_TARGET_LOCATION}'")

    pcd_origin_init = [0.0, 0.0, 0.0]
    print(f"pcd_origin_init: {pcd_origin_init}")

    if KNOWN_TARGET_LOCATION:
        pcd_origin_target = [-3.4854, 1.5203, 1.00]
        print(f"pcd_origin_target: {pcd_origin_target}")

    ##########################################################################
    # ## Apply a random transformation to the PCD.
    # This is useful for testing purposes if the input data is already aligned.
    # Set APPLY_RANDOM_TRANSFORM to False to disable.

    APPLY_RANDOM_TRANSFORM = False
    print(f"APPLY_RANDOM_TRANSFORM is set to '{APPLY_RANDOM_TRANSFORM}'")

    random_transform = np.eye(4)

    if APPLY_RANDOM_TRANSFORM == True:
        random_euler_deg = [0.0, 0.0, 33.0]
        random_translation = [32.0, 12.5, 4.4]

        # Build the 4x4 transformation matrix
        random_rotation = R.from_euler('xyz', random_euler_deg, degrees=True).as_matrix()
        random_transform[:3, :3] = random_rotation
        random_transform[:3, 3] = random_translation

        # Visualize before and after random transformation with IFC mesh
        original_pcd_viz = copy.deepcopy(downpcd_before_transform)
        transformed_pcd_viz = copy.deepcopy(downpcd_before_transform)
        transformed_pcd_viz.transform(random_transform)

        # Color the point clouds for clarity
        original_pcd_viz.paint_uniform_color([0, 0, 1])      # Blue: original
        transformed_pcd_viz.paint_uniform_color([1, 0, 0])   # Red: transformed

        o3d.visualization.draw_geometries([original_pcd_viz, transformed_pcd_viz])        
        
        print("Random Transform:")
        print(random_transform)

    downpcd = downpcd_before_transform.transform(random_transform)


    # ## Extract planar patches from the PCD.
    # This section detects planar patches in the point cloud using Open3D's `detect_planar_patches`.
    # The detected patches are then filtered based on minimum area and maximum thickness.
    # A subsampled version of each patch's point cloud is also stored.

    stime_pcd_patches = time.time()

    if PCDTYPE == "TLS":
        extract_patch_area_min = 0.25 #m^2
        point_area = pow(pcd_voxel_sampling_size, 2)
        my_min_plane_edge_length = 0.4
        normal_variance_threshold_deg = 10.0
        coplanarity_deg = 87.0
        outlier_ratio = 0.85
        min_num_points = int(round(extract_patch_area_min / point_area))
    else: # Other
        extract_patch_area_min = 0.25 #m^2
        point_area = pow(pcd_voxel_sampling_size, 2)
        my_min_plane_edge_length = 0.4
        normal_variance_threshold_deg = 10.0
        coplanarity_deg = 87.0
        outlier_ratio = 0.85
        min_num_points = int(round(extract_patch_area_min / point_area))

    oboxes = downpcd.detect_planar_patches(
        normal_variance_threshold_deg = normal_variance_threshold_deg,
        coplanarity_deg = coplanarity_deg,
        outlier_ratio = outlier_ratio,
        min_plane_edge_length = my_min_plane_edge_length,
        min_num_points = min_num_points,
        search_param = o3d.geometry.KDTreeSearchParamKNN(knn=30))

    print(f"Number of patches: {len(oboxes)}.")

    # Create a list of patches, filtering out those that are too small or too thick.
    print("Create patches list with bounding boxes and point clouds.")
    if PCDTYPE == "TLS":
        patch_thickness_max = 0.10
        patch_area_min = 1.0
    else: # Other
        patch_thickness_max = 0.10
        patch_area_min = 1.0

    print(f"Minimum patch area set to {patch_area_min} m^2, and minimum patch thickness set to {patch_thickness_max} m.")

    patch_sampling_count = 10
    patch_sampling_max = 1000
    print(f"Patch sampling: Count={patch_sampling_count}; Max={patch_sampling_max}")

    pcd_patches = []
    pcd_patches_too_small = []
    pcd_patches_too_thick = []

    for idb, patch_obox in enumerate(oboxes):
        mesh = o3d.geometry.TriangleMesh.create_from_oriented_bounding_box(patch_obox)
        mesh.paint_uniform_color(patch_obox.color)

        patch_indices = patch_obox.get_point_indices_within_bounding_box(downpcd.points)
        patch_pcd = downpcd.select_by_index(patch_indices)
        
        # Discard small patches
        patch_area = len(patch_pcd.points) * point_area
        if patch_area < patch_area_min:
            pcd_patches_too_small.append(patch_obox)  
            continue

        # Discard thick patches
        if np.min(patch_obox.extent) > patch_thickness_max:
            pcd_patches_too_thick.append(patch_obox)
            continue
        
        patch_pcd.paint_uniform_color(patch_obox.color)
        patch_center = np.asarray(patch_pcd.points).mean(axis=0)
        patch_normal = patch_obox.R[:,2]
        d = -np.dot(patch_normal, patch_center)
        a, b, c = patch_normal
        patch_model = [a, b, c, d]
        patch_downpcd = patch_pcd.uniform_down_sample(patch_sampling_count)
        patch_downpcd_size = len(patch_downpcd.points)
        if patch_downpcd_size > patch_sampling_max:
            indices = np.random.choice(patch_downpcd_size, patch_sampling_max, replace=False)
            patch_downpcd = patch_downpcd.select_by_index(indices)
        pcd_patches.append({"obox": patch_obox, 
                            "pcd": patch_pcd,
                            "downpcd": patch_downpcd,
                            "center": patch_center, 
                            "normal": patch_normal, 
                            "model": patch_model, 
                            "area": patch_area})

    # Sort all patches by area in descending order.
    pcd_patches = sorted(pcd_patches, key=lambda p: p['area'], reverse=True)

    etime_pcd_patches = time.time()
    duration_pcd_patches = etime_pcd_patches - stime_pcd_patches

    print(f"Number of patches extracted from pcd that are too small: {len(pcd_patches_too_small)}")
    print(f"Number of patches extracted from pcd that are too thick: {len(pcd_patches_too_thick)}")
    print(f"Number of patches extracted from pcd that are valid: {len(pcd_patches)}")
    print(f"Time to extract patches: {duration_pcd_patches:.2f} s")

    # Visualize the detected PCD patches.
    print("Draw detected patches, horizontal and/or no-horizontal, with meshes and/or point clouds")

    horizontal = False
    not_horizontal = False # Set to False to speed up visualization
    print(f"Horizontal planes displayed: {horizontal}; Non-horizontal patches displayed: {not_horizontal}")
    display_obox = False
    display_points = True

    geometries_pcd = []
    z_axis = [0, 0, 1.0]

    for patch in pcd_patches:
        patch_obox = patch["obox"]
        patch_pcd = patch["downpcd"]
        patch_normal = patch["normal"]

        abs_dotproduct = np.absolute(np.dot(patch_normal, z_axis))

        if not_horizontal == True and abs_dotproduct < 0.8 :
            if display_obox == True:
                geometries_pcd.append(patch_obox)

            if display_points == True:
                geometries_pcd.append(patch_pcd)

        elif horizontal == True and abs_dotproduct > 0.8 :
            if display_obox == True:
                geometries_pcd.append(patch_obox)

            if display_points == True:
                geometries_pcd.append(patch_pcd)

    # o3d.visualization.draw_geometries(geometries_pcd, window_name="PCD Patches")



    # ## Extract 4-Point Congruent Sets (4PlCSs) from the PCD patches.
    # A valid 4PlCS consists of four planar patches where:
    # - Two patches are parallel to each other.
    # - The other three patches are not pairwise parallel.
    # This section filters the combinations of patches to find valid and useful 4PlCSs.

    def angle_between_vectors_deg(n1: np.ndarray, n2: np.ndarray) -> float:
        """Returns angle in degrees between two normals."""
        cos_theta = np.clip(np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2)), -1.0, 1.0)
        angle = np.arccos(cos_theta)
        return np.degrees(angle)

    def are_normals_parallel(n1: np.ndarray, n2: np.ndarray, angle_thresh=5.0) -> bool:
        angle = angle_between_vectors_deg(n1, n2)
        return angle < angle_thresh or 180 - angle < angle_thresh

    def are_normals_parallel(angle_deg: float, angle_thresh=5.0) -> bool:
        return angle_deg < angle_thresh or (180 - angle_deg) < angle_thresh

    def are_coplanar(plane1: np.ndarray, plane2: np.ndarray, angle_thresh=5.0, dist_thresh=0.05) -> bool:
        """
        Checks if two planes are (nearly) coplanar: same normal direction and same offset from origin.
        """   
        n1 = plane1[:3]
        n2 = plane2[:3]

        if are_normals_parallel(n1, n2, angle_thresh):
            # Use normalized offset to compare distances
            d1 = plane1[3] / np.linalg.norm(n1)
            d2 = plane2[3] / np.linalg.norm(n2)

            dist_diff = abs(d1 - d2)
            return dist_diff < dist_thresh

    def intersect_planes(plane1: np.ndarray, plane2: np.ndarray, plane3: np.ndarray) -> np.ndarray:
        """
        Solves for the intersection point of 3 planes using linear algebra.
        """
        A = np.vstack([plane1[:3], plane2[:3], plane3[:3]])
        b = -np.array([plane1[3], plane2[3], plane3[3]])
        try:
            return np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            return None  # Planes do not intersect at a single point (e.g., parallel)

    def find_valid_patch_combinations(patches: list, parallel_angle_thresh: float = 5.0,
        coplanar_dist_thresh: float = 0.05, plane_dist_thresh: float = 3.0) -> List[Tuple[int, int, int, int]]:
        """
        Finds all valid 4-plane combinations from the given list of patches.

        Parameters:
        - patches: List of patches.
        - coplanar_threshold: Max allowed distance between parallel patches to consider them coplanar.
        - parallel_threshold: Max angle difference (in radians) to consider normals parallel.

        Returns:
        - List of valid combinations (each as a tuple of indices).
        """
        
        print(f"There are {len(patches)} patches.")

        # First compute all pairwise angles and distances:
        pairwise_angles = np.full((len(patches), len(patches)), -1.0, dtype=float)
        pairwise_distances = np.full((len(patches), len(patches)), -1.0, dtype=float)

        for pair in itertools.combinations(range(len(patches)), 2):
            i, j = pair
            patch_normal1 = patches[i]["normal"]
            patch_center1 = patches[i]["center"]

            patch_normal2 = patches[j]["normal"]
            patch_center2 = patches[j]["center"]
            
            pairwise_angles[i, j] = angle_between_vectors_deg(patch_normal1, patch_normal2)

            pairwise_distances[i, j] = np.absolute(np.dot(patch_center2 - patch_center1, patch_normal1))
        
        print(f"Pairwise angles and distances between all patches are pre-calculated.")
        
        # Now, create all sets of valid combinations
        valid_triplet = 0
        valid_combinations = []
        notideal_combinations = []
        notvalid_combinations = []

        combinations = itertools.combinations(range(len(patches)), 4)

        print(f"There are {len(list(combinations))} combinations of patches.")

        for comb in itertools.combinations(range(len(patches)), 4):
            
            nonparallel_triplet_found = False
            
            for triplet in itertools.combinations(comb, 3):
                parallel_patch_found = False
                for a, b in itertools.combinations(triplet, 2):
                    pairwise_angle = pairwise_angles[a, b]
                    if are_normals_parallel(pairwise_angle, parallel_angle_thresh):
                        parallel_patch_found = True
                        break
                if parallel_patch_found == False:
                    nonparallel_triplet_found = True
                    break
            if nonparallel_triplet_found == False:
                notvalid_combinations.append(comb)
                continue

            valid_triplet += 1

            # Check if two patches are nearly parallel and sufficiently distant or co-planar
            coplanar_found = False
            valid_parallel_pair_found = False

            for pair in itertools.combinations(comb, 2):
                i, j = pair
                pairwise_angle_ij = pairwise_angles[i, j]

                # Check if parallel:
                if are_normals_parallel(pairwise_angle_ij, parallel_angle_thresh):
                    distance_ij = pairwise_distances[i, j]
                    # Check if co-planar:
                    if distance_ij < coplanar_dist_thresh:
                        coplanar_found = True
                        break
                    # Check if parallel planes are far enough from each other.
                    elif  distance_ij > plane_dist_thresh: 
                        valid_parallel_pair_found = True
                        break  # No need to check other pairs
                                
            if coplanar_found:
                notvalid_combinations.append(comb)
                continue
            if not valid_parallel_pair_found:
                notideal_combinations.append(comb)
                continue

            # Calculate the all relative patch angles and patch orthogonal distance between i and j
            sublist = [x for x in comb if x != i and x != j]
            k, l = sublist
            angle_ik = pairwise_angles[i, k]
            angle_il = pairwise_angles[i, l]
            angle_kl = pairwise_angles[k, l]
            description = [angle_ik, angle_il, angle_kl, distance_ij]

            area = patches[i]["area"] + patches[j]["area"] + patches[k]["area"] + patches[l]["area"]
            
            comb_ordered = [i, k ,l, j]

            # Add plane combination to list of valid plane combinations, i.e. is a 4PlCS
            valid_combinations.append({"patch_ids": comb_ordered, "patch_geometry": description, "area": area})
        
        print(f"There are {valid_triplet} valid patch triplets.")
        print(f"There are {len(notvalid_combinations)} not valid 4PlCSs.")
        print(f"There are {len(notideal_combinations)} not ideal 4PlCSs.")
        print(f"There are {len(valid_combinations)} valid 4PlCSs.")
        
        return valid_combinations


    def find_largest_4PlCSs (a4PlCSs: list, patches: list, max_patches: int):
        
        # Select the max_patches patches with the largest area (patches is already sorted according to "area")
        patches_filtered = patches[:max_patches]
        
        # Select the 4PlCSs that contain at least one of the largest patches
        patches_filtered_ids = range(len(patches_filtered))
        filtered_4PlCSs = [a4PlCS for a4PlCS in a4PlCSs if any(id_ in patches_filtered_ids for id_ in a4PlCS['patch_ids'])] # type: ignore
        print (f"There are {len(filtered_4PlCSs)} filtered 4PlCSs.")

        return filtered_4PlCSs

    # Extract PCD 4PlCSs
    stime_pcd_4PlCSs = time.time()

    if PCDTYPE == "TLS":
        coplanar_angle_thresh = 5.0
        coplanar_dist_thresh = 0.05
        plane_dist_thresh = 1.5
    else:
        coplanar_angle_thresh = 5.0
        coplanar_dist_thresh = 0.05 
        plane_dist_thresh = 1.5

    pcd_4PlCSs = find_valid_patch_combinations(pcd_patches, coplanar_angle_thresh, coplanar_dist_thresh, plane_dist_thresh)


    # Select the top largest PCD 4PlCSs for matching.
    # NOTE: This could be improved by using a minimum area threshold instead of a fixed number.
    # --------------------------------
    pcd_patches_max = 20
    pcd_4PlCSs_largest = find_largest_4PlCSs(pcd_4PlCSs, pcd_patches, pcd_patches_max)
    pcd_4PlCSs_largest  = sorted(pcd_4PlCSs_largest , key=lambda p: p['area'], reverse=True)

    # Filter to the top N 4PlCSs
    pcd_4PlCSs_max = 100
    pcd_4PlCSs_filtered = pcd_4PlCSs_largest[:pcd_4PlCSs_max]
    print (f"There are {len(pcd_4PlCSs_filtered)} filtered 4PlCSs in the pcd.")

    etime_pcd_4PlCSs = time.time()
    duration_pcd_4PlCSs = etime_pcd_4PlCSs - stime_pcd_4PlCSs
    print(f"Time to extract 4PlCSs: {duration_pcd_4PlCSs:.2f} s")



    ###############################################################################
    # # 2. Loading and Processing IFC Model
    ###############################################################################

    # ## Load the IFC file and convert its geometry to Open3D meshes.
    # We load an IFC file and keep only the elements of certain classes (currently "IfcWall", "IfcSlab" and "IfcBeam").
    def ifc_shape_to_open3d_mesh(shape):
        vertices = np.array(shape.geometry.verts, dtype=np.float64).reshape(-1, 3)
        triangles = np.array(shape.geometry.faces, dtype=np.int32).reshape(-1, 3)

        mesh = o3d.geometry.TriangleMesh()
        mesh.vertices = o3d.utility.Vector3dVector(vertices)
        mesh.triangles = o3d.utility.Vector3iVector(triangles)
        mesh.compute_triangle_normals()
        return mesh

    # Define the element types for the analysis
    filtered_element_types = ["IfcWall","IfcSlab","IfcBeam"]

    # Load ifc file
    print(f"Input IFC model path: {ifc_path}")

    if not os.path.exists(ifc_path):
        print(f"Error: File not found at {ifc_path}")
    else:
        print(f"Opening file: {ifc_path}")

    # Load the IFC model and filter objects by specified element type
    ifc_file = ifcopenshell.open(ifc_path)

    # Initialize geometry processing settings and iterator
    tree = ifcopenshell.geom.tree()
    settings = ifcopenshell.geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)
    #settings.set(settings.SEW_SHELLS, True)

    # Initialize Open3D structures
    elements = []
    for type_name in filtered_element_types:
        elements.extend(ifc_file.by_type(type_name))

    # Convert elements to triangle meshes
    all_triangles = []
    all_meshes = []

    for element in elements:
        shape = ifcopenshell.geom.create_shape(settings, element)
        element_mesh = ifc_shape_to_open3d_mesh(shape)
        element_mesh.compute_triangle_normals()
        all_meshes.append(element_mesh)

    # Plot loaded IFC file
    # o3d.visualization.draw_geometries(all_meshes, mesh_show_back_face=True, window_name="IFC Model")

    # ## Extract planar patches from the IFC mesh.
    # This section extracts planar patches from the IFC mesh using a region-growing approach based on triangle normal similarity.
    # The detected patches are then filtered based on a minimum area.
    def build_triangle_adjacency(mesh):
        triangles = np.asarray(mesh.triangles)
        edge_to_triangles = defaultdict(list)

        # For each triangle, map its edges to the triangle index
        for idx, tri in enumerate(triangles):
            edges = [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])]
            for edge in edges:
                edge_sorted = tuple(sorted(edge))  # order doesn't matter for undirected edges
                edge_to_triangles[edge_sorted].append(idx)

        # Build adjacency list
        adjacency = [[] for _ in range(len(triangles))]
        for edge, tris in edge_to_triangles.items():
            if len(tris) == 2:
                t1, t2 = tris
                adjacency[t1].append(t2)
                adjacency[t2].append(t1)

        return adjacency

    def extract_planar_patches(mesh, angle_threshold_deg=10):
        triangles = np.asarray(mesh.triangles)
        triangle_normals = np.asarray(mesh.triangle_normals)
        
        # Build adjacency list
        adjacency = build_triangle_adjacency(mesh)
        
        angle_threshold_rad = math.radians(angle_threshold_deg)
        
        visited = np.zeros(len(triangles), dtype=bool)
        patches = []
        
        for i in range(len(triangles)):
            if visited[i]:
                continue
            
            # Start a new patch
            patch = [i]
            visited[i] = True
            queue = [i]
            
            while queue:
                current = queue.pop()
                current_normal = triangle_normals[current]
                
                for neighbor in adjacency[current]:
                    if visited[neighbor]:
                        continue
                    
                    neighbor_normal = triangle_normals[neighbor]
                    
                    angle = np.arccos(
                        np.clip(np.dot(current_normal, neighbor_normal), -1.0, 1.0)
                    )
                    
                    if angle < angle_threshold_rad:
                        visited[neighbor] = True
                        patch.append(neighbor)
                        queue.append(neighbor)
            
            patches.append(patch)
        
        return patches  # List of lists of triangle indices

    def triangles_to_mesh(triangles):
        # Flatten triangle points into a list of points
        all_points = np.vstack(triangles)

        # Use numpy unique to remove duplicates and build index mapping
        unique_points, inverse_indices = np.unique(all_points.round(decimals=8), axis=0, return_inverse=True)

        # Rebuild triangles as indices
        triangle_indices = []
        for i in range(len(triangles)):
            idx = inverse_indices[i * 3 : (i + 1) * 3]
            triangle_indices.append(idx)

        # Create Open3D mesh
        mesh = o3d.geometry.TriangleMesh()
        mesh.vertices = o3d.utility.Vector3dVector(unique_points)
        mesh.triangles = o3d.utility.Vector3iVector(triangle_indices)

        return mesh

    def extract_submesh(mesh, triangle_indices):
        triangles = np.asarray(mesh.triangles)
        vertices = np.asarray(mesh.vertices)

        selected_triangles = triangles[triangle_indices]

        # Find the unique vertices used in the selected triangles
        unique_vertex_indices, new_indices = np.unique(selected_triangles.flatten(), return_inverse=True)

        # Map triangles to new indices
        new_triangles = new_indices.reshape(-1, 3)

        # Create new mesh
        new_vertices = vertices[unique_vertex_indices]

        submesh = o3d.geometry.TriangleMesh()
        submesh.vertices = o3d.utility.Vector3dVector(new_vertices)
        submesh.triangles = o3d.utility.Vector3iVector(new_triangles)
        submesh.compute_vertex_normals()
        submesh.compute_triangle_normals()

        return submesh

    def compute_patch_mesh_normal_center_model(mesh, patch_indices):
        # Model
        patch_mesh = extract_submesh(mesh, patch_indices)

        # Normal and Center 
        triangles = np.asarray(mesh.triangles)
        vertices = np.asarray(mesh.vertices)
        triangle_normals = np.asarray(mesh.triangle_normals)

        area_weighted_normal = np.zeros(3)
        patch_centroid = np.zeros(3)
        total_area = 0.0

        for idx in patch_indices:
            tri = triangles[idx]
            v0, v1, v2 = vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]

            # Triangle area
            area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
            total_area += area

            # Weighted normal
            area_weighted_normal += triangle_normals[idx] * area

            # Centroid contribution
            tri_centroid = (v0 + v1 + v2) / 3.0
            patch_centroid += tri_centroid * area

        if total_area > 0:
            patch_normal = area_weighted_normal / np.linalg.norm(area_weighted_normal)
            patch_center = patch_centroid / total_area
        else:
            patch_normal = np.array([0.0, 0.0, 1.0])  # Fallback normal
            patch_center = np.zeros(3)

        # Plane model: ax + by + cz + d = 0
        a, b, c = patch_normal
        d = -np.dot(patch_normal, patch_center)

        plane_model = (a, b, c, d)

        return patch_mesh, patch_normal, patch_center, plane_model

    def compute_patch_obox(mesh, thickness):
        # Get vertices
        vertices = np.asarray(mesh.vertices)

        # Compute centroid
        centroid = vertices.mean(axis=0)

        # Subtract mean
        centered = vertices - centroid

        # PCA to get normal
        cov = np.cov(centered.T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        normal = eigvecs[:, 0]  # Smallest eigenvector → plane normal

        # In-plane basis vectors
        v1 = eigvecs[:, 1]
        v2 = eigvecs[:, 2]

        # Project to 2D plane coordinate system
        projections = np.dot(centered, np.stack([v1, v2], axis=1))

        # 2D bounding box
        min_proj = projections.min(axis=0)
        max_proj = projections.max(axis=0)

        # Reconstruct corners in 3D and add thickness along normal
        corners = []
        for dx in [0, 1]:
            for dy in [0, 1]:
                for dz in [-thickness / 2, thickness / 2]:  # Add thickness
                    point2d = min_proj + (max_proj - min_proj) * [dx, dy]
                    corner3d = centroid + point2d[0] * v1 + point2d[1] * v2 + dz * normal
                    corners.append(corner3d)

        # Create OBB from corners
        obox = o3d.geometry.OrientedBoundingBox.create_from_points(o3d.utility.Vector3dVector(corners))
        return obox


    def compute_patch_area(mesh):
        triangles = np.asarray(mesh.triangles)
        vertices = np.asarray(mesh.vertices)

        total_area = 0.0

        for tri in triangles:
            v0, v1, v2 = vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]
            area = 0.5 * np.linalg.norm(np.cross(v1 - v0, v2 - v0))
            total_area += area

        return total_area

    # Extract IFC Patches
    stime_ifc_patches = time.time()

    if PCDTYPE == "TLS":
        ifc_patch_obox_thickness = 0.10
    else:
        ifc_patch_obox_thickness = 0.10


    ifc_patches = []
    ifc_patches_too_small = []

    ifc_patches_total = 0

    for mesh in all_meshes:
        patches = extract_planar_patches(mesh, angle_threshold_deg=normal_variance_threshold_deg)
        ifc_patches_total += len(patches)

        for patch in patches:
            patch_mesh, patch_normal, patch_center, patch_model = compute_patch_mesh_normal_center_model(mesh, patch)
            
            patch_area = compute_patch_area(patch_mesh)
            
            # Discard patches that are too small:
            if patch_area < patch_area_min:
                ifc_patches_too_small.append(patch_mesh)
                continue

            # Discard wrongly oriented patches (if known location):
            if KNOWN_TARGET_LOCATION:
                if np.dot(pcd_origin_target - patch_center, patch_normal) < 0.0:
                    ifc_patches_wrong_orientation.append(patch_mesh)
                    continue

            patch_obox = compute_patch_obox(patch_mesh, thickness=ifc_patch_obox_thickness)
            patch_obox.color = np.array([0.3, 0.3, 0.3], dtype=np.float32)

            ifc_patches.append({"mesh": patch_mesh, 
                                "center": patch_center, 
                                "normal": patch_normal, 
                                "model": patch_model, 
                                "obox": patch_obox,
                                "area": patch_area})

    # Sort all patches by area
    ifc_patches = sorted(ifc_patches, key=lambda p: p['area'], reverse=True)

    print(f"Found {ifc_patches_total} planar patches.")
    print(f"Found {len(ifc_patches)} valid planar patches.")
    if KNOWN_TARGET_LOCATION:
        print(f"Found {len(ifc_patches_wrong_orientation)} wrongly oriented planar patches (because known pcd origin target location).")
    print(f"Found {len(ifc_patches_too_small)} too small planar patches.")

    etime_ifc_patches = time.time()
    duration_ifc_patches = etime_ifc_patches - stime_ifc_patches
    print(f"Time to extract ifc patches: {duration_ifc_patches:.2f} s")



    # ## Extract 4PlCSs from the IFC patches.
    # This section extracts valid and useful 4PlCSs from the IFC patches using the same method as for the point cloud patches.

    # Extract IFC 4PlCSs
    stime_ifc_4PlCSs = time.time()

    ifc_4PlCSs = find_valid_patch_combinations(ifc_patches, coplanar_angle_thresh, coplanar_dist_thresh, plane_dist_thresh)

    # Select Largest IFC 4PlCSs
    ifc_patches_max = 50
    ifc_4PlCSs_largest = find_largest_4PlCSs(ifc_4PlCSs, ifc_patches, ifc_patches_max)
    ifc_4PlCSs_largest = sorted(ifc_4PlCSs_largest, key=lambda p: p['area'], reverse=True)

    match_4PlCSs_max = 100000
    ifc_4PlCSs_max = round(match_4PlCSs_max / pcd_4PlCSs_max)
    ifc_4PlCSs_filtered = ifc_4PlCSs_largest[:ifc_4PlCSs_max]
    print (f"There are {len(ifc_4PlCSs_filtered)} largest 4PlCSs in the ifc.")

    etime_ifc_4PlCSs = time.time()
    duration_ifc_4PlCSs = etime_ifc_4PlCSs - stime_ifc_4PlCSs
    print(f"Time to extract ifc 4PlCSs: {duration_ifc_4PlCSs:.2f} s")



    ###########################################################################
    # # 3. Match 4PlCSs

    # ## Find the best transformation by matching 4PlCSs.
    # For each pair of PCD and IFC 4PlCSs, this section:
    # 1. Checks if their internal geometries match.
    # 2. Calculates the possible rigid transformations between them.
    # 3. Computes the support for each transformation (how many PCD points align with IFC patches).
    # 4. Filters transformations based on minimum support and removes duplicates.

    # Plot all pcd and ifc patches together for visualization
    geom = []
    for ifc_patch in ifc_patches:
        ifc_patch_mesh = ifc_patch["obox"]
        geom.append(ifc_patch_mesh)

    for pcd_patch in pcd_patches:
        pcd_patch_pcd = pcd_patch["downpcd"]
        geom.append(pcd_patch_pcd)

    # o3d.visualization.draw_geometries(geom, mesh_show_back_face=True)

    def intersect_three_planes(plane1, plane2, plane3):
        """
        Computes the intersection point of three planes given in (a, b, c, d) form.
        """
        # Normals of the planes
        n1 = np.array(plane1[:3])
        n2 = np.array(plane2[:3])
        n3 = np.array(plane3[:3])
        
        # Coefficient matrix (normals stacked as rows)
        A = np.vstack([n1, n2, n3])
        
        # Right-hand side: -d terms
        d = np.array([-plane1[3], -plane2[3], -plane3[3]])
        
        # Solve the linear system
        try:
            point = np.linalg.solve(A, d)
        except np.linalg.LinAlgError as e:
            raise ValueError("Planes do not intersect at a single point (they may be parallel or coplanar).") from e
        
        return point

    def align_vectors(A, B):
        """
        Aligns vector A to vector B using a rigid transformation (rotation + translation)
        to minimize the MSE between corresponding points.
        """
        assert A.shape == B.shape, "Vectors must have the same shape"
        
        # Compute centroids
        centroid_A = A.mean(axis=0)
        centroid_B = B.mean(axis=0)

        # Center the vectors
        AA = A - centroid_A
        BB = B - centroid_B

        # Compute optimal rotation using SVD
        H = AA.T @ BB
        U, _, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T

        # Fix improper rotation (reflection)
        if np.linalg.det(R) < 0:
            Vt[2, :] *= -1
            R = Vt.T @ U.T

        # Compute translation
        t = centroid_B - R @ centroid_A

        # Apply transformation
        A_aligned = (R @ A.T).T + t

        # Compute mean squared error
        mse = np.mean(np.linalg.norm(A_aligned - B, axis=1)**2)

        return {"R": R, "t": t, "mes": mse}

    def transform_plane(plane, transformation):
        """
        Applies rotation R and translation t to a plane equation (a, b, c, d).
        """
        n = np.array(plane[:3])  # normal
        d = plane[3]

        R = transformation[:3, :3]
        t = transformation[:3, 3]

        # Rotate the normal
        n_rot = R @ n

        # Get a point on the original plane
        p = -d * n / np.dot(n, n)

        # Transform the point
        p_transformed = R @ p + t

        # Compute new d
        d_new = -np.dot(n_rot, p_transformed)

        return np.append(n_rot, d_new)

    def find_rotation_angles_around_line(p0, v, planeA1, planeA2, planeB1, planeB2, angle_tol=0.05):
        """
        Finds all rotation angles (in radians) around line (p0, v) that align planeA1,A2 with planeB1,B2.
        Returns a list of (angle_rad, rotation_matrix_4x4) tuples.
        """
        v = v / np.linalg.norm(v)
        nA1 = np.array(planeA1[:3])
        nA2 = np.array(planeA2[:3])
        nB1 = np.array(planeB1[:3])
        nB2 = np.array(planeB2[:3])
        
        # Remove axial components (only perpendicular rotation matters)
        def orth_proj(n): return n - np.dot(n, v) * v
        nA1_proj = orth_proj(nA1)
        nA2_proj = orth_proj(nA2)
        nB1_proj = orth_proj(nB1)
        nB2_proj = orth_proj(nB2)

        # Normalize projected normals
        nA1_proj /= np.linalg.norm(nA1_proj)
        nA2_proj /= np.linalg.norm(nA2_proj)
        nB1_proj /= np.linalg.norm(nB1_proj)
        nB2_proj /= np.linalg.norm(nB2_proj)

        def angle_between(u, v):
            cross = np.cross(u, v)
            dot = np.dot(u, v)
            angle = np.arctan2(np.linalg.norm(cross), dot)
            return angle

        def rotate_around_line(point, axis, angle_rad):
            """
            Returns a 4x4 transformation matrix rotating around a line defined by (point, axis).
            """
            axis = axis / np.linalg.norm(axis)
            rot = R.from_rotvec(angle_rad * axis).as_matrix()
            Transformation = np.eye(4)
            Transformation[:3, :3] = rot
            Transformation[:3, 3] = point - rot @ point
            return Transformation

        # Try two angle candidates: θ1 from A1→B1, θ2 from A1→-B1
        candidate_angles = []
        for sign1 in [+1, -1]:
            θ1 = angle_between(sign1 * nA1_proj, nB1_proj)
            R1 = R.from_rotvec(θ1 * v).as_matrix()
            nA2_rotated = R1 @ nA2_proj

            for sign2 in [+1, -1]:
                angle_diff = angle_between(nA2_rotated, sign2 * nB2_proj)
                if angle_diff < angle_tol:
                    # Valid alignment found
                    T1 = rotate_around_line(p0, v, θ1)
                    if not any(np.isclose(θ1, v["angle"], atol=0.02) for v in candidate_angles):
                        candidate_angles.append({"angle": θ1, "transformation": T1, "angle_diff": angle_diff})
                    θ2 = θ1 + np.pi
                    T2 = rotate_around_line(p0, v, θ2)
                    if not any(np.isclose(θ2, v["angle"], atol=0.02) for v in candidate_angles):
                        candidate_angles.append({"angle": θ2, "transformation": T2, "angle_diff": angle_diff})
        
        return candidate_angles

    def angle_sets_equal(set1: np.ndarray, set2: np.ndarray, tol:float=1.0) -> bool:
        unmatched_set2 = list(set2)
        for val_1 in set1:
            found_match = False
            for i, val_2 in enumerate(unmatched_set2):
                if math.isclose(val_1, val_2, abs_tol=tol):
                    unmatched_set2.pop(i)
                    found_match = True
                    break
            if not found_match:
                return False
        return True

    def match_geometry(pcd_geometry, ifc_geometry, diff_distance_max = 0.2, diff_angle_max = 4.0):
        # Match geometry (distance):
        diff_distance = abs(ifc_geometry[3] - pcd_geometry[3])

        if diff_distance > diff_distance_max:
            return False

        # Match geometry (angles):
        angles_equal = angle_sets_equal(pcd_geometry[:3], ifc_geometry[:3], diff_angle_max)
        if angles_equal == False:
            return False
        
        return True

    def compute_custom_obb_for_planar_mesh(mesh, tol):
        # Get vertices
        vertices = np.asarray(mesh.vertices)

        # Compute centroid
        centroid = vertices.mean(axis=0)

        # Subtract mean
        centered = vertices - centroid

        # PCA to get normal
        cov = np.cov(centered.T)
        eigvals, eigvecs = np.linalg.eigh(cov)
        normal = eigvecs[:, 0]  # Smallest eigenvector → plane normal

        # In-plane basis vectors
        v1 = eigvecs[:, 1]
        v2 = eigvecs[:, 2]

        # Project to 2D plane coordinate system
        projections = np.dot(centered, np.stack([v1, v2], axis=1))

        # 2D bounding box
        min_proj = projections.min(axis=0)
        max_proj = projections.max(axis=0)

        # Reconstruct corners in 3D and add thickness along normal
        corners = []
        for dx in [0, 1]:
            for dy in [0, 1]:
                for dz in [-tol / 2, tol / 2]:  # Add thickness
                    point2d = min_proj + (max_proj - min_proj) * [dx, dy]
                    corner3d = centroid + point2d[0] * v1 + point2d[1] * v2 + dz * normal
                    corners.append(corner3d)

        # Create OBB from corners
        obox = o3d.geometry.OrientedBoundingBox.create_from_points(o3d.utility.Vector3dVector(corners))
        return obox

    def count_support(pcd, mesh, max_dist):
        # Compute oriented bounding box (OBB)
        obox = compute_custom_obb_for_planar_mesh(mesh,max_dist)

        # Crop the point cloud using the OBB
        cropped_pcd = pcd.crop(obox)

        # Count points
        return len(cropped_pcd.points)

    def count_support(pcd, obox):
        # Crop the point cloud using the OBB
        cropped_pcd = pcd.crop(obox)

        # Count points
        return len(cropped_pcd.points)

    def point_in_triangle(p, a, b, c):
        # Using Barycentric coordinate check
        v0 = c - a
        v1 = b - a
        v2 = p - a

        dot00 = np.dot(v0, v0)
        dot01 = np.dot(v0, v1)
        dot02 = np.dot(v0, v2)
        dot11 = np.dot(v1, v1)
        dot12 = np.dot(v1, v2)

        denom = dot00 * dot11 - dot01 * dot01
        if denom == 0:
            return False  # Degenerate triangle

        u = (dot11 * dot02 - dot01 * dot12) / denom
        v = (dot00 * dot12 - dot01 * dot02) / denom

        return (u >= 0) and (v >= 0) and (u + v <= 1)

    def project_point_to_plane(p, a, n):
        # Orthogonal projection of point `p` onto the plane with point `a` and normal `n`
        return p - np.dot(p - a, n) * n

    def get_points_projecting_inside_mesh(pcd, mesh):
        points = np.asarray(pcd.points)
        vertices = np.asarray(mesh.vertices)
        triangles = np.asarray(mesh.triangles)

        selected_points = []

        for p in points:
            for tri in triangles:
                a, b, c = vertices[tri[0]], vertices[tri[1]], vertices[tri[2]]

                # Compute triangle normal
                normal = np.cross(b - a, c - a)
                norm = np.linalg.norm(normal)
                normal = normal / norm

                # Project point onto triangle plane
                proj = project_point_to_plane(p, a, normal)

                # Check if projected point lies within triangle
                if point_in_triangle(proj, a, b, c):
                    selected_points.append(p)
                    break  # No need to check other triangles
        
        if selected_points:
            return o3d.geometry.PointCloud(points=o3d.utility.Vector3dVector(np.array(selected_points)))
        else:
            return o3d.geometry.PointCloud()

    def find_best_transformations_for_patchset_pair(pcd_4PlCS, ifc_4PlCS, min_support_ratio, no_vertical_inversion=False):

        pcd_planes = []
        for i in pcd_4PlCS["patch_ids"]:
            pcd_patch_plane = pcd_patches[i]["model"]
            pcd_planes.append(pcd_patch_plane)

        ifc_planes = []
        for i in ifc_4PlCS["patch_ids"]:
            ifc_patch_plane = ifc_patches[i]["model"]
            ifc_planes.append(ifc_patch_plane)

        # Calculate pcd 4PlCS vectors between intersection points
        pcd_intersection1 = intersect_three_planes(pcd_planes[0], pcd_planes[1], pcd_planes[2])
        pcd_intersection2 = intersect_three_planes(pcd_planes[3], pcd_planes[1], pcd_planes[2])
        pcd_intersection_vectors = [np.array([pcd_intersection1, pcd_intersection2]), np.array([pcd_intersection2, pcd_intersection1])]

        # Calculate ifc 4PlCS vector between intersection points
        ifc_intersection1 = intersect_three_planes(ifc_planes[0], ifc_planes[1], ifc_planes[2])
        ifc_intersection2 = intersect_three_planes(ifc_planes[3], ifc_planes[1], ifc_planes[2])
        ifc_intersection_vector = np.array([ifc_intersection1, ifc_intersection2])


        # Find best transformations between pcd 4PlCS and ifc 4PlCS
        candidates = []

        for option_i, pcd_intersection_vector in enumerate(pcd_intersection_vectors):
            result = align_vectors(pcd_intersection_vector, ifc_intersection_vector)
            Transformation12 = np.eye(4)
            Transformation12[:3, :3] = result["R"]
            Transformation12[:3, 3] = result["t"]

            # Transform pcd_planes to pcd_planes2
            pcd_planes2 = [transform_plane(pcd_plane, Transformation12) for pcd_plane in pcd_planes]

            # Find rotations around ifc_intersection_vector
            # NOTE: angle_tol could be passed as a parameter
            point = ifc_intersection_vector[0]
            vect = ifc_intersection_vector[1] - ifc_intersection_vector[0]
            planeA1 = pcd_planes2[1]
            planeA2 = pcd_planes2[2]
            planeB1 = ifc_planes[1]
            planeB2 = ifc_planes[2]
            candidate_angles = find_rotation_angles_around_line(point, vect, planeA1, planeA2, planeB1, planeB2, angle_tol=0.05)
            
            if not candidate_angles:
                continue
            
            for candidate_angle in candidate_angles:
                Transformation23 = candidate_angle["transformation"]
            
                # Calculate combined transformation
                Transformation13 = Transformation23 @ Transformation12

                # Test if passes transformation constraint
                if no_vertical_inversion == True and Transformation13[2, 2] < -0.9:
                    continue

                # Test if passes known location constraint (if relevant):
                if KNOWN_TARGET_LOCATION == True:
                    v1_h = np.append(pcd_origin_init, 1.0)
                    v2_h = Transformation13 @ v1_h  # matrix multiplication
                    pcd_origin_transformed = v2_h[:3] / v2_h[3] # Convert back to 3D coordinates
                    pcd_origin_diff = pcd_origin_transformed - pcd_origin_target
                    #print(f"{v1_h}, {pcd_origin_transformed}, {pcd_origin_target}, {pcd_origin_diff}")
                    if np.linalg.norm(pcd_origin_diff) > 1.0:
                        continue
                    # NOTE: In addition, this should check if the matched model patches point towards the `pcd_origin_transformed`. 

                # Calculate support from patch points and ifc oboxes
                point_support_4PlCS = 0.0           
                for i in pcd_4PlCS["patch_ids"]:
                    pcd_patch_downpcd =  copy.deepcopy(pcd_patches[i]["downpcd"])
                    pcd_patch_downpcd.transform(Transformation13)

                    for j in ifc_4PlCS["patch_ids"]:
                        ifc_patch_obox = ifc_patches[j]["obox"]
                        
                        pcd_patch_downpcd_cropped = pcd_patch_downpcd.crop(ifc_patch_obox)
                        point_support_4PlCS += len(pcd_patch_downpcd_cropped.points)
                       
                # Calculate minimum support
                a4PlCS_point_total = 0.0  
                for i in pcd_4PlCS["patch_ids"]:
                    a4PlCS_point_total += len(pcd_patches[i]["downpcd"].points)
                point_support_4PlCS_min = int(round(a4PlCS_point_total * min_support_ratio))

                if point_support_4PlCS < point_support_4PlCS_min:
                    continue

                candidate = {"transformation12": Transformation12, 
                            "transformation23": Transformation23,
                            "transformation13": Transformation13, 
                            "support-4PlCS": point_support_4PlCS}
                candidates.append(candidate)

        sorted_candidates = sorted(candidates, key=lambda x: x["support-4PlCS"], reverse=True)

        return sorted_candidates


    def is_similar_matrix(T1, T2, tol=1e-5):
        return np.allclose(T1, T2, atol=tol)

    def remove_duplicate_transformations1(best_4PlCS_pairs, tol=1e-5):
        unique_transforms = []
        for a4PlCS_pair in best_4PlCS_pairs:
            T = a4PlCS_pair["transformation"]
            if not any(is_similar_transform(T, other["transformation"], tol) for other in unique_transforms):
                unique_transforms.append(a4PlCS_pair)
        return unique_transforms

    def remove_duplicate_transformations2(best_4PlCS_pairs, angle_tol=0.1, dist_tol=0.1):
        unique_transforms = []
        for a4PlCS_pair in best_4PlCS_pairs:
            T = a4PlCS_pair["transformation"]
            R = T[:3,:3]
            t = T[:3,3]
            found_similar = False
            for unique_transf in unique_transforms:
                Tu = unique_transf["transformation"]
                Ru = Tu[:3,:3]
                tu = Tu[:3,3]
                if is_similar_matrix(R, Ru, angle_tol) and is_similar_matrix(t, tu, dist_tol):
                    found_similar = True
            if found_similar == False:
                unique_transforms.append(a4PlCS_pair)
        return unique_transforms

    def rotation_matrix_to_euler_angles(R_mat, order='xyz', degrees=True):
        """
        Convert a 3x3 rotation matrix to Euler angles.

        Parameters:
            R_mat (numpy.ndarray): 3x3 rotation matrix
            order (str): Axis order for Euler angles (e.g., 'xyz', 'zyx')
            degrees (bool): If True, return angles in degrees; otherwise in radians

        Returns:
            numpy.ndarray: Array of 3 Euler angles
        """
        r = R.from_matrix(R_mat)
        return r.as_euler(order, degrees=degrees)

    # Find 4PlCS pairs with best support
    stime_match_4PlCSs = time.time()

    total = len(pcd_4PlCSs_filtered) * len(ifc_4PlCSs_filtered)
    print(f"Number of 4PlCS pairs tested: {len(pcd_4PlCSs_filtered)} x {len(ifc_4PlCSs_filtered)} = {total}")

    # Threshold for geometry similarity check
    if PCDTYPE == "TLS":
        diff_distance_max = 0.2
        diff_angle_max = 4.0
        min_support_ratio = 0.3
    else:
        diff_distance_max = 0.2
        diff_angle_max = 4.0
        min_support_ratio = 0.3


    # Threshold for point count support when matching pcd and ifc patches
    # min_support_downpcd_count = int(round(patch_sampling_max / 2))
    best_4PlCS_pairs = []

    for pcd_i, pcd_4PlCS in enumerate(pcd_4PlCSs_filtered):
        pcd_geometry = pcd_4PlCS["patch_geometry"]

        for ifc_i, ifc_4PlCS in enumerate(ifc_4PlCSs_filtered):
            ifc_geometry = ifc_4PlCS["patch_geometry"]

            # Check if pcd_geometry and ifc_geometry match
            geometry_match = match_geometry(pcd_geometry, ifc_geometry, diff_distance_max, diff_angle_max)

            if not geometry_match:
                continue

            # Find best transformations for the patch set pair
            top_transformations = find_best_transformations_for_patchset_pair(pcd_4PlCS, ifc_4PlCS, min_support_ratio, no_vertical_inversion=True)
            
            if not top_transformations:
                continue

            for transf in top_transformations:
                good_4PlCS_pair = {"pcd_4PlCS": pcd_4PlCS, 
                                   "ifc_4PlCS": ifc_4PlCS, 
                                   "transformation": transf["transformation13"], 
                                   "support-4PlCS": transf["support-4PlCS"]}
                
                best_4PlCS_pairs.append(good_4PlCS_pair)

    print (f"Number of 4PlCS pairs retained with min support ratio = {min_support_ratio}: {len(best_4PlCS_pairs)}")

    # Keep 500 transformations with largest support and remove duplicates
    max_trans_count = 500
    sorted_best_4PlCS_pairs = sorted(best_4PlCS_pairs, key=lambda c: c["support-4PlCS"], reverse=True)
    top_4PlCS_pairs = sorted_best_4PlCS_pairs[:max_trans_count]
    stime_duplicate = time.time()
    unique_top_4PlCS_pairs = remove_duplicate_transformations2(top_4PlCS_pairs, angle_tol=0.1, dist_tol=0.1)
    print (f"Number of top unique 4PlCS pairs: {len(unique_top_4PlCS_pairs)}")

    # Print Top Pairs
    # Calculate inverse of initial random transformation for comparison
    random_transform_inv = np.linalg.inv(random_transform)
    random_rotation_inv = random_transform_inv[:3, :3]
    random_euler_deg_inv = rotation_matrix_to_euler_angles(random_rotation_inv, order='xyz', degrees=True)
    random_translation_inv = random_transform_inv[:3, 3]
    print(random_transform_inv)

    print_count = 5
    print(f"Printing the top {print_count} transformations in terms of point support.")

    for a4PlCS_pair in unique_top_4PlCS_pairs[:print_count]:
        support_4PlCS = a4PlCS_pair["support-4PlCS"]
        print(support_4PlCS)

        transformation = a4PlCS_pair["transformation"]
        print(transformation)

        R_est = transformation[:3, :3]
        euler_deg = rotation_matrix_to_euler_angles(R_est, order='xyz', degrees=True)
        diff_euler = euler_deg - random_euler_deg_inv
        print(f"Euler angles (deg): {euler_deg} with expected {random_euler_deg_inv} and diff: {diff_euler}")
        
        t_est = transformation[:3, 3]
        diff_translation = t_est - random_translation_inv
        print(f"Translation: {t_est} with expected {random_translation_inv} diff: {diff_translation}")

    etime_match_4PlCSs = time.time()
    duration_match_4PlCSs = etime_match_4PlCSs - stime_match_4PlCSs
    print(f"Time to match 4PlCSs: {duration_match_4PlCSs:.2f} s")
    duration_total = duration_pcd_patches + duration_pcd_4PlCSs + duration_ifc_patches + duration_ifc_4PlCSs + duration_match_4PlCSs
    print(f"Time total: {duration_total:.2f} s")

    # ## Refine the order of the best transformations.
    # This section re-evaluates the top transformations by calculating support based on projecting PCD points onto the actual IFC patch meshes, rather than just their bounding boxes. This provides a more accurate ranking.

    def calculate_support_patches(pcd_patches_in, transformation_in, ifc_patches_in):
        total_count = 0
        
        for pcd_patch in pcd_patches_in:
            pcd_patch_pcd = copy.deepcopy(pcd_patch["downpcd"])
            pcd_patch_pcd.transform(transformation_in)

            for ifc_patch in ifc_patches_in:
                ifc_patch_obox = ifc_patch["obox"]
                cropped_pcd = pcd_patch_pcd.crop(ifc_patch_obox)
                total_count += len(cropped_pcd.points)
        
        return total_count

    # Refine selection based on support from all patches
    stime_order_4PlCSs = time.time()

    order_max_count = 200
    final_4PlCS_pairs = []

    for a4PlCS_pair in unique_top_4PlCS_pairs[:order_max_count]:
        transformation = a4PlCS_pair["transformation"]

        # total_count = 0
        # for pcd_patch in pcd_patches:
        #     pcd_patch_pcd = copy.deepcopy(pcd_patch["downpcd"])
        #     pcd_patch_pcd.transform(transformation)
        #     for ifc_patch in ifc_patches:
        #         ifc_patch_mesh = ifc_patch["mesh"]
        #         ifc_patch_obox = ifc_patch["obox"]
        #         cropped_pcd = pcd_patch_pcd.crop(ifc_patch_obox)
        #         if float(len(cropped_pcd.points)) / float(len(pcd_patch_pcd.points)) > 0.8:
        #             pcd_matching_triangles = get_points_projecting_inside_mesh(cropped_pcd, ifc_patch_mesh)
        #             total_count += len(pcd_matching_triangles.points)

        total_count = calculate_support_patches(pcd_patches, transformation, ifc_patches)
        
        new_4PlCS_pair = copy.deepcopy(a4PlCS_pair)
        new_4PlCS_pair["support-patches"] = total_count
        final_4PlCS_pairs.append(new_4PlCS_pair)

    final_4PlCS_pairs = sorted(final_4PlCS_pairs, key=lambda c: c["support-patches"], reverse=True)
    print (f"Number of 4PlCS pairs final: {len(final_4PlCS_pairs)}")

    # Print Final Pairs
    if APPLY_RANDOM_TRANSFORM:
        print(random_transform_inv)

    print_count = 10
    print(f"Printing the top {print_count} transformations in terms of point support.")

    for a4PlCS_pair in final_4PlCS_pairs[:print_count]:
        support = a4PlCS_pair["support-4PlCS"]
        print(f"4PlCS Support: {support}")
        support_mesh = a4PlCS_pair["support-patches"]
        print(support_mesh)

        transformation = a4PlCS_pair["transformation"]
        print(transformation)
        
        R_est = transformation[:3, :3]
        euler_deg = rotation_matrix_to_euler_angles(R_est, order='xyz', degrees=True)
        diff_euler = euler_deg - random_euler_deg_inv
        print(f"Euler angles (deg): {euler_deg} with expected {random_euler_deg_inv} and diff: {diff_euler}")

        t_est = transformation[:3, 3]
        diff_translation = t_est - random_translation_inv
        print(f"Translation: {t_est} with expected {random_translation_inv} diff: {diff_translation}")


    etime_order_4PlCSs = time.time()
    duration_order_4PlCSs = etime_order_4PlCSs - stime_order_4PlCSs
    print(f"Time to order top {order_count} 4PlCSs: {duration_order_4PlCSs:.2f} s")
    duration_total = duration_pcd_patches + duration_pcd_4PlCSs + duration_ifc_patches + duration_ifc_4PlCSs + duration_match_4PlCSs + duration_order_4PlCSs
    print(f"Time total: {duration_total:.2f} s")


    # Plot Final Top Pairs
    plot_count = min(5, len(final_4PlCS_pairs))

    if plot_count == 0:
        print("No suitable transformations found.")
        return

    selected_index = -1
    for i, a4PlCS_pair in enumerate(final_4PlCS_pairs[:plot_count]):
        print(f"Displaying candidate {i+1}")
        geom = []

        transformation = a4PlCS_pair["transformation"]
        pcd_viz = copy.deepcopy(downpcd)
        pcd_viz.transform(transformation)
        geom.append(pcd_viz)

        #for pcd_patch in pcd_patches[:30]:
        #    pcd_patch_pcd = copy.deepcopy(pcd_patch["downpcd"])
        #    pcd_patch_pcd.transform(transformation)
        #    geom.append(pcd_patch_pcd)
        
        geom.extend(all_meshes)

        o3d.visualization.draw_geometries(geom, window_name=f"Candidate {i+1}", mesh_show_back_face=True)

        if selected_index == -1 and input(f"Accept this alignment? (y/n): ").lower() == 'y':
            selected_index = i
            break

    # # 4. Final ICP on selected best match
    # The user selects the preferred transformation from the refined list, and ICP is applied to fine-tune the alignment.
    if selected_index == -1:
        print("No candidate selected. Displaying all candidates again.")
        
        while selected_index < 0 or selected_index >= plot_count:
            try:
                user_input = input(f"Enter the number of the best alignment (1 to {plot_count}): ")
                selected_index = int(user_input) - 1
                if selected_index < 0 or selected_index >= plot_count:
                    print("Invalid selection. Please try again.")
            except ValueError:
                print("Invalid input. Please enter a number.")

    # Create combined mesh from all IFC elements
    combined_mesh = o3d.geometry.TriangleMesh()
    for mesh in all_meshes:
        combined_mesh += mesh

    combined_mesh.compute_vertex_normals()
    mesh_sampling_count = 100000
    mesh_pcd = combined_mesh.sample_points_uniformly(number_of_points=mesh_sampling_count)

    # ICP registration
    init_transform = final_4PlCS_pairs[selected_index]["transformation"]
    result = o3d.pipelines.registration.registration_icp(
        downpcd,
        mesh_pcd,
        max_correspondence_distance=0.05,
        init=init_transform,
        estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPlane()
    )

    # Print final transformation result
    final_transformation = copy.deepcopy(result.transformation)
    print(final_transformation)

    R_est = final_transformation[:3, :3]
    euler_deg = rotation_matrix_to_euler_angles(R_est, order='xyz', degrees=True)
    diff_euler = euler_deg - random_euler_deg_inv
    print(f"Euler angles (deg): {euler_deg} with expected {random_euler_deg_inv} and diff: {diff_euler}")

    t_est = final_transformation[:3, 3]
    diff_translation = t_est - random_translation_inv
    print(f"Translation: {t_est} with expected {random_translation_inv} diff: {diff_translation}")

    # Save transformed point cloud
    output_dir = os.path.join("output")
    os.makedirs(output_dir, exist_ok=True)
    pcd_transformed = copy.deepcopy(downpcd_before_transform)
    pcd_transformed.transform(final_transformation)
    output_path = os.path.join(output_dir, f"{pcd_filename}_aligned.ply")
    o3d.io.write_point_cloud(output_path, pcd_transformed)
    print(f"Transformed point cloud saved to {output_path}")

    # Plot the final alignment
    print("Displaying the final alignment of the full point cloud and the full IFC mesh.")
    o3d.visualization.draw_geometries([pcd_transformed, combined_mesh],
                                      window_name="Final Alignment: Point Cloud and IFC Mesh", mesh_show_back_face=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="4PlCS registration of a point cloud to an IFC model.")
    parser.add_argument("ifc_path", nargs='?', default=None, help="Path to the IFC file.")
    parser.add_argument("pcd_path", nargs='?', default=None, help="Path to the point cloud file.")
    args = parser.parse_args()

    ifc_path = args.ifc_path
    pcd_path = args.pcd_path

    if not ifc_path or not pcd_path:
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()

        if not ifc_path:
            ifc_path = filedialog.askopenfilename(title="Select IFC file", filetypes=[("IFC files", "*.ifc")])
        
        if not pcd_path:
            pcd_path = filedialog.askopenfilename(title="Select PCD file", filetypes=[("Point Cloud Data", "*.pcd *.ply")])

    if ifc_path and pcd_path:
        run_4plcs(ifc_path, pcd_path)
    else:
        print("IFC or PCD file not selected. Exiting.")
