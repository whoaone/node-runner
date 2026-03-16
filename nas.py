# nas.py

import math
import datetime
import numpy as np
from io import StringIO
from pyNastran.bdf.bdf import BDF
from scipy.spatial import cKDTree
from pyNastran.bdf.cards.elements.shell import CQUAD4, CTRIA3

class NastranModelGenerator:
    def __init__(self, params=None):
        self.params = params if params else {}
        self.model = BDF(debug=False)


    def parse_from_text(self, text):
        bdf_string_io = StringIO(text)
        try:
            self.model.read_bdf(bdf_string_io, punch=False)
    
        except Exception as e:
            raise RuntimeError(f"pyNastran failed to parse BDF: {e}")

    def _write_bdf(self, output_path):
        self.model.write_bdf(output_path, size=8, is_double=False)
    
    def _generate_nodes(self):
        radius = self.params['fuselage_radius']; num_bays = self.params['num_bays']; frame_spacing = self.params['frame_spacing']; num_stringers = self.params['num_stringers']
        num_frames = num_bays + 1; angle_increment = 360.0 / num_stringers
        node_id = 1
        for i in range(num_frames):
            x_coord = i * frame_spacing
            for j in range(num_stringers):
                angle_rad = math.radians(j * angle_increment)
                y_coord, z_coord = radius * math.cos(angle_rad), radius * math.sin(angle_rad)
                self.model.add_grid(node_id, [x_coord, y_coord, z_coord])
                node_id += 1

    def _generate_elements(self):
        num_stringers = self.params['num_stringers']; num_frames = self.params['num_bays'] + 1
        pid_skin, pid_frame, pid_stringer = 1, 2, 3
        mid_skin, mid_frame, mid_stringer = self.params.get('skin_material_id', 1), self.params.get('frame_material_id', 1), self.params.get('stringer_material_id', 1)
        skin_elem_type, frame_elem_type, stringer_elem_type = self.params.get('skin_element_type', 'CQUAD4'), self.params.get('frame_element_type', 'CBEAM'), self.params.get('stringer_element_type', 'CBEAM')
        
        pshell = self.model.add_pshell(pid_skin, t=self.params.get('skin_thickness', 0.01), comment='Skin'); pshell.mid = mid_skin
        frame_props = self._calculate_beam_props(self.params.get('frame_section_type', 'L'), self.params.get('frame_dims', {}))
        stringer_props = self._calculate_beam_props(self.params.get('stringer_section_type', 'L'), self.params.get('stringer_dims', {}))

        if frame_elem_type == 'CROD': self.model.add_prod(pid_frame, mid_frame, A=frame_props['A'], comment='Frames')
        else: self.model.add_pbeam(pid_frame, mid_frame, xxb=[0.0], so=['C'], area=[frame_props['A']], i1=[frame_props['I1']], i2=[frame_props['I2']], i12=[0.0], j=[frame_props['J']], comment='Frames')
        if stringer_elem_type == 'CROD': self.model.add_prod(pid_stringer, mid_stringer, A=stringer_props['A'], comment='Stringers')
        else: self.model.add_pbeam(pid_stringer, mid_stringer, xxb=[0.0], so=['C'], area=[stringer_props['A']], i1=[stringer_props['I1']], i2=[stringer_props['I2']], i12=[0.0], j=[stringer_props['J']], comment='Stringers')

        elem_id = 1
        for i in range(num_frames - 1):
            for j in range(num_stringers):
                n1, n2 = (i * num_stringers) + j + 1, (i * num_stringers) + 1 if (j + 1) >= num_stringers else (i * num_stringers) + j + 2
                n3, n4 = n2 + num_stringers, n1 + num_stringers
                if skin_elem_type == 'CQUAD4': self.model.add_cquad4(elem_id, pid_skin, [n1, n2, n3, n4])
                elif skin_elem_type == 'CMEMBRAN': self.model.add_card(['CMEMBRAN', elem_id, pid_skin, n1, n2, n3, n4], 'CMEMBRAN')
                elem_id += 1
        for i in range(num_frames):
            for j in range(num_stringers):
                n1, n2 = (i * num_stringers) + j + 1, (i * num_stringers) + 1 if (j + 1) >= num_stringers else (i * num_stringers) + j + 2
                if frame_elem_type == 'CBEAM': self.model.add_cbeam(elem_id, pid_frame, [n1, n2], [0., 0., 1.], g0=None)
                elif frame_elem_type == 'CBAR': self.model.add_cbar(elem_id, pid_frame, [n1, n2], [0., 0., 1.], g0=None)
                elif frame_elem_type == 'CROD': self.model.add_crod(elem_id, pid_frame, [n1, n2])
                elem_id += 1
        for i in range(num_frames - 1):
            for j in range(num_stringers):
                n1, n2 = (i * num_stringers) + j + 1, (i * num_stringers) + j + 1 + num_stringers
                if stringer_elem_type == 'CBEAM': self.model.add_cbeam(elem_id, pid_stringer, [n1, n2], [0., 0., 1.], g0=None)
                elif stringer_elem_type == 'CBAR': self.model.add_cbar(elem_id, pid_stringer, [n1, n2], [0., 0., 1.], g0=None)
                elif stringer_elem_type == 'CROD': self.model.add_crod(elem_id, pid_stringer, [n1, n2])
                elem_id += 1
    
    def generate_floor_structure(self):
        if not self.params.get("add_floor", False): return
        
        method = self.params.get("connection_method")
        if method == 'modify':
            self._create_floor_by_modifying()
        elif method == 'snap':
            self._create_floor_by_snapping()
            
    def _create_floor_by_snapping(self):
        print("--- Running Snap to Node ---")
        radius, frame_spacing = self.params['fuselage_radius'], self.params['frame_spacing']
        num_frames, z_floor = self.params['num_bays'] + 1, self.params['floor_z']
        force_vertical = self.params.get('force_vertical', False)

        pid_floor_beam, pid_stanchion = 4, 5 
        self.model.add_pbeam(pid_floor_beam, self.params.get('frame_material_id', 1), xxb=[0.0], so=['C'], area=[1.0], i1=[1.0], i2=[1.0], i12=[0.0], j=[1.0], comment='FloorBeams')
        self.model.add_pbeam(pid_stanchion, self.params.get('frame_material_id', 1), xxb=[0.0], so=['C'], area=[0.5], i1=[0.5], i2=[0.5], i12=[0.0], j=[0.5], comment='Stanchions')

        for i in range(num_frames):
            x_frame = i * frame_spacing
            frame_nodes = {nid: node.xyz for nid, node in self.model.nodes.items() if np.isclose(node.xyz[0], x_frame)}
            if not frame_nodes: continue
            if abs(z_floor) >= radius: continue
            
            y_end = math.sqrt(radius**2 - z_floor**2)
            left_nid = self._find_closest_node(np.array([x_frame, -y_end, z_floor]), frame_nodes)
            right_nid = self._find_closest_node(np.array([x_frame, y_end, z_floor]), frame_nodes)
            if not left_nid or not right_nid: continue

            avg_z = (self.model.nodes[left_nid].xyz[2] + self.model.nodes[right_nid].xyz[2]) / 2.0
            
            num_beam_elems = self.params['floor_elem_count']
            interior_beam_nids = self.add_nodes_between(left_nid, right_nid, num_beam_elems - 1)
            all_beam_nids = [left_nid] + interior_beam_nids + [right_nid]
            
            for nid in all_beam_nids:
                self.model.nodes[nid].xyz[2] = avg_z

            floor_beam_elements_to_create = []
            for j in range(len(all_beam_nids) - 1):
                floor_beam_elements_to_create.append({'n1': all_beam_nids[j], 'n2': all_beam_nids[j+1]})

            all_beam_nodes = {nid: self.model.nodes[nid].xyz for nid in all_beam_nids}
            for y_stanchion in self.params.get('stanchion_y_coords', []):
                if abs(y_stanchion) >= radius: continue
                
                z_skin = -math.sqrt(radius**2 - y_stanchion**2)
                skin_nid = self._find_closest_node(np.array([x_frame, y_stanchion, z_skin]), frame_nodes)
                
                beam_nid = None
                if not force_vertical:
                    beam_nid = self._find_closest_node(np.array([x_frame, y_stanchion, avg_z]), all_beam_nodes)
                else:
                    y_vertical = self.model.nodes[skin_nid].xyz[1]
                    target_beam_coord = np.array([x_frame, y_vertical, avg_z])
                    
                    parent_beam_elem = None
                    for beam_elem in floor_beam_elements_to_create:
                        y1 = self.model.nodes[beam_elem['n1']].xyz[1]
                        y2 = self.model.nodes[beam_elem['n2']].xyz[1]
                        if min(y1, y2) - 1e-9 <= y_vertical <= max(y1, y2) + 1e-9:
                            parent_beam_elem = beam_elem
                            break
                    
                    if parent_beam_elem:
                        beam_nid = self.add_node(None, *target_beam_coord)
                        all_beam_nodes[beam_nid] = target_beam_coord
                        floor_beam_elements_to_create.remove(parent_beam_elem)
                        floor_beam_elements_to_create.append({'n1': parent_beam_elem['n1'], 'n2': beam_nid})
                        floor_beam_elements_to_create.append({'n1': beam_nid, 'n2': parent_beam_elem['n2']})

                if beam_nid is None or beam_nid == skin_nid: continue

                num_stanchion_elems = self.params['stanchion_elem_count']
                interior_stanchion_nids = self.add_nodes_between(beam_nid, skin_nid, num_stanchion_elems - 1)
                all_stanchion_nids = [beam_nid] + interior_stanchion_nids + [skin_nid]
                for k in range(len(all_stanchion_nids) - 1):
                    self.add_line_element(None, pid_stanchion, all_stanchion_nids[k], all_stanchion_nids[k+1], 'CBEAM', {'method': 'vector', 'values': [1., 0., 0.]})
            
            for elem_conn in floor_beam_elements_to_create:
                self.add_line_element(None, pid_floor_beam, elem_conn['n1'], elem_conn['n2'], 'CBEAM', {'method': 'vector', 'values': [0., 0., 1.]})

    def _create_floor_by_modifying(self):
        print("--- Running Modify Fuselage Mesh (pyNastran) ---")
        radius, frame_spacing, num_frames, z_floor = self.params['fuselage_radius'], self.params['frame_spacing'], self.params['num_bays'] + 1, self.params['floor_z']
        
        pid_floor_beam, pid_stanchion = 4, 5
        self.model.add_pbeam(pid_floor_beam, self.params.get('frame_material_id', 1), xxb=[0.0], so=['C'], area=[1.0], i1=[1.0], i2=[1.0], i12=[0.0], j=[1.0], comment='FloorBeams')
        self.model.add_pbeam(pid_stanchion, self.params.get('frame_material_id', 1), xxb=[0.0], so=['C'], area=[0.5], i1=[0.5], i2=[0.5], i12=[0.0], j=[0.5], comment='Stanchions')

        print("--- PHASE 1: Creating new nodes and splitting frames... ---")
        new_nodes_by_frame, frame_mods, split_frame_map = {i: {} for i in range(num_frames)}, {}, {}
        for i in range(num_frames):
            x_frame = i * frame_spacing
            frame_nodes = {nid: node.xyz for nid, node in self.model.nodes.items() if np.isclose(node.xyz[0], x_frame)}
            if not frame_nodes: continue
            current_frame_elems = [elem for elem in self.model.elements.values() if elem.type in ['CBEAM', 'CBAR'] and elem.nodes[0] in frame_nodes and elem.nodes[1] in frame_nodes]

            targets = []
            if abs(z_floor) < radius:
                y_end = math.sqrt(radius**2 - z_floor**2)
                targets.extend([{'label': 'floor_beam_left', 'coord': (-y_end, z_floor)}, {'label': 'floor_beam_right', 'coord': (y_end, z_floor)}])
            for idx, y_stanchion in enumerate(self.params.get('stanchion_y_coords', [])):
                if abs(y_stanchion) < radius:
                    z_skin = -math.sqrt(radius**2 - y_stanchion**2)
                    targets.append({'label': f'stanchion_{idx}_skin', 'coord': (y_stanchion, z_skin)})
            
            for target in sorted(targets, key=lambda t: math.atan2(t['coord'][1], t['coord'][0])):
                y_target, z_target = target['coord']
                parent_elem = self._find_parent_frame_element_pynastran((y_target, z_target), current_frame_elems)
                if parent_elem:
                    new_nid = self.add_node(None, x_frame, y_target, z_target)
                    new_nodes_by_frame[i][target['label']] = new_nid
                    if parent_elem.eid not in frame_mods: frame_mods[parent_elem.eid] = {'elem': parent_elem, 'new_nodes': []}
                    frame_mods[parent_elem.eid]['new_nodes'].append(new_nid)

        for eid, mod_info in frame_mods.items():
            elem, n1_coord = mod_info['elem'], self.model.nodes[mod_info['elem'].nodes[0]].xyz
            sorted_new_nodes = sorted(mod_info['new_nodes'], key=lambda nid: np.linalg.norm(self.model.nodes[nid].xyz - n1_coord))
            all_split_nodes = [elem.nodes[0]] + sorted_new_nodes + [elem.nodes[1]]
            split_frame_map[tuple(sorted(elem.nodes))] = all_split_nodes
            print(f"  - Splitting Frame EID {eid} ({elem.nodes[0]}-{elem.nodes[1]}) with nodes {sorted_new_nodes}")
            del self.model.elements[eid]
            for i in range(len(all_split_nodes) - 1):
                self.add_line_element(None, elem.pid, all_split_nodes[i], all_split_nodes[i+1], elem.type, {'method':'vector', 'values':[0.,0.,1.]})
        print("--- PHASE 1: Complete. ---")

        print("\n--- PHASE 2: Re-meshing affected skin panels... ---")
        all_shell_eids = [eid for eid, elem in self.model.elements.items() if elem.type in ['CQUAD4', 'CTRIA3']]
        affected_count = 0
        for eid in all_shell_eids:
            shell = self.model.elements.get(eid)
            if not shell: continue 
            pid, n_orig = shell.pid, shell.nodes
            edge1_key = tuple(sorted((n_orig[0], n_orig[1])))
            edge2_key = tuple(sorted((n_orig[3], n_orig[2])))
            if edge1_key in split_frame_map or edge2_key in split_frame_map:
                affected_count += 1
                s1 = set(split_frame_map.get(edge1_key, [])) - set(edge1_key)
                s2 = set(split_frame_map.get(edge2_key, [])) - set(edge2_key)
                int_node1, int_node2 = s1.pop() if s1 else None, s2.pop() if s2 else None
                print(f"  - Re-meshing affected Shell EID {shell.eid} ({shell.nodes_ref}) -> New Nodes: {int_node1}, {int_node2}")
                del self.model.elements[shell.eid]
                if int_node1 and int_node2:
                    self.add_plate_element(None, pid, [n_orig[0], int_node1, int_node2, n_orig[3]], 'CQUAD4')
                    self.add_plate_element(None, pid, [int_node1, n_orig[1], n_orig[2], int_node2], 'CQUAD4')
                elif int_node1:
                    self.add_plate_element(None, pid, [n_orig[0], int_node1, n_orig[3]], 'CTRIA3')
                    self.add_plate_element(None, pid, [int_node1, n_orig[1], n_orig[2]], 'CTRIA3')
                    self.add_plate_element(None, pid, [int_node1, n_orig[2], n_orig[3]], 'CTRIA3')
                elif int_node2:
                    self.add_plate_element(None, pid, [n_orig[0], n_orig[1], int_node2], 'CTRIA3')
                    self.add_plate_element(None, pid, [n_orig[0], int_node2, n_orig[3]], 'CTRIA3')
                    self.add_plate_element(None, pid, [n_orig[1], n_orig[2], int_node2], 'CTRIA3')
        print(f"  - Found and re-meshed {affected_count} shells.")
        print("--- PHASE 2: Complete. ---")

        print("\n--- PHASE 3: Creating floor beams and stanchions... ---")
        for i in range(num_frames):
            left_nid, right_nid = new_nodes_by_frame.get(i, {}).get('floor_beam_left'), new_nodes_by_frame.get(i, {}).get('floor_beam_right')
            if left_nid and right_nid:
                all_beam_nids = [left_nid] + self.add_nodes_between(left_nid, right_nid, self.params['floor_elem_count'] - 1) + [right_nid]
                for j in range(len(all_beam_nids) - 1): self.add_line_element(None, pid_floor_beam, all_beam_nids[j], all_beam_nids[j+1], 'CBEAM', {'method':'vector', 'values':[0.,0.,1.]})
                
                beam_nodes = {nid: self.model.nodes[nid].xyz for nid in all_beam_nids}
                for idx, y_stanchion in enumerate(self.params.get('stanchion_y_coords', [])):
                    ideal_node = np.array([i * frame_spacing, y_stanchion, z_floor])
                    beam_nid = self._find_closest_node(ideal_node, beam_nodes)
                    skin_nid = new_nodes_by_frame.get(i, {}).get(f'stanchion_{idx}_skin')
                    if beam_nid and skin_nid and beam_nid != skin_nid:
                        all_stanchion_nids = [beam_nid] + self.add_nodes_between(beam_nid, skin_nid, self.params['stanchion_elem_count'] - 1) + [skin_nid]
                        for k in range(len(all_stanchion_nids) - 1): self.add_line_element(None, pid_stanchion, all_stanchion_nids[k], all_stanchion_nids[k+1], 'CBEAM', {'method':'vector', 'values':[1.,0.,0.]})
        print("--- PHASE 3: Complete. ---")
        
    def _find_parent_frame_element_pynastran(self, target_yz, frame_elements):
        target_angle = (math.atan2(target_yz[1], target_yz[0]) + 2 * math.pi) % (2 * math.pi)
        for elem in frame_elements:
            p1_xyz, p2_xyz = self.model.nodes[elem.nodes[0]].xyz, self.model.nodes[elem.nodes[1]].xyz
            angle1 = (math.atan2(p1_xyz[2], p1_xyz[1]) + 2 * math.pi) % (2 * math.pi)
            angle2 = (math.atan2(p2_xyz[2], p2_xyz[1]) + 2 * math.pi) % (2 * math.pi)
            if abs(angle1 - angle2) > math.pi:
                if target_angle >= max(angle1, angle2) or target_angle <= min(angle1, angle2): return elem
            elif min(angle1, angle2) <= target_angle <= max(angle1, angle2): return elem
        return None
        
    def _find_closest_node(self, target_coord, node_dict):
        if not node_dict: return None
        node_ids = list(node_dict.keys())
        node_coords = np.array(list(node_dict.values()))
        distances = np.linalg.norm(node_coords - target_coord, axis=1)
        closest_index = np.argmin(distances)
        return node_ids[closest_index]


    def _calculate_beam_props(self, section_type, dims):
        """Calculates beam section properties based on type and dimensions."""
        props = {'A': 0., 'I1': 0., 'I2': 0., 'J': 0.}
        try:
            stype = section_type.upper()
            if stype == 'L':
                props.update(self._calc_l_section_props(dims))
            elif stype == 'T':
                props.update(self._calc_t_section_props(dims))
            elif stype == 'Z':
                props.update(self._calc_z_section_props(dims))
            elif stype == 'J':
                props.update(self._calc_j_section_props(dims))
        except (KeyError, ZeroDivisionError):
            # Return default properties if dims are missing or invalid
            pass
        return props

    def _calc_l_section_props(self, dims):
        """Calculates properties for an L-section."""
        H, W, t = dims['H'], dims['W'], dims['t']
        A = t * (H + W - t)
        # Centroid calculation
        yc = (W*t*(W/2) + (H-t)*t*(t/2)) / A
        zc = (W*t*(t/2) + (H-t)*t*((H+t)/2)) / A
        # Moments of Inertia
        I1 = (W*t**3)/12 + W*t*(zc-t/2)**2 + (t*(H-t)**3)/12 + t*(H-t)*(zc-((H+t)/2))**2
        I2 = (t*W**3)/12 + W*t*(yc-W/2)**2 + ((H-t)*t**3)/12 + (H-t)*t*(yc-t/2)**2
        # Torsional constant
        J = (1/3) * (H*t**3 + (W-t)*t**3)
        return {'A': A, 'I1': I1, 'I2': I2, 'J': J}
    
    def _calc_t_section_props(self, dims):
        """Calculates properties for a T-section."""
        W, tf, H, tw = dims['W'], dims['tf'], dims['H'], dims['tw']
        A = W*tf + H*tw
        # Centroid calculation
        zc = ((W*tf)*(tf/2) + (H*tw)*(tf + H/2)) / A
        # Moments of Inertia
        I1 = (W*tf**3)/12 + W*tf*(zc-tf/2)**2 + (tw*H**3)/12 + H*tw*(tf+H/2-zc)**2
        I2 = (tf*W**3)/12 + (tw*H**3)/12 # Note: Original formula seems simplified; this matches original code.
        # Torsional constant
        J = (1/3)*(W*tf**3 + H*tw**3)
        return {'A':A,'I1':I1,'I2':I2,'J':J}

    def _calc_z_section_props(self, dims):
        """Calculates properties for a Z-section."""
        H, W, t = dims['H'], dims['W'], dims['t']
        A = H*t + 2*W*t
        I1 = (t*H**3)/12 + 2*((W*t**3)/12 + W*t*(H/2)**2)
        I2 = (H*t**3)/12 + 2*((t*W**3)/12 + W*t*(W/2)**2)
        J = (1/3) * (H*t**3 + 2*W*t**3)
        return {'A': A, 'I1': I1, 'I2': I2, 'J': J}

    def _calc_j_section_props(self, dims):
        """Calculates properties for a J-section."""
        H,W1,W2,t = dims['H'],dims['W1'],dims['W2'],dims['t']
        A = H*t + W1*t + W2*t
        # Centroid calculation
        yc=((H*t)*(t/2)+(W1*t)*(t+W1/2)+(W2*t)*(t+W2/2))/A
        zc=((H*t)*(H/2)+(W1*t)*(H-t/2)+(W2*t)*(t/2))/A
        # Moments of Inertia
        I1_web=(t*H**3)/12+H*t*(zc-H/2)**2
        I1_fl1=(W1*t**3)/12+W1*t*(zc-(H-t/2))**2
        I1_fl2=(W2*t**3)/12+W2*t*(zc-t/2)**2
        I2_web=(H*t**3)/12+H*t*(yc-t/2)**2
        I2_fl1=(t*W1**3)/12+W1*t*(yc-(t+W1/2))**2
        I2_fl2=(t*W2**3)/12+W2*t*(yc-(t+W2/2))**2
        # Torsional constant
        J=(1/3)*(H*t**3+W1*t**3+W2*t**3)
        return {'A':A,'I1':I1_web+I1_fl1+I1_fl2,'I2':I2_web+I2_fl1+I2_fl2,'J':J}
    

    # --- NEW: Method to find duplicate nodes (Phase 1) ---
    def find_duplicate_nodes(self, tolerance):
        """
        Finds groups of duplicate nodes within a given tolerance using a k-d tree.

        Args:
            tolerance (float): The distance within which nodes are considered duplicates.

        Returns:
            list[list[int]]: A list of groups, where each inner list contains the IDs
                             of duplicate nodes. e.g., [[101, 201], [105, 304, 501]]
        """
        if not self.model.nodes or tolerance <= 0:
            return []

        # Extract node IDs and their corresponding coordinates into numpy arrays
        node_ids = np.array(list(self.model.nodes.keys()))
        coords = np.array([node.get_position() for node in self.model.nodes.values()])

        # Build the k-d tree for efficient spatial searching
        tree = cKDTree(coords)

        # Find all pairs of nodes within the given tolerance radius.
        # This returns a list where each element is a list of indices of neighboring points.
        pairs = tree.query_ball_tree(tree, r=tolerance)

        # Process the raw pairs into unique, sorted groups of duplicate node IDs.
        # A set is used to automatically handle uniqueness.
        duplicate_sets = set()
        for i, neighbors_indices in enumerate(pairs):
            # A group must have more than one node to be considered a duplicate set.
            if len(neighbors_indices) > 1:
                # Convert indices back to actual node IDs and store as a sorted tuple.
                group_ids = tuple(sorted([node_ids[j] for j in neighbors_indices]))
                duplicate_sets.add(group_ids)

        # Convert the set of tuples back to a list of lists for the final output.
        return [list(group) for group in sorted(list(duplicate_sets))]
    

    # --- NEW: Method to merge duplicate nodes (Phase 2) ---
    def merge_duplicate_nodes(self, tolerance):
        """
        Finds and merges duplicate nodes, updating all element and entity references.

        Args:
            tolerance (float): The distance within which nodes are considered duplicates.

        Returns:
            dict: A summary of the operation, e.g.,
                  {'groups_found': 5, 'nodes_merged': 10}
        """
        # Step 1: Find the groups of duplicate nodes using our Phase 1 function.
        duplicate_groups = self.find_duplicate_nodes(tolerance)
        if not duplicate_groups:
            return {'groups_found': 0, 'nodes_merged': 0}

        # Step 2: For each group, determine the "master" node (the one with the lowest ID)
        # and create a mapping from the other nodes in the group to the master.
        remap_dict = {}
        nodes_to_delete = set()
        for group in duplicate_groups:
            master_nid = min(group)
            for nid in group:
                if nid != master_nid:
                    remap_dict[nid] = master_nid
                    nodes_to_delete.add(nid)

        # Step 3: Iterate through all entities in the model and update their node references.
        # Combine standard elements and rigid elements for updating.
        all_elements = list(self.model.elements.values()) + \
                       list(self.model.rigid_elements.values())

        for element in all_elements:
            # Re-create the element's node list, replacing any old IDs with master IDs.
            element.nodes = [remap_dict.get(nid, nid) for nid in element.nodes]

        # Update node references in constraints (SPCs).
        for sid, spc_cards in self.model.spcs.items():
            for card in spc_cards:
                new_node_list = [remap_dict.get(nid, nid) for nid in card.nodes]
                # Use set() to remove potential duplicate node IDs in the list after remapping.
                card.nodes = list(set(new_node_list))

        # Update node references in loads (FORCE, MOMENT, etc.).
        for sid, load_cards in self.model.loads.items():
            for card in load_cards:
                # These cards have a single node ID attribute.
                if hasattr(card, 'node_id') and card.node_id in remap_dict:
                    card.node_id = remap_dict[card.node_id]
                # Accommodate for legacy or different card types that might use '.node'.
                elif hasattr(card, 'node') and card.node in remap_dict:
                    card.node = remap_dict[card.node]

        # Step 4: After all references are updated, delete the old nodes from the model.
        for nid in nodes_to_delete:
            if nid in self.model.nodes:
                del self.model.nodes[nid]

        # Step 5: Return a summary of the operation.
        return {
            'groups_found': len(duplicate_groups),
            'nodes_merged': len(nodes_to_delete)
        }


     # --- NEW: Method to merge a specific list of node groups (Revised Phase 4) ---
    def merge_node_groups(self, groups_to_merge):
        """
        Merges a specific list of node groups, updating all entity references.

        Args:
            groups_to_merge (list[list[int]]): A list of groups to merge,
                e.g., [[101, 201], [105, 304, 501]].

        Returns:
            dict: A summary of the operation.
        """
        if not groups_to_merge:
            return {'groups_found': 0, 'nodes_merged': 0}

        remap_dict = {}
        nodes_to_delete = set()
        for group in groups_to_merge:
            master_nid = min(group)
            for nid in group:
                if nid != master_nid:
                    remap_dict[nid] = master_nid
                    nodes_to_delete.add(nid)

        all_elements = list(self.model.elements.values()) + \
                       list(self.model.rigid_elements.values())
        for element in all_elements:
            element.nodes = [remap_dict.get(nid, nid) for nid in element.nodes]

        for sid, spc_cards in self.model.spcs.items():
            for card in spc_cards:
                card.nodes = list(set([remap_dict.get(nid, nid) for nid in card.nodes]))

        for sid, load_cards in self.model.loads.items():
            for card in load_cards:
                if hasattr(card, 'node_id') and card.node_id in remap_dict:
                    card.node_id = remap_dict[card.node_id]
                elif hasattr(card, 'node') and card.node in remap_dict:
                    card.node = remap_dict[card.node]

        for nid in nodes_to_delete:
            if nid in self.model.nodes:
                del self.model.nodes[nid]

        return {
            'groups_found': len(groups_to_merge),
            'nodes_merged': len(nodes_to_delete)
        }

    




    def flip_shell_element_normals(self, eids_to_flip):
        """
        (DEBUGGING VERSION) Reverses node connectivity for shell elements.
        Includes extensive print statements to diagnose the filtering logic.
        """
        # --- Start of Debug Block ---
        print("\n--- DEBUG: Inside flip_shell_element_normals ---")
        print(f"Received {len(eids_to_flip)} EIDs to potentially flip: {eids_to_flip}")
        print("-------------------------------------------------")
        print("Checking element types for filtering...")

        valid_shells_to_flip = []
        for eid in eids_to_flip:
            if eid in self.model.elements:
                element = self.model.elements[eid]
                elem_type = element.type
                # This line prints the type of every element passed into the function
                print(f"  - EID {eid}: Type is '{elem_type}'")
                
                if elem_type in ['CQUAD4', 'CTRIA3', 'CMEMBRAN']:
                    valid_shells_to_flip.append(element)
                    print(f"    -> VALID. Added to flip list.")
                else:
                    print(f"    -> INVALID. Skipped.")
            else:
                print(f"  - EID {eid}: Not found in model.elements. Skipped.")
        
        print("-------------------------------------------------")
        print(f"Filter complete. Final list to be flipped contains {len(valid_shells_to_flip)} elements.")
        # --- End of Debug Block ---

        # The original flipping logic remains the same
        if not valid_shells_to_flip:
            return 0

        for element in valid_shells_to_flip:
            nodes = element.nodes
            if element.type in ['CQUAD4', 'CMEMBRAN'] and len(nodes) == 4:
                element.nodes = [nodes[0], nodes[3], nodes[2], nodes[1]]
            elif element.type == 'CTRIA3' and len(nodes) == 3:
                element.nodes = [nodes[0], nodes[2], nodes[1]]
        
        final_count = len(valid_shells_to_flip)
        print(f"Function finished. Returning count: {final_count}")
        print("--- END DEBUG ---\n")
        return final_count
    
    
    
    
    def transform_nodes(self, node_ids, params):
        if not (transform_type := params.get('type')) or not node_ids: return
        node_coords = np.array([self.model.nodes[nid].get_position() for nid in node_ids])
        centroid = np.mean(node_coords, axis=0)
        
        if transform_type == 'translate':
            delta = np.array(params.get('delta', [0,0,0]), dtype=float)
            for nid in node_ids: self.model.nodes[nid].xyz += delta
        elif transform_type == 'move_to_point':
            delta = np.array(params.get('target', [0,0,0]), dtype=float) - centroid
            for nid in node_ids: self.model.nodes[nid].xyz += delta
        elif transform_type == 'scale':
            factors = np.array(params.get('factors', [1,1,1]), dtype=float)
            center = np.array([0,0,0])
            if params.get('center_type') == 'centroid': center = centroid
            elif params.get('center_type') == 'custom': center = np.array(params.get('custom_center', [0,0,0]), dtype=float)
            for nid in node_ids: 
                original_pos = self.model.nodes[nid].get_position()
                new_pos = center + (original_pos - center) * factors
                self.model.nodes[nid].set_position(new_pos)
        elif transform_type == 'rotate':
            angle_rad, axis = np.radians(params.get('angle', 0.0)), params.get('axis', 'z')
            center = np.array([0,0,0])
            if params.get('center_type') == 'centroid': center = centroid
            elif params.get('center_type') == 'custom': center = np.array(params.get('custom_center', [0,0,0]), dtype=float)
            c, s = np.cos(angle_rad), np.sin(angle_rad)
            R = {'x': np.array([[1,0,0], [0,c,-s], [0,s,c]]), 'y': np.array([[c,0,s], [0,1,0], [-s,0,c]]), 'z': np.array([[c,-s,0], [s,c,0], [0,0,1]])}[axis]
            for nid in node_ids:
                original_pos = self.model.nodes[nid].get_position()
                new_pos = center + R @ (original_pos - center)
                self.model.nodes[nid].set_position(new_pos)

    def get_next_available_id(self, entity_type='node'):
        if entity_type == 'node':
            ids = self.model.nodes.keys()
        elif entity_type == 'element':
            ids = list(self.model.elements.keys()) + list(self.model.rigid_elements.keys())
        else:
            ids = []
        return max(ids) + 1 if ids else 1
        
    def add_node(self, node_id, x, y, z):
        if not node_id: node_id = self.get_next_available_id('node')
        return self.model.add_grid(node_id, [x, y, z]).nid
        
    def add_nodes_between(self, start_nid, end_nid, num_nodes_to_add):
        start_coord, end_coord = self.model.nodes[start_nid].get_position(), self.model.nodes[end_nid].get_position()
        new_nids = []
        for i in range(1, num_nodes_to_add + 1):
            fraction = i / (num_nodes_to_add + 1.0)
            new_coord = start_coord + fraction * (end_coord - start_coord)
            new_nids.append(self.add_node(None, *new_coord))
        return new_nids
        
    def add_line_element(self, eid, pid, n1, n2, type, orientation):
        if not eid: eid = self.get_next_available_id('element')
        g0, x = None, [0.,0.,1.] 
        if orientation and orientation.get('method') == 'node': g0 = orientation['values'][0]
        elif orientation and orientation.get('method') == 'vector': x = orientation['values']
        
        if type == 'CBEAM': return self.model.add_cbeam(eid, pid, [n1, n2], x, g0=g0).eid
        elif type == 'CBAR': return self.model.add_cbar(eid, pid, [n1, n2], x, g0=g0).eid
        elif type == 'CROD': return self.model.add_crod(eid, pid, [n1, n2]).eid
        
    def add_plate_element(self, eid, pid, nodes, type):
        if not eid: eid = self.get_next_available_id('element')
        if type == 'CQUAD4': return self.model.add_cquad4(eid, pid, nodes).eid
        elif type == 'CTRIA3': return self.model.add_ctria3(eid, pid, nodes).eid
        
   

    def add_rbe_element(self, eid, elem_type, indep_nodes, dep_nodes, dof):
        if not eid: eid = self.get_next_available_id('element')
        if elem_type == 'RBE2': return self.model.add_rbe2(eid, indep_nodes[0], dof, dep_nodes).eid
        elif elem_type == 'RBE3': return self.model.add_rbe3(eid, dep_nodes[0], dof, indep_nodes).eid



    def add_nodal_load(self, sid, node_ids, load_type, components):
        """Adds FORCE or MOMENT cards to a list of nodes.

        Args:
            sid (int): The load set ID.
            node_ids (list[int]): A list of node IDs to apply the load to.
            load_type (str): Either 'FORCE' or 'MOMENT'.
            components (list[float]): The [X, Y, Z] components of the load.

        Returns:
            int: The number of cards added.
        """
        mag = np.linalg.norm(components)
        if np.isclose(mag, 0.0):
            return 0  # Do not add a zero-magnitude load

        xyz = np.array(components) / mag
        added_count = 0
        
        if not isinstance(node_ids, list):
            node_ids = [node_ids]

        for nid in node_ids:
            if load_type.upper() == 'FORCE':
                self.model.add_force(sid, nid, mag, xyz)
                added_count += 1
            elif load_type.upper() == 'MOMENT':
                self.model.add_moment(sid, nid, mag, xyz)
                added_count += 1
        return added_count


    def add_pressure_load(self, sid, element_ids, pressure):
        """Adds a PLOAD4 card to a list of shell elements.

        Args:
            sid (int): The load set ID.
            element_ids (list[int]): List of CQUAD4 or CTRIA3 element IDs.
            pressure (float): The pressure value to apply.

        Returns:
            int: The number of elements the load was applied to.
        """
        if not element_ids:
            return 0
        
        for eid in element_ids:
            # --- FIX: Removed the non-existent 'is_blank_if_default' keyword argument. ---
            card_fields = ['PLOAD4', sid, eid, pressure]
            self.model.add_card(card_fields, 'PLOAD4')
            
        return len(element_ids)

    def add_default_temperature(self, sid, temperature):
        """Adds a TEMPD card for a default temperature.

        Args:
            sid (int): The temperature set ID.
            temperature (float): The temperature value.

        Returns:
            int: The number of cards added (always 1).
        """
        # --- FIX: Use the low-level add_card method to ensure correct formatting. ---
        # This prevents the writer from incorrectly converting a Python list to a string.
        card_fields = ['TEMPD', sid, temperature]
        self.model.add_card(card_fields, 'TEMPD')
        return 1

    def add_nodal_constraint(self, sid, node_ids, dof):
        """Adds an SPC1 card to constrain a list of nodes.

        Args:
            sid (int): The constraint set ID.
            node_ids (list[int]): A list of node IDs to constrain.
            dof (str): A string of constrained degrees of freedom (e.g., "123").

        Returns:
            int: The number of nodes constrained.
        """
        if not node_ids or not dof:
            return 0
        
        self.model.add_spc1(sid, str(dof), node_ids)
        return len(node_ids)
        

    def delete_elements(self, eids_to_delete):
        deleted_count = 0
        for eid in eids_to_delete:
            if self.model.elements.pop(eid, None):
                deleted_count += 1
            elif self.model.rigid_elements.pop(eid, None):
                deleted_count += 1
        return deleted_count


    def delete_nodes(self, nids_to_delete):
        nids_to_delete_set = set(nids_to_delete)
        
        # Find all elements connected to the nodes being deleted
        eids_to_delete = set()
        all_elements = {**self.model.elements, **self.model.rigid_elements}
        for eid, element in all_elements.items():
            if not nids_to_delete_set.isdisjoint(element.nodes):
                eids_to_delete.add(eid)
        
        # Delete the connected elements
        deleted_elements_count = self.delete_elements(list(eids_to_delete))

        # Delete the nodes
        deleted_nodes_count = 0
        for nid in nids_to_delete:
            if self.model.nodes.pop(nid, None):
                deleted_nodes_count += 1
        
        return {'nodes': deleted_nodes_count, 'elements': deleted_elements_count}
    

    def delete_loads_by_sid(self, sids_to_delete):
        """Deletes load sets from the model by their SID."""
        if not isinstance(sids_to_delete, list):
            sids_to_delete = [sids_to_delete]
        
        deleted_count = 0
        for sid in sids_to_delete:
            if self.model.loads.pop(sid, None):
                deleted_count += 1
        return deleted_count
    

    def delete_constraints_by_sid(self, sids_to_delete):
        """Deletes constraint sets from the model by their SID."""
        if not isinstance(sids_to_delete, list):
            sids_to_delete = [sids_to_delete]

        deleted_count = 0
        for sid in sids_to_delete:
            if self.model.spcs.pop(sid, None):
                # Also need to clear the spcadds/spc1s if they exist
                if sid in self.model.spcadds: self.model.spcadds.pop(sid, None)
                if sid in self.model.spcs: self.model.spcs.pop(sid, None) # legacy
                deleted_count += 1
        return deleted_count
    
      # --- NEW: Element Quality Calculation Methods (Phase 1) ---

    def calculate_element_quality(self):
        """
        Calculates quality metrics for all shell elements in the model.

        Returns:
            dict: A dictionary mapping element IDs to their quality metrics.
                  e.g., {101: {'warp': 0.05, 'aspect': 1.5, 'skew': 12.0}, ...}
        """
        quality_data = {}
        # Get all node coordinates once to avoid repeated lookups inside the loop
        all_node_coords = {nid: node.xyz for nid, node in self.model.nodes.items()}

        for eid, element in self.model.elements.items():
            try:
                # Get the actual 3D coordinates for each of the element's nodes
                node_coords = [all_node_coords[nid] for nid in element.nodes]
            except KeyError:
                continue # Skip elements with missing nodes

            metrics = {}
            if element.type in ['CQUAD4', 'CMEMBRAN'] and len(node_coords) == 4:
                metrics['warp'] = self._calculate_quad_warp(node_coords)
                metrics['aspect'] = self._calculate_aspect_ratio(node_coords)
                metrics['skew'] = self._calculate_quad_skew(node_coords)
            elif element.type == 'CTRIA3' and len(node_coords) == 3:
                # Triangles do not have warping
                metrics['aspect'] = self._calculate_aspect_ratio(node_coords)
                metrics['skew'] = self._calculate_tria_skew(node_coords)

            if metrics:
                quality_data[eid] = metrics

        return quality_data

    def _calculate_quad_warp(self, points):
        """Calculates the warping angle of a quad element in degrees."""
        p0, p1, p2, p3 = points
        # Create normals for two triangles formed by splitting the quad along a diagonal
        n1 = np.cross(p1 - p0, p2 - p0)
        n2 = np.cross(p2 - p0, p3 - p0)
        
        norm_n1 = np.linalg.norm(n1)
        norm_n2 = np.linalg.norm(n2)
        if norm_n1 < 1e-9 or norm_n2 < 1e-9: return 0.0

        # Find the angle between the two normals
        dot_product = np.dot(n1 / norm_n1, n2 / norm_n2)
        angle = np.rad2deg(np.arccos(np.clip(dot_product, -1.0, 1.0)))
        return angle

    def _calculate_aspect_ratio(self, points):
        """Calculates the ratio of the longest edge to the shortest edge."""
        num_points = len(points)
        edge_lengths = []
        for i in range(num_points):
            p1 = points[i]
            p2 = points[(i + 1) % num_points] # Wraps around for the last edge
            edge_lengths.append(np.linalg.norm(p1 - p2))
        
        min_edge = min(edge_lengths)
        if not edge_lengths or min_edge < 1e-9: return 1.0
        
        return max(edge_lengths) / min_edge

    def _calculate_quad_skew(self, points):
        """Calculates the skew of a quad as the deviation from 90 degrees between diagonals."""
        p0, p1, p2, p3 = points
        # Diagonals
        v1 = p2 - p0
        v2 = p3 - p1
        
        dot_product = np.dot(v1, v2)
        mag_prod = np.linalg.norm(v1) * np.linalg.norm(v2)
        if mag_prod < 1e-9: return 0.0

        angle = np.rad2deg(np.arccos(np.clip(dot_product / mag_prod, -1.0, 1.0)))
        return abs(90.0 - angle)

    def _calculate_tria_skew(self, points):
        """Calculates the skew of a triangle as the max deviation from an ideal angle (60 deg)."""
        p0, p1, p2 = points
        # Get edge lengths
        a = np.linalg.norm(p1 - p2)
        b = np.linalg.norm(p0 - p2)
        c = np.linalg.norm(p0 - p1)
        
        try: # Use the Law of Cosines to find the angles
            alpha = np.rad2deg(np.arccos((b**2 + c**2 - a**2) / (2 * b * c)))
            beta = np.rad2deg(np.arccos((a**2 + c**2 - b**2) / (2 * a * c)))
            gamma = 180.0 - alpha - beta
        except (ValueError, ZeroDivisionError):
            return 90.0 # Max possible skew for a degenerate triangle

        # Find the angle that deviates the most from the ideal 60 degrees
        min_angle = min(alpha, beta, gamma)
        max_angle = max(alpha, beta, gamma)
        return max(max_angle - 60, 60 - min_angle)
    



    def add_coordinate_system(self, cid, origin, z_axis_point, xz_plane_point):
        """
        Adds a CORD2R (rectangular) coordinate system to the model.

        Args:
            cid (int): The coordinate system ID (must be > 0).
            origin (list[float]): The [x, y, z] coordinates of the origin point.
            z_axis_point (list[float]): The [x, y, z] of a point on the new z-axis.
            xz_plane_point (list[float]): The [x, y, z] of a point in the new xz-plane.

        Returns:
            The created pyNastran CORD2R card object.
        """
        # --- FIX: Corrected keyword arguments from z_axis/xz_plane to zaxis/xzplane ---
        return self.model.add_cord2r(cid, rid=0, origin=origin, zaxis=z_axis_point, xzplane=xz_plane_point)



    def add_cylindrical_coord_system(self, cid, origin, z_axis_point, xz_plane_point):
        """Adds a CORD2C (cylindrical) coordinate system to the model."""
        return self.model.add_cord2c(cid, rid=0, origin=origin, zaxis=z_axis_point, xzplane=xz_plane_point, comment=f'$ Global Cylindrical')

    def add_spherical_coord_system(self, cid, origin, z_axis_point, xz_plane_point):
        """Adds a CORD2S (spherical) coordinate system to the model."""
        return self.model.add_cord2s(cid, rid=0, origin=origin, zaxis=z_axis_point, xzplane=xz_plane_point, comment=f'$ Global Spherical')


    def _build_rotation_matrix(self, rotations_deg):
        """Builds a 3x3 ZYX Euler rotation matrix from degrees."""
        rx, ry, rz = np.deg2rad(rotations_deg)
        Rx = np.array([[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]])
        Ry = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
        Rz = np.array([[np.cos(rz), -np.sin(rz), 0], [np.sin(rz), np.cos(rz), 0], [0, 0, 1]])
        # Combine in ZYX order
        return Rz @ Ry @ Rx

    def add_coord_by_translate_rotate(self, cid, ref_cid, translations, rotations, coord_type):
        """Creates a new coordinate system by transforming a reference system."""
        ref_coord = self.model.coords[ref_cid]
        
        # --- FIX: Changed ref_coord.Origin() to ref_coord.origin ---
        ref_origin = ref_coord.origin
        
        ref_matrix = ref_coord.beta() # Gets the 3x3 orientation matrix of the reference Csys

        # 1. Create the new orientation matrix by applying the new rotations
        new_orientation_matrix = self._build_rotation_matrix(rotations) @ ref_matrix

        # 2. Calculate the new origin by applying translation in the reference frame
        new_origin = ref_origin + (ref_matrix @ np.array(translations))

        # 3. Derive the 3 "defining points" required by Nastran from the new matrix
        z_axis_vector = new_orientation_matrix[:, 2] # The 3rd column is the Z-axis
        z_axis_point = new_origin + z_axis_vector

        x_axis_vector = new_orientation_matrix[:, 0] # The 1st column is the X-axis
        xz_plane_point = new_origin + x_axis_vector
        
        # 4. Call the correct low-level card-adding method based on the desired type
        if coord_type == 'rectangular':
            self.model.add_cord2r(cid, rid=ref_cid, origin=new_origin, zaxis=z_axis_point, xzplane=xz_plane_point)
        elif coord_type == 'cylindrical':
            self.model.add_cord2c(cid, rid=ref_cid, origin=new_origin, zaxis=z_axis_point, xzplane=xz_plane_point)
        elif coord_type == 'spherical':
            self.model.add_cord2s(cid, rid=ref_cid, origin=new_origin, zaxis=z_axis_point, xzplane=xz_plane_point)

        return self.model.coords[cid]
    

    def import_and_append_bdf(self, filepath):
        """
        Reads a BDF file and merges its contents into the current model.
        This method can handle both full decks and partial (punch) files.
        Any entities in the imported file with IDs that already exist in
        the current model will be overwritten.

        Args:
            filepath (str): The path to the BDF file to import.

        Returns:
            dict: A summary of the number of each entity type imported.
        """
        summary = {
            'nodes': 0, 'elements': 0, 'properties': 0, 'materials': 0,
            'loads': 0, 'constraints': 0, 'coords': 0
        }
        
        temp_model = BDF(debug=False)
        
        try:
            temp_model.read_bdf(filepath, punch=False)
        except Exception:
            try:
                temp_model.read_bdf(filepath, punch=True)
            except Exception as e:
                raise RuntimeError(f"Failed to parse BDF as full or partial deck: {e}")

        try:
            # Merge Nodes (GRID)
            for nid, node in temp_model.nodes.items():
                self.model.add_grid(nid, node.xyz, node.cp, node.cd, node.ps, node.seid, comment=node.comment)
                summary['nodes'] += 1

            # Merge Elements (CQUAD4, CTRIA3, CBEAM, etc.)
            for eid, element in temp_model.elements.items():
                # --- FIX: Changed write_card_list() to repr_fields() ---
                self.model.add_card(element.repr_fields(), element.type, is_list=True, comment=element.comment)
                summary['elements'] += 1
            for eid, element in temp_model.rigid_elements.items():
                # --- FIX: Changed write_card_list() to repr_fields() ---
                self.model.add_card(element.repr_fields(), element.type, is_list=True, comment=element.comment)
                summary['elements'] += 1

            # Merge Properties (PSHELL, PCOMP, etc.)
            for pid, prop in temp_model.properties.items():
                # --- FIX: Changed write_card_list() to repr_fields() ---
                self.model.add_card(prop.repr_fields(), prop.type, is_list=True, comment=prop.comment)
                summary['properties'] += 1

            # Merge Materials (MAT1, MAT8, etc.)
            for mid, mat in temp_model.materials.items():
                # --- FIX: Changed write_card_list() to repr_fields() ---
                self.model.add_card(mat.repr_fields(), mat.type, is_list=True, comment=mat.comment)
                summary['materials'] += 1
                
            # Merge Coordinate Systems (CORD2R, etc.)
            for cid, coord in temp_model.coords.items():
                if cid > 2:
                    # --- FIX: Changed write_card_list() to repr_fields() ---
                    self.model.add_card(coord.repr_fields(), coord.type, is_list=True, comment=coord.comment)
                    summary['coords'] += 1

            # Merge Loads (FORCE, PLOAD4, etc.)
            for sid, load_cards in temp_model.loads.items():
                if sid not in self.model.loads:
                    self.model.loads[sid] = []
                # For loads/spcs, we append the card objects directly
                for card in load_cards:
                    self.model.loads[sid].append(card)
                    summary['loads'] += 1

            # Merge Constraints (SPC, SPC1, etc.)
            for sid, spc_cards in temp_model.spcs.items():
                if sid not in self.model.spcs:
                    self.model.spcs[sid] = []
                for card in spc_cards:
                    self.model.spcs[sid].append(card)
                    summary['constraints'] += 1

        except Exception as e:
            raise RuntimeError(f"Failed to merge BDF data into current model: {e}")

        return summary


# START: New method in nas.py
    def add_cbush(self, eid, pid, n1, n2, orientation):
        """
        Adds a CBUSH element to the model.

        Args:
            eid (int): Element ID.
            pid (int): Property ID (must be a PBUSH).
            n1 (int): Node ID 1.
            n2 (int): Node ID 2 (can be None for grounding).
            orientation (dict): Defines the element orientation.
                {'method': 'vector', 'values': [x,y,z]}
                {'method': 'node', 'values': [g0]}
                {'method': 'cid', 'values': [cid]}
                {'method': 'default'}
        """
        if not eid:
            eid = self.get_next_available_id('element')

        nids = [n1, n2] if n2 is not None else [n1]
        x, g0, cid = None, None, None
        
        method = orientation.get('method')
        if method == 'vector':
            x = orientation['values']
        elif method == 'node':
            g0 = orientation['values'][0]
        elif method == 'cid':
            cid = orientation['values'][0]
        
        return self.model.add_cbush(eid, pid, nids, x=x, g0=g0, cid=cid).eid
# END: New method in nas.py

# START: Final replacement for get_cbush_orientation_matrix in nas.py
    def get_cbush_orientation_matrix(self, eid):
        """
        Calculates the 3x3 orientation matrix for a CBUSH element.
        
        Returns:
            np.ndarray: A 3x3 numpy array representing the rotation matrix.
                        Columns are the element's x, y, and z axes.
        """
        if eid not in self.model.elements or self.model.elements[eid].type != 'CBUSH':
            return np.identity(3)

        elem = self.model.elements[eid]
        n1 = elem.nodes[0]
        n2 = elem.nodes[1] if len(elem.nodes) > 1 else None
        
        p1 = self.model.nodes[n1].get_position()
        p2 = self.model.nodes[n2].get_position() if n2 is not None else p1

        # --- THIS IS THE DEFINITIVE FIX ---
        # A blank vector from a BDF is parsed as [None, None, None]. 
        # This check correctly identifies if a real vector has been provided.
        has_orientation_vector = elem.x is not None and not all(v is None for v in elem.x)

        # If no orientation method is specified, it aligns with the global system.
        if elem.g0 is None and not has_orientation_vector and elem.cid is None:
            return np.identity(3)

        # For all other cases, the element x-axis is from node 1 to node 2.
        x_axis = p2 - p1
        if np.linalg.norm(x_axis) < 1e-9: # Zero-length element case
            x_axis = np.array([1., 0., 0.])
        else:
            x_axis = x_axis / np.linalg.norm(x_axis)

        # Determine the vector that defines the XY plane
        v_plane = None
        if has_orientation_vector:
            v_plane = np.array(elem.x, dtype=float)
        elif elem.g0 is not None:
            p0 = self.model.nodes[elem.g0].get_position()
            v_plane = p0 - p1
        elif elem.cid is not None:
            v_plane = self.model.coords[elem.cid].beta()[:, 1]

        # Fallback logic for invalid or collinear orientation vectors
        if v_plane is None or np.linalg.norm(v_plane) < 1e-9:
            v_plane = np.array([0., 0., 1.])

        if np.linalg.norm(np.cross(x_axis, v_plane)) < 1e-9:
            v_plane = np.array([0., 1., 0.])
        
        v_plane = v_plane / np.linalg.norm(v_plane)

        # Calculate the final axes using cross products
        z_axis = np.cross(x_axis, v_plane)
        z_axis = z_axis / np.linalg.norm(z_axis)
        y_axis = np.cross(z_axis, x_axis)
        
        return np.column_stack([x_axis, y_axis, z_axis])
# END: Final replacement for get_cbush_orientation_matrix in nas.py


# START: New helper method in nas.py
    def get_beam_orientation_matrix(self, eid):
        """
        Calculates the 3x3 orientation matrix for a CBEAM or CBAR element.
        """
        if eid not in self.model.elements:
            return np.identity(3)
        
        elem = self.model.elements[eid]
        if elem.type not in ['CBEAM', 'CBAR']:
            return np.identity(3)

        n1, n2 = elem.nodes
        p1 = self.model.nodes[n1].get_position()
        p2 = self.model.nodes[n2].get_position()

        # The element x-axis is from node 1 to node 2.
        x_axis = p2 - p1
        if np.linalg.norm(x_axis) < 1e-9: # Zero-length element
            return np.identity(3)
        x_axis = x_axis / np.linalg.norm(x_axis)

        # Determine the vector that defines the XY plane
        v_plane = None
        if elem.x is not None and not all(v is None for v in elem.x):
            v_plane = np.array(elem.x, dtype=float)
        elif elem.g0 is not None:
            p0 = self.model.nodes[elem.g0].get_position()
            v_plane = p0 - p1
        
        if v_plane is None or np.linalg.norm(v_plane) < 1e-9 or np.linalg.norm(np.cross(x_axis, v_plane)) < 1e-9:
            # Fallback for collinear or undefined vectors
            v_plane = np.array([0., 0., 1.])
            if np.linalg.norm(np.cross(x_axis, v_plane)) < 1e-9:
                v_plane = np.array([0., 1., 0.])
        
        v_plane = v_plane / np.linalg.norm(v_plane)

        z_axis = np.cross(x_axis, v_plane)
        z_axis = z_axis / np.linalg.norm(z_axis)
        y_axis = np.cross(z_axis, x_axis)
        
        return np.column_stack([x_axis, y_axis, z_axis])
# END: New helper method in nas.py