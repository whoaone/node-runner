# nas.py

import math
import datetime
import numpy as np
from io import StringIO
from collections import namedtuple
from pyNastran.bdf.bdf import BDF
from scipy.spatial import cKDTree
from pyNastran.bdf.cards.elements.shell import CQUAD4, CTRIA3

# Data structures for lenient BDF import
SkippedCard = namedtuple('SkippedCard', ['line', 'card', 'raw', 'error'])
LenientResult = namedtuple('LenientResult', ['model', 'skipped', 'counts'])

# Result wrapper carrying field-format detection alongside the parsed model.
# A `None` LenientResult means strict parse succeeded; format is still set.
ImportResult = namedtuple('ImportResult', ['model', 'lenient_result', 'detected_format'])

# Valid Nastran field formats. Used by both the writer (size selection) and
# the reader (format detection result).
FIELD_FORMATS = ("short", "long", "free")


def _convert_fixed_line_to_free(line):
    """Convert one fixed-field BDF line into comma-separated free format.

    Handles both short (8-char) and long (16-char) field cards.  Keeps the
    line ending intact.  Trailing empty fields are dropped.
    """
    # Strip the line ending so we can put it back unchanged.
    if line.endswith('\r\n'):
        eol = '\r\n'
        body = line[:-2]
    elif line.endswith('\n'):
        eol = '\n'
        body = line[:-1]
    else:
        eol = ''
        body = line

    if not body.strip():
        return line

    # Card name is always the first 8 chars in fixed-field BDF.
    card_name = body[:8].strip()
    is_long = card_name.endswith('*')

    # Continuation cards: first field is blank or starts with '+' / '*'.
    # We treat them like normal cards for chunking purposes; the leading
    # token preserves the continuation marker.
    rest = body[8:]
    field_width = 16 if is_long else 8

    fields = [card_name]
    for i in range(0, len(rest), field_width):
        chunk = rest[i:i + field_width]
        fields.append(chunk.strip())

    # Drop trailing empty fields (keep at least the card name).
    while len(fields) > 1 and fields[-1] == '':
        fields.pop()

    return ','.join(fields) + eol


_EXEC_CASE_CONTROL_KEYWORDS = (
    'SOL ', 'SOL,', 'CEND', 'TITLE', 'SUBCASE', 'LABEL', 'ECHO', 'SUBTITLE',
    'TIME ', 'DIAG', 'ANALYSIS', 'DISP', 'STRESS', 'STRAIN', 'FORCE', 'MPCFORCE',
    'SPCFORCE', 'OLOAD', 'GPFORCE', 'ESE', 'METHOD', 'LOAD ', 'SPC ', 'MPC ',
    'INCLUDE', 'OUTPUT', 'SET ', 'PARAM,', 'PARAM ',
)


def _convert_bdf_to_free(path):
    """Rewrite a fixed-field BDF in place as comma-separated free format.

    Handles both full decks (with executive/case control + ``BEGIN BULK``)
    and punch decks (bulk data only). For full decks the executive- and
    case-control sections pass through untouched; for punch decks every
    non-comment line is treated as a bulk-data card.

    Comments (``$``-lines) are preserved verbatim. Lines that already
    contain a comma are assumed to be free-format already and skipped.
    """
    with open(path, 'r', encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()

    # Detect whether this is a full deck (has BEGIN BULK) or a punch file.
    has_begin_bulk = any(line.lstrip().upper().startswith('BEGIN BULK') for line in lines)

    out = []
    in_bulk = not has_begin_bulk  # punch -> entire file is bulk data
    for line in lines:
        upper = line.lstrip().upper()
        if upper.startswith('BEGIN BULK'):
            in_bulk = True
            out.append(line)
            continue
        if upper.startswith('ENDDATA'):
            in_bulk = False
            out.append(line)
            continue
        stripped = line.strip()
        if not stripped or stripped.startswith('$'):
            out.append(line)
            continue
        if not in_bulk:
            # Executive or case-control region of a full deck - leave intact.
            out.append(line)
            continue
        # Defensive: if a recognized exec/case keyword leaked into the bulk
        # region (mis-classified file), pass it through unchanged.
        if any(upper.startswith(kw) for kw in _EXEC_CASE_CONTROL_KEYWORDS):
            out.append(line)
            continue
        if ',' in line:
            out.append(line)
            continue
        out.append(_convert_fixed_line_to_free(line))

    with open(path, 'w', encoding='utf-8') as fh:
        fh.writelines(out)


def detect_bdf_field_format(filepath):
    """Inspect a BDF file and return one of 'short', 'long', 'free'.

    Heuristic, looks at the first few non-blank, non-comment, non-control
    bulk-data lines:
      - any comma in a card line  -> 'free'
      - card name ends with '*'   -> 'long' (e.g. 'GRID*')
      - else                      -> 'short' (default 8-char fixed)

    Returns 'short' as the safe default if the file can't be sniffed.
    """
    try:
        with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
            checked = 0
            for raw in fh:
                line = raw.rstrip('\n').rstrip('\r')
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith('$'):
                    continue
                upper = stripped.upper()
                # Skip executive/case-control keywords commonly above bulk data.
                if upper.startswith(('SOL ', 'CEND', 'BEGIN BULK', 'ENDDATA',
                                     'TITLE', 'SUBCASE', 'LABEL', 'ECHO',
                                     'SUBTITLE', 'INCLUDE', 'PARAM', 'TIME ',
                                     'DIAG', 'ANALYSIS')):
                    continue
                # Bulk-data line: classify.
                if ',' in line:
                    return 'free'
                # First field is the card name; if it ends with '*' it's long.
                head = line.lstrip().split(None, 1)[0] if line.strip() else ''
                if head.endswith('*'):
                    return 'long'
                # Confidence builds across a few cards; we just need one
                # decisive observation. Default at the end is 'short'.
                checked += 1
                if checked >= 5:
                    break
    except OSError:
        pass
    return 'short'


def _humanize_bdf_error(card_name, raw_msg):
    """Translate cryptic pyNastran exceptions into user-friendly messages."""
    import re
    msg = raw_msg.strip()

    # --- Duplicate ID errors ---
    # "self.elements IDs are not unique=[2862 3938]"
    m = re.match(
        r"self\.(elements|properties|materials|coords|thermal_materials"
        r"|massses|masses)\s+IDs are not unique\s*=\s*\[([^\]]*)\]", msg)
    if m:
        kind_map = {
            'elements': 'Element', 'properties': 'Property',
            'materials': 'Material', 'coords': 'Coordinate system',
            'thermal_materials': 'Thermal material',
            'massses': 'Mass element', 'masses': 'Mass element',
        }
        kind = kind_map.get(m.group(1), m.group(1).title())
        ids = m.group(2).strip()
        return f"Duplicate {kind} ID(s): [{ids}] - already defined earlier in the file"

    # --- Field type errors ---
    # "eid = '3.0' (field #1) on card must be an integer (not float)."
    m = re.search(
        r"(\w+)\s*=\s*'?([^']*?)'?\s*\(field #(\d+)\).*must be (an? \w+)",
        msg)
    if m:
        field_name = m.group(1)
        value = m.group(2)
        field_num = m.group(3)
        expected = m.group(4)
        return (f"Field '{field_name}' (#{field_num}) has value '{value}' "
                f"- expected {expected}")

    # --- Positive ID assertion ---
    # "eid=0 must be positive; elem=\n..."
    m = re.match(r"(\w+)=(-?\d+)\s+must be positive", msg)
    if m:
        field = m.group(1).upper()
        val = m.group(2)
        return f"Invalid {field} = {val} (must be a positive integer)"

    # --- Card too short / missing fields ---
    if 'len(card)' in msg or 'nfields' in msg:
        m = re.search(r'len\(card\)\s*=\s*(\d+)', msg)
        n = m.group(1) if m else '?'
        return f"Card has too few fields ({n}) - data may be truncated or malformed"

    # --- Unknown card type ---
    if 'not a valid key' in msg.lower() or 'not found' in msg.lower():
        return f"Unsupported card type '{card_name}' - not recognized by the parser"

    # --- General cleanup: strip "self." prefixes, trim multiline ---
    msg = msg.replace('self.', '')
    # Take only first line if multiline
    first_line = msg.split('\n')[0].strip()
    if len(first_line) > 160:
        first_line = first_line[:157] + '...'
    return first_line

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

    def _write_bdf(self, output_path, field_format='short'):
        """Write the model as a Nastran BDF in the requested field format.

        field_format: 'short' (size=8), 'long' (size=16), or 'free' (comma).
        For 'free', pyNastran's BDF.write_bdf is invoked at size=8 and the
        resulting file is post-processed in place by `_convert_bdf_to_free`.
        """
        if field_format not in FIELD_FORMATS:
            field_format = 'short'

        if field_format == 'long':
            self.model.write_bdf(output_path, size=16, is_double=False)
        else:
            self.model.write_bdf(output_path, size=8, is_double=False)
            if field_format == 'free':
                _convert_bdf_to_free(output_path)
    
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

    @staticmethod
    def calc_i_section_props(dims):
        """I/H beam: H (total height), W_top, W_bot, t_web, t_ftop, t_fbot."""
        H = dims['H']; Wt = dims['W_top']; Wb = dims['W_bot']
        tw = dims['t_web']; tft = dims['t_ftop']; tfb = dims['t_fbot']
        hw = H - tft - tfb
        A_top = Wt * tft; A_web = hw * tw; A_bot = Wb * tfb
        A = A_top + A_web + A_bot
        zc = (A_bot * tfb / 2 + A_web * (tfb + hw / 2) + A_top * (H - tft / 2)) / A
        I1 = ((Wb * tfb**3) / 12 + A_bot * (zc - tfb / 2)**2
              + (tw * hw**3) / 12 + A_web * (zc - tfb - hw / 2)**2
              + (Wt * tft**3) / 12 + A_top * (zc - H + tft / 2)**2)
        I2 = (tfb * Wb**3 + hw * tw**3 + tft * Wt**3) / 12
        J = (Wb * tfb**3 + hw * tw**3 + Wt * tft**3) / 3
        return {'A': A, 'I1': I1, 'I2': I2, 'J': J}

    @staticmethod
    def calc_c_section_props(dims):
        """C-channel: H (total height), W (flange width), t_web, t_flange."""
        H = dims['H']; W = dims['W']; tw = dims['t_web']; tf = dims['t_flange']
        hw = H - 2 * tf
        A_web = hw * tw; A_fl = W * tf
        A = A_web + 2 * A_fl
        yc = (A_web * tw / 2 + 2 * A_fl * W / 2) / A
        I1 = (tw * hw**3) / 12 + 2 * ((W * tf**3) / 12 + A_fl * (H / 2 - tf / 2)**2)
        I2 = ((hw * tw**3) / 12 + A_web * (yc - tw / 2)**2
              + 2 * ((tf * W**3) / 12 + A_fl * (yc - W / 2)**2))
        J = (hw * tw**3 + 2 * W * tf**3) / 3
        return {'A': A, 'I1': I1, 'I2': I2, 'J': J}

    @staticmethod
    def calc_box_section_props(dims):
        """Rectangular box/tube: H, W, t_web, t_flange."""
        H = dims['H']; W = dims['W']; tw = dims['t_web']; tf = dims['t_flange']
        A = 2 * W * tf + 2 * (H - 2 * tf) * tw
        hw = H - 2 * tf
        I1 = (W * H**3 - (W - 2 * tw) * hw**3) / 12
        I2 = (H * W**3 - hw * (W - 2 * tw)**3) / 12
        Am = (H - tf) * (W - tw)
        perimeter = 2 * ((H - tf) / tw + (W - tw) / tf)
        J = 4 * Am**2 / perimeter if perimeter > 0 else 0.
        return {'A': A, 'I1': I1, 'I2': I2, 'J': J}

    @staticmethod
    def calc_tube_section_props(dims):
        """Round tube: R (outer radius), t (wall thickness)."""
        import math
        R = dims['R']; t = dims['t']
        Ri = R - t
        A = math.pi * (R**2 - Ri**2)
        I = math.pi / 4 * (R**4 - Ri**4)
        J = math.pi / 2 * (R**4 - Ri**4)
        return {'A': A, 'I1': I, 'I2': I, 'J': J}

    @staticmethod
    def calc_solid_rect_props(dims):
        """Solid rectangle: H (height), W (width)."""
        H = dims['H']; W = dims['W']
        A = H * W
        I1 = W * H**3 / 12
        I2 = H * W**3 / 12
        a, b = max(H, W) / 2, min(H, W) / 2
        J = a * b**3 * (16 / 3 - 3.36 * b / a * (1 - b**4 / (12 * a**4)))
        return {'A': A, 'I1': I1, 'I2': I2, 'J': J}

    @staticmethod
    def calc_solid_circle_props(dims):
        """Solid circle: R (radius)."""
        import math
        R = dims['R']
        A = math.pi * R**2
        I = math.pi / 4 * R**4
        J = math.pi / 2 * R**4
        return {'A': A, 'I1': I, 'I2': I, 'J': J}


    # Section library: {pid: {'type': str, 'dims': dict}}
    # Populated when BeamSectionLibraryDialog is used.
    section_library = {}

    @staticmethod
    def _section_polygon(section_type, dims, n_circle=24):
        """Return a closed 2D polygon [(y, z), ...] for a cross-section type.

        The polygon is centred on the section centroid in y-z space.
        """
        if section_type in ('I / H Beam', 'I'):
            H = dims['H']; Wt = dims.get('W_top', dims.get('W', 1))
            Wb = dims.get('W_bot', dims.get('W', 1))
            tft = dims.get('t_ftop', dims.get('tf', 0.1))
            tfb = dims.get('t_fbot', dims.get('tf', 0.1))
            tw = dims.get('t_web', dims.get('tw', 0.1))
            # Outline going clockwise from bottom-left
            return [
                (-Wb/2, 0), (Wb/2, 0), (Wb/2, tfb),
                (tw/2, tfb), (tw/2, H-tft),
                (Wt/2, H-tft), (Wt/2, H),
                (-Wt/2, H), (-Wt/2, H-tft),
                (-tw/2, H-tft), (-tw/2, tfb),
                (-Wb/2, tfb), (-Wb/2, 0),
            ]
        elif section_type in ('C Channel', 'C'):
            H = dims['H']; W = dims['W']
            tf = dims.get('t_flange', dims.get('tf', 0.1))
            tw = dims.get('t_web', dims.get('tw', 0.1))
            return [
                (0, 0), (W, 0), (W, tf), (tw, tf),
                (tw, H-tf), (W, H-tf), (W, H),
                (0, H), (0, 0),
            ]
        elif section_type in ('Box (Rect Tube)', 'Box'):
            H = dims['H']; W = dims['W']
            tf = dims.get('t_flange', dims.get('tf', 0.1))
            tw = dims.get('t_web', dims.get('tw', 0.1))
            outer = [(-W/2, -H/2), (W/2, -H/2), (W/2, H/2), (-W/2, H/2), (-W/2, -H/2)]
            inner = [(-W/2+tw, -H/2+tf), (W/2-tw, -H/2+tf),
                     (W/2-tw, H/2-tf), (-W/2+tw, H/2-tf), (-W/2+tw, -H/2+tf)]
            return outer + [(None, None)] + inner  # hole marker
        elif section_type in ('Round Tube', 'Tube'):
            R = dims['R']; t = dims['t']
            import math
            angles = [2*math.pi*i/n_circle for i in range(n_circle+1)]
            outer = [(R*math.cos(a), R*math.sin(a)) for a in angles]
            Ri = R - t
            inner = [(Ri*math.cos(a), Ri*math.sin(a)) for a in angles]
            return outer + [(None, None)] + inner
        elif section_type in ('Solid Rectangle', 'Rect'):
            H = dims['H']; W = dims['W']
            return [(-W/2, -H/2), (W/2, -H/2), (W/2, H/2), (-W/2, H/2), (-W/2, -H/2)]
        elif section_type in ('Solid Circle', 'Circle'):
            import math
            R = dims['R']
            angles = [2*math.pi*i/n_circle for i in range(n_circle+1)]
            return [(R*math.cos(a), R*math.sin(a)) for a in angles]
        elif section_type in ('L',):
            H = dims.get('H', 1); W = dims.get('W', 1)
            tf = dims.get('tf', 0.1); tw = dims.get('tw', 0.1)
            return [(0,0), (W,0), (W,tf), (tw,tf), (tw,H), (0,H), (0,0)]
        elif section_type in ('T',):
            H = dims.get('H', 1); W = dims.get('W', 1)
            tf = dims.get('tf', 0.1); tw = dims.get('tw', 0.1)
            return [(-tw/2,0), (tw/2,0), (tw/2,H-tf),
                    (W/2,H-tf), (W/2,H), (-W/2,H),
                    (-W/2,H-tf), (-tw/2,H-tf), (-tw/2,0)]
        elif section_type in ('Z',):
            H = dims.get('H', 1); W = dims.get('W', 1)
            tf = dims.get('tf', 0.1); tw = dims.get('tw', 0.1)
            return [(0,0), (W,0), (W,tf), (tw,tf), (tw,H-tf),
                    (W,H-tf), (W,H), (0,H), (0,H-tf),
                    (-W+tw,H-tf), (-W+tw,tf), (0,tf), (0,0)]
        else:
            # Fallback: small circle
            import math
            R = max(dims.get('R', 0.5), 0.5)
            angles = [2*math.pi*i/12 for i in range(13)]
            return [(R*math.cos(a), R*math.sin(a)) for a in angles]

    def get_beam_section_polygon(self, pid):
        """Return cross-section polygon for a beam property.

        Uses section_library if available, otherwise approximates from
        property cross-section area.

        Returns:
            list of (y, z) or None if no section info available.
        """
        if pid in self.section_library:
            info = self.section_library[pid]
            return self._section_polygon(info['type'], info['dims'])

        # Try to infer from property
        prop = self.model.properties.get(pid)
        if prop is None:
            return None
        A = getattr(prop, 'A', None) or getattr(prop, 'area', None)
        if A and A > 0:
            import math
            R = math.sqrt(A / math.pi)
            return self._section_polygon('Circle', {'R': R})
        return None

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
            ids = list(self.model.elements.keys()) + list(self.model.rigid_elements.keys()) + list(self.model.masses.keys()) + list(self.model.plotels.keys())
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

    def add_solid_element(self, eid, pid, nodes, type):
        """Add a CHEXA, CTETRA, or CPENTA solid element."""
        if not eid: eid = self.get_next_available_id('element')
        if type == 'CHEXA': return self.model.add_chexa(eid, pid, nodes).eid
        elif type == 'CTETRA': return self.model.add_ctetra(eid, pid, nodes).eid
        elif type == 'CPENTA': return self.model.add_cpenta(eid, pid, nodes).eid

    def add_cshear(self, eid, pid, nodes):
        """Add a CSHEAR shear panel element."""
        if not eid: eid = self.get_next_available_id('element')
        return self.model.add_cshear(eid, pid, nodes).eid

    def add_cgap(self, eid, pid, n1, n2, orientation):
        """Add a CGAP gap element."""
        if not eid: eid = self.get_next_available_id('element')
        g0, x = None, [0., 0., 1.]
        if orientation and orientation.get('method') == 'node': g0 = orientation['values'][0]
        elif orientation and orientation.get('method') == 'vector': x = orientation['values']
        return self.model.add_cgap(eid, pid, [n1, n2], x, g0).eid

    def add_plotel(self, eid, nodes):
        """Add a PLOTEL visualization element."""
        if not eid: eid = self.get_next_available_id('element')
        return self.model.add_plotel(eid, nodes).eid

    def add_rbe_element(self, eid, elem_type, indep_nodes, dep_nodes, dof):
        if not eid: eid = self.get_next_available_id('element')
        if elem_type == 'RBE2': return self.model.add_rbe2(eid, indep_nodes[0], dof, dep_nodes).eid
        elif elem_type == 'RBE3': return self.model.add_rbe3(eid, dep_nodes[0], dof, indep_nodes).eid

    def add_conm2(self, eid, nid, mass, cid=0, X=None, I=None):
        """Add a CONM2 concentrated mass element.

        Args:
            eid (int): Element ID (0 for auto).
            nid (int): Node ID where the mass is applied.
            mass (float): The mass value.
            cid (int): Coordinate system ID for offset (default 0).
            X (list[float]): Offset [X1, X2, X3] (default [0,0,0]).
            I (list[float]): Inertia [I11, I21, I22, I31, I32, I33] (default all zeros).

        Returns:
            int: The created element ID.
        """
        if not eid: eid = self.get_next_available_id('element')
        return self.model.add_conm2(eid, nid, mass, cid=cid, X=X, I=I).eid

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

    def add_gravity_load(self, sid, scale, N, cid=0):
        """Adds a GRAV card for a gravity load.

        Args:
            sid (int): The load set ID.
            scale (float): Acceleration magnitude.
            N (list[float]): Direction vector [Nx, Ny, Nz].
            cid (int): Coordinate system ID (default 0).

        Returns:
            int: The number of cards added (always 1).
        """
        self.model.add_grav(sid, scale, N, cid=cid)
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
                metrics['jacobian'] = self._calculate_quad_jacobian_ratio(node_coords)
                metrics['taper'] = self._calculate_quad_taper(node_coords)
                min_a, max_a = self._calculate_interior_angles(node_coords)
                metrics['min_angle'] = min_a
                metrics['max_angle'] = max_a
            elif element.type == 'CTRIA3' and len(node_coords) == 3:
                # Triangles do not have warping, jacobian, or taper
                metrics['aspect'] = self._calculate_aspect_ratio(node_coords)
                metrics['skew'] = self._calculate_tria_skew(node_coords)
                min_a, max_a = self._calculate_interior_angles(node_coords)
                metrics['min_angle'] = min_a
                metrics['max_angle'] = max_a
            elif element.type in ('CTETRA', 'CTETRA4', 'CTETRA10') and len(node_coords) >= 4:
                pts = node_coords[:4]
                metrics['aspect'] = self._calculate_solid_aspect_ratio(pts, edges=(
                    (0, 1), (1, 2), (2, 0), (0, 3), (1, 3), (2, 3),
                ))
                metrics['volume'] = self._calculate_tet_volume(pts)
                metrics['jacobian'] = self._calculate_tet_jacobian(pts)
            elif element.type in ('CHEXA', 'CHEXA8', 'CHEXA20') and len(node_coords) >= 8:
                pts = node_coords[:8]
                metrics['aspect'] = self._calculate_solid_aspect_ratio(pts, edges=(
                    (0, 1), (1, 2), (2, 3), (3, 0),
                    (4, 5), (5, 6), (6, 7), (7, 4),
                    (0, 4), (1, 5), (2, 6), (3, 7),
                ))
                metrics['volume'] = self._calculate_hex_volume(pts)
                metrics['jacobian'] = self._calculate_hex_jacobian(pts)
            elif element.type in ('CPENTA', 'CPENTA6', 'CPENTA15') and len(node_coords) >= 6:
                pts = node_coords[:6]
                metrics['aspect'] = self._calculate_solid_aspect_ratio(pts, edges=(
                    (0, 1), (1, 2), (2, 0),
                    (3, 4), (4, 5), (5, 3),
                    (0, 3), (1, 4), (2, 5),
                ))
                metrics['volume'] = self._calculate_penta_volume(pts)

            if metrics:
                quality_data[eid] = metrics

        return quality_data

    @staticmethod
    def _calculate_solid_aspect_ratio(points, edges):
        """Generic solid aspect ratio: max edge / min edge."""
        ls = []
        for a, b in edges:
            ls.append(float(np.linalg.norm(points[a] - points[b])))
        ls = [l for l in ls if l > 1e-12]
        if not ls:
            return 0.0
        return max(ls) / min(ls)

    @staticmethod
    def _calculate_tet_volume(points):
        p0, p1, p2, p3 = points
        return float(abs(np.dot(p1 - p0, np.cross(p2 - p0, p3 - p0))) / 6.0)

    @staticmethod
    def _calculate_tet_jacobian(points):
        """Min Jacobian determinant ratio across the 4 corners of a tet."""
        p0, p1, p2, p3 = points
        # Tet has constant Jacobian for linear shape functions, but we
        # report the determinant magnitude normalized by ideal-tet volume
        # of the same edge length.
        edge_lens = [
            float(np.linalg.norm(p1 - p0)),
            float(np.linalg.norm(p2 - p0)),
            float(np.linalg.norm(p3 - p0)),
        ]
        if min(edge_lens) < 1e-12:
            return 0.0
        v = abs(np.dot(p1 - p0, np.cross(p2 - p0, p3 - p0)))
        # Ideal tet volume scaled to mean edge length
        mean_edge = sum(edge_lens) / 3.0
        ideal_v = (mean_edge ** 3) / (6.0 * np.sqrt(2.0))
        if ideal_v < 1e-12:
            return 0.0
        return float(v / 6.0 / ideal_v)

    @staticmethod
    def _calculate_hex_volume(points):
        """Hex volume via 6-tet decomposition sharing diagonal 0-6.

        For a unit cube this returns 1.0 exactly. The decomposition uses
        the body diagonal from corner 0 to corner 6 as the common edge
        across all six tets.
        """
        p = points
        tets = (
            (0, 1, 2, 6),
            (0, 2, 3, 6),
            (0, 3, 7, 6),
            (0, 7, 4, 6),
            (0, 4, 5, 6),
            (0, 5, 1, 6),
        )
        v = 0.0
        for a, b, c, d in tets:
            v += abs(np.dot(p[b] - p[a], np.cross(p[c] - p[a], p[d] - p[a]))) / 6.0
        return float(v)

    @staticmethod
    def _calculate_hex_jacobian(points):
        """Min/max Jacobian determinant ratio at the 8 corners of a hex.

        At each corner we form the local edge frame from the three edges
        leaving that corner and compute its determinant. Returns
        min(det) / max(det) - 1.0 is a perfect cube; <0 indicates inverted.
        """
        p = points
        # Corner -> three edge-neighbors (along x, y, z faces of unit hex)
        corners = (
            (0, 1, 3, 4),
            (1, 2, 0, 5),
            (2, 3, 1, 6),
            (3, 0, 2, 7),
            (4, 5, 7, 0),
            (5, 6, 4, 1),
            (6, 7, 5, 2),
            (7, 4, 6, 3),
        )
        dets = []
        for c, a, b, d in corners:
            v1 = points[a] - points[c]
            v2 = points[b] - points[c]
            v3 = points[d] - points[c]
            dets.append(float(np.dot(v1, np.cross(v2, v3))))
        dets = [abs(x) for x in dets if abs(x) > 1e-12]
        if not dets:
            return 0.0
        return min(dets) / max(dets)

    @staticmethod
    def _calculate_penta_volume(points):
        """Wedge / penta volume via 3-tet decomposition."""
        p = points
        tets = (
            (0, 1, 2, 5),
            (0, 1, 5, 4),
            (0, 4, 5, 3),
        )
        v = 0.0
        for a, b, c, d in tets:
            v += abs(np.dot(p[b] - p[a], np.cross(p[c] - p[a], p[d] - p[a]))) / 6.0
        return float(v)

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

    def _calculate_quad_jacobian_ratio(self, points):
        """Calculates the Jacobian ratio of a quad (min/max determinant at corners).
        Projects the 3D quad onto its local 2D plane, then evaluates the Jacobian
        at each isoparametric corner. Returns min(J)/max(J). 1.0 = perfect."""
        p0, p1, p2, p3 = [np.array(p) for p in points]

        # Build local 2D coordinate system on the element plane
        center = (p0 + p1 + p2 + p3) / 4.0
        e1 = (p1 - p0) + (p2 - p3)
        e1 = e1 / (np.linalg.norm(e1) + 1e-30)
        normal = np.cross(p2 - p0, p3 - p1)
        normal = normal / (np.linalg.norm(normal) + 1e-30)
        e2 = np.cross(normal, e1)
        e2 = e2 / (np.linalg.norm(e2) + 1e-30)

        # Project nodes to 2D
        pts_2d = []
        for p in [p0, p1, p2, p3]:
            d = p - center
            pts_2d.append([np.dot(d, e1), np.dot(d, e2)])
        x = [pt[0] for pt in pts_2d]
        y = [pt[1] for pt in pts_2d]

        # Evaluate Jacobian determinant at 4 isoparametric corners
        # Shape function derivatives: dN/dxi, dN/deta at corners (xi, eta)
        corners = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        jacobians = []
        for xi, eta in corners:
            # dN/dxi = [-(1-eta), (1-eta), (1+eta), -(1+eta)] / 4
            # dN/deta = [-(1-xi), -(1+xi), (1+xi), (1-xi)] / 4
            dNdxi = [-(1 - eta) / 4, (1 - eta) / 4, (1 + eta) / 4, -(1 + eta) / 4]
            dNdeta = [-(1 - xi) / 4, -(1 + xi) / 4, (1 + xi) / 4, (1 - xi) / 4]

            dxdxi = sum(dNdxi[i] * x[i] for i in range(4))
            dydxi = sum(dNdxi[i] * y[i] for i in range(4))
            dxdeta = sum(dNdeta[i] * x[i] for i in range(4))
            dydeta = sum(dNdeta[i] * y[i] for i in range(4))

            det_j = dxdxi * dydeta - dydxi * dxdeta
            jacobians.append(abs(det_j))

        max_j = max(jacobians)
        if max_j < 1e-30:
            return 0.0
        return min(jacobians) / max_j

    def _calculate_quad_taper(self, points):
        """Calculates the taper ratio of a quad element.
        Splits quad along both diagonals, compares area ratios. 1.0 = perfect."""
        p0, p1, p2, p3 = [np.array(p) for p in points]

        def tri_area(a, b, c):
            return np.linalg.norm(np.cross(b - a, c - a)) / 2.0

        # Diagonal 0-2: triangles (0,1,2) and (0,2,3)
        a1 = tri_area(p0, p1, p2)
        a2 = tri_area(p0, p2, p3)
        # Diagonal 1-3: triangles (0,1,3) and (1,2,3)
        a3 = tri_area(p0, p1, p3)
        a4 = tri_area(p1, p2, p3)

        ratios = []
        if min(a1, a2) > 1e-30:
            ratios.append(max(a1, a2) / min(a1, a2))
        if min(a3, a4) > 1e-30:
            ratios.append(max(a3, a4) / min(a3, a4))

        return max(ratios) if ratios else 1.0

    def _calculate_interior_angles(self, points):
        """Calculates min and max interior angles of a polygon (quad or tria).
        Returns (min_angle_degrees, max_angle_degrees)."""
        n = len(points)
        angles = []
        for i in range(n):
            p_prev = np.array(points[(i - 1) % n])
            p_curr = np.array(points[i])
            p_next = np.array(points[(i + 1) % n])

            v1 = p_prev - p_curr
            v2 = p_next - p_curr
            len1 = np.linalg.norm(v1)
            len2 = np.linalg.norm(v2)
            if len1 < 1e-12 or len2 < 1e-12:
                angles.append(0.0)
                continue
            cos_a = np.clip(np.dot(v1, v2) / (len1 * len2), -1.0, 1.0)
            angles.append(np.degrees(np.arccos(cos_a)))

        return (min(angles), max(angles)) if angles else (0.0, 0.0)

    def calculate_mass_summary(self):
        """Calculates total mass, center of gravity, and per-property mass breakdown.

        Returns:
            dict with keys: total_mass, cg (3-element list), conm2_mass,
            structural_mass, by_pid ({pid: mass}), notes (list of strings).
        """
        model = self.model
        all_node_coords = {nid: node.xyz for nid, node in model.nodes.items()}
        by_pid = {}
        total_weighted_pos = np.zeros(3)
        total_mass = 0.0
        conm2_mass = 0.0
        notes = []

        # --- CONM2 masses ---
        for eid, mass_elem in model.masses.items():
            m = mass_elem.mass
            conm2_mass += m
            nid = mass_elem.node_ids[0] if hasattr(mass_elem, 'node_ids') else mass_elem.nid
            if nid in all_node_coords:
                total_weighted_pos += m * np.array(all_node_coords[nid])
            total_mass += m

        # --- Shell elements (CQUAD4, CTRIA3, CMEMBRAN) ---
        shell_types = {'CQUAD4', 'CTRIA3', 'CMEMBRAN'}
        for eid, elem in model.elements.items():
            if elem.type not in shell_types:
                continue
            try:
                pid = elem.pid
                prop = model.properties.get(pid)
                if prop is None:
                    continue

                # Get thickness
                thickness = getattr(prop, 't', None)
                if thickness is None or thickness <= 0:
                    continue

                # Get material density
                mid = getattr(prop, 'mid', None) or getattr(prop, 'mid1', None)
                if mid is None:
                    continue
                mat = model.materials.get(mid)
                if mat is None:
                    continue
                rho = getattr(mat, 'rho', 0.0)
                if rho <= 0:
                    continue

                # Compute area
                coords = [np.array(all_node_coords[nid]) for nid in elem.nodes]
                if len(coords) == 4:
                    area = np.linalg.norm(np.cross(coords[2] - coords[0], coords[3] - coords[1])) / 2.0
                elif len(coords) == 3:
                    area = np.linalg.norm(np.cross(coords[1] - coords[0], coords[2] - coords[0])) / 2.0
                else:
                    continue

                m = area * thickness * rho
                by_pid[pid] = by_pid.get(pid, 0.0) + m
                centroid = np.mean(coords, axis=0)
                total_weighted_pos += m * centroid
                total_mass += m
            except (KeyError, TypeError, AttributeError):
                continue

        # --- Beam/Bar elements (CBEAM, CBAR) ---
        beam_types = {'CBEAM', 'CBAR'}
        for eid, elem in model.elements.items():
            if elem.type not in beam_types:
                continue
            try:
                pid = elem.pid
                prop = model.properties.get(pid)
                if prop is None:
                    continue

                area = getattr(prop, 'A', None) or getattr(prop, 'area', None)
                if area is None or area <= 0:
                    continue

                mid = getattr(prop, 'mid', None) or getattr(prop, 'mid1', None)
                if mid is None:
                    continue
                mat = model.materials.get(mid)
                if mat is None:
                    continue
                rho = getattr(mat, 'rho', 0.0)
                if rho <= 0:
                    continue

                coords = [np.array(all_node_coords[nid]) for nid in elem.nodes[:2]]
                length = np.linalg.norm(coords[1] - coords[0])

                m = area * length * rho
                by_pid[pid] = by_pid.get(pid, 0.0) + m
                centroid = (coords[0] + coords[1]) / 2.0
                total_weighted_pos += m * centroid
                total_mass += m
            except (KeyError, TypeError, AttributeError):
                continue

        # --- Rod elements (CROD) ---
        for eid, elem in model.elements.items():
            if elem.type != 'CROD':
                continue
            try:
                pid = elem.pid
                prop = model.properties.get(pid)
                if prop is None:
                    continue

                area = getattr(prop, 'A', None)
                if area is None or area <= 0:
                    continue

                mid = getattr(prop, 'mid', None)
                if mid is None:
                    continue
                mat = model.materials.get(mid)
                if mat is None:
                    continue
                rho = getattr(mat, 'rho', 0.0)
                if rho <= 0:
                    continue

                coords = [np.array(all_node_coords[nid]) for nid in elem.nodes[:2]]
                length = np.linalg.norm(coords[1] - coords[0])

                m = area * length * rho
                by_pid[pid] = by_pid.get(pid, 0.0) + m
                centroid = (coords[0] + coords[1]) / 2.0
                total_weighted_pos += m * centroid
                total_mass += m
            except (KeyError, TypeError, AttributeError):
                continue

        # --- Compute CG ---
        structural_mass = total_mass - conm2_mass
        cg = (total_weighted_pos / total_mass).tolist() if total_mass > 1e-30 else [0.0, 0.0, 0.0]

        # Check for solids
        solid_types = {'CHEXA', 'CTETRA', 'CPENTA'}
        has_solids = any(e.type in solid_types for e in model.elements.values())
        if has_solids:
            notes.append("Solid element mass not included (volume calculation not supported).")

        # --- Inertia tensor about CG ---
        # Collect (mass, centroid) pairs for a second pass
        mass_centroids = []

        # Rebuild centroid list (reuse logic above but store pairs)
        for eid, mass_elem in model.masses.items():
            m = mass_elem.mass
            nid = mass_elem.node_ids[0] if hasattr(mass_elem, 'node_ids') else mass_elem.nid
            if nid in all_node_coords and m > 0:
                mass_centroids.append((m, np.array(all_node_coords[nid])))

        for eid, elem in model.elements.items():
            if elem.type not in shell_types | beam_types | {'CROD'}:
                continue
            try:
                pid = elem.pid
                prop = model.properties.get(pid)
                if prop is None:
                    continue
                if elem.type in shell_types:
                    t = getattr(prop, 't', None)
                    if not t or t <= 0:
                        continue
                    mid = getattr(prop, 'mid', None) or getattr(prop, 'mid1', None)
                else:
                    t = None
                    area = getattr(prop, 'A', None) or getattr(prop, 'area', None)
                    if not area or area <= 0:
                        continue
                    mid = getattr(prop, 'mid', None) or getattr(prop, 'mid1', None)
                if mid is None:
                    continue
                mat = model.materials.get(mid)
                if mat is None:
                    continue
                rho = getattr(mat, 'rho', 0.0)
                if rho <= 0:
                    continue

                coords = [np.array(all_node_coords[nid]) for nid in elem.nodes[:4 if elem.type in shell_types else 2]]
                if elem.type in shell_types:
                    if len(coords) == 4:
                        a = np.linalg.norm(np.cross(coords[2] - coords[0], coords[3] - coords[1])) / 2.0
                    elif len(coords) == 3:
                        a = np.linalg.norm(np.cross(coords[1] - coords[0], coords[2] - coords[0])) / 2.0
                    else:
                        continue
                    m = a * t * rho
                else:
                    length = np.linalg.norm(coords[1] - coords[0])
                    m = area * length * rho
                centroid = np.mean(coords, axis=0)
                if m > 0:
                    mass_centroids.append((m, centroid))
            except (KeyError, TypeError, AttributeError):
                continue

        inertia = np.zeros(6)  # [Ixx, Iyy, Izz, Ixy, Ixz, Iyz]
        cg_arr = np.array(cg)
        for m, c in mass_centroids:
            r = c - cg_arr
            inertia[0] += m * (r[1]**2 + r[2]**2)  # Ixx
            inertia[1] += m * (r[0]**2 + r[2]**2)  # Iyy
            inertia[2] += m * (r[0]**2 + r[1]**2)  # Izz
            inertia[3] -= m * r[0] * r[1]           # Ixy
            inertia[4] -= m * r[0] * r[2]           # Ixz
            inertia[5] -= m * r[1] * r[2]           # Iyz

        return {
            'total_mass': total_mass,
            'cg': cg,
            'conm2_mass': conm2_mass,
            'structural_mass': structural_mass,
            'by_pid': by_pid,
            'notes': notes,
            'inertia': inertia.tolist(),
        }

    def find_orphans(self):
        """Find unused properties, unused materials, and unreferenced nodes.

        Returns:
            dict with keys:
                unused_pids: list of property IDs not referenced by any element
                unused_mids: list of material IDs not referenced by any property
                orphan_nids: list of node IDs not used by elements, masses, or SPCs
        """
        model = self.model

        # PIDs referenced by elements
        used_pids = set()
        for elem in model.elements.values():
            pid = getattr(elem, 'pid', None)
            if pid is not None:
                used_pids.add(pid)
        unused_pids = sorted(set(model.properties.keys()) - used_pids)

        # MIDs referenced by properties
        used_mids = set()
        for prop in model.properties.values():
            for attr in ('mid', 'mid1', 'mid2', 'mid3', 'mid4'):
                mid = getattr(prop, attr, None)
                if mid is not None and mid > 0:
                    used_mids.add(mid)
        unused_mids = sorted(set(model.materials.keys()) - used_mids)

        # Nodes used by elements, masses, rigid elements, SPCs
        used_nids = set()
        for elem in model.elements.values():
            used_nids.update(elem.nodes)
        for elem in model.rigid_elements.values():
            if hasattr(elem, 'independent_nodes'):
                used_nids.update(elem.independent_nodes)
            if hasattr(elem, 'dependent_nodes'):
                used_nids.update(elem.dependent_nodes)
            if hasattr(elem, 'nodes'):
                used_nids.update(elem.nodes)
        for mass_elem in model.masses.values():
            nid = mass_elem.node_ids[0] if hasattr(mass_elem, 'node_ids') else getattr(mass_elem, 'nid', None)
            if nid is not None:
                used_nids.add(nid)
        for sid, spc_list in model.spcs.items():
            for spc in spc_list:
                if hasattr(spc, 'nodes'):
                    used_nids.update(spc.nodes)
                elif hasattr(spc, 'node_ids'):
                    used_nids.update(spc.node_ids)
        orphan_nids = sorted(set(model.nodes.keys()) - used_nids)

        return {
            'unused_pids': unused_pids,
            'unused_mids': unused_mids,
            'orphan_nids': orphan_nids,
        }

    def import_cad_file(self, filepath, mesh_size, elem_preference, pid):
        """Import STEP/IGES file via gmsh, creating meshed shell elements.

        Returns:
            tuple: (created_nids, created_eids)
        """
        try:
            import gmsh
        except ImportError:
            raise RuntimeError("gmsh is not installed. Run: pip install gmsh")

        created_nids, created_eids = [], []

        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("cad_import")

        try:
            gmsh.model.occ.importShapes(filepath)
            gmsh.model.occ.synchronize()

            gmsh.option.setNumber("Mesh.MeshSizeMin", mesh_size * 0.5)
            gmsh.option.setNumber("Mesh.MeshSizeMax", mesh_size * 1.5)
            gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 12)
            if elem_preference == 'quad':
                gmsh.option.setNumber("Mesh.RecombineAll", 1)
                gmsh.option.setNumber("Mesh.Algorithm", 8)

            gmsh.model.mesh.generate(2)

            node_tags, coords, _ = gmsh.model.mesh.getNodes()
            coords = coords.reshape(-1, 3)

            gmsh_to_nid = {}
            for i, gtag in enumerate(node_tags):
                nid = self.add_node(
                    None, float(coords[i, 0]), float(coords[i, 1]),
                    float(coords[i, 2]))
                gmsh_to_nid[int(gtag)] = nid
                created_nids.append(nid)

            elem_types, elem_tags_list, elem_node_tags_list = \
                gmsh.model.mesh.getElements(dim=2)

            for et, e_tags, en_tags in zip(
                    elem_types, elem_tags_list, elem_node_tags_list):
                n_per_elem = len(en_tags) // len(e_tags)
                for j in range(len(e_tags)):
                    gnodes = [int(en_tags[j * n_per_elem + k])
                              for k in range(n_per_elem)]
                    fem_nodes = [gmsh_to_nid[g] for g in gnodes]
                    if n_per_elem == 4:
                        eid = self.add_plate_element(
                            None, pid, fem_nodes, 'CQUAD4')
                    elif n_per_elem == 3:
                        eid = self.add_plate_element(
                            None, pid, fem_nodes, 'CTRIA3')
                    else:
                        continue
                    created_eids.append(eid)
        finally:
            gmsh.finalize()

        return created_nids, created_eids

    def import_stl_file(self, filepath, pid):
        """Import STL file, deduplicating vertices, creating CTRIA3 elements.

        Returns:
            tuple: (created_nids, created_eids)
        """
        import struct

        created_nids, created_eids = [], []

        # Try binary first, fall back to ASCII
        vertices, faces = [], []
        try:
            with open(filepath, 'rb') as f:
                header = f.read(80)
                n_tris = struct.unpack('<I', f.read(4))[0]
                for _ in range(n_tris):
                    data = f.read(50)  # 12 floats + 2 bytes attribute
                    vals = struct.unpack('<12fH', data)
                    # Skip normal (vals[0:3]), get 3 vertices
                    v1 = vals[3:6]
                    v2 = vals[6:9]
                    v3 = vals[9:12]
                    base = len(vertices)
                    vertices.extend([v1, v2, v3])
                    faces.append((base, base + 1, base + 2))
        except (struct.error, ValueError):
            # ASCII STL
            vertices, faces = [], []
            with open(filepath, 'r') as f:
                current_tri = []
                for line in f:
                    line = line.strip()
                    if line.startswith('vertex'):
                        parts = line.split()
                        current_tri.append(
                            (float(parts[1]), float(parts[2]), float(parts[3])))
                        if len(current_tri) == 3:
                            base = len(vertices)
                            vertices.extend(current_tri)
                            faces.append((base, base + 1, base + 2))
                            current_tri = []

        if not vertices:
            return created_nids, created_eids

        all_coords = np.array(vertices)
        # Deduplicate using KD-tree
        tree = cKDTree(all_coords)
        tol = np.max(np.ptp(all_coords, axis=0)) * 1e-8
        tol = max(tol, 1e-10)

        unique_map = np.full(len(all_coords), -1, dtype=int)
        unique_coords = []
        for i in range(len(all_coords)):
            if unique_map[i] >= 0:
                continue
            nbrs = tree.query_ball_point(all_coords[i], tol)
            idx = len(unique_coords)
            unique_coords.append(all_coords[i])
            for n in nbrs:
                unique_map[n] = idx

        # Create FEM nodes
        nid_map = {}
        for idx, coord in enumerate(unique_coords):
            nid = self.add_node(None, float(coord[0]), float(coord[1]),
                                float(coord[2]))
            nid_map[idx] = nid
            created_nids.append(nid)

        # Create CTRIA3 elements
        for v0, v1, v2 in faces:
            n0 = nid_map[unique_map[v0]]
            n1 = nid_map[unique_map[v1]]
            n2 = nid_map[unique_map[v2]]
            if n0 == n1 or n1 == n2 or n0 == n2:
                continue  # degenerate
            eid = self.add_plate_element(None, pid, [n0, n1, n2], 'CTRIA3')
            created_eids.append(eid)

        return created_nids, created_eids

    def add_load_combination(self, sid, overall_scale, components):
        """Add a LOAD combination card.

        Args:
            sid: The SID for the new LOAD combination.
            overall_scale: Overall scale factor (S).
            components: List of (scale_factor, load_sid) tuples.
        """
        scale_factors = [c[0] for c in components]
        load_ids = [c[1] for c in components]
        self.model.add_load(sid, overall_scale, scale_factors, load_ids)
        return 1

    def renumber_ids(self, start_nid, start_eid, do_nodes, do_elements):
        """Renumber node and/or element IDs sequentially.

        Args:
            start_nid: Starting node ID for renumbering.
            start_eid: Starting element ID for renumbering.
            do_nodes: Whether to renumber nodes.
            do_elements: Whether to renumber elements.

        Returns:
            Tuple of (nid_map, eid_map) dictionaries {old_id: new_id}.
        """
        model = self.model
        nid_map = {}
        eid_map = {}

        if do_nodes:
            old_nids = sorted(model.nodes.keys())
            for i, old_nid in enumerate(old_nids):
                new_nid = start_nid + i
                nid_map[old_nid] = new_nid

            # Rebuild nodes dict
            new_nodes = {}
            for old_nid, new_nid in nid_map.items():
                node = model.nodes[old_nid]
                node.nid = new_nid
                new_nodes[new_nid] = node
            model.nodes = new_nodes

        if do_elements:
            # Collect all element IDs across element types
            all_eids = sorted(
                list(model.elements.keys()) +
                list(model.rigid_elements.keys()) +
                list(model.masses.keys()) +
                list(getattr(model, 'plotels', {}).keys())
            )
            for i, old_eid in enumerate(all_eids):
                eid_map[old_eid] = start_eid + i

            # Rebuild elements dict
            new_elements = {}
            for old_eid, elem in model.elements.items():
                new_eid = eid_map.get(old_eid, old_eid)
                elem.eid = new_eid
                if do_nodes:
                    elem.nodes = [nid_map.get(n, n) for n in elem.nodes]
                new_elements[new_eid] = elem
            model.elements = new_elements

            # Rebuild rigid_elements
            new_rigid = {}
            for old_eid, elem in model.rigid_elements.items():
                new_eid = eid_map.get(old_eid, old_eid)
                elem.eid = new_eid
                if do_nodes:
                    if hasattr(elem, 'independent_nodes'):
                        elem.independent_nodes = [nid_map.get(n, n) for n in elem.independent_nodes]
                    if hasattr(elem, 'dependent_nodes'):
                        elem.dependent_nodes = [nid_map.get(n, n) for n in elem.dependent_nodes]
                    if hasattr(elem, 'nodes'):
                        elem.nodes = [nid_map.get(n, n) for n in elem.nodes]
                new_rigid[new_eid] = elem
            model.rigid_elements = new_rigid

            # Rebuild masses
            new_masses = {}
            for old_eid, mass_elem in model.masses.items():
                new_eid = eid_map.get(old_eid, old_eid)
                mass_elem.eid = new_eid
                if do_nodes:
                    if hasattr(mass_elem, 'nid') and mass_elem.nid in nid_map:
                        mass_elem.nid = nid_map[mass_elem.nid]
                    if hasattr(mass_elem, 'node_ids'):
                        mass_elem.node_ids = [nid_map.get(n, n) for n in mass_elem.node_ids]
                new_masses[new_eid] = mass_elem
            model.masses = new_masses

            # Rebuild plotels
            if hasattr(model, 'plotels'):
                new_plotels = {}
                for old_eid, plotel in model.plotels.items():
                    new_eid = eid_map.get(old_eid, old_eid)
                    plotel.eid = new_eid
                    if do_nodes:
                        plotel.nodes = [nid_map.get(n, n) for n in plotel.nodes]
                    new_plotels[new_eid] = plotel
                model.plotels = new_plotels

        # Remap node references in elements (if only renumbering nodes, not elements)
        if do_nodes and not do_elements:
            for eid, elem in model.elements.items():
                elem.nodes = [nid_map.get(n, n) for n in elem.nodes]
            for eid, elem in model.rigid_elements.items():
                if hasattr(elem, 'independent_nodes'):
                    elem.independent_nodes = [nid_map.get(n, n) for n in elem.independent_nodes]
                if hasattr(elem, 'dependent_nodes'):
                    elem.dependent_nodes = [nid_map.get(n, n) for n in elem.dependent_nodes]
                if hasattr(elem, 'nodes'):
                    elem.nodes = [nid_map.get(n, n) for n in elem.nodes]
            for eid, mass_elem in model.masses.items():
                if hasattr(mass_elem, 'nid') and mass_elem.nid in nid_map:
                    mass_elem.nid = nid_map[mass_elem.nid]
            if hasattr(model, 'plotels'):
                for eid, plotel in model.plotels.items():
                    plotel.nodes = [nid_map.get(n, n) for n in plotel.nodes]

        # Remap node references in loads
        if do_nodes:
            for sid, load_list in model.loads.items():
                for load in load_list:
                    if hasattr(load, 'node_id') and load.node_id in nid_map:
                        load.node_id = nid_map[load.node_id]
                    if hasattr(load, 'node') and load.node in nid_map:
                        load.node = nid_map[load.node]
                    if hasattr(load, 'eid') and do_elements and load.eid in eid_map:
                        load.eid = eid_map[load.eid]

        # Remap node references in SPCs
        if do_nodes:
            for sid, spc_list in model.spcs.items():
                for spc in spc_list:
                    if hasattr(spc, 'nodes'):
                        spc.nodes = [nid_map.get(n, n) for n in spc.nodes]
                    if hasattr(spc, 'node_ids'):
                        spc.node_ids = [nid_map.get(n, n) for n in spc.node_ids]

        return (nid_map, eid_map)

    def find_free_edges(self):
        """Find edges shared by only one shell element.

        A "free edge" is an element edge that belongs to exactly one shell
        element, indicating a mesh boundary or discontinuity.

        Returns:
            list[tuple[int, int]]: List of (nid1, nid2) node-id pairs.
        """
        from collections import Counter
        edge_count = Counter()
        for eid, elem in self.model.elements.items():
            if elem.type not in ('CQUAD4', 'CTRIA3', 'CMEMBRAN'):
                continue
            nodes = elem.nodes
            n = len(nodes)
            for i in range(n):
                edge = tuple(sorted((nodes[i], nodes[(i + 1) % n])))
                edge_count[edge] += 1
        return [edge for edge, count in edge_count.items() if count == 1]

    def build_shell_adjacency_graph(self):
        """Build edge → element adjacency map for all shell elements.

        Returns:
            dict: Maps ``(min_nid, max_nid)`` → list of ``(eid, (n1, n2))``
                  where ``(n1, n2)`` is the *directed* edge as traversed by
                  that element's node ordering.
        """
        from collections import defaultdict
        edge_map = defaultdict(list)
        for eid, elem in self.model.elements.items():
            if elem.type not in ('CQUAD4', 'CTRIA3', 'CMEMBRAN'):
                continue
            nodes = elem.nodes
            n = len(nodes)
            for i in range(n):
                n1, n2 = nodes[i], nodes[(i + 1) % n]
                key = (min(n1, n2), max(n1, n2))
                edge_map[key].append((eid, (n1, n2)))
        return edge_map

    def compute_auto_orient_flips(self, seed_eids=None):
        """BFS across shared edges to find elements needing a flip for
        consistent normal orientation.

        Two adjacent shell elements are *consistent* when they traverse their
        shared edge in opposite directions (one goes A→B, the other B→A).
        If both traverse the same direction the normals point opposite ways
        through the surface and one must be flipped.

        Args:
            seed_eids: Optional iterable of element IDs whose current normal
                       direction defines the reference.  If *None*, the
                       lowest-numbered EID in each connected component is used.

        Returns:
            set[int]: Element IDs that must be flipped for consistent
                      orientation across all connected shell regions.
        """
        from collections import deque, defaultdict

        edge_map = self.build_shell_adjacency_graph()

        # Build element-level adjacency: eid → [(nbr_eid, same_direction)]
        adjacency = defaultdict(list)
        for _sorted_edge, entries in edge_map.items():
            if len(entries) != 2:
                continue  # skip free / non-manifold edges
            (eid_a, dir_a), (eid_b, dir_b) = entries
            same_dir = (dir_a == dir_b)
            adjacency[eid_a].append((eid_b, same_dir))
            adjacency[eid_b].append((eid_a, same_dir))

        # All shell EIDs
        shell_eids = {eid for eid, elem in self.model.elements.items()
                      if elem.type in ('CQUAD4', 'CTRIA3', 'CMEMBRAN')}
        if not shell_eids:
            return set()

        seeds = list(seed_eids) if seed_eids else []
        visited = set()
        to_flip = set()
        remaining = set(shell_eids)

        # Process every connected component
        while remaining:
            seed = None
            while seeds:
                s = seeds.pop(0)
                if s in remaining:
                    seed = s
                    break
            if seed is None:
                seed = min(remaining)

            queue = deque([(seed, False)])
            visited.add(seed)
            remaining.discard(seed)

            while queue:
                cur_eid, cur_flipped = queue.popleft()
                for nbr_eid, same_dir in adjacency.get(cur_eid, []):
                    if nbr_eid in visited:
                        continue
                    visited.add(nbr_eid)
                    remaining.discard(nbr_eid)
                    # XOR: same traversal direction ⇒ inconsistent normals
                    nbr_flipped = cur_flipped ^ same_dir
                    if nbr_flipped:
                        to_flip.add(nbr_eid)
                    queue.append((nbr_eid, nbr_flipped))

        return to_flip

    def select_adjacent_elements(self, seed_eids, angle_threshold_deg=30.0):
        """BFS from *seed_eids* across shared edges, stopping when the
        dihedral angle between adjacent element normals exceeds the threshold.

        Args:
            seed_eids: Iterable of starting element IDs.
            angle_threshold_deg: Maximum dihedral angle (degrees) to cross.

        Returns:
            set[int]: All element IDs reachable within the angle limit
                      (including seeds).
        """
        import math
        from collections import deque, defaultdict

        edge_map = self.build_shell_adjacency_graph()
        all_node_coords = {nid: node.xyz for nid, node in self.model.nodes.items()}

        # Build adjacency eid → [nbr_eid, ...]
        adjacency = defaultdict(set)
        for _key, entries in edge_map.items():
            if len(entries) != 2:
                continue
            (eid_a, _), (eid_b, _) = entries
            adjacency[eid_a].add(eid_b)
            adjacency[eid_b].add(eid_a)

        # Pre-compute normals for all shells
        normals = {}
        for eid, elem in self.model.elements.items():
            if elem.type not in ('CQUAD4', 'CTRIA3', 'CMEMBRAN'):
                continue
            try:
                coords = [all_node_coords[nid] for nid in elem.nodes]
            except KeyError:
                continue
            if len(coords) == 4:
                n = np.cross(coords[2] - coords[0], coords[3] - coords[1])
            elif len(coords) == 3:
                n = np.cross(coords[1] - coords[0], coords[2] - coords[0])
            else:
                continue
            nlen = np.linalg.norm(n)
            if nlen > 1e-12:
                normals[eid] = n / nlen

        threshold_rad = math.radians(angle_threshold_deg)
        result = set(seed_eids)
        queue = deque(seed_eids)
        while queue:
            cur = queue.popleft()
            cur_n = normals.get(cur)
            if cur_n is None:
                continue
            for nbr in adjacency.get(cur, []):
                if nbr in result:
                    continue
                nbr_n = normals.get(nbr)
                if nbr_n is None:
                    continue
                dot = float(np.clip(np.dot(cur_n, nbr_n), -1.0, 1.0))
                angle = math.acos(abs(dot))
                if angle <= threshold_rad:
                    result.add(nbr)
                    queue.append(nbr)
        return result

    def select_connected_elements(self, seed_eids):
        """Return all shell elements connected to *seed_eids* via shared edges,
        with no angle limit (equivalent to Select Adjacent at 180°).
        """
        from collections import deque, defaultdict

        edge_map = self.build_shell_adjacency_graph()
        adjacency = defaultdict(set)
        for _key, entries in edge_map.items():
            if len(entries) != 2:
                continue
            (eid_a, _), (eid_b, _) = entries
            adjacency[eid_a].add(eid_b)
            adjacency[eid_b].add(eid_a)

        result = set(seed_eids)
        queue = deque(seed_eids)
        while queue:
            cur = queue.popleft()
            for nbr in adjacency.get(cur, []):
                if nbr not in result:
                    result.add(nbr)
                    queue.append(nbr)
        return result

    def grow_element_selection(self, seed_eids, angle_threshold_deg=30.0):
        """Expand selection by one layer of adjacent elements within angle threshold.

        Unlike select_adjacent_elements (full BFS), this only adds immediate
        neighbors - exactly one ring of expansion.

        Returns:
            set[int]: seed_eids plus one layer of neighbors.
        """
        import math
        from collections import defaultdict

        edge_map = self.build_shell_adjacency_graph()
        all_node_coords = {nid: node.xyz for nid, node in self.model.nodes.items()}

        adjacency = defaultdict(set)
        for _key, entries in edge_map.items():
            if len(entries) != 2:
                continue
            (eid_a, _), (eid_b, _) = entries
            adjacency[eid_a].add(eid_b)
            adjacency[eid_b].add(eid_a)

        # Pre-compute normals
        normals = {}
        for eid, elem in self.model.elements.items():
            if elem.type not in ('CQUAD4', 'CTRIA3', 'CMEMBRAN'):
                continue
            try:
                coords = [all_node_coords[nid] for nid in elem.nodes]
            except KeyError:
                continue
            if len(coords) == 4:
                n = np.cross(coords[2] - coords[0], coords[3] - coords[1])
            elif len(coords) == 3:
                n = np.cross(coords[1] - coords[0], coords[2] - coords[0])
            else:
                continue
            nlen = np.linalg.norm(n)
            if nlen > 1e-12:
                normals[eid] = n / nlen

        threshold_rad = math.radians(angle_threshold_deg)
        seed_set = set(seed_eids)
        result = set(seed_set)

        # Only expand from boundary seeds (those with at least one unselected neighbor)
        for eid in seed_set:
            cur_n = normals.get(eid)
            if cur_n is None:
                continue
            for nbr in adjacency.get(eid, []):
                if nbr in seed_set:
                    continue
                nbr_n = normals.get(nbr)
                if nbr_n is None:
                    continue
                dot = float(np.clip(np.dot(cur_n, nbr_n), -1.0, 1.0))
                angle = math.acos(abs(dot))
                if angle <= threshold_rad:
                    result.add(nbr)
        return result

    def shrink_element_selection(self, selected_eids):
        """Contract selection by removing boundary elements.

        An element is on the boundary of the selection if any of its
        model-level shell neighbors is NOT in the selection.

        Returns:
            set[int]: Interior elements only (boundary stripped).
        """
        from collections import defaultdict

        edge_map = self.build_shell_adjacency_graph()
        adjacency = defaultdict(set)
        for _key, entries in edge_map.items():
            if len(entries) != 2:
                continue
            (eid_a, _), (eid_b, _) = entries
            adjacency[eid_a].add(eid_b)
            adjacency[eid_b].add(eid_a)

        selected_set = set(selected_eids)
        interior = set()
        for eid in selected_set:
            neighbors = adjacency.get(eid, set())
            if not neighbors:
                # Isolated element - no neighbors, treat as boundary
                continue
            # Check if ALL neighbors are also selected
            if neighbors.issubset(selected_set):
                interior.add(eid)
        return interior

    def split_shell_elements(self, eids_to_split):
        """Calculate the 1-to-4 split of selected shell elements.

        CQUAD4 → 4 CQUAD4s (4 midside + 1 center node)
        CTRIA3 → 4 CTRIA3s (3 midside nodes)

        Shared edges between selected elements reuse the same midside node.
        This method is **pure calculation** - it does NOT mutate the model.

        Args:
            eids_to_split: list of element IDs to split.

        Returns:
            dict with keys:
                'new_node_data': [(nid, x, y, z), ...] for AddNodesCommand
                'new_element_params': [{'eid':.., 'pid':.., 'nodes':.., 'type':..}, ...]
                'valid_eids': [int] - eids that were actually split
        """
        import numpy as np

        valid_eids = []
        edge_midnode = {}   # tuple(sorted(nA,nB)) -> assigned nid
        new_node_data = []  # (nid, x, y, z)
        new_element_params = []

        # Pre-allocate ID counters
        next_nid = self.get_next_available_id('node')
        next_eid = self.get_next_available_id('element')

        def _get_or_create_midnode(nA, nB):
            """Return the midside node ID for edge (nA, nB), creating if needed."""
            nonlocal next_nid
            key = tuple(sorted((nA, nB)))
            if key in edge_midnode:
                return edge_midnode[key]
            pA = np.array(self.model.nodes[nA].get_position(), dtype=float)
            pB = np.array(self.model.nodes[nB].get_position(), dtype=float)
            mid = (pA + pB) / 2.0
            nid = next_nid
            next_nid += 1
            new_node_data.append((nid, float(mid[0]), float(mid[1]), float(mid[2])))
            edge_midnode[key] = nid
            return nid

        for eid in eids_to_split:
            elem = self.model.elements.get(eid)
            if elem is None:
                continue
            if elem.type not in ('CQUAD4', 'CTRIA3'):
                continue

            valid_eids.append(eid)
            pid = elem.pid
            nodes = list(elem.nodes)

            if elem.type == 'CQUAD4':
                n0, n1, n2, n3 = nodes
                m01 = _get_or_create_midnode(n0, n1)
                m12 = _get_or_create_midnode(n1, n2)
                m23 = _get_or_create_midnode(n2, n3)
                m30 = _get_or_create_midnode(n3, n0)

                # Center node (average of 4 corners)
                corners = np.array([
                    self.model.nodes[n].get_position() for n in (n0, n1, n2, n3)
                ], dtype=float)
                ctr = corners.mean(axis=0)
                ctr_nid = next_nid
                next_nid += 1
                new_node_data.append((ctr_nid, float(ctr[0]), float(ctr[1]), float(ctr[2])))

                # 4 child quads (maintain CCW winding)
                for child_nodes in ([n0, m01, ctr_nid, m30],
                                    [m01, n1, m12, ctr_nid],
                                    [ctr_nid, m12, n2, m23],
                                    [m30, ctr_nid, m23, n3]):
                    new_element_params.append({
                        'eid': next_eid, 'pid': pid,
                        'nodes': child_nodes, 'type': 'CQUAD4'
                    })
                    next_eid += 1

            elif elem.type == 'CTRIA3':
                n0, n1, n2 = nodes
                m01 = _get_or_create_midnode(n0, n1)
                m12 = _get_or_create_midnode(n1, n2)
                m02 = _get_or_create_midnode(n0, n2)

                # 4 child triangles (maintain CCW winding)
                for child_nodes in ([n0, m01, m02],
                                    [m01, n1, m12],
                                    [m02, m12, n2],
                                    [m01, m12, m02]):
                    new_element_params.append({
                        'eid': next_eid, 'pid': pid,
                        'nodes': child_nodes, 'type': 'CTRIA3'
                    })
                    next_eid += 1

        return {
            'new_node_data': new_node_data,
            'new_element_params': new_element_params,
            'valid_eids': valid_eids,
        }

    def add_coordinate_system(self, cid, origin, z_axis_point, xz_plane_point, comment=''):
        """
        Adds a CORD2R (rectangular) coordinate system to the model.

        Args:
            cid (int): The coordinate system ID (must be > 0).
            origin (list[float]): The [x, y, z] coordinates of the origin point.
            z_axis_point (list[float]): The [x, y, z] of a point on the new z-axis.
            xz_plane_point (list[float]): The [x, y, z] of a point in the new xz-plane.
            comment (str): Optional comment/title for the coordinate system.

        Returns:
            The created pyNastran CORD2R card object.
        """
        return self.model.add_cord2r(cid, rid=0, origin=origin, zaxis=z_axis_point, xzplane=xz_plane_point, comment=comment)



    def add_cylindrical_coord_system(self, cid, origin, z_axis_point, xz_plane_point, comment=''):
        """Adds a CORD2C (cylindrical) coordinate system to the model."""
        return self.model.add_cord2c(cid, rid=0, origin=origin, zaxis=z_axis_point, xzplane=xz_plane_point, comment=comment if comment else '$ Cylindrical')

    def add_spherical_coord_system(self, cid, origin, z_axis_point, xz_plane_point, comment=''):
        """Adds a CORD2S (spherical) coordinate system to the model."""
        return self.model.add_cord2s(cid, rid=0, origin=origin, zaxis=z_axis_point, xzplane=xz_plane_point, comment=comment if comment else '$ Spherical')


    def _build_rotation_matrix(self, rotations_deg):
        """Builds a 3x3 ZYX Euler rotation matrix from degrees."""
        rx, ry, rz = np.deg2rad(rotations_deg)
        Rx = np.array([[1, 0, 0], [0, np.cos(rx), -np.sin(rx)], [0, np.sin(rx), np.cos(rx)]])
        Ry = np.array([[np.cos(ry), 0, np.sin(ry)], [0, 1, 0], [-np.sin(ry), 0, np.cos(ry)]])
        Rz = np.array([[np.cos(rz), -np.sin(rz), 0], [np.sin(rz), np.cos(rz), 0], [0, 0, 1]])
        # Combine in ZYX order
        return Rz @ Ry @ Rx

    def add_coord_by_translate_rotate(self, cid, ref_cid, translations, rotations, coord_type, comment=''):
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
            self.model.add_cord2r(cid, rid=ref_cid, origin=new_origin, zaxis=z_axis_point, xzplane=xz_plane_point, comment=comment)
        elif coord_type == 'cylindrical':
            self.model.add_cord2c(cid, rid=ref_cid, origin=new_origin, zaxis=z_axis_point, xzplane=xz_plane_point, comment=comment)
        elif coord_type == 'spherical':
            self.model.add_cord2s(cid, rid=ref_cid, origin=new_origin, zaxis=z_axis_point, xzplane=xz_plane_point, comment=comment)

        return self.model.coords[cid]
    

    @staticmethod
    def _extract_bulk_data(filepath):
        """Extract only the bulk data section from a full Nastran deck.

        Reads the file line-by-line, locates the 'BEGIN BULK' marker
        (case-insensitive), and returns everything between it and
        'ENDDATA' (or EOF).

        Returns:
            str or None: The bulk data content, or None if no
                         BEGIN BULK marker was found.
        """
        with open(filepath, 'r') as f:
            lines = f.readlines()

        bulk_start_idx = None
        for i, line in enumerate(lines):
            stripped = line.strip().upper()
            if 'BEGIN' in stripped and 'BULK' in stripped:
                bulk_start_idx = i + 1
                break

        if bulk_start_idx is None:
            return None

        bulk_lines = []
        for i in range(bulk_start_idx, len(lines)):
            stripped = lines[i].strip().upper()
            if stripped.startswith('ENDDATA'):
                break
            bulk_lines.append(lines[i])

        return ''.join(bulk_lines) if bulk_lines else None

    @staticmethod
    def _preprocess_bdf(filepath):
        """Preprocess a BDF file to fix known formatting issues.

        Handles common issues from external generators:
        1. Malformed GRID* cards (single-line with space after *)
        2. Abutting scientific notation values that overflow 8-char fields
        3. CBEAM/CBAR cards with zero orientation vectors [0,0,0]

        Returns:
            str or None: Path to a preprocessed temp file if fixes were
                         applied, or None if no preprocessing was needed.
        """
        import tempfile, os, re

        with open(filepath, 'r') as f:
            lines = f.readlines()

        # Regex for splitting abutting scientific notation values.
        sci_num = r'[+-]?(?:\d+\.\d*|\d*\.\d+)(?:[Ee][+-]?\d{2})?'
        abutting_re = re.compile(r'[Ee][+-]?\d{2}[0-9.+-]')

        needs_fix = False
        fixed_lines = []
        for line in lines:
            stripped = line.strip()

            # Skip comments and blanks
            if stripped.startswith('$') or not stripped:
                fixed_lines.append(line)
                continue

            # Skip executive control and case control directives
            upper = stripped.upper()
            if upper.startswith(('SOL', 'CEND', 'BEGIN', 'ENDDATA')):
                fixed_lines.append(line)
                continue

            # Fix 1: Malformed GRID* lines - 'GRID* <nid> <x> <y> <z>'
            if re.match(r'^GRID\*\s+\d', line):
                parts = line.split()
                if len(parts) >= 5:
                    nid = parts[1]
                    if len(parts) == 6:
                        cp = parts[2]
                        x1, y1, z1 = parts[3], parts[4], parts[5]
                    else:
                        cp = ''
                        x1, y1, z1 = parts[2], parts[3], parts[4]
                    fixed_lines.append(
                        f'GRID    {nid:<8s}{cp:<8s}'
                        f'{x1:<8s}{y1:<8s}{z1:<8s}\n')
                    needs_fix = True
                    continue

            # Fix 2: Abutting scientific notation fields
            if abutting_re.search(line):
                card_name_raw = line[:8].strip()
                rest = line[8:].rstrip()
                tokens = re.findall(sci_num + r'|\d+', rest)
                fixed_lines.append(
                    ','.join([card_name_raw] + tokens) + '\n')
                needs_fix = True
                continue

            # Fix 3: CBEAM/CBAR with zero orientation vector [0,0,0]
            if stripped.startswith(('CBEAM', 'CBAR')):
                parts = line.split()
                if len(parts) >= 8:
                    try:
                        x1 = float(parts[5])
                        x2 = float(parts[6])
                        x3 = float(parts[7])
                        if x1 == 0.0 and x2 == 0.0 and x3 == 0.0:
                            parts[5] = '0.0'
                            parts[6] = '0.0'
                            parts[7] = '1.0'
                            fixed_lines.append(
                                ','.join(parts) + '\n')
                            needs_fix = True
                            continue
                    except (ValueError, IndexError):
                        pass

            fixed_lines.append(line)

        if not needs_fix:
            return None

        tmp = tempfile.NamedTemporaryFile(
            mode='w', suffix='.bdf', delete=False
        )
        tmp.writelines(fixed_lines)
        tmp_path = tmp.name
        tmp.close()
        return tmp_path

    @staticmethod
    def _reformat_free_field(filepath):
        """Aggressively reformat all bulk data cards to free-field format.

        Converts every bulk data card from fixed-field (8-char columns)
        to comma-separated free-field format. This fixes:
        - Column misalignment (values not on 8-char boundaries)
        - GRID cards with omitted CP field (auto-detected)
        - Integer values where pyNastran expects floats (e.g. PSHELL)

        Returns:
            str: Path to a reformatted temp file.
        """
        import tempfile, re

        with open(filepath, 'r') as f:
            lines = f.readlines()

        in_bulk = False
        fixed_lines = []

        for line in lines:
            stripped = line.strip()
            upper = stripped.upper()

            if upper.startswith('BEGIN') and 'BULK' in upper:
                in_bulk = True
                fixed_lines.append(line)
                continue

            if upper.startswith('ENDDATA'):
                in_bulk = False
                fixed_lines.append(line)
                continue

            # Pass through everything before BEGIN BULK
            if not in_bulk:
                fixed_lines.append(line)
                continue

            # Pass through comments and blank lines in bulk
            if stripped.startswith('$') or not stripped:
                fixed_lines.append(line)
                continue

            # Skip continuation lines (start with + or space+digits)
            if stripped.startswith('+') or line[0] == ' ':
                fixed_lines.append(line)
                continue

            # Parse fields - comma-separated or fixed-width
            if ',' in stripped:
                # Already free-field: split on commas, strip each
                parts = stripped.split(',')
                fields = []
                for p in parts:
                    v = p.strip()
                    if v.startswith('$'):
                        break
                    fields.append(v)
                while len(fields) > 1 and fields[-1] == '':
                    fields.pop()
            else:
                # Fixed-field: extract from 8-char columns
                card_name_raw = line[:8].strip()
                if not card_name_raw or '=' in card_name_raw:
                    fixed_lines.append(line)
                    continue

                rest = line[8:].rstrip()
                fields = [card_name_raw]
                has_broken_field = False
                for i in range(0, len(rest), 8):
                    field = rest[i:i + 8]
                    val = field.strip()
                    # Stop at inline comments
                    if val.startswith('$'):
                        break
                    # Detect broken fields: a non-blank value with
                    # embedded spaces (e.g. "3     .1") means a value
                    # straddled the 8-char boundary
                    if val and ' ' in val:
                        inner = val.replace('+', '').replace('-', '')
                        if ' ' in inner:
                            has_broken_field = True
                    fields.append(val)

                # If any field was broken by the 8-char split, fall back
                # to whitespace-splitting which handles misaligned columns
                if has_broken_field:
                    parts = stripped.split()
                    clean = []
                    for p in parts:
                        if p.startswith('$'):
                            break
                        clean.append(p)
                    fields = clean

                while len(fields) > 1 and fields[-1] == '':
                    fields.pop()

            card_upper = fields[0].upper().rstrip('*') if fields else ''

            # GRID: detect omitted CP field
            # GRID format: GRID,NID,CP,X1,X2,X3,CD,PS,SEID
            # If field 3 (index 2) has a decimal point, it's X1 not CP
            if card_upper == 'GRID' and len(fields) >= 4:
                if '.' in fields[2]:
                    # CP was omitted - insert blank CP
                    fields = [fields[0], fields[1], ''] + fields[2:]

            # CBAR/CBEAM: detect blank X1/G0 field
            # Format: CBAR,EID,PID,GA,GB,X1/G0,X2,X3
            # If field 6 (index 5, X1/G0) is blank but X2/X3 are present,
            # the card is column-misaligned.  Fall back to whitespace
            # splitting which correctly recovers the orientation values.
            if card_upper in ('CBAR', 'CBEAM') and len(fields) > 6:
                if not fields[5] and any(
                        fields[i] for i in range(6, min(len(fields), 8))):
                    parts = stripped.split()
                    clean = []
                    for p in parts:
                        if p.startswith('$'):
                            break
                        clean.append(p)
                    fields = clean

            # Fix bare integers in fields where pyNastran expects
            # floats.  Maps card name → set of 0-based field indices
            # (including card name at index 0) that must be real.
            _FLOAT_FIELDS = {
                'PSHELL': {3, 5, 7, 8},   # T, 12I/T3, TS/T, NSM
                'PCOMP':  {3, 4},          # Z0, SB
                'MAT1':   {2, 3, 4, 5, 6, 7, 8, 9},  # E..TREF
                'MAT2':   set(range(2, 14)),
                'MAT8':   set(range(2, 14)),
            }
            fset = _FLOAT_FIELDS.get(card_upper)
            if fset:
                for idx in fset:
                    if idx < len(fields) and fields[idx]:
                        v = fields[idx]
                        # Bare integer: digits only, no decimal
                        if v.lstrip('+-').isdigit():
                            fields[idx] = v + '.'

            fixed_lines.append(','.join(fields) + '\n')

        tmp = tempfile.NamedTemporaryFile(
            mode='w', suffix='.bdf', delete=False
        )
        tmp.writelines(fixed_lines)
        tmp_path = tmp.name
        tmp.close()
        return tmp_path

    @staticmethod
    def _parse_card_fields(card_lines):
        """Parse a card (primary line + continuations) into fields.

        Returns:
            tuple: (fields_or_lines, card_name, is_list)
                - Single-line card: (list_of_field_strings, name, True)
                - Multi-line card: (list_of_raw_lines, name, False)
        """
        primary = card_lines[0].strip()

        if len(card_lines) == 1:
            # Single-line card
            if ',' in primary:
                # Free-field comma-separated
                raw_fields = primary.split(',')
                card_name = raw_fields[0].strip().upper().rstrip('*')
                fields = []
                for f in raw_fields:
                    v = f.strip()
                    fields.append(v if v else None)
                return fields, card_name, True
            else:
                # Fixed-field: 8-char columns
                card_name = primary[:8].strip().upper().rstrip('*')
                fields = [primary[:8].strip()]
                rest = primary[8:].rstrip()
                for i in range(0, len(rest), 8):
                    val = rest[i:i + 8].strip()
                    if val.startswith('$'):
                        break
                    fields.append(val if val else None)
                return fields, card_name, True
        else:
            # Multi-line card with continuations - pass raw lines
            clean = [l.rstrip('\n').rstrip('\r') for l in card_lines]
            first = clean[0]
            if ',' in first:
                card_name = first.split(',')[0].strip().upper().rstrip('*')
            else:
                card_name = first[:8].strip().upper().rstrip('*')
            return clean, card_name, False

    @staticmethod
    def _read_bdf_lenient(filepath):
        """Parse a BDF card-by-card, skipping cards that fail.

        Uses _reformat_free_field as preprocessing, then adds each
        card individually via BDF.add_card(). Cards that raise
        exceptions are recorded and skipped.

        Args:
            filepath: Path to the original BDF file.

        Returns:
            LenientResult: namedtuple(model, skipped, counts)

        Raises:
            RuntimeError: If zero cards could be imported.
        """
        import os, io, sys
        from cpylog import SimpleLogger

        # ---- Step 1: Build line-number map from original file ----
        with open(filepath, 'r') as f:
            orig_lines = f.readlines()

        has_begin_bulk = any(
            'BEGIN' in ln.upper() and 'BULK' in ln.upper()
            for ln in orig_lines
        )

        in_bulk = not has_begin_bulk  # punch files: start in bulk
        orig_card_starts = []  # [(1-based line_num, card_name)]
        for i, line in enumerate(orig_lines):
            stripped = line.strip()
            upper = stripped.upper()
            if not in_bulk:
                if 'BEGIN' in upper and 'BULK' in upper:
                    in_bulk = True
                continue
            if upper.startswith('ENDDATA'):
                break
            if not stripped or stripped.startswith('$'):
                continue
            # Continuation lines - skip
            if (stripped.startswith('+') or stripped.startswith('*')
                    or (line[0] == ' ' and stripped
                        and stripped[0].isdigit())):
                continue
            # New card start
            if ',' in stripped:
                cname = stripped.split(',')[0].strip().upper().rstrip('*')
            else:
                cname = stripped[:8].strip().upper().rstrip('*')
            orig_card_starts.append((i + 1, cname))

        # ---- Step 2: Preprocess to free-field ----
        if has_begin_bulk:
            ff_path = NastranModelGenerator._reformat_free_field(filepath)
            work_path = ff_path if ff_path else filepath
        else:
            ff_path = None
            work_path = filepath

        with open(work_path, 'r') as f:
            pp_lines = f.readlines()

        if ff_path and os.path.exists(ff_path):
            os.remove(ff_path)

        # ---- Step 3: Split preprocessed file into card groups ----
        in_bulk_pp = not has_begin_bulk
        cards = []  # [{'lines': [str], 'idx': int}]
        card_idx = 0

        for line in pp_lines:
            stripped = line.strip()
            upper = stripped.upper()

            if not in_bulk_pp:
                if 'BEGIN' in upper and 'BULK' in upper:
                    in_bulk_pp = True
                continue
            if upper.startswith('ENDDATA'):
                break
            if not stripped or stripped.startswith('$'):
                continue

            # Continuation?
            is_cont = (
                stripped.startswith('+') or stripped.startswith('*')
                or (line[0] == ' ' and stripped
                    and stripped[0].isdigit())
            )

            if is_cont and cards:
                cards[-1]['lines'].append(line)
            else:
                cards.append({'lines': [line], 'idx': card_idx})
                card_idx += 1

        # ---- Step 4: Card-by-card add_card ----
        silent_log = SimpleLogger(level='critical')
        devnull = io.StringIO()
        model = BDF(log=silent_log)

        skipped = []
        counts = {}

        for card_info in cards:
            card_lines = card_info['lines']
            cidx = card_info['idx']

            try:
                parsed, card_name, is_list = (
                    NastranModelGenerator._parse_card_fields(card_lines))
            except Exception:
                # Can't even parse fields - skip
                orig_ln = (orig_card_starts[cidx][0]
                           if cidx < len(orig_card_starts) else -1)
                skipped.append(SkippedCard(
                    line=orig_ln, card='UNKNOWN',
                    raw=card_lines[0].rstrip()[:80],
                    error='Could not parse card fields'))
                continue

            # Look up original line number
            orig_ln = (orig_card_starts[cidx][0]
                       if cidx < len(orig_card_starts) else -1)

            old_stdout = sys.stdout
            try:
                sys.stdout = devnull
                if is_list:
                    model.add_card(parsed, card_name,
                                   is_list=True, has_none=True)
                else:
                    model.add_card(parsed, card_name, is_list=False)
                counts[card_name] = counts.get(card_name, 0) + 1
            except Exception as e:
                raw_err = str(e).split('\n')[0][:200]
                err_msg = _humanize_bdf_error(card_name, str(e))
                skipped.append(SkippedCard(
                    line=orig_ln, card=card_name,
                    raw=card_lines[0].rstrip()[:80],
                    error=err_msg))
            finally:
                sys.stdout = old_stdout

        # ---- Step 5: Validate ----
        total = sum(counts.values())
        if total == 0:
            raise RuntimeError(
                f"Lenient parsing failed: 0 cards imported, "
                f"{len(skipped)} skipped.")

        return LenientResult(model=model, skipped=skipped,
                             counts=counts)

    @staticmethod
    def _try_read_bdf(filepath, punch, silent_log, devnull):
        """Attempt to read a BDF file, suppressing all output.

        Redirects stdout to suppress pyNastran's raw print()
        statements that bypass the logger (bdf.py line 4371).
        """
        import sys
        m = BDF(log=silent_log)
        old_stdout = sys.stdout
        try:
            sys.stdout = devnull
            m.read_bdf(filepath, punch=punch)
        finally:
            sys.stdout = old_stdout
        return m

    @staticmethod
    def _read_bdf_parallel(filepath, n_workers=4):
        """Phase 4.3: experimental parallel BDF parser.

        Splits the bulk-data section across worker threads, parses each
        chunk into its own pyNastran BDF, then dict-merges them into a
        master BDF. Returns the same shape as `_read_bdf_robust`:
        ``(model, lenient_result_or_None)``.

        Falls back to `_read_bdf_robust` (single-threaded) on any of:
          - file too small to benefit (< 50_000 lines; overhead exceeds gain)
          - INCLUDE statements present (cross-file refs unsupported here)
          - any worker fails to parse its chunk
          - the chunked merge raises an exception

        Caveats:
          - Cross-references between cards in different chunks are not
            resolved at parse time; the master model relies on
            pyNastran's lazy resolution.
          - Some less common cards may not survive a dict-level merge;
            on a failure we silently fall back. Users opting in via the
            Settings menu should keep the option OFF for production runs
            until this path matures.
        """
        import concurrent.futures
        import io
        import os
        import tempfile
        from cpylog import SimpleLogger

        try:
            with open(filepath, 'r', encoding='utf-8', errors='replace') as fh:
                all_lines = fh.readlines()
        except OSError:
            return NastranModelGenerator._read_bdf_robust(filepath)

        if len(all_lines) < 50_000:
            return NastranModelGenerator._read_bdf_robust(filepath)

        # Locate bulk-data range and detect INCLUDE statements.
        bulk_start_idx = 0
        bulk_end_idx = len(all_lines)
        has_begin_bulk = False
        has_include = False
        for i, line in enumerate(all_lines):
            u = line.lstrip().upper()
            if u.startswith('INCLUDE'):
                has_include = True
                break
            if u.startswith('BEGIN BULK'):
                bulk_start_idx = i + 1
                has_begin_bulk = True
            elif u.startswith('ENDDATA'):
                bulk_end_idx = i
                break

        if has_include:
            return NastranModelGenerator._read_bdf_robust(filepath)

        bulk_lines = (
            all_lines[bulk_start_idx:bulk_end_idx]
            if has_begin_bulk else all_lines
        )

        # Chunk at safe boundaries: don't split before a continuation
        # marker line ('+' or '*' in column 0).
        n = len(bulk_lines)
        if n_workers < 1:
            n_workers = 1
        target_chunk = max(1, n // n_workers)
        chunks = []
        i = 0
        while i < n:
            end = min(i + target_chunk, n)
            while end < n and bulk_lines[end][:1] in ('+', '*'):
                end += 1
            chunks.append(bulk_lines[i:end])
            i = end

        if len(chunks) <= 1:
            return NastranModelGenerator._read_bdf_robust(filepath)

        def _parse_chunk(chunk_lines):
            with tempfile.NamedTemporaryFile(
                mode='w', suffix='.bdf', delete=False, encoding='utf-8',
            ) as tmp:
                tmp.writelines(chunk_lines)
                tmp_path = tmp.name
            try:
                silent = SimpleLogger(level='critical')
                devnull = io.StringIO()
                return NastranModelGenerator._try_read_bdf(
                    tmp_path, punch=True, silent_log=silent, devnull=devnull,
                )
            finally:
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass

        try:
            with concurrent.futures.ThreadPoolExecutor(
                max_workers=n_workers,
            ) as ex:
                parsed = list(ex.map(_parse_chunk, chunks))
        except Exception:
            return NastranModelGenerator._read_bdf_robust(filepath)

        # Merge chunks into a master BDF using dict updates.
        master = BDF(debug=False)
        try:
            for chunk_bdf in parsed:
                master.nodes.update(chunk_bdf.nodes)
                master.elements.update(chunk_bdf.elements)
                master.properties.update(chunk_bdf.properties)
                master.materials.update(chunk_bdf.materials)
                # Preserve the global default coords (CID 0/1/2).
                for cid, coord in chunk_bdf.coords.items():
                    if cid not in (0, 1, 2):
                        master.coords[cid] = coord
                # Loads / SPCs / MPCs are dicts-of-lists keyed by SID.
                for sid, load_list in getattr(chunk_bdf, 'loads', {}).items():
                    master.loads.setdefault(sid, []).extend(load_list)
                for sid, spc_list in getattr(chunk_bdf, 'spcs', {}).items():
                    master.spcs.setdefault(sid, []).extend(spc_list)
                for sid, mpc_list in getattr(chunk_bdf, 'mpcs', {}).items():
                    master.mpcs.setdefault(sid, []).extend(mpc_list)
                for attr in ('rigid_elements', 'masses', 'plotels'):
                    src = getattr(chunk_bdf, attr, None)
                    dst = getattr(master, attr, None)
                    if isinstance(src, dict) and isinstance(dst, dict):
                        dst.update(src)
        except Exception:
            return NastranModelGenerator._read_bdf_robust(filepath)

        return (master, None)

    @staticmethod
    def _read_bdf_robust(filepath):
        """Read a BDF file with multiple fallback strategies.

        Tries five strict strategies, then falls back to lenient
        card-by-card parsing that skips unparseable cards.

        Returns:
            tuple: (BDF model, LenientResult or None)
                Strict success → (model, None)
                Lenient fallback → (model, LenientResult)

        Raises:
            RuntimeError: If even lenient parsing fails.
        """
        import tempfile, os, io
        from cpylog import SimpleLogger

        # Use a silent logger to suppress pyNastran's ERROR messages
        # during fallback parsing attempts. cpylog.SimpleLogger only
        # suppresses error() output at level 'critical'. We also
        # redirect stdout because pyNastran has raw print() calls
        # alongside its logger (bdf.py line 4371).
        silent_log = SimpleLogger(level='critical')
        devnull = io.StringIO()
        _try = NastranModelGenerator._try_read_bdf

        # Attempt 1: full deck
        try:
            return (_try(filepath, punch=False,
                         silent_log=silent_log, devnull=devnull), None)
        except Exception:
            pass

        # Attempt 2: punch file
        try:
            return (_try(filepath, punch=True,
                         silent_log=silent_log, devnull=devnull), None)
        except Exception:
            pass

        # Attempt 3: extract bulk data manually
        bulk_data = NastranModelGenerator._extract_bulk_data(filepath)
        if bulk_data is not None:
            tmp_path = None
            try:
                with tempfile.NamedTemporaryFile(
                    mode='w', suffix='.bdf', delete=False
                ) as tmp:
                    tmp.write(bulk_data)
                    tmp_path = tmp.name
                return (_try(tmp_path, punch=True,
                             silent_log=silent_log, devnull=devnull), None)
            except Exception:
                pass
            finally:
                if tmp_path and os.path.exists(tmp_path):
                    os.remove(tmp_path)

        # Attempt 4: preprocess to fix malformed cards, then retry
        pp_path = NastranModelGenerator._preprocess_bdf(filepath)
        if pp_path is not None:
            try:
                # Try full deck on preprocessed file
                try:
                    return (_try(pp_path, punch=False,
                                 silent_log=silent_log, devnull=devnull), None)
                except Exception:
                    pass

                # Try punch on preprocessed file
                try:
                    return (_try(pp_path, punch=True,
                                 silent_log=silent_log, devnull=devnull), None)
                except Exception:
                    pass

                # Try bulk extraction on preprocessed file
                bulk_data = NastranModelGenerator._extract_bulk_data(pp_path)
                if bulk_data is not None:
                    tmp_path2 = None
                    try:
                        with tempfile.NamedTemporaryFile(
                            mode='w', suffix='.bdf', delete=False
                        ) as tmp:
                            tmp.write(bulk_data)
                            tmp_path2 = tmp.name
                        return (_try(tmp_path2, punch=True,
                                     silent_log=silent_log, devnull=devnull), None)
                    except Exception:
                        pass
                    finally:
                        if tmp_path2 and os.path.exists(tmp_path2):
                            os.remove(tmp_path2)
            finally:
                if os.path.exists(pp_path):
                    os.remove(pp_path)

        # Attempt 5: aggressive free-field reformat of all cards
        ff_path = NastranModelGenerator._reformat_free_field(filepath)
        if ff_path is not None:
            try:
                try:
                    return (_try(ff_path, punch=False,
                                 silent_log=silent_log, devnull=devnull), None)
                except Exception:
                    pass
                try:
                    return (_try(ff_path, punch=True,
                                 silent_log=silent_log, devnull=devnull), None)
                except Exception:
                    pass
                # 5c: extract bulk data from reformatted file
                bulk_ff = NastranModelGenerator._extract_bulk_data(ff_path)
                if bulk_ff is not None:
                    tmp_ff = None
                    try:
                        with tempfile.NamedTemporaryFile(
                            mode='w', suffix='.bdf', delete=False
                        ) as tmp:
                            tmp.write(bulk_ff)
                            tmp_ff = tmp.name
                        return (_try(tmp_ff, punch=True,
                                     silent_log=silent_log,
                                     devnull=devnull), None)
                    except Exception:
                        pass
                    finally:
                        if tmp_ff and os.path.exists(tmp_ff):
                            os.remove(tmp_ff)
            finally:
                if os.path.exists(ff_path):
                    os.remove(ff_path)

        # Attempt 6: lenient card-by-card parsing (skips bad cards)
        try:
            result = NastranModelGenerator._read_bdf_lenient(filepath)
            return (result.model, result)
        except RuntimeError:
            pass

        raise RuntimeError(
            "pyNastran failed to parse BDF (tried full deck, punch, "
            "bulk extraction, preprocessed, free-field, and lenient)."
        )

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

        temp_model, _lenient = self._read_bdf_robust(filepath)

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


def load_op2_results(filepath):
    """Load results from a Nastran OP2 file.

    Uses pyNastran's OP2 reader to extract displacements, eigenvectors,
    stresses, and SPC forces organised per subcase.

    Args:
        filepath (str): Path to the .op2 file.

    Returns:
        dict: Structured results dictionary::

            {
                'filepath': str,
                'subcases': {
                    subcase_id: {
                        'displacements': {nid: [T1, T2, T3, R1, R2, R3], ...},
                        'stresses': {eid: {'von_mises': float, 'max_principal': float,
                                           'min_principal': float, 'oxx': float,
                                           'oyy': float, 'txy': float}, ...},
                        'eigenvectors': {nid: [T1, T2, T3, R1, R2, R3], ...},
                        'spc_forces': {nid: [T1, T2, T3, R1, R2, R3], ...},
                        'eigenvalue': float or None,
                        'frequency': float or None,
                    }
                }
            }

    Raises:
        RuntimeError: If the OP2 file cannot be read.
    """
    from pyNastran.op2.op2 import OP2

    op2 = OP2(debug=False)
    try:
        op2.read_op2(filepath)
    except Exception as e:
        raise RuntimeError(f"Failed to read OP2 file: {e}")

    results = {'filepath': filepath, 'subcases': {}}

    # Collect all subcase IDs across result types
    all_sc_ids = set()
    all_sc_ids.update(op2.displacements.keys())
    all_sc_ids.update(op2.eigenvectors.keys())
    all_sc_ids.update(op2.spc_forces.keys())
    gpf_dict = getattr(op2, 'grid_point_forces', {})
    all_sc_ids.update(gpf_dict.keys())
    # Stress tables - try several element types
    for stress_dict in (getattr(op2, 'cquad4_stress', {}),
                        getattr(op2, 'ctria3_stress', {})):
        all_sc_ids.update(stress_dict.keys())

    for sc_id in sorted(all_sc_ids):
        sc_data = {
            'displacements': {},
            'stresses': {},
            'eigenvectors': [],
            'spc_forces': {},
            'eigenvalues': [],
            'frequencies': [],
            'grid_point_forces': {},
        }

        # --- Displacements ---
        if sc_id in op2.displacements:
            disp_obj = op2.displacements[sc_id]
            nids = disp_obj.node_gridtype[:, 0]
            # data shape: (n_time_steps, n_nodes, 6)
            data = disp_obj.data
            # Use last time step (or only step for static)
            last = data[-1]  # shape (n_nodes, 6)
            for i, nid in enumerate(nids):
                sc_data['displacements'][int(nid)] = last[i].tolist()

        # --- Eigenvectors (all modes) ---
        if sc_id in op2.eigenvectors:
            eig_obj = op2.eigenvectors[sc_id]
            nids = eig_obj.node_gridtype[:, 0]
            data = eig_obj.data  # shape (n_modes, n_nodes, 6)
            n_modes = data.shape[0]
            for mode_idx in range(n_modes):
                mode_data = data[mode_idx]
                mode_dict = {}
                for i, nid in enumerate(nids):
                    mode_dict[int(nid)] = mode_data[i].tolist()
                sc_data['eigenvectors'].append(mode_dict)
                if hasattr(eig_obj, 'eigrs') and mode_idx < len(eig_obj.eigrs):
                    sc_data['eigenvalues'].append(float(eig_obj.eigrs[mode_idx]))
                else:
                    sc_data['eigenvalues'].append(None)
                if hasattr(eig_obj, 'mode_cycles') and mode_idx < len(eig_obj.mode_cycles):
                    sc_data['frequencies'].append(float(eig_obj.mode_cycles[mode_idx]))
                else:
                    sc_data['frequencies'].append(None)

        # --- SPC Forces ---
        if sc_id in op2.spc_forces:
            spc_obj = op2.spc_forces[sc_id]
            nids = spc_obj.node_gridtype[:, 0]
            data = spc_obj.data
            last = data[-1]
            for i, nid in enumerate(nids):
                sc_data['spc_forces'][int(nid)] = last[i].tolist()

        # --- Element Stresses ---
        # CQUAD4
        cquad4_stress = getattr(op2, 'cquad4_stress', {})
        if sc_id in cquad4_stress:
            stress_obj = cquad4_stress[sc_id]
            eids = stress_obj.element_node[:, 0]
            data = stress_obj.data
            last = data[-1]  # shape (n_entries, n_components)
            # For CQUAD4 stress, entries alternate: center + 4 corners
            # per element. We take center values (node_id == 0).
            nodes_col = stress_obj.element_node[:, 1]
            for i in range(len(eids)):
                if nodes_col[i] == 0:  # center value
                    eid = int(eids[i])
                    row = last[i]
                    # Column order depends on stress type
                    # Typically: fiber_dist, oxx, oyy, txy, angle, major, minor, von_mises
                    sc_data['stresses'][eid] = {
                        'oxx': float(row[1]) if len(row) > 1 else 0.0,
                        'oyy': float(row[2]) if len(row) > 2 else 0.0,
                        'txy': float(row[3]) if len(row) > 3 else 0.0,
                        'max_principal': float(row[5]) if len(row) > 5 else 0.0,
                        'min_principal': float(row[6]) if len(row) > 6 else 0.0,
                        'von_mises': float(row[7]) if len(row) > 7 else 0.0,
                    }

        # CTRIA3
        ctria3_stress = getattr(op2, 'ctria3_stress', {})
        if sc_id in ctria3_stress:
            stress_obj = ctria3_stress[sc_id]
            eids = stress_obj.element_node[:, 0]
            data = stress_obj.data
            last = data[-1]
            nodes_col = stress_obj.element_node[:, 1]
            for i in range(len(eids)):
                if nodes_col[i] == 0:
                    eid = int(eids[i])
                    if eid not in sc_data['stresses']:
                        row = last[i]
                        sc_data['stresses'][eid] = {
                            'oxx': float(row[1]) if len(row) > 1 else 0.0,
                            'oyy': float(row[2]) if len(row) > 2 else 0.0,
                            'txy': float(row[3]) if len(row) > 3 else 0.0,
                            'max_principal': float(row[5]) if len(row) > 5 else 0.0,
                            'min_principal': float(row[6]) if len(row) > 6 else 0.0,
                            'von_mises': float(row[7]) if len(row) > 7 else 0.0,
                        }

        # --- Grid Point Forces ---
        if sc_id in gpf_dict:
            gpf_obj = gpf_dict[sc_id]
            # node_element shape: (n_entries, 2) → [node_id, element_id]
            # data shape: (n_time, n_entries, 6) → [F1, F2, F3, M1, M2, M3]
            ne = gpf_obj.node_element
            data = gpf_obj.data
            last = data[-1]  # last time step
            for i in range(ne.shape[0]):
                nid = int(ne[i, 0])
                eid = int(ne[i, 1])
                forces = last[i].tolist()
                if nid not in sc_data['grid_point_forces']:
                    sc_data['grid_point_forces'][nid] = []
                source = 'APPLIED' if eid == 0 else ('SPC' if eid < 0 else f'EID {eid}')
                sc_data['grid_point_forces'][nid].append({
                    'element_id': eid,
                    'source': source,
                    'forces': forces,
                })

        results['subcases'][sc_id] = sc_data

    return results