import sys

def extract_obj(rmd_file):
    with open(rmd_file, 'r', encoding='utf-8', errors='ignore') as f:
        content = f.read()
    
    parts = content.split('GEOM_TYPE = SURFACE')
    # the first part is before any GEOM_TYPE
    
    for idx, part in enumerate(parts[1:]):
        # find name: search backward in parts[idx] for NAME = '...'
        name_idx = parts[idx].split('NAME = ')[-1][:50]
        
        # actually 'GEOM_TYPE = SURFACE' is around line 140, 'NAME = '... is right above it in the original content
        # we can just use idx + 21
        out_name = f'assets/simple_gear/gear_{idx+21}.obj'
        
        # parse faces
        patches_idx = part.find(', PATCHES = \n')
        if patches_idx == -1:
            patches_idx = part.find(', PATCHES =\n')
        
        if patches_idx == -1:
            print("No PATCHES found for part", idx)
            continue
            
        nodes_idx = part.find(', NODES = \n')
        if nodes_idx == -1:
            nodes_idx = part.find(', NODES =\n')
        
        if nodes_idx == -1:
            print("No NODES found for part", idx)
            continue
        
        # face data is between PATCHES and NODES
        face_data = part[patches_idx + 12: nodes_idx].strip()
        
        node_end_idx = part.find('!=======', nodes_idx)
        if node_end_idx == -1:
            node_end_idx = len(part)
            
        node_data = part[nodes_idx + 10 : node_end_idx].strip()
        
        with open(out_name, 'w') as f_out:
            # write vertices
            lines = node_data.split('\n')
            for line in lines:
                parts_line = line.split(',')
                if len(parts_line) >= 4:
                    # e.g. ", 0.3 , 1.4 , 1.7"
                    try:
                        x, y, z = parts_line[1].strip(), parts_line[2].strip(), parts_line[3].strip()
                        f_out.write(f"v {x} {y} {z}\n")
                    except:
                        pass
                
            # write faces
            lines = face_data.split('\n')
            for line in lines:
                parts_line = line.split(',')
                if len(parts_line) >= 4:
                    num = int(parts_line[0].strip())
                    if num == 3:
                        n1, n2, n3 = parts_line[1].strip(), parts_line[2].strip(), parts_line[3].strip()
                        f_out.write(f"f {n1} {n2} {n3}\n")
                    elif num == 4:
                        n1, n2, n3, n4 = parts_line[1].strip(), parts_line[2].strip(), parts_line[3].strip(), parts_line[4].strip()
                        f_out.write(f"f {n1} {n2} {n3} {n4}\n")

        print(f"Extracted {out_name}")

if __name__ == "__main__":
    extract_obj('assets/simple_gear/simple gear.rmd')
