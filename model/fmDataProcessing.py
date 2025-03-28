import numpy as np

# --- Re–Named Constants ---
OFFSETS = {
    'A': 0x3D8,
    'B': 0x3D4,
    'C': 0x25C,
    'C_alt': 0x3CC,
    'D': 0x258
}

# --- Utility Functions ---
def diff_bits(a, b):
    return bin(a ^ b).count("1")

def matches_syndrome(synd, target, tol=1):
    return diff_bits(synd, target) <= tol

def syndrome_calc(bit_seq):
    poly = 0b11111001110
    reg = 0
    # Process only the first 16 bits
    for bit in bit_seq[:16]:
        reg = (reg << 1) | int(bit)
        if reg & (1 << 16):
            reg ^= poly << 6
    return reg & 0x3FF

def sync_groups(bit_stream):
    block_len = 26
    labels = ['A', 'B', 'C', 'D']
    groups_found = []
    # Try start offsets from 2 to block_len-1
    for start_offset in range(2, block_len):
        seg = bit_stream[start_offset:]
        pos = 0
        while pos + block_len * 4 <= len(seg):
            grp = {}
            valid = True
            for idx, lab in enumerate(labels):
                blk = seg[pos + idx * block_len : pos + (idx + 1) * block_len]
                if len(blk) < block_len:
                    valid = False
                    break
                synd = syndrome_calc(blk)
                # For block C, accept either possible syndrome
                if lab == 'C':
                    if not (matches_syndrome(synd, OFFSETS['C']) or matches_syndrome(synd, OFFSETS['C_alt'])):
                        valid = False
                        break
                else:
                    if not matches_syndrome(synd, OFFSETS[lab]):
                        valid = False
                        break
                grp[lab] = blk[:16]
            if valid:
                print(f"[SYNCED] Group synchronized at offset {start_offset}, position {pos}")
                groups_found.append(grp)
                return groups_found
            pos += 1
    if not groups_found:
        print("No groups synchronized.")
    return groups_found

def bits_to_hex_str(bits_list):
    mapping = {
        (0,0,0,0): "0", (0,0,0,1): "1", (0,0,1,0): "2", (0,0,1,1): "3",
        (0,1,0,0): "4", (0,1,0,1): "5", (0,1,1,0): "6", (0,1,1,1): "7",
        (1,0,0,0): "8", (1,0,0,1): "9", (1,0,1,0): "A", (1,0,1,1): "B",
        (1,1,0,0): "C", (1,1,0,1): "D", (1,1,1,0): "E", (1,1,1,1): "F"
    }
    hex_str = ""
    for i in range(0, len(bits_list), 4):
        nibble = tuple(bits_list[i:i+4])
        hex_str += mapping.get(nibble, "?")
    return hex_str

def extract_pty(bits_B):
    try:
        return int("".join(map(str, bits_B[6:11])), 2)
    except Exception:
        return -1

# --- Helper Function to Extract Text ---
def extract_text_from_block(bits):
    text = ""
    for pos in range(0, len(bits), 8):
        chunk = bits[pos:pos+8]
        try:
            text += chr(int("".join(map(str, chunk)), 2))
        except Exception:
            text += "?"
    return text

# --- Decoding Functions for Each Group Type --- LETS MANHANDLE THIS SHHHHHHH
def decode_grp0(grp, version):
    pty = extract_pty(grp['B'])
    if version == 0:
        # For version A, get text from both blocks C and D
        ps = extract_text_from_block(grp['C']) + extract_text_from_block(grp['D'])
    else:
        # For version B, only decode block C
        ps = extract_text_from_block(grp['C'])
    return f"PTY: {pty}, PS: {ps}"

def decode_grp1(grp, version):
    pty = extract_pty(grp['B'])
    try:
        pin = int("".join(map(str, grp['C'])), 2)
    except Exception:
        pin = -1
    return f"PTY: {pty}, PIN: {pin}"

def decode_grp2(grp, version):
    pty = extract_pty(grp['B'])
    rt = extract_text_from_block(grp['C']) + extract_text_from_block(grp['D'])
    return f"PTY: {pty}, RT: {rt}"

def decode_grp4(grp, version):
    pty = extract_pty(grp['B'])
    try:
        bits_D = grp['D']
        mins = int("".join(map(str, bits_D[:6])), 2)
        hrs = int("".join(map(str, bits_D[6:12])), 2)
        offs = int("".join(map(str, bits_D[12:])), 2)
        time_str = f"{hrs:02d}:{mins:02d} (offset {offs})"
    except Exception:
        time_str = "??"
    return f"PTY: {pty}, Clock: {time_str}"

def decode_grp10(grp, version):
    pty = extract_pty(grp['B'])
    pt_name = extract_text_from_block(grp['C']) + extract_text_from_block(grp['D'])
    return f"PTY: {pty}, PT Name: {pt_name}"

def raw_decode(grp, grp_label):
    pty = extract_pty(grp['B'])
    return f"PTY: {pty}, {grp_label}: B={bits_to_hex_str(grp['B'])} " \
           f"C={bits_to_hex_str(grp['C'])} D={bits_to_hex_str(grp['D'])}"

def decode_group(grp):
    PI = bits_to_hex_str(grp['A'])
    try:
        grp_type = int("".join(map(str, grp['B'][:4])), 2)
    except Exception:
        grp_type = -1
    version = grp['B'][4]
    msg = f"PI: {PI}, Group {grp_type}{'A' if version == 0 else 'B'}: "
    if grp_type == 0:
        msg += decode_grp0(grp, version)
    elif grp_type == 1:
        msg += decode_grp1(grp, version)
    elif grp_type == 2:
        msg += decode_grp2(grp, version)
    elif grp_type == 4:
        msg += decode_grp4(grp, version)
    elif grp_type == 10:
        msg += decode_grp10(grp, version)
    elif grp_type in [3, 5, 6, 7, 8, 9, 11, 12, 13, 14, 15]:
        type_names = {
            3: "Group 3", 5: "Group 5", 6: "Group 6", 7: "Group 7",
            8: "Group 8", 9: "Group 9", 11: "Group 11", 12: "Group 12",
            13: "Group 13", 14: "Group 14", 15: "Group 15"
        }
        msg += raw_decode(grp, type_names.get(grp_type, f"Group {grp_type}"))
    else:
        msg += f"(Unknown group {grp_type})"
    return msg

def decode_all_groups(grp_list):
    if not grp_list:
        return "No valid groups"
    return "\n".join(decode_group(g) for g in grp_list)

def process_rds_from_constellation(diff_bits):
    grps = sync_groups(diff_bits)
    result = decode_all_groups(grps)
    print("[Decoded RDS]:", result)
    return result
