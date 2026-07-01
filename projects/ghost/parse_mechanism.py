"""
Parse Chemkin THERM.DAT and chem.inp to generate C++ mechanism data.
Filters to C/H/O species only + N2 as inert.
"""
import re
import os

# Paths
base_dir = r"C:\Users\yanch\OneDrive - Concordia University - Canada\Documents\论文\wang\NS"
therm_file = os.path.join(base_dir, "detailed-CH4-Air", "THERM.DAT")
cheminp_file = os.path.join(base_dir, "detailed-CH4-Air", "chem.inp")
chemout_file = os.path.join(base_dir, "detailed-CH4-Air", "chem.out")
output_file = os.path.join(base_dir, "projects", "ghost", "mechanism_data.ipp")

# ============================================================
# Step 1: Identify C/H/O species from chem.out
# ============================================================
cho_species = []  # list of (name, mol_wt)
n2_found = False

with open(chemout_file, 'r') as f:
    in_species = False
    for line in f:
        if 'SPECIES' in line and 'CONSIDERED' in line:
            in_species = True
            continue
        if in_species and '----' in line:
            continue
        if in_species and line.strip() == '':
            break
        if in_species:
            # Parse species line: "  1. H                G 0   1.00797  200.0 6000.0  1  0  0  0  0"
            parts = line.split()
            if len(parts) >= 10:
                try:
                    idx = int(parts[0].rstrip('.'))
                    name = parts[1]
                    mol_wt = float(parts[3])
                    # Element counts: H C O N AR
                    elem_h = int(parts[6])
                    elem_c = int(parts[7])
                    elem_o = int(parts[8])
                    elem_n = int(parts[9])
                    elem_ar = int(parts[10]) if len(parts) > 10 else 0
                    
                    if elem_n == 0 and elem_ar == 0:
                        cho_species.append((name, mol_wt))
                    elif name == 'N2':
                        n2_found = True
                except:
                    pass

print(f"Found {len(cho_species)} C/H/O species")
for s in cho_species:
    print(f"  {s[0]} ({s[1]})")

# ============================================================
# Step 2: Parse THERM.DAT for NASA coefficients
# ============================================================
# Build a set of species names we need
target_species = set(name for name, _ in cho_species)
target_species.add('N2')

therm_data = {}  # name -> (low_T, high_T, mid_T, coeffs_low[7], coeffs_high[7])

with open(therm_file, 'r') as f:
    lines = [line.rstrip() for line in f.readlines()]

i = 0
while i < len(lines):
    line = lines[i].strip()
    if not line or line.startswith('THERMO') or line.startswith('END'):
        i += 1
        continue
    
    # First line of species entry
    # Format: NAME              notes....  T_low  T_high  T_mid  ...
    # Fields are at fixed positions in standard NASA format
    if len(line) >= 80:
        name = line[0:18].strip()
        # Check if N2 with REF ELEMENT
        if 'REF ELEMENT' in line:
            name = line.split('REF')[0].strip()
        
        if name in target_species and name not in therm_data:
            try:
                # Temperature ranges from line 1 positions 45-75
                t_range = line[45:75]
                t_parts = t_range.split()
                if len(t_parts) >= 3:
                    t_low = float(t_parts[0])
                    t_high = float(t_parts[1])
                    t_mid = float(t_parts[2])
                elif len(t_parts) >= 2:
                    t_low = float(t_parts[0])
                    t_high = float(t_parts[1])
                    t_mid = 1000.0
                else:
                    t_low = 300.0
                    t_high = 5000.0
                    t_mid = 1000.0
                
                # Lines 2, 3, 4 contain the coefficients
                coeff_text = ""
                for k in range(1, 4):
                    if i + k < len(lines):
                        coeff_text += lines[i + k]
                
                # Parse coefficients: 5 per line, format E notation
                coeffs = []
                for match in re.finditer(r'([+-]?\d+\.\d+E[+-]\d+)', coeff_text.replace(' ', '')):
                    coeffs.append(float(match.group(1)))
                
                # NASA format: 14 coefficients total
                # Line 2: a1(high), a2(high), a3(high), a4(high), a5(high)
                # Line 3: a6(high), a7(high), a1(low), a2(low), a3(low)
                # Line 4: a4(low), a5(low), a6(low), a7(low)
                if len(coeffs) >= 14:
                    coeff_high = coeffs[0:7]  # a1-a7 high
                    coeff_low = coeffs[7:14]   # a1-a7 low
                    therm_data[name] = (t_low, t_high, t_mid, coeff_low, coeff_high)
                    print(f"  Parsed THERM for {name}: T=[{t_low}, {t_high}], Tmid={t_mid}")
                else:
                    print(f"  WARNING: {name} has only {len(coeffs)} coefficients")
            except Exception as e:
                print(f"  ERROR parsing {name}: {e}")
    
    i += 1

print(f"\nParsed {len(therm_data)} species from THERM.DAT")

# ============================================================
# Step 3: Parse chem.inp for reactions
# ============================================================
# Build species name -> index mapping
species_map = {}
for idx, (name, mw) in enumerate(cho_species):
    species_map[name] = idx
species_map['N2'] = len(cho_species)  # N2 is last

reactions = []
current_species_block = []

with open(cheminp_file, 'r') as f:
    content = f.read()

# Parse reaction section
# Reactions are in format:
# REACTANTS => PRODUCTS  A  b  Ea
# Lines with / are third-body efficiencies
# LOW /  A  b  Ea  for low-pressure limit
# TROE /  alpha  T***  T*  T**  for Troe parameters

# Find reaction section
reaction_match = re.search(r'REACTIONS\s+.*?\n', content)
if not reaction_match:
    reaction_match = re.search(r'REAC.*?\n', content)

reaction_start = reaction_match.end() if reaction_match else 0
reaction_text = content[reaction_start:]

# Split into individual reactions
# Each reaction ends with a line that has the rate parameters
lines = reaction_text.split('\n')

def parse_reaction_line(line):
    """Parse a reaction line like: H+O2=OH+O  3.52E+16  -0.7  17070.0"""
    # Remove comments
    if '!' in line:
        line = line[:line.index('!')]
    
    line = line.strip()
    if not line:
        return None
    
    # Split by = to get reactants and products
    if '=' not in line and '=>' not in line:
        return None
    
    # Handle different arrow formats
    if '=>' in line:
        parts = line.split('=>')
    elif '=' in line:
        parts = line.split('=')
    else:
        return None
    
    if len(parts) < 2:
        return None
    
    reactants_str = parts[0].strip()
    products_and_rate = parts[1].strip()
    
    # Parse reactants
    reactants = []
    for r in reactants_str.split('+'):
        r = r.strip()
        if not r:
            continue
        # Handle stoichiometric coefficient (e.g., 2H2)
        stoich = 1
        match = re.match(r'(\d+)(.*)', r)
        if match:
            stoich = int(match.group(1))
            r = match.group(2).strip()
        reactants.append((r, stoich))
    
    # Parse products and rate parameters
    # Split by whitespace - last 3 values are A, b, Ea
    rate_parts = products_and_rate.split()
    if len(rate_parts) < 3:
        return None
    
    try:
        Ea = float(rate_parts[-1])
        b = float(rate_parts[-2])
        A = float(rate_parts[-3])
    except:
        return None
    
    # Products are everything before the rate parameters
    products_str = ' '.join(rate_parts[:-3])
    
    products = []
    for p in products_str.split('+'):
        p = p.strip()
        if not p:
            continue
        stoich = 1
        match = re.match(r'(\d+)(.*)', p)
        if match:
            stoich = int(match.group(1))
            p = match.group(2).strip()
        products.append((p, stoich))
    
    return {
        'reactants': reactants,
        'products': products,
        'A': A,
        'b': b,
        'Ea': Ea
    }

i = 0
while i < len(lines):
    line = lines[i].strip()
    
    if not line or line.startswith('!'):
        i += 1
        continue
    
    if line.upper() == 'END':
        break
    
    # Check if this is a reaction line
    rxn = parse_reaction_line(line)
    if rxn is None:
        i += 1
        continue
    
    # Check if all species are in our target set
    all_valid = True
    for name, _ in rxn['reactants']:
        if name not in species_map:
            if '+M' in line and name == 'M':
                continue  # M is third body
            all_valid = False
            break
    
    if not all_valid:
        i += 1
        continue
    
    for name, _ in rxn['products']:
        if name not in species_map:
            if '+M' in line and name == 'M':
                continue
            all_valid = False
            break
    
    if not all_valid:
        i += 1
        continue
    
    # Check for third-body efficiencies (lines starting with species/efficiency pairs)
    third_body_eff = {}
    has_third_body = '+M' in line
    low_params = None
    troe_params = None
    is_duplicate = 'DUPLICATE' in line.upper() or 'DUP' in line.upper()
    
    j = i + 1
    while j < len(lines) and j < i + 10:
        next_line = lines[j].strip()
        if not next_line or next_line.startswith('!'):
            j += 1
            continue
        
        # Check for third-body efficiencies
        if '/' in next_line and not next_line.upper().startswith('LOW') and not next_line.upper().startswith('TROE'):
            # Parse efficiency pairs: species/efficiency/
            for match in re.finditer(r'(\w+)\s*/\s*([\d.]+)\s*/', next_line):
                sp_name = match.group(1)
                eff = float(match.group(2))
                if sp_name in species_map:
                    third_body_eff[sp_name] = eff
            j += 1
            continue
        
        # Check for LOW pressure limit
        if next_line.upper().startswith('LOW'):
            low_match = re.match(r'LOW\s*/\s*([\d.]+E[+-]\d+)\s+([+-]?\d+\.?\d*)\s+([\d.]+E[+-]\d+)', next_line.upper())
            if low_match:
                low_params = {
                    'A': float(low_match.group(1)),
                    'b': float(low_match.group(2)),
                    'Ea': float(low_match.group(3))
                }
            j += 1
            continue
        
        # Check for TROE parameters
        if next_line.upper().startswith('TROE'):
            troe_match = re.match(r'TROE\s*/\s*([\d.]+E[+-]\d+)\s+([\d.]+E[+-]\d+)\s+([\d.]+E[+-]\d+)\s*([\d.]+E[+-]\d+)?', next_line.upper())
            if troe_match:
                troe_params = {
                    'alpha': float(troe_match.group(1)),
                    'T3': float(troe_match.group(2)),
                    'T1': float(troe_match.group(3)),
                    'T2': float(troe_match.group(4)) if troe_match.group(4) else None
                }
            j += 1
            continue
        
        break
    
    i = j
    
    # Determine reaction type
    if is_duplicate:
        rxn_type = 'DUPLICATE'
    elif troe_params is not None:
        rxn_type = 'TROE'
    elif low_params is not None:
        rxn_type = 'LINDEMANN'
    elif has_third_body:
        rxn_type = 'THIRD_BODY'
    else:
        rxn_type = 'ARRHENIUS'
    
    reactions.append({
        'reactants': rxn['reactants'],
        'products': rxn['products'],
        'A': rxn['A'],
        'b': rxn['b'],
        'Ea': rxn['Ea'],
        'type': rxn_type,
        'third_body': third_body_eff,
        'low': low_params,
        'troe': troe_params
    })

print(f"\nParsed {len(reactions)} reactions")

# ============================================================
# Step 4: Generate C++ output
# ============================================================
with open(output_file, 'w') as f:
    f.write("/*\n")
    f.write("  Auto-generated mechanism data from Chemkin files\n")
    f.write(f"  Species: {len(cho_species)} C/H/O + N2 = {len(cho_species)+1}\n")
    f.write(f"  Reactions: {len(reactions)}\n")
    f.write("*/\n\n")
    
    # Species names
    f.write("// Species name strings\n")
    f.write("__device__ const char* SPECIES_NAMES[] = {\n")
    for name, _ in cho_species:
        f.write(f'    "{name}",\n')
    f.write('    "N2"\n')
    f.write("};\n\n")
    
    # Molecular weights
    f.write("// Molecular weights [kg/mol]\n")
    f.write("__device__ __constant__ real MOL_WT[NUM_SPECIES] = {\n")
    for _, mw in cho_species:
        f.write(f"    {mw}e-3f,  // kg/mol\n")
    f.write("    28.01340e-3f\n")
    f.write("};\n\n")
    
    # NASA Tmid
    f.write("// NASA temperature midpoints [K]\n")
    f.write("__device__ __constant__ real NASA_TMID[NUM_SPECIES] = {\n")
    for name, _ in cho_species:
        if name in therm_data:
            f.write(f"    {therm_data[name][2]}f,\n")
        else:
            f.write(f"    1000.0f,  // default for {name}\n")
    if 'N2' in therm_data:
        f.write(f"    {therm_data['N2'][2]}f\n")
    else:
        f.write("    1000.0f\n")
    f.write("};\n\n")
    
    # NASA coefficients
    f.write("// NASA polynomial coefficients [14 per species]\n")
    f.write("// Layout: low(a1..a7), high(a1..a7)\n")
    f.write("__device__ __constant__ real NASA_COEFFS[NUM_SPECIES][14] = {\n")
    for name, _ in cho_species:
        if name in therm_data:
            _, _, _, coeff_low, coeff_high = therm_data[name]
            all_coeffs = coeff_low + coeff_high
            coeff_str = ", ".join(f"{c:.8e}f" for c in all_coeffs)
            f.write(f"    {{ {coeff_str} }},  // {name}\n")
        else:
            f.write(f"    {{ 0.0f }},  // {name} - MISSING THERM DATA\n")
    if 'N2' in therm_data:
        _, _, _, coeff_low, coeff_high = therm_data['N2']
        all_coeffs = coeff_low + coeff_high
        coeff_str = ", ".join(f"{c:.8e}f" for c in all_coeffs)
        f.write(f"    {{ {coeff_str} }}   // N2\n")
    else:
        f.write("    { 0.0f }   // N2 - MISSING\n")
    f.write("};\n\n")
    
    # Reactions
    f.write(f"// Reaction data: {len(reactions)} reactions\n")
    f.write("__device__ __constant__ ReactionData REACTIONS[MAX_REACTIONS] = {\n")
    
    for idx, rxn in enumerate(reactions):
        # Map reactants
        reactants = [-1, -1, -1]
        stoich_r = [0, 0, 0]
        for k, (name, s) in enumerate(rxn['reactants'][:3]):
            if name == 'M':
                continue
            if name in species_map:
                reactants[k] = species_map[name]
                stoich_r[k] = s
        
        products = [-1, -1, -1, -1]
        stoich_p = [0, 0, 0, 0]
        for k, (name, s) in enumerate(rxn['products'][:4]):
            if name == 'M':
                continue
            if name in species_map:
                products[k] = species_map[name]
                stoich_p[k] = s
        
        # Third body efficiencies
        tb_species = [-1]*8
        tb_eff = [0.0]*8
        for k, (name, eff) in enumerate(rxn['third_body'].items()):
            if k < 8 and name in species_map:
                tb_species[k] = species_map[name]
                tb_eff[k] = eff
        
        has_tb = len(rxn['third_body']) > 0
        
        low_A = rxn['low']['A'] if rxn['low'] else 0.0
        low_b = rxn['low']['b'] if rxn['low'] else 0.0
        low_Ea = rxn['low']['Ea'] if rxn['low'] else 0.0
        
        troe_alpha = rxn['troe']['alpha'] if rxn['troe'] else 0.0
        troe_T3 = rxn['troe']['T3'] if rxn['troe'] else 0.0
        troe_T1 = rxn['troe']['T1'] if rxn['troe'] else 0.0
        troe_T2 = rxn['troe']['T2'] if rxn['troe'] and rxn['troe']['T2'] else 0.0
        
        f.write(f"    {{  // Reaction {idx+1}: ")
        f.write("+".join(f"{s}{n}" if s > 1 else n for n, s in rxn['reactants']))
        f.write(" => ")
        f.write("+".join(f"{s}{n}" if s > 1 else n for n, s in rxn['products']))
        f.write(f"\n")
        f.write(f"        {rxn['A']:.6e}f, {rxn['b']:.6e}f, {rxn['Ea']:.6e}f,\n")
        f.write(f"        {rxn['type']}, {len(rxn['reactants'])}, {len(rxn['products'])},\n")
        f.write(f"        {{ {reactants[0]}, {reactants[1]}, {reactants[2]} }},\n")
        f.write(f"        {{ {products[0]}, {products[1]}, {products[2]}, {products[3]} }},\n")
        f.write(f"        {{ {stoich_r[0]}, {stoich_r[1]}, {stoich_r[2]} }},\n")
        f.write(f"        {{ {stoich_p[0]}, {stoich_p[1]}, {stoich_p[2]}, {stoich_p[3]} }},\n")
        f.write(f"        {'true' if has_tb else 'false'},\n")
        f.write(f"        {{ {tb_species[0]}, {tb_species[1]}, {tb_species[2]}, {tb_species[3]}, {tb_species[4]}, {tb_species[5]}, {tb_species[6]}, {tb_species[7]} }},\n")
        f.write(f"        {{ {tb_eff[0]}f, {tb_eff[1]}f, {tb_eff[2]}f, {tb_eff[3]}f, {tb_eff[4]}f, {tb_eff[5]}f, {tb_eff[6]}f, {tb_eff[7]}f }},\n")
        f.write(f"        {low_A:.6e}f, {low_b:.6e}f, {low_Ea:.6e}f,\n")
        f.write(f"        {troe_alpha:.6e}f, {troe_T3:.6e}f, {troe_T1:.6e}f, {troe_T2:.6e}f\n")
        f.write(f"    }},\n")
    
    f.write("};\n\n")
    f.write(f"__device__ __constant__ int NUM_REACTIONS = {len(reactions)};\n")

print(f"\nGenerated {output_file}")
print(f"  Species: {len(cho_species) + 1}")
print(f"  Reactions: {len(reactions)}")