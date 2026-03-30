'''
This code is used to calculate the buried surface area (BSA) and relative buried surface area (rBSA) for each CDR loop—or for each individual amino acid within the loops—of the wild-type antibody in complex with human CD98hc.
'''
from Bio.PDB import PDBParser, ShrakeRupley

# File paths
complex_pdb = "complex.pdb"        # Complex structure PDB file
chains_pdb = {"G+K": "GK_free.pdb"}  # Isolated chain PDB file (combined chains)
residues_of_interest = {
    'G': [31, 32, 33, 34, 35, 50, 51, 52, '52a', 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 95, 96, 97, 98, 99, 100, 
          '100a', '100b', '100c', '100d', '100e', '100f', '100g', '100h', 101, 102],
    'K': [24, 25, 26, 27, '27a', '27b', '27c', '27d', '27e', '27f', 28, 29, 30, 31, 32, 33, 34, 89, 90, 91, 92, 93, 94, 95, 96, 97]
}

# =========================
# Calculate residue-level SASA
# =========================
def calc_residue_sasa(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('X', pdb_file)
    
    sr = ShrakeRupley()
    sr.compute(structure, level="R")  
    
    residue_sasa = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != ' ':  
                    continue
                chain_id = chain.id
                resseq = residue.id[1]
                icode = residue.id[2].strip()  
                key = (chain_id, resseq, icode)
                sasa = getattr(residue, 'sasa', 0.0)
                residue_sasa[key] = sasa
    return residue_sasa

# =========================
# Calculate BSA and rBSA
# =========================
def calc_bsa_and_rBSA(complex_pdb, chains_pdb, residues_of_interest):
    print("计算 complex_pdb 的 SASA...")
    complex_sasa = calc_residue_sasa(complex_pdb)
    
    chain_sasa = {}
    for chain_name, pdb_file in chains_pdb.items():
        print(f"计算 {pdb_file} 的 SASA...")
        chain_sasa[chain_name] = calc_residue_sasa(pdb_file)
    
    print("\n计算 BSA 和 rBSA...")
    for chain in residues_of_interest:
        for res_id in residues_of_interest[chain]:
            if isinstance(res_id, str) and res_id[-1].isalpha():
                base_res_id = int(res_id[:-1])
                icode = res_id[-1]
            else:
                base_res_id = int(res_id)
                icode = ''
            
            sasa_complex = complex_sasa.get((chain, base_res_id, icode), None)
            if sasa_complex is None:
                print(f"{chain} {res_id}: WARNING: Residue not found in the complex structure")
                continue
            
            sasa_free = chain_sasa['G+K'].get((chain, base_res_id, icode), None)
            if sasa_free is None:
                print(f"{chain} {res_id}: WARNING: Residue not found in the isolated chain structure")
                continue
            
            bsa = sasa_free - sasa_complex
            rBSA = bsa / sasa_free if sasa_free > 0 else 0.0
            
            print(f"{chain} {res_id}: SASA_free = {sasa_free:.2f}, "
                  f"SASA_complex = {sasa_complex:.2f}, BSA = {bsa:.2f}, rBSA = {rBSA:.2f}")

# =========================
# running
# =========================
calc_bsa_and_rBSA(complex_pdb, chains_pdb, residues_of_interest)
