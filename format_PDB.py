#!/usr/bin/env python

from Bio import PDB, pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_sequence(structure):
    sequences = {}
    for model in structure:
        for chain in model:
            seq = Seq(''.join(residue.resname for residue in chain))
            sequences[chain.id] = SeqRecord(seq, id=chain.id)
    return sequences

def match_chains(template_sequences, model_sequences):
    matches = {}
    for template_id, template_seq in template_sequences.items():
        best_match = None
        best_score = float('-inf')
        for model_id, model_seq in model_sequences.items():
            alignments = pairwise2.align.globalxx(str(template_seq.seq), str(model_seq.seq))
            if alignments:
                score = alignments[0][2] 
                if score > best_score:
                    best_score = score
                    best_match = model_id
        if best_match:
            matches[template_id] = best_match
    return matches

def process_pdb(template_file, model_file, output_file):
    parser = PDB.PDBParser(QUIET=True)
    template_structure = parser.get_structure("template", template_file)
    model_structure = parser.get_structure("model", model_file)

    template_sequences = extract_sequence(template_structure)
    model_sequences = extract_sequence(model_structure)
    chain_matches = match_chains(template_sequences, model_sequences)
    #chain_matches = {'a': 'a', 'b': 'b', 'c': 'c', 'd': 'd', 'f': 'f', 'g': 'g', 'h': 'h', 'e': 'e'}

    print("Chain mapping:", chain_matches)

    with open(template_file, 'r') as template, open(output_file, 'w') as output:
        for line in template:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                chain_id = line[21]
                if chain_id in chain_matches:
                    residue_number = int(line[22:26])
                    atom_name = line[12:16].strip()
                    model_chain = model_structure[0][chain_matches[chain_id]]
                    if residue_number in model_chain:
                        model_residue = model_chain[residue_number]
                        if atom_name in model_residue:
                            model_atom = model_residue[atom_name]
                            new_line = (
                                line[:30] +
                                f"{model_atom.coord[0]:8.3f}{model_atom.coord[1]:8.3f}{model_atom.coord[2]:8.3f}" +
                                line[54:]
                            )
                            output.write(new_line)
                        else:
                            output.write(line)  # If atom not found in model, keep original line
                            print (f"Warning: Atom not found: {line}")
                    else:
                        output.write(line)  # If residue not found in model, keep original line
                        print (f"Warning: Residue not found: {line}")
                else:
                    output.write(line)
                    print (f"Warning: Chain not found: {line}")
            else:
                output.write(line)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Process PDB files")
    parser.add_argument("--template", required=True, help="Path to template PDB file")
    parser.add_argument("--model", required=True, help="Path to model PDB file")
    args = parser.parse_args()

    output_file = f"casp_{args.model}"
    process_pdb(args.template, args.model, output_file)
    print(f"Output written to {output_file}")
