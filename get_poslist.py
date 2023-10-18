# coding: utf-8
from Bio import PDB
from Bio.PDB.Polypeptide import three_to_one
import sys

parser = PDB.PDBParser()
s = parser.get_structure('test', sys.argv[1])

m = s[0]

with open('poslist.txt', 'w') as fh:
    for c in m.get_chains():
        for r in c.get_residues():
            fh.write(f"{c.id}.{three_to_one(r.resname)}.{str(r.id[1])} {c.id}\n")

