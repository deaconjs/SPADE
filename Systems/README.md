# Demo systems

Three protein–protein complexes, ready to open in the viewer and run through the
Calculate menu, including conservation and interface statistics:

| System | Complex | Chains |
|---|---|---|
| `1brs` | barnase–barstar | A, D |
| `1a0o` | CheY bound to a CheA fragment | A, B |
| `1tx4` | RhoGAP–RhoA | A, B |

```
python MolecularViewer.py Systems/1brs/1brs.pdb
```

Each folder holds the structure (`.pdb`), a ClustalW multiple sequence alignment per
chain (`<pdb><chain>.aln`, built from NCBI BLAST searches in 2005–06), and caches SPADE
would otherwise compute on first use: solvent accessibility (`.bsa`), molecular
surfaces (`.sms`), shielding contacts (`<pdb><chain>_p<placement>_d<cut-off>.ctc`) and a
thumbnail (`.jpg`). Deleting any cache file just makes SPADE recompute it.
