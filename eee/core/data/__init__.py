"""
Data used in calculations. 
"""

AA = [("A","ALA"),
      ("C","CYS"),
      ("D","ASP"),
      ("E","GLU"),
      ("F","PHE"),
      ("G","GLY"),
      ("H","HIS"),
      ("I","ILE"),
      ("K","LYS"),
      ("L","LEU"),
      ("M","MET"),
      ("N","ASN"),
      ("P","PRO"),
      ("Q","GLN"),
      ("R","ARG"),
      ("S","SER"),
      ("T","THR"),
      ("V","VAL"),
      ("W","TRP"),
      ("Y","TYR")]

AA_3TO1 = dict([(a[1],a[0]) for a in AA])
AA_1TO3 = dict([(a[0],a[1]) for a in AA])

AA_TO_INT = dict([(a,i) for i, a in enumerate(AA_1TO3)])
AA_TO_INT["-"] = len(AA_TO_INT) + 1
INT_TO_AA = list(AA_TO_INT.keys())

GAS_CONSTANT = 0.001987
