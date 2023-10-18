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


AA_3TO1.update({'H1S':'H','H2S':'H','H3S':'H'})

GAS_CONSTANT = 0.001987
