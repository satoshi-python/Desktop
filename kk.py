
import tool
trr = "/lustre7/home/lustre3/satoshi/MED/aff4/test_all.trr"
pdb = "/lustre7/home/lustre3/satoshi/MED/aff4/HEN.pdb"
List = [262440]
List = [int(i) for i in List]
tool.MD_to_pdb_chain(trr, pdb, List)
