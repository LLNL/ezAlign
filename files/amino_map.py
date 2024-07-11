# amino_map[AA_resname]= atname_map
# atname_map[AA_atname] = CG_atname
amino_map = {
	'ALA':{'CA':'BB'},
	'ARG':{'CA':'BB','CG':'SC1','CZ':'SC2'},
	'ASN':{'CA':'BB','CG':'SC1'},
	'ASP':{'CA':'BB','CG':'SC1'},
	'CYM':{'CA':'BB','SG':'SC1'},
	'CYS':{'CA':'BB','SG':'SC1'},
	'CYX':{'CA':'BB','SG':'SC1'},
	'GLN':{'CA':'BB','CG':'SC1'},
	'GLU':{'CA':'BB','CG':'SC1'},
	'GLY':{'CA':'BB'},
	'HID':{'CA':'BB','CG':'SC1','CE1':'SC2','CD2':'SC3'},
	'HIE':{'CA':'BB','CG':'SC1','CE1':'SC2','CD2':'SC3'},
	'HSE':{'CA':'BB','CG':'SC1','CE1':'SC2','CD2':'SC3'},
	'HSD':{'CA':'BB','CG':'SC1','CE1':'SC2','CD2':'SC3'},
	'HSP':{'CA':'BB','CG':'SC1','CE1':'SC2','CD2':'SC3'},
	'HIP':{'CA':'BB','CG':'SC1','CE1':'SC2','CD2':'SC3'},
	'HIS':{'CA':'BB','CG':'SC1','CE1':'SC2','CD2':'SC3'},
	'ILE':{'CA':'BB','CG1':'SC1'},
	'LEU':{'CA':'BB','CG':'SC1'},
	'LYS':{'CA':'BB','CG':'SC1','NZ':'SC2'},
	'MET':{'CA':'BB','CG':'SC1'},
	'PHE':{'CA':'BB','CG':'SC1','CE1':'SC2','CE2':'SC3'},
	'PRO':{'CA':'BB','CG':'SC1'},
	'SER':{'CA':'BB','OG':'SC1'},
	'SEP':{'CA':'BB','P':'SC1'},
	'THR':{'CA':'BB','CB':'SC1'},
	'TRP':{'CA':'BB','CG':'SC1','NE1':'SC2','CE3':'SC3','CH2':'SC4'},
	'TYR':{'CA':'BB','CD1':'SC1','CD2':'SC2','CZ':'SC3'},
	'VAL':{'CA':'BB','CB':'SC1'},
#	'CYF':{'C':'BB','SG':'SC1','CL':'F1','CO':'F2','CU':'F3','CW':'F4'}
	'CYF':{'C':'BB','SG':'C1','CK':'C2','CP':'C3','CU':'C4','CM':'C5','CR':'C6','CW':'C7'}
}

# Protein-associated residue map
# CG and AA resname must be the same here
passoc_map = {
	'CGW':{'OW':'CGW','OH2':'CGW'},
	'ZN2':{'ZN':'ZN'},
	'MG':{'MG':'MG'},
	'ATP':{'PEnd':'PG','PB':'PB','PA':'PA','C5\'':'RB1','C3\'':'RB2','C1\'':'SC1','C8':'SC2','C5':'SC3','C6':'SC4','N3':'SC5'},
	'GTP':{'PEnd':'PG','PB':'PB','PA':'PA','C5\'':'RB1','C3\'':'RB2','C1\'':'SC1','C8':'SC2','C5':'SC3','C6':'SC4','N3':'SC5'}

}
