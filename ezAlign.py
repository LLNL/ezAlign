#!/usr/bin/env python3

import subprocess
import os
import time
import shutil
import math
import re
import sys
import warnings
from sys import argv
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import argparse
from collections import OrderedDict

PDBFMTSTR='ATOM  {!s:>5.5} {:<4s} {:<4s}X{!s:>4.4}    {:>8.3f}{:>8.3f}{:>8.3f}                         \n'

# calls the gromacs executable and errors if it failed
# nt allows custom arg.cex -n for debugging
def call_gromacs(args,line,prefix='',nt=None, nosuff=False):
	if line[0:5] == 'mdrun':
		if args.cex != '':
			pre = args.cex
			if nt is not None:
				pre += ' -n ' + str(nt)
			line = pre + ' ' + args.gex + ' ' + line
		else:
			line = args.gex + ' ' + line
		if args.suff != '' and not nosuff:
			line = line + args.suff
		if args.threads != '0':
			line += ' -nt ' + args.threads
	else:
		line = args.gex + ' ' + line
	line = prefix + line
	retval = subprocess.call(line,shell=True)
	if retval != 0:
		raise RuntimeError(str(retval))

# returns the residue map dictionary keyed by
# the CG resname read from residues.map and cmdline
def get_residue_maps(BaseDir,args):
	maps = dict()
	if args.m is not None:
		file_to_resmap(maps,args.m)
	map_file = BaseDir+"/files/residues.map"
	file_to_resmap(maps,map_file)
	return maps

# Returns the mapping index for x
def map_ind(x):
	if x == '-':
		return float('nan')
	else:
		return int(x)

# Updates maps keyed by
# the CG resname read from map_file
# earliest entry takes precedence
def file_to_resmap(maps,map_file):
	with open(map_file,'r') as myfile:
		lines = myfile.readlines()
	j = 0
	N = len(lines)
	while j < N:
		line = lines[j]
		if line[0] == '#':
			j += 1
			continue
		line = re.split('[ ,\t\n]+',' ' + line)
		if len(line) != 4:
			j += 1
			continue
		my_map = [line[2]]
		CG_resname = line[1]
		j += 1
		line = lines[j]
		line = re.split('[ ,\t\n]+',' ' + line)
		indices = tuple(map_ind(x) for x in line[1:-1])
		my_map.append(indices)
		if CG_resname not in maps:
			maps[CG_resname] = my_map
	return maps

# replaces the CG resnames with AA resnames by executing sed
# for file filename
def sed_resnames(res_maps,filename):
	for CG_resname in res_maps.keys():
		AT_resname = res_maps[CG_resname][0]
		if CG_resname == AT_resname:
			continue
		line = "sed -i 's/{:s}/{:s}/g' {:s}".format(
			CG_resname,AT_resname,filename)
		os.system(line)

# writes ezAlign position restraints to itp files.
# overwrites if already present and first restraints
# are different
def write_pos_restraints(resmap,BaseDir,args):
	if args.t is not None:
		cmdline_res = re.split('[/.]',args.t)[-2]
		if cmdline_res == resmap[0]:
			filename = args.t
		else:
			filename = BaseDir + '/files/{:s}.itp'.format(resmap[0])
	else:
		filename = BaseDir + '/files/{:s}.itp'.format(resmap[0])
	try:
		with open(filename,'r') as myfile:
			buff = myfile.read()
	except:
		return
	if len(buff) == 0:
		return
	posres_bufflist = ['\n#ifdef POSRES\n',
		'[ position_restraints ]\n',
		'; ai  funct  fcx    fcy    fcz\n']
	softres_bufflist = ['\n#ifdef SOFTRES\n',
                '[ position_restraints ]\n',
                '; ai  funct  fcx    fcy    fcz\n']
	fpos = 1000 # kJ/mol/nm^2
	fpos_soft = 10 # kJ/mol/nm^2
	posstr = "{:>6d}     1     {:d}   {:d}   {:d}\n"
	for j in resmap[1]:
		if not math.isnan(j):
			posres_line = posstr.format(
				j,fpos,fpos,fpos)
			softres_line = posstr.format(
				j,fpos_soft,fpos_soft,fpos_soft)
			posres_bufflist.append(posres_line)
			softres_bufflist.append(softres_line)
	posres_bufflist.append('#endif\n')
	softres_bufflist.append('#endif\n')
	overwrite_posres = args.overwrite_posres
	try:
		j = buff.index("\n#ifdef POSRES")
		k = buff.index('\n',j+1)+1
		k = buff.index('\n',k)+1
		k = buff.index('\n',k)+1
		k = buff.index('\n',k)+1
		pos_buff = "{:d}   {:d}   {:d}\n".format(fpos,fpos,fpos)
		n = len(pos_buff)
		if buff[k-n:k] != pos_buff:
			overwrite_posres = True
		k = buff.index("\n#ifdef SOFTRES",k)
		k = buff.index('\n',k+1)+1
		k = buff.index('\n',k)+1
		k = buff.index('\n',k)+1
		k = buff.index('\n',k)+1
		pos_buff = "{:d}   {:d}   {:d}\n".format(fpos_soft,fpos_soft,fpos_soft)
		n = len(pos_buff)
		if buff[k-n:k] != pos_buff:
			overwrite_posres = True
	except:
		j = None
		overwrite_posres = True
	if overwrite_posres == False:
		return

	#CHOL, SAPI, SAPC, POPE, POPC, SAPE, SAPS, PSM
	# G, G, G, G, G, G 
	#test_res = {'CHOL','SAPI','SAPC','POPE','POPC','SAPE','SAPS','PSM'}
	#test_res = {'CHOL','SAPI','SAPC','POPE','POPC','SAPE','PSM'}
	#if resmap[0] in test_res:
	#	bufflist = [buff[:j]] + posres_bufflist + softres_bufflist
	#else:
	#	bufflist = [buff[:j]]
	#print(resmap[0])
	bufflist = [buff[:j]] + posres_bufflist + softres_bufflist

	buff = ''.join(bufflist)
	with open(filename,'w') as myfile:
		myfile.write(buff)		

# writes required itp files to include using aa_resnames
def include_itps(filename,aa_resnames,BaseDir,args):
	with open(filename,'r') as myfile:
		buff = myfile.read()
	itp_bufflist = []
	if args.t is not None:
		cmdline_res = re.split('[/.]',args.t)[-2]
	for aa_resname in aa_resnames:
		if args.t is not None:
			if aa_resname == cmdline_res:
				line = "#include \"{:s}\"\n".format(args.t)
			else:
				line = "#include \"{:s}/files/{:s}.itp\"\n".format(BaseDir,aa_resname)
		else:
			line = "#include \"{:s}/files/{:s}.itp\"\n".format(BaseDir,aa_resname)
		itp_bufflist.append(line)
	j = 0
	while buff[j].startswith('#') or buff[j].startswith(';'):
		j = buff.index('\n',j) + 1
	bufflist = ([buff[:j]] + itp_bufflist +
		[buff[j:]])
	buff = ''.join(bufflist)
	with open(filename,'w') as myfile:
		myfile.write(buff)

# Returns a list of molecule names from the listed
# .itp files in filename. Also returns the list of
# .itp files.
def get_molnames(filename):
	molnames = []
	prot_itps = []
	with open(filename,'r') as myfile:
		buff = myfile.read()
	filelist = re.split('\n+',buff)
	N = len(filelist)
	for j in range(N):
		file = filelist[j]
		if file == '':
			continue
		prot_itp = re.split('/',file)[-1]
		prot_itps.append(prot_itp)
		molnames.append(get_molname(file))
	return molnames,prot_itps

# Returns the molecule name of .itp filename
def get_molname(filename):
	with open(filename,'r') as myfile:
		buff = myfile.read()
	j = buff.index("[ moleculetype ]")
	j = buff.index('\n',j) + 1
	while(buff[j][0] == ';'):
		j = buff.index('\n',j) + 1
	k = buff.index('\n',j)
	line = buff[j:k]
	fields = re.split('[ \t]+',line)
	return fields[0]
	
# Writes the full topology to aa.top from aa_lipid.top
def write_full_top(args,prot_molnames,prot_itps,
	num_water,num_NA,num_CL):
	with open('aa_lipid.top','r') as myfile:
		buff = myfile.read()
	j = buff.index('[ system ]')
	bufflist = [buff[:j]]
	prot_itps_listed = set()
	for prot_itp in prot_itps:
		if prot_itp not in prot_itps_listed:
			bufflist.append('#include "{:s}"\n'.format(prot_itp))
			prot_itps_listed.add(prot_itp)
	bufflist.append('\n')
	bufflist.append(buff[j:])
	if args.pi is not None:
		for molname in prot_molnames:
			bufflist.append("{:s}         1\n".format(molname))
	bufflist.append("SOL    %d\n NA    %d\n CL     %d\n" % (num_water, num_NA, num_CL))
	with open('aa.top','w') as myfile:
		myfile.write(''.join(bufflist))

# Returns an mdanalysis CA atomgroup that is protein-like:
# only residues that have atom names CA, N, and C
# NOTE: if resids are not monotonic and you have residues that
# contain CA, N, or C and aren't protein, this function could
# fail to produce the correct CA atomgroup in some edge cases
def get_aa_CA(aa_protein_template):
	sels = []
	sels.append(aa_protein_template.select_atoms('name CA'))
	sels.append(aa_protein_template.select_atoms('name N'))
	sels.append(aa_protein_template.select_atoms('name C'))
	nsels = len(sels)
	resids = np.zeros(nsels,dtype=int)
	numsels = np.zeros(nsels,dtype=int)
	for j in range(nsels):
		resids[j] = sels[j].resids[0]
		numsels[j] = sels[j].positions.shape[0]
	agslice = []
	selinds = np.zeros((nsels),dtype=int)
	while True:
		if np.all(resids == resids[0]):
			agslice.append(selinds[0])
			selinds += 1
		else:
			selinds[np.argmin(resids)] += 1
		if not (np.all(selinds < numsels)):
			break
		for j in range(nsels):
			resids[j] = sels[j].resids[selinds[j]]
	return sels[0][agslice]

# Returns aa and cg indices specified by amino_map
def map_amino(aa_ca_sel,cg_bb_sel,aa_template,cg_template):
	from files.amino_map import amino_map
	aa_ref_ixs = []
	cg_ref_ixs = []
	for j in range(aa_ca_sel.n_atoms):
		resname = aa_ca_sel[j].resname
		atname_map = amino_map[resname]
		for aa_name in atname_map:
			cg_name = atname_map[aa_name]
			aa_index = get_index(aa_name,j,aa_ca_sel,
				aa_template)
			cg_index = get_index(cg_name,j,cg_bb_sel,
				cg_template)
			aa_ref_ixs.append(aa_index)
			cg_ref_ixs.append(cg_index)
	return aa_ref_ixs,cg_ref_ixs

# Returns aa and cg indices specified by passoc_map
# starting from cg index k. Also returns new k.
def map_passoc(aa_template,cg_template,k):
	from files.amino_map import passoc_map
	aa_ixs = []
	cg_ixs = []
	n_aa = len(aa_template.atoms)
	n_cg = len(cg_template.atoms)
	j = 0
	while j < n_aa:
		resname = aa_template.atoms[j].resname
		if (resname in passoc_map):
			j,k = map_passoc_residue(aa_template,cg_template,
				aa_ixs,cg_ixs,n_aa,n_cg,passoc_map,j,k)
		else:
			j += 1
	return aa_ixs,cg_ixs,k

# Appends aa_ixs and cg_ixs for an aa residue
# starting at j, from cg index starting from k0.
# Returns j (new aa-index starting point).
# Returns k (new cg-index starting point).
def map_passoc_residue(aa_template,cg_template,aa_ixs,
	cg_ixs,n_aa,n_cg,passoc_map,j,k0):
	resname = aa_template.atoms[j].resname
	resid = aa_template.atoms[j].resid
	cur_resid = resid
	cur_aname = aa_template.atoms[j].name
	cg_resname = cg_template.atoms[k0].resname
	while cg_resname != resname:
		k0 += 1
		cg_resname = cg_template.atoms[k0].resname
	while cur_resid == resid:
		if cur_aname in passoc_map[resname]:
			aa_ixs.append(aa_template.atoms[j].ix)
			cg_ixs.append(map_passoc_atom(cg_template,resname,
				passoc_map[resname][cur_aname],k0))
		j += 1
		if j == n_aa:
			break
		cur_resid = aa_template.atoms[j].resid
		cur_aname = aa_template.atoms[j].name
	k = cg_ixs[-1] + 1
	cg_resid = cg_template.atoms[cg_ixs[-1]].resid
	while k < n_cg:
		cur_cg_resid = cg_template.atoms[k].resid
		if cur_cg_resid != cg_resid:
			break
		k += 1
	return j,k	
	
# returns the index corresponding to resname and
# aname starting from k
def map_passoc_atom(template,resname,aname,k):
	cur_aname = template.atoms[k].name
	cur_resname = template.atoms[k].resname
	resid = template.atoms[k].resid
	cur_resid = resid
	while cur_resid == resid and cur_resname == resname:
		if cur_aname == aname:
			return k
		k += 1
		cur_aname = template.atoms[k].name
		cur_resname = template.atoms[k].resname
		cur_resid = template.atoms[k].resid
	assert(0)

# Writes protein position restraints for relaxation according 
# to AA and CG names from amino_map.py.  Also writes .itp and pdbs.
def write_prot_posres(aa_ca_sel,cg_bb_sel,aa_template,cg_template,
	prot_itp,j_in,k_in):
	posres_bufflist = ['\n#ifdef POSRES_PROT\n',
		'[ position_restraints ]\n',
		'; ai  funct  fcx    fcy    fcz\n']
	posres2_bufflist = ['\n#ifdef POSRES_PROT_SOFT\n',
		'[ position_restraints ]\n',
		'; ai  funct  fcx    fcy    fcz\n']
	pdb_bufflist = []
	aa_amino_ixs,cg_amino_ixs = map_amino(aa_ca_sel,
		cg_bb_sel,aa_template,cg_template)
	aa_passoc_ixs,cg_passoc_ixs,k_out = map_passoc(
		aa_template,cg_template,k_in)
	aa_ref_ixs = aa_amino_ixs + aa_passoc_ixs
	cg_ref_ixs = cg_amino_ixs + cg_passoc_ixs
	zipped = list(zip(aa_ref_ixs,cg_ref_ixs))
	sortedzip = sorted(zipped,key=lambda x: x[0])
	aa_ref_ixs, cg_ref_ixs = zip(*sortedzip)
	k = 0
	N_k = len(aa_ref_ixs)
	f_const = 1000 #kj/mol/nm**2
	f_soft  = 100
	n_atoms = len(aa_template.atoms)
	for j in range(n_atoms):
		cur_ix = aa_template.atoms[j].ix
		resname = aa_template.atoms[j].resname
		name = aa_template.atoms[j].name
		resid = aa_template.atoms[j].resid
		crd = aa_template.atoms[j].position
		if k < N_k and cur_ix == aa_ref_ixs[k]:
			cg_pos = cg_template.atoms[cg_ref_ixs[k]].position
			pdb_bufflist.append(PDBFMTSTR.format(
				cur_ix+1,name,resname,
				resid,*cg_pos))
			posres_bufflist.append('{:>10d} 1 {:>10d} {:>10d} {:>10d}\n'.format(
				cur_ix+1,f_const,f_const,f_const))
			posres2_bufflist.append('{:>10d} 1 {:>10d} {:>10d} {:>10d}\n'.format(
				cur_ix+1,f_soft,f_soft,f_soft))
			k += 1
		else:
			pdb_bufflist.append(PDBFMTSTR.format(
				cur_ix+1,name,resname,
				resid,*crd))
	posres_bufflist.append('#endif\n')
	posres2_bufflist.append('#endif\n')
	with open(prot_itp,'a') as myfile:
		myfile.write(''.join(posres_bufflist
			+posres2_bufflist))
	if j_in == 0:
		mode = 'w'
	else:
		mode = 'a'
	with open('posres_prot.pdb',mode) as myfile:
		myfile.write(''.join(pdb_bufflist))
	return k_out

# Returns index for name using ref_sel[j]
# such that index is within the same residue
def get_index(name,j,ref_sel,template):
	resid = ref_sel.resids[j]
	ref_ix = ref_sel[j].ix
	cur_ix = ref_ix
	N_atoms = len(template.atoms)
	while cur_ix < N_atoms:
		cur_name = template.atoms[cur_ix].name
		cur_resid = template.atoms[cur_ix].resid
		if cur_resid != resid:
			break
		if cur_name == name:
			return cur_ix
		cur_ix += 1
	cur_ix = ref_ix-1
	while True:
		cur_name = template.atoms[cur_ix].name
		cur_resid = template.atoms[cur_ix].resid
		if cur_resid != resid:
			print('\nERROR: unmappable atom name "{:s}" in residue "{:d}, resname {:s}".\n'.format(
				name,resid,
				template.atoms[ref_ix].resname))
			assert(0)
		if cur_name == name:
			return cur_ix
		cur_ix -= 1

# Appends prot_molnames to aa_prot.top
# NOTE: currently only supports one protein
# per moleculename at a time
def write_prot_top(prot_molnames,prot_itps):
	with open('aa_lipid.top','r') as myfile:
		buff = myfile.read()
	j = buff.index('[ system ]')
	bufflist = [buff[:j]]
	for prot_itp in prot_itps:
		bufflist.append('#include "{:s}"\n'.format(prot_itp))
	bufflist.append('\n')
		# NOTE: this doesn't support duplicate .itps
	k = buff.index('[ molecules ]')
	k = buff.index('\n',k)+1
	bufflist.append(buff[j:k])
	for molname in prot_molnames:
		bufflist.append("{:s}      1\n".format(molname))
	with open('aa_prot.top','w') as myfile:
		myfile.write(''.join(bufflist))
	
# Returns the number of steric clashes
# in coordinate file filename
def get_contacts(filename):
	u = mda.Universe(filename,in_memory=True)
	badContacts = mda.lib.distances.self_capped_distance(
		u.atoms.positions, max_cutoff=0.7, min_cutoff=0,
		 box=u.dimensions, return_distances=False)
	return len(badContacts)

# initialization
def init():
	parser=argparse.ArgumentParser(prog="ezAlign",
		description="This program builds all-atom models",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-f","--pdbfile",default="input_CG.pdb",help="the CG input pdb file")
	parser.add_argument("-p","--topfile",default="cg.top",help="the CG input top file")
	parser.add_argument("-m",default=None,help="Additional residue mappings file.  Must follow same format as $EZALIGN_BASE/files/residues.map.  Does not support relative paths outside the working directory.")
	parser.add_argument("-d",default=None,help="PDB file of single molecule for molecules not included in $EZALIGN_BASE/files. Should be named \"one_RESNAME.pdb\".  Does not support relative paths outside the working directory.")
	parser.add_argument("-t",default=None,help="ITP file of single molecule for molecules not included in $EZALIGN_BASE/files. Should be named \"RESNAME.itp\".  Does not support relative paths outside the working directory.")
	parser.add_argument("-pp",default=None,help="If file type (.pdb): PDB file of a protein and associated residues contained in the system. If file type (.txt): list of PDB files of proteins and associated residues. Protein must be listed after lipids and small molecules but before water/ions.")
	parser.add_argument("-pi",default=None,help="If file type (.itp): ITP file of a protein and associated residues contained in the system. If file type (.txt): list of ITP files of proteins and associated residues.")
	parser.add_argument("-tt",default=None,help="Custom template topology file for atomistic simulations. If specified, this is used in place of \"$EZALIGN_BASE/files/aa_lipid.top.template\". Useful when working with custom forcefields.  This will still automatically include relevant .itps from \"$EZALIGN_BASE/files\".")
	parser.add_argument("-o", default='ezAligned',help="Base name of equilibrated output files (.pdb, .cpt, .top).")
	parser.add_argument("-nt","--threads",default="0",help="Number of thread-MPI threads for gromacs mdrun. NOTE: do not specify if using regular MPI. If zero then the parameter is unspecified, and GROMACS uses its default configuration.")
	parser.add_argument("-gex",default="gmx",help="GROMACS executable")
	parser.add_argument("-cex",default='',help="Cluster scheduler executable (eg. srun). Prefixes gromacs mdrun commands.")
	parser.add_argument("-suff",default='',help="Suffix to provide gmx mdrun commands (eg. \"-nb gpu\").")
	parser.add_argument("--overwrite_posres",action='store_true',help="Force overwrite position restraints for molecules in $EZALIGN_BASE/files. Useful if a residue's mapping has been modified.")
	parser.add_argument("--keep",action='store_true',help="For debugging purposes.  Retains the ezAlign working directory '.ezAlign'.")
	args=parser.parse_args()

	if "EZALIGN_BASE" not in os.environ:
		print('ERROR: Environment variable $EZALIGN_BASE must specify the ezAlign base directory. Hint: on a system running bash, execute ./configure and follow the instructions.')
		exit(2)

	if args.suff != '' and not args.suff.startswith(' '):
		args.suff = ' ' + args.suff

	if(args.pdbfile=="input_CG.pdb"):
		print("Using default CG pdb...")
	if(args.topfile=="cg.top"):
		print("Using default cg.top...")
	if(args.threads=="1"):
		print("Using only 1 thread...")
	return args

# generates strings from mapping indices
def maps_to_strs(aa_inds,cg_inds):
	N = len(aa_inds)
	cg_list = []
	aa_list = []
	for j in range(N):
		if not math.isnan(aa_inds[j]):
			aa_list.append(str(aa_inds[j]))
			cg_list.append(str(cg_inds[j] + 1))
	aa_str = ' '.join(aa_list)
	cg_str = ' '.join(cg_list)
	return (aa_str, cg_str)

#converts aa_map to returned aa_dict
#such that aa_dict[aa_ind] = cg_ind
#Note: aa_ind is 1-indexed, while
# cg_ind is 0-index
#Note: aa_map is a list, and
#aa_map[cg_ind] = aa_ind, while
#aa_dict[aa_ind] = cg_ind
def aa_map_to_dict(aa_map):
	aa_dict = dict()
	N = len(aa_map)
	offset = 0
	for j in range(N):
		if not math.isnan(aa_map[j]):
			aa_dict[aa_map[j]] = j + offset
		else:
			offset -= 1
	return aa_dict

# copies the listed files in filename to target
def copyfilelist(filename,target):
	with open(filename,'r') as myfile:
		buff = myfile.read()
	filelist = re.split('\n+',buff)
	for file in filelist:
		if file != '':
			shutil.copy(file,target)

# returns a list of protein .pdbs
# from filename
def get_prot_pdbs(filename):
	prot_pdbs = []
	if filename[-4:] == '.txt':
		with open(filename,'r') as myfile:
			buff = myfile.read()
		prot_pdbs = re.split('\n+',buff)
		prot_pdbs = [j for j in prot_pdbs if '.pdb' in j]
	else:
		prot_pdbs = [filename]
	return prot_pdbs

# combines .pdbs in pdb_list and writes to
# filename. Ignores headers and such.
# Can't start with HETATM
def combine_pdbs(pdb_list,filename):
	bufflist = []
	for pdb in pdb_list:
		with open(pdb,'r') as myfile:
			buff = myfile.read()
		j = buff.index('ATOM')
		k = buff.index('\nEND')
		bufflist.append(buff[j:k])
		bufflist.append('\nEND\n')
	with open(filename,'w') as myfile:
		myfile.write(''.join(bufflist))

# rigid alignment followed by flexible relaxation
#  of proteins specified in prot_molnames
def align_prots(prot_molnames,prot_itps,all_CG,args,BaseDir):
	prot_pdbs = get_prot_pdbs(args.pp)
	n_prot = len(prot_pdbs)
	aligned_pdbs = []
	cg_CA_all = all_CG.select_atoms('name BB')
	cg_start = 0
	k = 0
	for j in range(n_prot):
		aa_protein_template = mda.Universe(prot_pdbs[j],in_memory=True)
		aa_CA = get_aa_CA(aa_protein_template)
		n_CA = len(aa_CA)
		cg_end = cg_start + n_CA
		cg_CA = cg_CA_all[cg_start:cg_end]
		cg_start = cg_end
		cg_prot_com = cg_CA.center(weights=None)
		aa_prot_com = aa_CA.center(weights=None)
		cg_prot_coord = cg_CA.positions - cg_prot_com
		aa_prot_coord = aa_CA.positions - aa_prot_com
		R, min_rmsd = align.rotation_matrix(aa_prot_coord, cg_prot_coord)
		aa_protein_template.atoms.translate(-aa_prot_com)
		aa_protein_template.atoms.rotate(R)
		aa_protein_template.atoms.translate(cg_prot_com)
		aligned_pdb = "prot{:d}.pdb".format(j)
		aa_protein_template.atoms.write(aligned_pdb)
		aligned_pdbs.append(aligned_pdb)
		k = write_prot_posres(aa_CA,cg_CA,aa_protein_template,all_CG,prot_itps[j],j,k)
	aligned_filename = 'EZALIGNED_TEMPLATE_AA_prot.pdb'
	combine_pdbs(aligned_pdbs,aligned_filename)
	write_prot_top(prot_molnames,prot_itps)
	box = all_CG.dimensions[:3] * 0.1
	call_gromacs(args,'editconf -f {:s} -o aa_prot1.gro -box {:.3f} {:.3f} {:.3f} -noc'.format(
		aligned_filename,*box))
	call_gromacs(args,'grompp -f '+BaseDir+'/files/em1_prot.mdp -c aa_prot1.gro -r posres_prot.pdb -p aa_prot.top -o em1_prot -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm em1_prot -c em1_prot.pdb',nt=1)
	call_gromacs(args,'grompp -f '+BaseDir+'/files/md1_prot.mdp -c em1_prot.pdb -r posres_prot.pdb -p aa_prot.top -o md1_prot -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm md1_prot -c md1_prot.pdb')
	call_gromacs(args,'grompp -f '+BaseDir+'/files/md2_prot.mdp -c md1_prot.pdb -r posres_prot.pdb -p aa_prot.top -o md2_prot -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm md2_prot -c md2_prot.pdb')
	subprocess.call('cat md2_prot.pdb | grep ATOM >> full_system.pdb', shell=True)


# rigid alignment followed by flexible relaxation
# of lipids and small molecules
def align_lipids(BaseDir,args,lipids,all_lipids,all_CG):
	sed_resnames(lipids,"aa_lipid.top")
	subprocess.call("cat "+args.pdbfile+" | grep CRYST > lipid.pdb",shell=True)
	subprocess.call("cat "+args.pdbfile+" | grep CRYST > lipid_pos.pdb",shell=True)
	
	###########################################################
	###### Align lipids and write positions and restraints ####
	###########################################################
	k=1
	bufflist_lipid = []
	bufflist_lipidpos = []
	if args.d is not None:
		cmdline_res = re.split('[/.]',args.d)[-2][4:]
	for j in range(len(all_lipids.residues)):
		aa_lipid = lipids[all_lipids.residues[j].resname][0]
		aa_map = lipids[all_lipids.residues[j].resname][1]
		aa_str, cg_str = maps_to_strs(aa_map,
			all_lipids.residues[j].atoms.ix)
		if args.d is not None:
			if aa_lipid == cmdline_res:
				template_AA= mda.Universe("one_"+aa_lipid+".pdb",in_memory=True)
			else:
				template_AA= mda.Universe(BaseDir+"/files/one_"+aa_lipid+".pdb",in_memory=True)
		else:
			template_AA= mda.Universe(BaseDir+"/files/one_"+aa_lipid+".pdb",in_memory=True)
		template_AA.dimensions =  all_CG.dimensions
		template_AA.atoms.masses  += 72
		sub_AA = template_AA.atoms[[]]
		for j in range(len(aa_map)): # to preserve order
			if not math.isnan(aa_map[j]):
				sub_AA += template_AA.select_atoms('bynum '
					+ str(aa_map[j]))
		sub_CG = all_CG.select_atoms('bynum ' + cg_str)
		cg_com = sub_CG.center(weights=None)
		aa_com = sub_AA.center(weights=None)
		cg_coord = sub_CG.positions - cg_com
		aa_coord = sub_AA.positions - aa_com
		R, min_rmsd = align.rotation_matrix(aa_coord, cg_coord)
		template_AA.atoms.translate(-aa_com)
		template_AA.atoms.rotate(R)
		template_AA.atoms.translate(cg_com)
		aa_pos = template_AA.atoms.positions
	#################################### 
	###         PRINT REF CG         ###  
	####################################
		aa_dict = aa_map_to_dict(aa_map)
		for i in range(len(aa_pos)):
			aa_name = template_AA.atoms[i].name
			aa_resid = template_AA.atoms[i].resid
			bufflist_lipid.append(PDBFMTSTR.format(
				k,aa_name,aa_lipid,aa_resid, 
				*aa_pos[i]))
			if i + 1 in aa_dict:
				cg_name = sub_CG.atoms[aa_dict[i+1]].name
				bufflist_lipidpos.append(PDBFMTSTR.format(
					k,cg_name,aa_lipid,aa_resid,
					*sub_CG.atoms.positions[aa_dict[i+1]]))
			else:
				bufflist_lipidpos.append(PDBFMTSTR.format(
					k,aa_name,aa_lipid,aa_resid,
					*aa_pos[i]))
			k += 1
	bufflist_lipid.append("TER\nENDMDL\n")
	bufflist_lipidpos.append("TER\nENDMDL\n")
	with open("lipid.pdb",'a') as f:
		f.write(''.join(bufflist_lipid))
	with open("lipid_pos.pdb",'a') as f:
		f.write(''.join(bufflist_lipidpos))
	
	#####################################
	##### EM ALL LIPIDS WITH POS RES ####
	#####################################
	
	call_gromacs(args,'grompp -f '+BaseDir+'/files/em_all_no_inter.mdp -c lipid.pdb -r lipid_pos.pdb -p aa_lipid.top -o t_EM_ALL -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm t_EM_ALL -c t_EM_ALL.gro -v')
	call_gromacs(args,'trjconv -f t_EM_ALL.gro -s t_EM_ALL.tpr -pbc nojump -o t_EM_pbc.pdb','echo 0 | ')
	call_gromacs(args,'grompp -f '+BaseDir+'/files/sd_no_inter.mdp -c t_EM_pbc.pdb  -r lipid_pos.pdb -p aa_lipid.top -o t_SD_ALL -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm t_SD_ALL -c t_SD_ALL.gro -v')
	call_gromacs(args,'trjconv -f t_SD_ALL.gro -s t_SD_ALL.tpr -pbc nojump -o t_SD_pbc.pdb','echo 0 | ')
	
	subprocess.call('cat t_SD_pbc.pdb | grep CRYST1 > full_system.pdb', shell=True)
	subprocess.call('cat t_SD_pbc.pdb | grep ATOM >> full_system.pdb', shell=True)

# Core ezAlign function, runs in RunDir.
# Outputs args.o (.pdb, .cpt, .top) files.
def ezAlign(args):
	RunDir = ".ezAlign"
	prot_molnames = []
	prot_itps = []
	
	if os.path.isdir(RunDir):
		shutil.rmtree(RunDir)
	os.mkdir(RunDir)
	
	shutil.copy(args.pdbfile,RunDir + '/' + "input_CG.pdb")
	shutil.copy(args.topfile,RunDir + '/' + "cg.top")
	if args.m is not None:
		shutil.copy(args.m,RunDir)
	if args.d is not None:
		shutil.copy(args.d,RunDir)
	if args.t is not None:
		shutil.copy(args.t,RunDir)
	if args.pi is not None and args.pi[-4:] == ".txt":
		copyfilelist(args.pi,RunDir)
		copyfilelist(args.pp,RunDir)
		prot_molnames,prot_itps = get_molnames(args.pi)
	elif args.pi is not None:
		shutil.copy(args.pi,RunDir)
		prot_molnames = {get_molname(args.pi)}
		prot_itps = [re.split('/',args.pi)[-1]]
	if args.pp is not None:
		shutil.copy(args.pp,RunDir)
	if args.tt is not None:
		shutil.copy(args.tt,RunDir+'/'+"aa_lipid.top")
	
	os.chdir(RunDir)
	args.pdbfile = "input_CG.pdb"
	args.topfile = "cg.top"

	BaseDir = os.environ["EZALIGN_BASE"]
	nt=args.threads
	input_top=args.topfile
	
	all_CG = mda.Universe(args.pdbfile,in_memory=True)
	
	lipids = get_residue_maps(BaseDir,args)
	for resmap in lipids.values():
		write_pos_restraints(resmap,BaseDir,args)
	lipids_str = (" ".join(str(x) for x in lipids.keys()))
	all_lipids = all_CG.select_atoms("resname "+lipids_str)

	cg_resnames = set(OrderedDict.fromkeys(
		all_CG.residues.resnames))

	cg_resnames = cg_resnames & set(lipids.keys())
	aa_resnames = set()
	for cg_resname in cg_resnames:
		aa_resnames.add(lipids[cg_resname][0])

	####################################################
	###### Write top files and other files set-up  #####
	####################################################
	
	if args.tt is None:
		os.system("cp "+BaseDir+"/files/aa_lipid.top.template aa_lipid.top")

	include_itps("aa_lipid.top",aa_resnames,BaseDir,args)
	os.system("sed -i 's|\"./files|\""+BaseDir+"/files|g' aa_lipid.top")
	
	fin_top=open(input_top,"r")
	fout=open("aa_lipid.top","a")
	start_read = 0
	n_moltypes = 0
	for line in fin_top:
		if "[ molecules ]" in line:
			start_read=1
			continue
		if start_read:
			data=line.strip().split()
			if len(data) and not line.startswith(";"):
				if (data[0]!="W" and data[0]!="WF" and 
					data[0]!="NA" and data[0]!="Na" and 
					data[0]!="NA+" and data[0]!="Na+" and 
					data[0]!="CL" and data[0]!="Cl" and 
					data[0]!="CL-" and data[0]!="Cl-"):
					if args.pi is None:
						fout.write(line)
						n_moltypes += 1
					elif data[0] not in prot_molnames:
						fout.write(line)
						n_moltypes += 1
	if start_read != 1:
		raise RuntimeError('"[ molecules ]" directive in cg topology not found.')
	fout.write("\n")
	fout.close()
	fin_top.close()

	if n_moltypes > 0:
		align_lipids(BaseDir,args,lipids,all_lipids,all_CG)
	else:
		subprocess.call("cat "+args.pdbfile+" | grep CRYST > full_system.pdb",shell=True)

	##############################
	##     PROTEIN ALIGNMENT    ##
	##############################

	if args.pi is not None and args.pp is not None:
		align_prots(prot_molnames,prot_itps,all_CG,args,BaseDir)

	###############################################
	### For Water:print four waters at beads COM ##
	###############################################
	water_CG = all_CG.select_atoms("resname W or resname WF")
	i = 0
	k=1
	num_water=0
	water_AA = mda.Universe(BaseDir+"/files/aa_water.pdb",in_memory=True)
	n_CGW = water_CG.positions.shape[0]
	aa_cog = np.average(water_AA.atoms.positions)
	bufflist_waters = []
	# could maybe speed below up via vectorization somehow
	for j in range(n_CGW):
		water_new = water_AA.atoms.positions + (water_CG.positions[j] - aa_cog) 
		for j in range(12):
			bufflist_waters.append(PDBFMTSTR.format(
				k,'OH','SOL',4, water_new[j][0], 
				water_new[j][1], water_new[j][2]))
			k += 1
		i += 1
		num_water += 4
	with open("full_system.pdb",'a') as f:
		f.write(''.join(bufflist_waters))

	##################################
	###          Add IONS          ###
	##################################
	ions_CG = all_CG.select_atoms("resname ION or resname NA or resname CL")
	i = 0
	k=1
	num_NA=0
	num_CL=0
	ion_AA = mda.Universe(BaseDir+"/files/aa_water_ion.pdb",in_memory=True)
	ion_AA_cog = np.average(ion_AA.atoms.positions)
	n_CGI = ions_CG.positions.shape[0]
	pnames = {'Na','Na+','NA+','SOD','NA'}
	nnames = {'Cl','Cl-','CL-','CLA','CL'}
	bufflist_ions = []
	for j in range(n_CGI):
		ion_new = ion_AA.atoms.positions + (ions_CG.positions[j] - ion_AA_cog)
		for l in range(12):
			bufflist_ions.append(PDBFMTSTR.format(
				k,'OH','SOL',4,ion_new[l][0], 
				ion_new[l][1],ion_new[l][2]))
			k += 1
		num_water += 4
		bufflist_ions.append(PDBFMTSTR.format(
			k,'OH','ION',4,ion_new[12][0], 
			ion_new[12][1],ion_new[12][2]))
		i += 1
		if (ions_CG[j].name in pnames):
			num_NA += 1
		elif (ions_CG[j].name in nnames):
			num_CL +=1
		else:
			print("NOT RIGHT ION\n")
			assert(0)	
	with open("ions.pdb",'w') as f:
		f.write(''.join(bufflist_ions))
	
	subprocess.call('cat ions.pdb | grep SOL >> full_system.pdb', shell=True)
	subprocess.call('cat ions.pdb | grep ION >> full_system.pdb', shell=True)
	subprocess.call('echo TER >> full_system.pdb', shell=True)
	subprocess.call('echo ENDMDL >> full_system.pdb', shell=True)

	write_full_top(args,prot_molnames,prot_itps,num_water,num_NA,num_CL)
	##########################################
	##  Minimization and equilibrate        ##
	##########################################
	#os.environ['GMX_SUPPRESS_DUMP'] = '1'

	#print('\n Clashes in full_system.pdb: {:d}\n'.format(
	#	get_contacts('full_system.pdb')))
	# above is broken, probably needs trjconv first

	call_gromacs(args,'grompp -f '+BaseDir+'/files/em1_sc.mdp -c full_system.pdb -r full_system.pdb -p aa.top -o EM_SC1 -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm EM_SC1 -v')
	call_gromacs(args,'grompp -f '+BaseDir+'/files/sd1_sc.mdp -c EM_SC1.gro -r EM_SC1.gro -p aa.top -o SD_SC1 -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm SD_SC1 -v')
	
	print('\n Clashes in SD_SC1.gro: {:d}\n'.format(
		get_contacts('SD_SC1.gro')))

	call_gromacs(args,'grompp -f '+BaseDir+'/files/em1.mdp -c SD_SC1.gro -r SD_SC1.gro -p aa.top -o EM -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm EM -c EM.gro -v')
	call_gromacs(args,'trjconv -f EM.gro -s EM.tpr -pbc nojump -o t_EM1.pdb','echo 0 | ')
	call_gromacs(args,'make_ndx -f EM.gro ','(echo del 0-100; echo r SOL NA CL; echo !0; echo name 0 SOL_ION; echo name 1 Bilayer; echo q) | ')
	
	call_gromacs(args,'grompp -f '+BaseDir+'/files/eq1.mdp -c t_EM1.pdb -r t_EM1.pdb -p aa.top -o eq_1 -n index.ndx -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm eq_1 -c eq_1.gro -v')
	call_gromacs(args,'trjconv -f eq_1.gro -s eq_1.tpr -pbc nojump -o t_EQ1.pdb','echo 0 | ')
	call_gromacs(args,'grompp -f '+BaseDir+'/files/eq2.mdp -c t_EQ1.pdb -r t_EQ1.pdb -p aa.top -o eq_2 -n index.ndx -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm eq_2 -c eq_2.pdb -v')
	call_gromacs(args,'grompp -f '+BaseDir+'/files/eq3.mdp -c eq_2.pdb  -p aa.top -o eq_3 -n index.ndx -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm eq_3 -c eq_3.pdb  -v')
	call_gromacs(args,'grompp -f '+BaseDir+'/files/eq4.mdp -c eq_3.pdb  -p aa.top -o eq_4 -n index.ndx -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm eq_4 -c eq_4.pdb -v')
	
	shutil.copyfile('eq_4.pdb','../'+args.o+'.pdb')
	shutil.copyfile('eq_4.cpt','../'+args.o+'.cpt')
	shutil.copyfile('aa.top','../'+args.o+'.top')

	if not args.keep:
		os.chdir('..')
		shutil.rmtree(RunDir,ignore_errors=True)

def main():
	args = init()
	ezAlign(args)
main()
