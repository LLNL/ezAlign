#!/usr/bin/env python3

import subprocess
import os
import time
import shutil
import re
import sys
import warnings
from sys import argv
import numpy
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import argparse
from collections import OrderedDict

# calls the gromacs executable and errors if it failed
def call_gromacs(args,line,prefix=''):
	if line[0:5] == 'mdrun':
		if args.cex != '':
			line = args.cex + ' ' + args.gex + ' ' + line
		else:
			line = args.gex + ' ' + line
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
		indices = tuple(int(x) for x in line[1:-1])
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

# writes ezAlign position restraints to itp files
# if they are not already present
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
	if "\n#ifdef POSRES" in buff:
		return
	posres_bufflist = ['\n#ifdef POSRES\n',
		'[ position_restraints ]\n',
		'; ai  funct  fcx    fcy    fcz\n']
	softres_bufflist = ['\n#ifdef SOFTRES\n',
                '[ position_restraints ]\n',
                '; ai  funct  fcx    fcy    fcz\n']
	for j in resmap[1]:
		posres_line = "{:>6d}     1     100000   100000   100000\n".format(j)
		softres_line = "{:>6d}     1     10   10   10\n".format(j)
		posres_bufflist.append(posres_line)
		softres_bufflist.append(softres_line)
	posres_bufflist.append('#endif\n')
	softres_bufflist.append('#endif\n')
	bufflist = [buff] + posres_bufflist + softres_bufflist
	buff = ''.join(bufflist)
	with open(filename,'w') as myfile:
		myfile.write(buff)		

# writes required itp files to include using aa_resnames
# into topology file filename
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
	j = buff.index("charmm_other.prm")
	j = buff.index('\n',j)
	bufflist = ([buff[:j+1]] + itp_bufflist +
		[buff[j+2:]])
	buff = ''.join(bufflist)
	with open(filename,'w') as myfile:
		myfile.write(buff)

def init():
	parser=argparse.ArgumentParser(prog="ezAlign",
		description="This program builds all-atom models",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)

	parser.add_argument("-f","--pdbfile",default="input_CG.pdb",help="the CG input pdb file")
	parser.add_argument("-p","--topfile",default="cg.top",help="the CG input top file")
	parser.add_argument("-m",default=None,help="Additional residue mappings file.  Must follow same format as $EZALIGN_BASE/files/residues.map.  Does not support relative paths outside the working directory.")
	parser.add_argument("-d",default=None,help="PDB file of single molecule for molecules not included in $EZALIGN_BASE/files. Should be named \"one_RESNAME.pdb\".  Does not support relative paths outside the working directory.")
	parser.add_argument("-t",default=None,help="ITP file of single molecule for molecules not included in $EZALIGN_BASE/files. Should be named \"RESNAME.itp\".  Does not support relative paths outside the working directory.")
	parser.add_argument("-o", default='ezAligned',help="Base name of equilibrated output files (.pdb, .cpt, .top).")
	parser.add_argument("-nt","--threads",default="0",help="Number of thread-MPI threads for gromacs mdrun. NOTE: do not specify if using regular MPI. If zero then the parameter is unspecified, and GROMACS uses its default configuration.")
	parser.add_argument("-gex",default="gmx",help="GROMACS executable")
	parser.add_argument("-cex",default='',help="Cluster scheduler executable (eg. srun). Prefixes gromacs mdrun commands.")
	parser.add_argument("--keep",action='store_true',help="For debugging purposes.  Retains the ezAlign working directory '.ezAlign'.")
	args=parser.parse_args()

	if "EZALIGN_BASE" not in os.environ:
		print('ERROR: Environment variable $EZALIGN_BASE must specify the ezAlign base directory. Hint: on a system running bash, execute ./configure and follow the instructions.')
		exit(2)

	if(args.pdbfile=="input_CG.pdb"):
		print("Using default CG pdb...")
	if(args.topfile=="cg.top"):
		print("Using default cg.top...")
	if(args.threads=="1"):
		print("Using only 1 thread...")
	return args

def ezAlign(args):
	RunDir = ".ezAlign"
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
	os.chdir(RunDir)
	args.pdbfile = "input_CG.pdb"
	args.topfile = "cg.top"

	all_CG = mda.Universe(args.pdbfile,in_memory=True)
	BaseDir = os.environ["EZALIGN_BASE"]
	nt=args.threads
	input_top=args.topfile
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

	os.system("cp "+BaseDir+"/files/aa_lipid.top.template aa_lipid.top")

	include_itps("aa_lipid.top",aa_resnames,BaseDir,args)
	fin_top=open(input_top,"r")
	fout=open("aa_lipid.top","a")
	start_read = 0
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
					fout.write(line)
	fout.write("\n")
	fout.close()
	fin_top.close()

	os.system("sed -i 's|\"./files|\""+BaseDir+"/files|g' aa_lipid.top")

	sed_resnames(lipids,"aa_lipid.top")

	subprocess.call("cat "+args.pdbfile+" | grep CRYST > lipid.pdb",shell=True)
	subprocess.call("cat "+args.pdbfile+" | grep CRYST > lipid_pos.pdb",shell=True)

	try:
		OF = open("log.dat", 'w')
		OFlipid = open("lipid.pdb", 'a')
		OFlipidPos = open("lipid_pos.pdb", 'a')
	except IOError:
		exit(2)
	###########################################################
	###### Align lipids and write positions and restraints ####
	###########################################################
	k=1
	if args.d is not None:
		cmdline_res = re.split('[/.]',args.d)[-2][4:]
	for j in range(len(all_lipids.residues)):
		aa_lipid = lipids[all_lipids.residues[j].resname][0]
		aa_map = lipids[all_lipids.residues[j].resname][1]
		if args.d is not None:
			if aa_lipid == cmdline_res:
				template_AA= mda.Universe("one_"+aa_lipid+".pdb",in_memory=True)
			else:
				template_AA= mda.Universe(BaseDir+"/files/one_"+aa_lipid+".pdb",in_memory=True)
		else:
			template_AA= mda.Universe(BaseDir+"/files/one_"+aa_lipid+".pdb",in_memory=True)
		template_AA.dimensions =  all_CG.dimensions
		template_AA.atoms.masses  += 72
		cg_str = (" ".join(str(x+1) for x in all_lipids.residues[j].atoms.ix))
		aa_str = (" ".join(str(x) for x in aa_map))
		for i in range(len(aa_map)):
			if i == 0:
				sub_template_AA = template_AA.select_atoms('bynum ' + str(aa_map[i]))
			else:
				sub_template_AA += template_AA.select_atoms('bynum ' + str(aa_map[i]))
			sub_CG = all_CG.select_atoms('bynum ' + str(cg_str))
		cg_com = sub_CG.center(weights=None)
		aa_com = sub_template_AA.center(weights=None)
		cg_coord = sub_CG.positions - cg_com
		aa_coord = sub_template_AA.positions - aa_com
		R, min_rmsd = align.rotation_matrix(aa_coord, cg_coord)
		template_AA.atoms.translate(-aa_com)
		template_AA.atoms.rotate(R)
		template_AA.atoms.translate(cg_com)
		old_pos = template_AA.atoms.positions
	#################################### 
	###         PRINT REF CG         ###  
	####################################
		sel2 = all_CG.select_atoms('bynum '+cg_str)
		template_POS = template_AA
		cross_AA = template_POS.atoms.positions
		for k in range(len(aa_map)):
			cross_AA[((aa_map[k])-1),:] = sel2.atoms.positions[(k),:]	
		template_POS.atoms.positions = cross_AA
		i = 0
		for i in range(len(template_POS.atoms.positions)):
			OFlipid.write("ATOM  %6d CH  LIP     4    %8.3f%8.3f%8.3f  1.00  0.00\n" % (k, old_pos[i][0], old_pos[i][1], old_pos[i][2]))
			OFlipidPos.write("ATOM  %6d CH  LIP     4    %8.3f%8.3f%8.3f  1.00  0.00\n" % (k, template_POS.atoms.positions[i][0], template_POS.atoms.positions[i][1], template_POS.atoms.positions[i][2]))
			k += 1

	OFlipid.write("TER\n ENDMDL\n")
	OFlipidPos.write("TER\n ENDMDL\n")
	OFlipid.close()
	OFlipidPos.close()
	
	#####################################
	##### EM ALL LIPIDS WITH POS RES ####
	#####################################
	call_gromacs(args,'grompp -f '+BaseDir+'/files/em_all_no_inter.mdp -c lipid.pdb -r lipid_pos.pdb -p aa_lipid.top -o t_EM_ALL -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm t_EM_ALL -c t_EM_ALL.gro -v')
	call_gromacs(args,'trjconv -f t_EM_ALL.gro -s t_EM_ALL.tpr -pbc nojump -o t_EM_pbc.pdb','echo 0 | ')
	call_gromacs(args,'grompp -f '+BaseDir+'/files/sd_no_inter.mdp -c t_EM_pbc.pdb  -r lipid_pos.pdb -p aa_lipid.top -o t_SD_ALL -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm t_SD_ALL -c t_SD_ALL.gro -v')
	call_gromacs(args,'trjconv -f t_SD_ALL.gro -s t_SD_ALL.tpr -pbc nojump -o t_SD_pbc.pdb','echo 0 | ')

	####################################
	#### CHECK FOR BAD LIPID CONTACTS ##
	####################################
	l_all = mda.Universe('t_SD_pbc.pdb',in_memory=True)
	badContacts = mda.lib.distances.self_capped_distance(l_all.atoms.positions, max_cutoff=0.1, min_cutoff=0, box=l_all.dimensions, return_distances=False)
	if len(badContacts) > 0:
		call_gromacs(args,'grompp -f '+BaseDir+'/files/g_eq.mdp -c t_SD_pbc.pdb  -r t_SD_pbc.pdb -p aa_lipid.top -o g_1 -maxwarn 34')
		call_gromacs(args,'mdrun -deffnm g_1 -c g_1.gro -v')
		call_gromacs(args,'trjconv -f g_1.gro -s g_1.tpr -pbc nojump -o t_g1_pbc.pdb','echo 0 | ')
		call_gromacs(args,'grompp -f '+BaseDir+'/files/g_eq2.mdp -c t_g1_pbc.pdb  -r t_g1_pbc.pdb -p aa_lipid.top -o g_2 -maxwarn 34')
		call_gromacs(args,'mdrun -deffnm g_2 -c g_2.gro -v')
		call_gromacs(args,'trjconv -f g_2.gro -s g_2.tpr -pbc nojump -o grow.pdb','echo 0 | ')
	else:
		subprocess.call('cp t_SD_pbc.pdb grow.pdb', shell=True)
	l_all2 = mda.Universe("grow.pdb",in_memory=True)
	badContacts2 = mda.lib.distances.self_capped_distance(l_all2.atoms.positions, max_cutoff=0.1, min_cutoff=0, box=l_all.dimensions, return_distances=False)
	subprocess.call('cat grow.pdb | grep -v ENDM | grep -v TER > full_system.pdb', shell=True)	
	OF.write("Bad contacts are %d before, and  %d after" % (len(badContacts), len(badContacts2)))
	OF.close()

	try:
		OF2 = open("full_system.pdb", 'a')
		OF3 = open("ions.pdb", 'w')
	except IOError:
		exit(2)

	###############################################
	### For Water:print four waters at beads COM ##
	###############################################
	water_CG = all_CG.select_atoms("resname W or resname WF")
	i = 0
	k=1
	num_water=0
	water_AA = mda.Universe(BaseDir+"/files/aa_water.pdb",in_memory=True)
	for residue in water_CG.residues:
		water_new = water_AA.atoms.positions + (residue.atoms.positions -  numpy.average(water_AA.atoms.positions)) 
		for j in range(12):
			OF2.write("ATOM  %6d OH  SOL     4    %8.3f%8.3f%8.3f  1.00  0.00\n" % (k, water_new[j][0], water_new[j][1], water_new[j][2]))
			k += 1
		i += 1
		num_water += 4
	OF2.close()

	##################################
	###          Add IONS          ###
	##################################
	ions_CG = all_CG.select_atoms("resname ION")
	i = 0
	k=1
	num_NA=0
	num_CL=0
	ion_AA = mda.Universe(BaseDir+"/files/aa_water_ion.pdb",in_memory=True)
	for residue in ions_CG.residues:
		ion_new = ion_AA.atoms.positions + (residue.atoms.positions - numpy.average(ion_AA.atoms.positions))
		for j in range(12):
			OF3.write("ATOM  %6d OH  SOL     4    %8.3f%8.3f%8.3f  1.00  0.00\n" % (k, ion_new[j][0], ion_new[j][1], ion_new[j][2]))
			k += 1
		num_water += 4
		OF3.write("ATOM  %6d OH  ION     4    %8.3f%8.3f%8.3f  1.00  0.00\n" % (k, ion_new[12][0], ion_new[12][1], ion_new[12][2]))
		i += 1
		if residue.atoms.names == 'Na' or residue.atoms.names  == 'Na+' or residue.atoms.names == 'NA+' or residue.atoms.names == 'SOD'or residue.atoms.names == 'NA':
			num_NA += 1
		elif residue.atoms.names == 'Cl' or residue.atoms.names == 'Cl-' or residue.atoms.names == 'CL-' or residue.atoms.names == 'CLA'or residue.atoms.names == 'CL':
			num_CL +=1
		else:
			print("NOT RIGHT ION\n")

	OF3.close()
	subprocess.call('cat ions.pdb | grep SOL >> full_system.pdb', shell=True)
	subprocess.call('cat ions.pdb | grep ION >> full_system.pdb', shell=True)
	subprocess.call('echo TER >> full_system.pdb', shell=True)
	subprocess.call('echo ENDMDL >> full_system.pdb', shell=True)

	############################
	##   Write full topology  ##
	############################
	subprocess.call('cp aa_lipid.top aa.top', shell=True)
	fout_aa=open("aa.top","a")
	fout_aa.write("SOL    %d\n NA    %d\n CL     %d\n" % (num_water, num_NA, num_CL))
	fout_aa.close()
	##########################################
	##  Minimization and equilibrate        ##
	##########################################
	call_gromacs(args,'grompp -f '+BaseDir+'/files/em1_sc.mdp -c full_system.pdb -r full_system.pdb -p aa.top -o EM_SC1 -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm EM_SC1 -v')
	call_gromacs(args,'grompp -f '+BaseDir+'/files/sd1_sc.mdp -c EM_SC1.gro -r EM_SC1.gro -p aa.top -o SD_SC1 -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm SD_SC1 -v')

	call_gromacs(args,'grompp -f '+BaseDir+'/files/em1.mdp -c SD_SC1.gro -r SD_SC1.gro -p aa.top -o EM -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm EM -c EM.gro -v')
	call_gromacs(args,'trjconv -f EM.gro -s EM.tpr -pbc nojump -o t_EM1.pdb','echo 0 | ')

	call_gromacs(args,'grompp -f '+BaseDir+'/files/eq1.mdp -c t_EM1.pdb -r t_EM1.pdb -p aa.top -o eq_1 -maxwarn 34')
	call_gromacs(args,'mdrun -deffnm eq_1 -c eq_1.gro -v')
	call_gromacs(args,'make_ndx -f eq_1.gro ','(echo del 0-100; echo r SOL NA CL; echo !0; echo name 0 SOL_ION; echo name 1 Bilayer; echo q) | ')
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
		shutil.rmtree(RunDir)

def main():
	args = init()
	ezAlign(args)
main()
