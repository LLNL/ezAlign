# ezAlign
Coarse grain to atomistic molecular coordinate and topology converter for molecular dynamics simulations.
## Description
ezAlign takes coarse grain (CG) coordinate and topology files and converts and outputs their corresponding atomistic formats using an alignment and relaxation procedure outlined in reference [1].  ezAlign is designed to convert complex, solvated biological systems including lipid membranes, drug-like molecules, and proteins using GROMACS [2]. A GROMACS checkpoint (.cpt) file is also outputted to enable continuation simulations that retain the equilibrated atomic velocities.  Note that ezAlign specifically requires GROMACS 2022 or later.

Independent atomistic coordinates and topologies for every molecule must already be included in `ezAlign/files`.  Single molecule coordinate, topology, and mapping may also be specified during execution using the appropriate command-line options.  A number of commonly simulated biological molecules are currently provided.  For instructions on including additional molecules, see the parameterization section.  A protein complex can now be mapped by specifying a template coordinate (-pp) and topology file (-pi).  Amino acid mappings are stored in `files/amino_map.py`, which can be modified to include non-standard amino acid mappings.

## Requirements
* Linux OS
* Python 3.x
* MDAnalysis 2.x.x or greater
* numpy 1.20.x or greater
* GROMACS 2022.x or 2023.x

Note: Gromacs 2024 has changed the minimum allowed short-range cutoffs, so some of ezAlign's .mdp files need to be adjusted to use it.

### Topology requirements
Input GROMACS topology directives must be specified as `[ directive ]` (not `[directive]`).

## Installation
1. Clone or download ezAlign
2. Execute `./configure`
3. Follow the outputted sourcing instructions

## Arguments
Execute `ezAlign.py -h` for a brief description of arguments.

## Examples
A working example of CG oleate in a POPC bilayer is provided in the examples folder `examples/oleate`.  Within the folder, simply execute
```
ezAlign.py -f input_CG.pdb -p cg.top
```
### Protein examples
Two single-complex proteins in POPC bilayers, hERG and GABA_A are provided.  Since these systems are larger than the oleate example, you may want to run ezAlign with some parallelization (see arguments `-nt`, `-gex`, and/or `-cex`) for completion within 30 minutes.

For hERG, within `examples/hERG`, execute
```
ezAlign.py -f cg.pdb -p cg.top -pi aa_herg.itp -pp aa_herg.pdb
```
For GABA_A, within `examples/GABAA`, execute
```
ezAlign.py -f cg.pdb -p cg.top -pi aa_prot.itp -pp aa_prot.pdb
```
## Parameterization
Including parameters for new molecules in ezAlign is fairly straightforward using CHARMM-GUI [3].

1. Create and login to your account at https://www.charmm-gui.org/.
2. Select Input Generator -> Ligand Reader & Modeler, or click https://www.charmm-gui.org/?doc=input/ligandrm.
3. Input the molecular structure of the molecule to be parameterized using the various options.  For example, butanol can be inputted by pasting `CCCCO` into the "Load SMILES" box and then clicking "Load SMILES".
4. Check the box next to "Find similar residues in the CHARMM FF".
5. Click "Next Step: Search ligand" at the bottom right of the page.
6. Review the residues found, and select one appropriately.  In the case of butanol, it is not in the CHARMM forcefield, so select the CGenFF option.  You may rename the residue in the box with "LIG", in this case we use "BOL" (martini's residue name for butanol).
7. Click "Next Step: Generate PDB" in the bottom right of the page and wait.
8. Click "download.tgz" in the top right corner.
9. Extract charmm-gui.tgz.
10. Copy charmm-gui/ligandrm.pdb to $EZALIGN_BASE/files/one_BOL.pdb
11. Copy charmm-gui/gromacs/BOL.itp to $EZALIGN_BASE/files/BOL.itp
12. Modify $EZALIGN_BASE/files/residues.map to contain a CG bead to AA atom mapping for butanol.  For butanol this may look like
```
BOL BOL
3
```
which maps the single martini CG bead to the central carbon (C3) of atomistic butanol.

## Protein backmapping
The ezAlign program currently supports both single-complex protein and multi-complex protein systems via the `-pi` and `-pp` options.  To map a single-complex protein system, simply specify `-pi` with the .itp file of the protein complex and `-pp` with the template .pdb file of the protein complex, which is aligned and relaxed into the CG conformation.  Proteins within a single-complex are aligned together during the initial rigid alignment.  If the proteins' relative positions are expected to deviate significantly from the provided template structure specified with `-pp`, it is advised to treat each protein as an independent complex, following the multi-complex protein insctructions in the section below. See the "Protein examples" section above for provided single-complex protein systems.  Note the name of the protein complex in the `[ molecules ]` directive of the CG .top file provided to `-p` should match the name of the protein complex specified by the `[ moleculetype]` directive in the .itp file provided to `-pi`.

### Multi-complex proteins
To backmap multiple protein complexes with ezAlign, specify `-pi` and `-pp` with .txt files that list full paths of each complex's .itp and template .pdb, respectively, delimited by newlines.  Each protein complex is independently aligned during the initial rigid alignment.  The names of each protein complex in the `[ molecules ]` directive of the CG .top file provided to `-p` should match the names of the protein complexes specified by the `[ moleculetype]` directive in the .itp files provided to `-pi`.  The ordering of the protein complexes should be consistent across `-p`, `-pi`, and `-pp`.

### Non-standard amino acids and protein-associated residues
Both non-standard amino acids and protein-associated residues are supported via modification of `$EZALIGN_BASE/files/amino_map.py`.  Simply add a corresponding entry to the `amino_map` python dictionary, following the format specification at the top of the file.  To add a protein-associated residue, similarly add a corresponding entry to the `passoc_map` python dictionary.  Note that protein-associated residues are residues that are included in the protein structure and topology files provided to `-pp` and `-pi`, respectively.  During rigid alignment, protein-associated residues are aligned with their corresponding protein complexes.  If a protein-associated residue is expected to deviate significantly from its relative position specified in `-pp`, consider specifying the residue as an independent molecule instead (see the "Parameterization" section above).

## References
1. Bennett, W. F. D., Bernardi, A., Ozturk, T. N., Ingólfsson, H. I., Fox, S. J., Sun, D., & Maupin, C. M. (2024). ezAlign: A Tool for Converting Coarse-Grained Molecular Dynamics Structures to Atomistic Resolution for Multiscale Modeling. Molecules, 29(15), 3557, DOI: [10.3390/molecules29153557](https://doi.org/10.3390/molecules29153557).
2. Abraham, M. J., Murtola, T., Schulz, R., Páll, S., Smith, J. C., Hess, B., & Lindahl, E. (2015). GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. SoftwareX, 1, 19-25.
3. Jo, S., Cheng, X., Lee, J., Kim, S., Park, S. J., Patel, D. S., ... & Im, W. (2017). CHARMM‐GUI 10 years for biomolecular modeling and simulation. Journal of computational chemistry, 38(15), 1114-1124.

## License
ezAlign is distributed under the GNU General Public License (GPL) v2.0 license. 

LLNL-CODE-846696
