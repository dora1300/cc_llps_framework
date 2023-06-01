# Protein building scripts

The python scripts in this directory are used for making 3D CC protein structures for the CC LLPS Framework. 

Using the two scripts plus a custom `design_file.csv`, you can make any desired CC protein you wish, as long as it's made only of coil segments and disordered linker segments.




## Description of files

### `model_maker.py`
`usage: model_maker.py [-h] [-df DF] -n N [-name NAME]

Model Maker Tool -- use this to turn your desired CC protein into a starting PDB structure! Easy, and corresponds with the Topology Assistant
too. Spreads on easy, like Ubik!

optional arguments:
  -h, --help  show this help message and exit
  -df DF      Definition file: this contains the information on how to build the protein, with alternating coil and linker segments.
  -n N        The number of beads in the CC protein.
  -name NAME  The name you'd like to give the CC protein. No extensions, please!`
  
This script uses PeptideBuilder, ProDy, and BioPython to make a 3D structure for a protein of your choosing. You create a `design_file.csv` which details how you want to build a CC protein and supply it to this script. It outputs 3 different files: (1) `output.pdb`, which is an all-atom representation of your protein; (2) `output_ca.pdb` which is the C-alpha coarse-grained version of the protein, and (2) `output_aa_dihedrals.csv` which calculates the dihedral angles, and the corresponding C-alpha pseudo-angles, for each residue.

IMPORTANT! -- PeptideBuilder requires a protein sequence to make a protein. The sequence you choose DOES NOT impact the behavior of the CC protein in the CC LLPS framework, because the topology files handle the interaction terms. The sequence exists so I can make a protein and so that it's easier to select beads in viewing software. Feel free to change the residues/sequence as you wish.


### `topology_assistant.py`
`usage: topology_assistant.py [-h] [-df DF] -n N [-itp] [-itp_filename ITP_FILENAME]

Topology Writer Tool -- Use this to make topology files for CC proteins, easy as pie! Cheap as Ubik, too!

optional arguments:
  -h, --help            show this help message and exit
  -df DF                Design file: this contains the information on how to build the protein, with alternating coil and linker segments.
  -n N                  The number of beads in the CC protein.
  -itp                  Activate to turn on .itp file writer.
  -itp_filename ITP_FILENAME
                        File name for the .itp file. Please provide extension.`

This script automatically makes the `.itp` topology file you need to represent your protein of interest. It uses the exact same `design_file.csv` that you use with `model_maker.py`. The overall topology `.top` file must be made manually, but I've provided `.top` files that I used in this study for your own use and modification as well.

Again, just like the above script, the sequence DOES NOT MATTER in that it has no impact on the behavior of the coils or linkers. What matters is the `type` in the `[ atoms ]` directive, which are the actual beads (coil-coil, multimer-driving, etc.) controlled by the interaction parameters set in the `.top` file.


## Explanation of `design_file.csv`
Designing CC proteins is handled in a `design_file.csv`, where you specify all the details of the CC protein. **This file must be .csv.** You can control the number of coil and linker segments, the length of each segment, the type of multimer that can form, a special interaction ID (useful for enforcing the maximally specific interaction regime), and the start of the first *'a'* heptad position. You can also use comments in the file, specified by a `#` at the start of the line. 

You specify each segment of a protein on separate line. The overall structure of a line is as follows:
Segment type , start index of segment (1-index) , end index of segment (1-index) , start index of first 'a' heptad position, type of multimer, special ID for multimer.

For example, here are two segment entries for a 1-coil-1-linker protein:
`coil,1,32,1,dim,1
linker,33,57`

Here are the acceptable values for each of the columns of a line. Case and spelling matter:
- Segment type
    - `coil` or `linker`
- Start index of segment (1-index)
    - any number, but must be consecutive with the previous segment or else the protein building won't work
- End index of segment (1-index)
    - any number, but must come after the start index obviously
- Start index of first 'a' heptad position
    - any number within the start and end indices (inclusive) where you want to start the first heptad of the coil
- Type of multimer
    - `dim` = dimer, `tri` = dimer, `tet` = tetramer
- Special ID for multimer
    - this can be any character or text. This allows you to specify unique bead types for the coil-coil sticky beads so you can enforce certain types of specific binding. **You must manually edit the .top file** to reflect the special IDs in your protein, if any.

Both the `model_maker.py` and `topology_assistant.py` scripts are not smart. There are precise and specific ways to specify each of the segment parameters, described above. If you deviate from the file setup that I've described then I can't guarantee that the scripts will work.



## Example / how-to use
I've included an example design file, `exampleProtein.csv`, for you to practice making proteins and topology files with. Play around with this file to explore the ways to make different proteins.