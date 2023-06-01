"""
@Title:             model_maker.py
@Name:              Mando A Ramirez
@Date:              2022 05 27

@Description:       This script is the main model maker for my coil models -- this produces the 3D atomistic and CA
coarse-grained starting structures for each model. This is an adaptation of the script I used in the DAR3-25 series
to make linkers. Now, it will be used to make the whole model.

IMPORTANTLY, it will take the SAME model_design.csv type of file used in the topology writing script
"coil_topology_assistant.py" because this script will incorporate different amino acids based on the different segments.
Specifically, it will assign the following:
    Coil sticky beads = ILE
    Multimer beads = LYS (pos 'e') and GLU (pos 'g')
    Inert beads in coils = ALA
    Inert beads in linkers = GLY

These amino acid assignments do not mean anything "physical" in the simulations, but they will be useful in handling PDB
files. I'm specifically thinking that it will be useful when selecting specific regions of models in PyMOL.
    Also, the DAR3-25 version of this script relied on a separate file titled "tozzini_converter.py" for some of the
math analysis. I am putting all those functions in this file to make things easier.

@Updates:
2022 10 28 -- updated to fix the model file parsing, so that I can remove comments.

2022 06 01 -- updated to correctly handling comments, which needs file decoding. Also updated the protein sequence so
that positions 'e' are LYS and 'g' are GLU
"""

import math as m
import argparse
import numpy as np
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB
import prody as pd


"""
Tozzini converter functions
"""
#GLOBAL ANGLE CONSTANTS -- provided in DEGREES but then converted into radians
tau = 111
omega = 180
y1 = 20.7
y2 = 14.7            # y2 should not be 20.7, causes the pseudoangles to be too low
t_r = tau * (m.pi / 180)
o_r = omega * (m.pi / 180)
y1_r = y1 * (m.pi / 180)
y2_r = y2 * (m.pi / 180)

def THETA(T, Y1, Y2, Phi1, Psi1):
    """
    The structure of this equation is complex. I will shorten it into this so that I can keep better
    track of it
        cos(theta) = cos(T)*A - B + sin(T)*C
    where,
        A = [cos(Y1)cos(Y2) - sin(Y1)sin(Y2)cos(Psi1)cos(Phi1)]
        B = sin(Y1)sin(Y2)sin(Psi1)sin(Phi1)
        c = [cos(Phi1)sin(Y1)cos(Y2) + cos(Psi1)cos(Y1)sin(Y2)]

    All variables are provided in RADIANS.
    The answer is returned in DEGREES
    """
    A = m.cos(Y1)*m.cos(Y2) - m.sin(Y1)*m.sin(Y2)*m.cos(Psi1)*m.cos(Phi1)
    B = m.sin(Y1)*m.sin(Y2)*m.sin(Psi1)*m.sin(Phi1)
    C = m.cos(Phi1)*m.sin(Y1)*m.cos(Y2) + m.cos(Psi1)*m.cos(Y1)*m.sin(Y2)

    cos_theta = m.cos(T)*A - B + m.sin(T)*C
    theta = m.acos(cos_theta) * (180 / m.pi)
    return theta


def ALPHA(Y1, Y2, Phi1, Psi1, Phi2, Psi2):
    """
    All variables are provided as RADIANS. It is lastly converted into degrees after the computation
    step.
    """
    a = (m.pi + Psi1 + Phi2 + Y1*m.sin(Psi2) + Y2*m.sin(Phi1)) * (180 / m.pi)
    return a

def ALPHA_COMPLEX(T, Y1, Y2, Phi, Psi):
    A = 0.25 * Y1**2 * m.sin(2*Psi)
    B = 0.25 * Y2**2 * m.sin(2*Phi)
    C = Y1 * Y2 * m.sin(Phi + Psi)
    D = Y1 * (T - (m.pi / 2)) * m.sin(Psi)
    E = Y2 * (T - (m.pi / 2)) * m.sin(Phi)

    ac_r = Phi + Psi + m.pi + (Y1 * m.sin(Psi)) + (Y2 * m.sin(Phi)) + A + B + C - D - E
    ac = ac_r * (180 / m.pi)    # answer is returned in DEGREES

    if ac > 180:
        ac = ac - 360
    else:
        pass
    return ac


"""
Peptide Builder Functions
"""
def build_peptide(PeptideName, Sequence):
    """
    This function makes the starting structure of a given coil model. It only operates in a single directory,
    with no functionality to specify different directory output, so be careful where you use it.
    :param PeptideName: the name given to the built model, user input
    :param Sequence: the sequence that corresponds to the model to be built. This is parsed by the main code.
    :return: None -- all the necessary files are handled and saved. Nothing is returned.
    """
    aa_name = f"{PeptideName}.pdb"
    ca_name = f"{PeptideName}_ca.pdb"
    output_file = f"{PeptideName}_aa_dihedrals.csv"

    # Start by setting up the structure and defining the phi/psi angles beyond the defaults. Leave the defaults for the
    # bond angles and bond lengths
    geo = Geometry.geometry(Sequence[0])
    # These angles were from DAR3-15
    geo.phi = -64.088
    geo.psi_im1 = -41.260
    linker = PeptideBuilder.initialize_res(geo)
    for res in Sequence[1:]:
        linker = PeptideBuilder.add_residue(linker, res, geo.phi, geo.psi_im1)

    out_linker = Bio.PDB.PDBIO()
    out_linker.set_structure(linker)
    out_linker.save(aa_name)

    # Now this part of the code will reopen this fully atomistic structure, read the torsion angles of it and convert into
    # Tozzini space, then save the C-alpha coarse grained model. It will also save the torsions and the pseudotorsions into
    # and output file
    linker_pdb = pd.parsePDB(aa_name)
    hv = pd.HierView(linker_pdb)
    outfile = open(output_file, 'w')
    outfile.write("Residue No.,Residue Name,phi (deg),psi (deg),alpha (deg),theta (deg)\n")
    #    Calculate the torsion angles and the pseudo angles and then put into the output file
    #    Keep in mind that ProDy 1-INDEXES the residues
    for i in range(1, len(Sequence) + 1):
        if i == 1 or i == len(Sequence):
            outfile.write(f"{i},{hv.getResidue('A', i).getResname()}\n")
        else:
            phi = pd.calcPhi(hv.getResidue('A', i))
            psi = pd.calcPsi(hv.getResidue('A', i))
            phi_r = phi * (np.pi / 180)
            psi_r = psi * (np.pi / 180)
            theta = THETA(t_r, y1_r, y2_r, phi_r, psi_r)
            alpha = ALPHA_COMPLEX(t_r, y1_r, y2_r, phi_r, psi_r)
            outfile.write(f"{i},{hv.getResidue('A', i).getResname()},{phi:.2f},{psi:.2f},{alpha:.2f},{theta:.2f}\n")
    outfile.close()
    # Lastly, save the structure in CA form
    linker_ca = linker_pdb.select("name CA")
    pd.writePDB(ca_name, linker_ca)
    return None


if __name__ == "__main__":
    # I don't anticipate needing to use the functions in this script elsewhere, so the
    # name==main style is a little overkill, but I'm including it here for good pythonic practice

    # Set up the argument parser!
    parser = argparse.ArgumentParser(description="Model Maker Tool -- use this to turn your desired CC protein into"
                                                 " a starting PDB structure! Easy, and corresponds with the Topology"
                                                 " Assistant too. Spreads on easy, like Ubik!")
    parser.add_argument("-df", help="Design file: this contains the information on how to build the protein, with"
                                    " alternating coil and linker segments.")
    parser.add_argument("-n", help="The number of beads in the CC protein.", type=int, required=True)
    parser.add_argument("-name", help="The name you'd like to give the CC protein. No "
                                      "extensions, please!", type=str, default="coil_model")

    # parse the arguments
    args = parser.parse_args()

    # --> This code is taken from the "topology_assistant_nonspecific_modularModel.py" script that I have.
    # this handles (ignores) comments
    # set up the list of segments, in order, and their corresponding indices, in order
    # Segments contains the names of each of the segments as they appear in order in the provided file
    # Positions contains the start/stop (and heptad) indices for the segments that correspond to the same index in the
    # Segments list
    Segments = []; Positions = []
    with open(args.df, "r", encoding="utf-8-sig") as file:
        for line in file:
            if line[0] == "#":
                # I am now adding in comment support. Comments can exist in the model generation file and can be denoted
                # as "#"
                continue
            else:
                line_split = line.rstrip("\n").split(",")
                Segments.append(str(line_split[0]))
                if str(line_split[0]) == "coil":   # uses 0-indexing, this now handles the final column being oligo type
                    Positions.append([int(line_split[1]), int(line_split[2]), int(line_split[3]),
                                    str(line_split[4]), line_split[5]])  # notice that index 5 has no type command!!!
                else:
                    Positions.append([int(line_split[1]), int(line_split[2])])

    # Now parse the model_design data and interpret the sequence
    model_sequence = ""
    for i in range(len(Segments)):
        seg = Segments[i]
        segStart = Positions[i][0]
        segStop = Positions[i][1]
        if seg == "linker":     # linker handling
            model_sequence += "G" * (segStop - segStart + 1)
        else:                   # coil handling
            a_posits = []
            d_posits = []
            e_posits = []
            g_posits = []
            a_index = Positions[i][2]
            d_index = Positions[i][2] + 3
            e_index = Positions[i][2] + 4
            g_index = Positions[i][2] + 6
            while a_index <= segStop:
                a_posits.append(a_index)
                a_index += 7
            while d_index <= segStop:
                d_posits.append(d_index)
                d_index += 7
            for c_i in range(segStart, segStop+1):
                if c_i in a_posits:
                    model_sequence += "I"
                elif c_i in d_posits:
                    model_sequence += "I"
                elif c_i in e_posits:
                    model_sequence += "A"
                elif c_i in g_posits:
                    model_sequence += "A"
                else:
                    model_sequence += "A"

    #verify that the length provided in the arguments matches the last position in the model design parameters matches
    # the length of the sequence
    assert args.n == Positions[-1][1], f"The size provided in the arguments does not match the size of the model in " \
                                       f"the file."
    assert args.n == len(model_sequence), f"The size provided in the arguments does not match the size of the model " \
                                          f"sequence made."

    build_peptide(args.name, model_sequence)