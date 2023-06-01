"""
@Title:             topology_assistant.py
@Name:              Mando A Ramirez
@Date:              2022 10 28

@Description:       This script is an extension of "topology_assistant_nonspecific.py" to handle the modular multimer
proteins that are described in DAR3-36b. This allows me to specify oligomerization state for the coils.

This code goes back and forth between using 1 indexing and 0 indexing, which ever is most convenient for the application.
Be mindful of this! I will try to denote what indexing is being used.

Also be mindful that there is *NO* difference in .itp writer for specific vs non-specific, because specificity is
handled at the .top [nonbond_params] level.

As of 20221028 the topology writer is not updated so I strongly recommend against using it! Just make your own .top
file for now!


@Updates:
2022 11 10 -- fixed an error in the atom_helper function. Originally, it only printed "Di" beads and I fixed it to
print the correct type of oligo bead based on the specified oligomer type.

2023 05 10 -- I'm removing the .top topology writing feature since I don't use it. It might get added back at some
point.

2023 06 01 -- updated to actually correctly handle comments (which needs decoding if file comes from Excel).
"""

import argparse

""" Constants """
# Since I have now found parameters to use in the framework, I will make those parameters constant
DAMP = 0.01                         # From DAR3-25p, used for damping linker parameters

MASS = 109.0
BOND_LENGTH = 0.385
BOND_FORCE = 50000.000

COIL_P14_C6 = 0.069130          # 1-4 pairs C6 (attractive) term
COIL_P14_C12 = 5.97404E-04      # 1-4 pairs C12 (repulsive) term
COIL_P15_C6 = 0.249383          # 1-5 pairs C6 (attractive) term
COIL_P15_C12 = 7.773996E-03     # 1-5 pairs C6 (attractive) term
LINK_P14_C6 = 0.069130 * DAMP          # 1-4 pairs C6 (attractive) term FOR LINKERS
LINK_P14_C12 = 5.97404E-04 * DAMP      # 1-4 pairs C12 (repulsive) term FOR LINKERS
LINK_P15_C6 = 0.249383 * DAMP          # 1-5 pairs C6 (attractive) term FOR LINKERS
LINK_P15_C12 = 7.773996E-03 * DAMP     # 1-5 pairs C6 (attractive) term FOR LINKERS

THETA_HELIX = 93
THETA_HELIX_FORCE = 50.0
ALPHA1_HELIX = -130.0
ALPHA1_HELIX_FORCE = 20
ALPHA3_HELIX = -30.0
ALPHA3_HELIX_FORCE = 20

THETA_LINK_FORCE = 50 * DAMP
ALPHA1_LINK_FORCE = 20 * DAMP
ALPHA3_LINK_FORCE = 20 * DAMP




""" Function Definitions """
def write_itp(itp_filename, Length, segments_list, positions_list):
    """

    This is the order that the .itp file will take
        [ moleculetype ]
        [ atoms ]
        [ bonds ]
        [ pairs ]
        [ exclusions ]
        [ angles ]
        [ dihedrals ]
    :return: None
    """
    A_BEADS = []
    A_BEADS_BY_COIL = []
    EG_BEADS_BY_COIL = []

    opener = f""";
; Individual topology file (.itp) for a CC protein
; Title: {itp_filename}
; Using the modular multimer scheme
; This file is generated automatically by "toplogy_assistant.py"
; This file MUST BE MANUALLY INCLUDED IN THE OVERALL .top FILE!
;

[ moleculetype ]
;name          nrexcl
coil_model     4

"""

    # next up do the atoms
    atom_information = f"""

[ atoms ]
;nr  type  resnr  residue  atom  cgnr  charge  mass
"""
    for index in range(len(segments_list)): # this is 0-indexed
        if segments_list[index] == "coil":
            temp_string, temp_a_posits, temp_d_posits, temp_e_posits, temp_g_posits = atom_helper(
                positions_list[index][0],
                positions_list[index][1],
                segment_type="coil",
                heptad_start=positions_list[index][2],
                oligo=positions_list[index][3],
                oligoID=positions_list[index][4]
            )   # all these are 0-indexed too
            atom_information += temp_string
            A_BEADS.extend(temp_a_posits)           # all these positions are 1-indexed
            A_BEADS.extend(temp_d_posits)           # all these positions are 1-indexed
            merged_coil_list = temp_a_posits + temp_d_posits
            merged_coil_list.sort()
            merged_eg_list = temp_e_posits + temp_g_posits  # important for doing the right exclusions, below
            merged_eg_list.sort()
            # please keep in mind that .sort() does NOT return a sorted listed, it simply sorts the list in place
            # it actually returns "None". So do the sort on the list first, then pass the sorted list wherever you want.
            A_BEADS_BY_COIL.append(merged_coil_list)
            EG_BEADS_BY_COIL.append(merged_eg_list)
        else:
            temp_string = atom_helper(
                positions_list[index][0],
                positions_list[index][1]
            )   # all these are 0-indexed too
            atom_information += temp_string

    # now write the bond information
    bond_information = f"""

[ bonds ]
;i    j    func   r0 (nm)        kb (kJ/mol/nm^2)
"""
    for posb in range(1, Length):    # 1-indexed
        bond_information += f" {posb}    {posb+1}    1        {BOND_LENGTH}        {BOND_FORCE}\n"


    # Now write out ALL the pairs
    pair_information = f"""

[ pairs ]
;i   j   func     C6(attrac)    C12(repul)
;    1-4 pairs
"""
    for index1 in range(len(segments_list)): # this is 0-indexed
        if segments_list[index1] == "coil": # but the reference to positions_list values are 1-indexed
            temp_string = pairs_helper(
                positions_list[index1][0],
                positions_list[index1][1],
                segment_type="coil",
                pair_type="1-4"
            )
            pair_information += temp_string
        else:
            temp_string = pairs_helper(
                positions_list[index1][0],
                positions_list[index1][1],
                segment_type="linker",
                pair_type="1-4"
            )
            pair_information += temp_string

    pair_information += ";    1-5 pairs\n"

    for index2 in range(len(segments_list)): # this is 0-indexed
        if segments_list[index2] == "coil": # but the reference to positions_list values are 1-indexed
            temp_string = pairs_helper(
                positions_list[index2][0],
                positions_list[index2][1],
                segment_type="coil",
                pair_type="1-5"
            )
            pair_information += temp_string
        else:
            temp_string = pairs_helper(
                positions_list[index2][0],
                positions_list[index2][1],
                segment_type="linker",
                pair_type="1-5"
            )
            pair_information += temp_string


    # now handle the Exclusions
    exclusions_information = f"""

[ exclusions ]
; Exclude Atype-Atype interactions within individual coil segments to prevent collapse
; Also exclude oligo-oligo interactions within individual coil segments for same reason
; specificity is controlled through [nonbond_params] in the .top file
"""
    coil_num = 1
    for segment_i in range(len(A_BEADS_BY_COIL)):
        exclusion_string = f"; exclusions for Coil segment: {coil_num}\n"
        current_A_list = A_BEADS_BY_COIL[segment_i]
        current_EG_list = EG_BEADS_BY_COIL[segment_i]
        for A_i in range(len(current_A_list)):
            if A_i != len(current_A_list)-1:
                for a_i_i in current_A_list[A_i:]:
                    exclusion_string += f"{a_i_i} "
            else:
                pass
            exclusion_string += "\n"
        for EG_i in range(len(current_EG_list)):
            if EG_i != len(current_EG_list) - 1:
                for eg_i_i in current_EG_list[EG_i:]:
                    exclusion_string += f"{eg_i_i} "
            else:
                pass
            exclusion_string += "\n"
        exclusion_string += "\n"
        exclusions_information += exclusion_string
        coil_num += 1

# This code is for generating exclusions so that NO coil segment within a single molecule can interact. I am commenting
# this out as of 2022 06 10
#     A_BEADS.sort()
#     for ai in range(len(A_BEADS)):  # this is handled by 0-indexed
#         exclusion_string = f""
#         if ai != len(A_BEADS)-1:
#             for tai in range(ai, len(A_BEADS)):
#                 exclusion_string += f"{A_BEADS[tai]} "
#             exclusion_string += "\n"
#             exclusions_information += exclusion_string
#         else:
#             pass


    # Now time to handle the angles
    angle_information = f"""

[ angles ]
;i    j    k    func   theta0 (deg)    Ktheta (kJ/mol/rad^2)
"""
    for i3 in range(len(segments_list)):    # this is 0-indexed
        if segments_list[i3] == "coil":     # the values in positions_list are 1-indexed
            coil_start = positions_list[i3][0]
            coil_stop = positions_list[i3][1]
            temp_string = angle_helper(
                coil_start,
                coil_stop,
                segment_type="coil"
            )
            angle_information += temp_string
            # now I need to handle the manual angles (2 need to be specified)
            if i3 == (len(segments_list)-1):
                continue
            else:
                angle_information += f" {coil_stop - 1}    {coil_stop}    {coil_stop + 1}     1       {THETA_HELIX}        {THETA_LINK_FORCE}  ; angle theta^'_CL\n"
                angle_information += f" {coil_stop}    {coil_stop + 1}    {coil_stop + 2}     1       {THETA_HELIX}        {THETA_LINK_FORCE}  ; angle theta^''_CL\n"
        else:
            linker_start = positions_list[i3][0]
            linker_stop = positions_list[i3][1]
            temp_string = angle_helper(
                linker_start,
                linker_stop,
                segment_type="linker"
            )
            angle_information += temp_string
            # now I need to handle the manual angles (2 need to be specified)
            if i3 == (len(segments_list)-1):    # stop if I'm at the end of the model obviously
                continue
            else:
                angle_information += f" {linker_stop - 1}    {linker_stop}    {linker_stop + 1}     1       {THETA_HELIX}        {THETA_LINK_FORCE}  ; angle theta^'_LC\n"
                angle_information += f" {linker_stop}    {linker_stop + 1}    {linker_stop + 2}     1       {THETA_HELIX}        {THETA_LINK_FORCE}  ; angle theta^''_LC\n"


    # now it's time to handle the dihedrals!
    dihedral_information = f"""

[ dihedrals ]
;i   j   k   l   func   phi0(deg)   kb (kJ/mol)    mult.
;    multiplicity = 1
"""
    for i4 in range(len(segments_list)):
        if segments_list[i4] == "coil":
            coil_d_start = positions_list[i4][0]
            coil_d_stop = positions_list[i4][1]
            temp_string = dihedral_helper(
                coil_d_start,
                coil_d_stop,
                1,
                segment_type="coil",
            )
            dihedral_information += temp_string
            # now I need to handle the 3 manual torsions
            if i4 == (len(segments_list)-1):
                continue
            else:
                dihedral_information += f" {coil_d_stop - 2}   {coil_d_stop - 1}   {coil_d_stop}   {coil_d_stop + 1}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^0_C\n"
                dihedral_information += f" {coil_d_stop - 1}   {coil_d_stop}   {coil_d_stop + 1}   {coil_d_stop + 2}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^'_CL\n"
                dihedral_information += f" {coil_d_stop}   {coil_d_stop + 1}   {coil_d_stop + 2}   {coil_d_stop + 3}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^''_CL\n"
        else:
            link_d_start = positions_list[i4][0]
            link_d_stop = positions_list[i4][1]
            temp_string = dihedral_helper(
                link_d_start,
                link_d_stop,
                1,
                segment_type="linker",
            )
            dihedral_information += temp_string
            # now I need to handle the 3 manual torsions
            if i4 == (len(segments_list)-1):
                continue
            else:
                dihedral_information += f" {link_d_stop - 2}   {link_d_stop - 1}   {link_d_stop}   {link_d_stop + 1}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^0_C\n"
                dihedral_information += f" {link_d_stop - 1}   {link_d_stop}   {link_d_stop + 1}   {link_d_stop + 2}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^'_CL\n"
                dihedral_information += f" {link_d_stop}   {link_d_stop + 1}   {link_d_stop + 2}   {link_d_stop + 3}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      1  ; angle alpha^''_CL\n"

    dihedral_information += f";    multiplicity = 3\n"

    for i4 in range(len(segments_list)):
        if segments_list[i4] == "coil":
            coil_d_start = positions_list[i4][0]
            coil_d_stop = positions_list[i4][1]
            temp_string = dihedral_helper(
                coil_d_start,
                coil_d_stop,
                3,
                segment_type="coil",
            )
            dihedral_information += temp_string
            # now I need to handle the 3 manual torsions
            if i4 == (len(segments_list)-1):
                continue
            else:
                dihedral_information += f" {coil_d_stop - 2}   {coil_d_stop - 1}   {coil_d_stop}   {coil_d_stop + 1}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^0_C\n"
                dihedral_information += f" {coil_d_stop - 1}   {coil_d_stop}   {coil_d_stop + 1}   {coil_d_stop + 2}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^'_CL\n"
                dihedral_information += f" {coil_d_stop}   {coil_d_stop + 1}   {coil_d_stop + 2}   {coil_d_stop + 3}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^''_CL\n"
        else:
            link_d_start = positions_list[i4][0]
            link_d_stop = positions_list[i4][1]
            temp_string = dihedral_helper(
                link_d_start,
                link_d_stop,
                3,
                segment_type="linker",
            )
            dihedral_information += temp_string
            # now I need to handle the 3 manual torsions
            if i4 == (len(segments_list)-1):
                continue
            else:
                dihedral_information += f" {link_d_stop - 2}   {link_d_stop - 1}   {link_d_stop}   {link_d_stop + 1}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^0_C\n"
                dihedral_information += f" {link_d_stop - 1}   {link_d_stop}   {link_d_stop + 1}   {link_d_stop + 2}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^'_CL\n"
                dihedral_information += f" {link_d_stop}   {link_d_stop + 1}   {link_d_stop + 2}   {link_d_stop + 3}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      3  ; angle alpha^''_CL\n"


    # NOW -- BIG FINALE -- write all the parts to a file
    with open(itp_filename, 'w') as f:
        f.write(opener)
        f.write(atom_information)
        f.write(bond_information)
        f.write(pair_information)
        f.write(exclusions_information)
        f.write(angle_information)
        f.write(dihedral_information)

    return None


def atom_helper(aStart, aStop, segment_type="linker", heptad_start=0, oligo=None, oligoID=0):
    # arguments provided as 1-indexed
    """

    :param aStart: the start (1-indexed) index of the segment in question
    :param aStop: the end (1-indexed) index of the segment in question
    :param segment_type: the switch for either linker or coil segment, behavior changes based on segment
    :param heptad_start: the location of the FIRST A bead in the coil segment, only has meaning for segment_type="coil"
    :param oligo:
    :param oligoID:
    :return: a complete string of every atom in the segment provided, bounded inclusive by aStart and aStop.
    Also return the list a_posits and d_posits for the given segment, if a coil, which includes the indices of the A and
    D position A-beads of the segment (1-indexed)
    """
    atom_string = ""
    OLIGO = ""
    a_posits = []           # stores A position indices by 1-indexing
    d_posits = []           # stores D position indices by 1-indexing
    e_posits = []           # stores E position indices by 1-indexing
    g_posits = []           # stores G position indices by 1-indexing

    if oligo == "dim":
        OLIGO = "Di "

    if oligo == "tri":
        OLIGO = "Tri"

    if oligo == "tet":
        OLIGO = "Tet"

    if oligoID == "" or oligoID is None:
        oligoID = ""

    if segment_type == "linker":    # control for making linker information
        for posl in range(aStart, aStop+1):  # 1-indexing
            atom_string += f" {posl}   B     {posl}    GLY     CA    {posl}    0.000    {MASS}\n"
        return atom_string
    else:                           # control for making coil information
        a_posits.append(heptad_start)
        while len(a_posits) < 200:
            new_pos = a_posits[-1] + 7
            if new_pos <= aStop:        # uses 1-indexing
                a_posits.append(new_pos)
            else:
                break
        d_posits.append(heptad_start + 3)
        while len(d_posits) < 200:
            new_pos = d_posits[-1] + 7
            if new_pos <= aStop:        # uses 1-indexing
                d_posits.append(new_pos)
            else:
                break
        e_posits.append(heptad_start + 4)
        while len(e_posits) < 200:
            new_pos = e_posits[-1] + 7
            if new_pos <= aStop:        # uses 1-indexing
                e_posits.append(new_pos)
            else:
                break
        g_posits.append(heptad_start + 6)
        while len(g_posits) < 200:
            new_pos = g_posits[-1] + 7
            if new_pos <= aStop:        # uses 1-indexing
                g_posits.append(new_pos)
            else:
                break

        for posc in range(aStart, aStop+1):        # uses 1-indexing
            if posc in a_posits:
                atom_string += f" {posc}   A{oligo}{oligoID}    {posc}    ILE     CA    {posc}    0.000    {MASS}\n"
            elif posc in d_posits:
                atom_string += f" {posc}   A{oligo}{oligoID}    {posc}    ILE     CA    {posc}    0.000    {MASS}\n"
            elif posc in e_posits:
                atom_string += f" {posc}   {OLIGO}   {posc}    LYS     CA    {posc}    0.000    {MASS}\n"
            elif posc in g_posits:
                atom_string += f" {posc}   {OLIGO}   {posc}    GLU     CA    {posc}    0.000    {MASS}\n"
            else:
                atom_string += f" {posc}   B     {posc}    ALA     CA    {posc}    0.000    {MASS}\n"

        return atom_string, a_posits, d_posits, e_posits, g_posits

def pairs_helper(pStart, pStop, segment_type="linker", pair_type="1-4"):
    """

    :param pStart:
    :param pStop:
    :param segment_type:
    :param pair_type:
    :return:
    """
    pairs_string = ""
    pairs_list = []
    # This part handles if the pairs are 1-4 pairs or not
    if pair_type == "1-4":
        for i in range(pStart, pStop-2):        # this is 1-indexed!
            pairs_list.append([i, i+3])
        if segment_type == "coil":
            for i1 in range(len(pairs_list)):
                pairs_string += f" {pairs_list[i1][0]}   {pairs_list[i1][1]}   1     {COIL_P14_C6}     {COIL_P14_C12}\n"
        else:
            for i1 in range(len(pairs_list)):
                pairs_string += f" {pairs_list[i1][0]}   {pairs_list[i1][1]}   1     {LINK_P14_C6:.10f}     {LINK_P14_C12:.10f}\n"
        return pairs_string
    elif pair_type == "1-5":
        for i in range(pStart, pStop-3):        # this is 1-indexed!
            pairs_list.append([i, i+4])
        if segment_type == "coil":
            for i1 in range(len(pairs_list)):
                pairs_string += f" {pairs_list[i1][0]}   {pairs_list[i1][1]}   1     {COIL_P15_C6}     {COIL_P15_C12}\n"
        else:
            for i1 in range(len(pairs_list)):
                pairs_string += f" {pairs_list[i1][0]}   {pairs_list[i1][1]}   1     {LINK_P15_C6:.10f}     {LINK_P15_C12:.10f}\n"
        return pairs_string
    else:
        print("The pair_type specifer doesn't match any expected pairs for the coil model. There must be a problem.")
        return None

def angle_helper(anStart, anStop, segment_type="linker"):
    """

    :param anStart:
    :param anStop:
    :param segment_type:
    :return:
    """
    angle_string = ""
    angle_list = []
    for i in range(anStart, anStop-1):
        angle_list.append([i, i+1, i+2])
    if segment_type == "coil":
        for j in range(len(angle_list)):
            angle_string += f" {angle_list[j][0]}    {angle_list[j][1]}    {angle_list[j][2]}     1       {THETA_HELIX}        {THETA_HELIX_FORCE}\n"
        return angle_string
    else:
        for j in range(len(angle_list)):
            angle_string += f" {angle_list[j][0]}    {angle_list[j][1]}    {angle_list[j][2]}     1       {THETA_HELIX}        {THETA_LINK_FORCE}\n"
        return angle_string

def dihedral_helper(dStart, dStop, multiplicity, segment_type="linker"):
    """

    :param dStart:
    :param dStop:
    :param multiplicity:
    :param segment_type:
    :return:
    """
    dihedral_string = ""
    dihedral_list = []
    for i in range(dStart, dStop-2):        # 1-indexed
        dihedral_list.append([i, i+1, i+2, i+3])
    if segment_type == "coil":
        if multiplicity == 1:
            for k in range(len(dihedral_list)):
                dihedral_string += f" {dihedral_list[k][0]}   {dihedral_list[k][1]}   {dihedral_list[k][2]}   {dihedral_list[k][3]}    1     {ALPHA1_HELIX}        {ALPHA1_HELIX_FORCE}      {multiplicity}\n"
            return dihedral_string
        if multiplicity == 3:
            for k in range(len(dihedral_list)):
                dihedral_string += f" {dihedral_list[k][0]}   {dihedral_list[k][1]}   {dihedral_list[k][2]}   {dihedral_list[k][3]}    1     {ALPHA3_HELIX}        {ALPHA3_HELIX_FORCE}      {multiplicity}\n"
            return dihedral_string
    else:   # linker segments
        if multiplicity == 1:
            for k in range(len(dihedral_list)):
                dihedral_string += f" {dihedral_list[k][0]}   {dihedral_list[k][1]}   {dihedral_list[k][2]}   {dihedral_list[k][3]}    1     {ALPHA1_HELIX}        {ALPHA1_LINK_FORCE}      {multiplicity}\n"
            return dihedral_string
        if multiplicity == 3:
            for k in range(len(dihedral_list)):
                dihedral_string += f" {dihedral_list[k][0]}   {dihedral_list[k][1]}   {dihedral_list[k][2]}   {dihedral_list[k][3]}    1     {ALPHA3_HELIX}        {ALPHA3_LINK_FORCE}      {multiplicity}\n"
            return dihedral_string

if __name__ == "__main__":
    # I don't anticipate needing to use the functions in this script elsewhere, so the
    # name==main style is a little overkill, but I'm including it here for good pythonic practice

    # Set up the argument parser!
    parser = argparse.ArgumentParser(description="Topology Writer Tool -- Use this to make topology files for CC "
                                                 "proteins, easy as pie! Cheap as Ubik, too!")
    parser.add_argument("-df", help="Design file: this contains the information on how to build the protein, with"
                                    " alternating coil and linker segments.")
    parser.add_argument("-n", help="The number of beads in the CC protein.", type=int, required=True)
    parser.add_argument("-itp", help="Activate to turn on .itp file writer.",
                        action="store_true")
    parser.add_argument("-itp_filename", help="File name for the .itp file. Please provide extension.")

    # parse the arguments
    args = parser.parse_args()

    # set up the list of segments, in order, and their corresponding indices, in order
    # Segments contains the names of each of the segments as they appear in order in the provided file
    # Positions contains the start/stop (and heptad) indices for the segments that correspond to the same index in the
    # Segments list
    Segments = []; Positions = []
    with open(args.df, "r", encoding="utf-8-sig") as file:
        for line in file:
            if line[0] == "#":
                # I am now adding in comment support. Comments can exist in the design file and can be denoted
                # as "#"
                continue
            else:
                line_split = line.rstrip("\n").split(",")
                Segments.append(str(line_split[0]))
                if str(line_split[0]) == "coil":   # uses 0-indexing, this now handles the final column being oligo type
                    Positions.append([int(line_split[1]), int(line_split[2]), int(line_split[3]),
                                    str(line_split[4]), line_split[5]]) # notice that index 5 has no type command!!!
                else:
                    Positions.append([int(line_split[1]), int(line_split[2])])
    if args.itp:
        write_itp(args.itp_filename, args.n, Segments, Positions)
