"""
@Title:             multimerization_analysis.py
@Name:              Mando A Ramirez
@Date:              2023 04 06

@Description:       This script is an extension and improvement upon its related sister code,
"coil_multimer_analysis.py".

The difference between this and the other script is that this script calculates distances and 
neighbors using the **centers of mass** of each coil's A beads. Additionally, the distance 
calculations are "vectorized", meaning they happen all at once and I then figure out neighbors
after that.


This script analyzes the multimeric state of coil simulations to produce statistics about the
different types of multimers and higher order assemblies that form with coils in a box. I am writing this script with
slab simulations in mind so I can track the types of multimers that form through simulation time. I will be calcualte
the following:
    - the types of each multimers through time (per-frame)
    - the distribution of each multimers for the whole simulation
    - the number-per-frame and distribution of self-vs-other interactions
    - the number-per-frame and distribution of free coils for the simulation

This script will first be written to handle only homomeric simulations of models, but can easily be modified to handle
heteromeric simulations.

@Updates:
2023 04 28 -- added a new feature to save the contacts in every analyzed frame. This will output the indices of coils
that are determined to be interacting, and save them to a .txt file. For a typical 10 us simulation, with every single
frame analyzed, the expected output file is somewhere around 10-20 MB, which isn't terrible.
"""

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os
import pandas as pd
import time

"""
Function definitions
"""
def calc_COM(trajF, POSI, mass):
    """
    :param traj:    a *single frame* from a trajectory. DO NOT pass the entire trajectory
    :param POSI:    the list of lists of coil positions, where each entry is the beads for 
        a single coil that the COM will be calculated for
    :param mass:    the mass used for ALL beads for the COM calculation
    :return COMS:   array containing the COM coordinates for all of the coils 
        enumerated in POSI
    """
    COMS = []
    
    for coilI, coilP in enumerate(POSI):
        # coilI is the index of a given coil, same index for POSI and COMS
        # coilP are the 0-indexed indices for the coil beads
        coords = trajF.xyz[:,coilP,:][0]
        coordsT = np.transpose(coords)
        # coords is a Nx3 array, where each row corresponds to the xyz coords for a bead
        # in the coilP
        # I transpose this to 3xN so that each row becomes one of the coordinate dimensions
        # for ALL The beads in the coil
        
        xCOM = float( mass*(np.sum(coordsT[0])) / (mass * len(coordsT[0])) )
        yCOM = float( mass*(np.sum(coordsT[1])) / (mass * len(coordsT[1])) )
        zCOM = float( mass*(np.sum(coordsT[2])) / (mass * len(coordsT[2])) )
        
        COMS.append(np.array([xCOM, yCOM, zCOM]))
    
    return np.array(COMS)


""" Set up the arg parser """
parser = argparse.ArgumentParser(description="Coil multimerization tool - this determines the types of multimers and"
                                             " higher order assemblies present in your coil simulations!")
parser.add_argument("-df", help="Definition file: this contains the information on how to build the model, with"
                                " alternating coil and linker segments.")
parser.add_argument("-nmodels", help="The total number of models in the simulation", required=True, type=int)
parser.add_argument("-N", help="The total number of coil segments in the simulation. Necessary for the distribution"
                               " plotting.", required=True, type=int)
parser.add_argument("-l", help="No. of beads in an ENTIRE SINGLE PROTEIN. This is not the size of an individual coil. This is "
                               "the size of an entire individual protein. As of right now, ONLY proteins that have the"
                               " same no. of beads can be analyzed.",
                    required=True, type=int)
parser.add_argument("-mass", help="The mass of the CG beads. Currently, all beads need to be "
    " the same mass.", type=float)
parser.add_argument("-cutoff", help="The cut-off distance used to determine if two beads are close enough to be in "
                                    "an multimer. [nm]", default=1.3, type=float)
parser.add_argument("-t", help="trajectory file with extension (.xtc)", required=True)
parser.add_argument("-p", help="toplogy file i.e. PDB file with extension", required=True)
parser.add_argument("-start", help="the start TIME [ps] to begin the analysis. Default = 0 [ps]", default=0,
                    type=int)
parser.add_argument("-stop", help="the stop TIME [ps] to end the analysis. Default is the entire simulation", type=int)
parser.add_argument("-f", help="the stride of analysis i.e. analyze every fth frame, default = 1", default=1, type=int)
parser.add_argument("-single", help="a switch to turn on if the analysis is only for a single structure file",
    action="store_true", default=False)
parser.add_argument("-name", help="name that you'd like to add to the analysis output files",
                    default="multimerAnalysis")
parser.add_argument("-timeseries", help="pass the flag to turn on plotting of multimer data as a time series",
                    action="store_true", default=False)
parser.add_argument("-distro", help="pass this flag to turn on plotting of average multimer population, as a "
                                    "distribution", action="store_true", default=False)
parser.add_argument("-verbose", help="pass this flag to output extra information to stdout.", action="store_true",
                    default=False)
parser.add_argument("-save", help="pass this flag to save the 'contact map', which is just a list of coil indices"
                                  "of every multimer, in every analyzed frame. Output file is .txt.",
                    action="store_true", default=False)

args = parser.parse_args()

trajectory= args.t
topology = args.p

""" Print a warning statement """
print("**********")
print("PLEASE BE ADVISED: This code automatically converts the information provided in the -df file into 0-INDEXED "
      "positions! Your input is 1-indexed, but this code will automatically change it to 0-indexed."
      " "
      "Also note that this code only needs to run once, since it outputs the data and then you can modify plots of that"
      " data without doing the entire analysis all over again.")
print("**********\n")



""" 
Begin the timing procedure, added 20230221 
"""
time_start = time.perf_counter()



"""
Load the trajectory!
"""
if args.single:
    traj = md.load(trajectory, top=topology)
    args.f = 1
    args.timeseries = False     # can't do timeseries for just 1 frame!
        # this also adds a user sanity check to make sure that the code doesn't 
        
        # try to do something that doesn't make sense.
    remaining_frames = 1
    analyzed_frames = 1
else:
    traj_load = md.load(trajectory, top=topology)
    frame_start = 0
    while frame_start <= traj_load.n_frames:
        if int(frame_start * traj_load.timestep) >= args.start:
            break
        else:
            frame_start += 1
    frame_end = 0
    if args.stop is None:
        frame_end = traj_load.n_frames
    else:
        while frame_end <= traj_load.n_frames:
            if int(frame_end * traj_load.timestep) >= args.stop:
                break
            else:
                frame_end += 1
    if frame_end < frame_start:
        raise ValueError("The specified end point is before the start point! Exiting because this doesn't make sense")
        exit(1)
    traj = traj_load[frame_start:frame_end:]
    remaining_frames = traj_load.n_frames - frame_start
    analyzed_frames = remaining_frames / args.f



"""
Useful global constants
"""
CUTOFF = args.cutoff                 # in nanometers!
NMODELS = args.nmodels
PSIZE = args.l
MASS = args.mass

if args.verbose:
    print("Global constants used in analysis:")
    print(f"Number of proteins (models): {NMODELS}")
    print(f"Length of each protein: {PSIZE} beads")
    print(f"Mass of each CG bead: {MASS} amu")
    print(f"Cutoff length used to define a neighbor: {CUTOFF} nm")
    print()



"""
Step 1 - read the file about the model description to find where all the coils and A beads are.
code taken from other things like model_maker.py for file reading
"""
# Read the file and save only the information about the coils. Save the start/end indices for the coils and the start
# of the A-bead positions.
# POSITIONS ARE 0-INDEXED
Coil_positions_df = []
coil_counter = 0                # this is how many coils are present in a model
with open(args.df, "r") as file:
    for line in file:
        line_split = line.rstrip("\n").split(",")
        if str(line_split[0]) == "coil":
            # information in line_split = [coil_start, coil_end, a_bead_start]
            # information originally 1-indexed
            # this information needs to be 0-indexed, so I'm converting it!!
            coil_counter += 1
            Coil_positions_df.append([int(line_split[1])-1, int(line_split[2])-1, int(line_split[3])-1])
        else:
            pass

print(f"Coil model file: {args.df}")
print(f"Parsed coil positions for a single coil (0-indexed): {Coil_positions_df}")


# Part of the sanity check ==> check that the size of each protein provided in the args
# matches what is in the protein model file.
prot_size = Coil_positions_df[-1][1]
try:
    (prot_size + 1) == PSIZE
except:
    print("The protein size (no. beads) does not match the size in the file. Exiting now!")
    print(f"Provided protein size in args: {PSIZE}")
    print(f"Protein size from model file: {prot_size + 1}")
    exit(1)

# Now, for each coil, generate the A- and D- position A-beads FOR THE FIRST MODEL ONLY. I need to enumerate all the
# a-bead positions first before I can generate all the a-beads for all the models
# ALL POSITIONS ARE 0-INDEXED
Coil1_A_beads = []
for i in range(len(Coil_positions_df)):
    segStart = Coil_positions_df[i][0]
    segStop = Coil_positions_df[i][1]
    a_posits = []
    d_posits = []
    a_index = Coil_positions_df[i][2]
    d_index = Coil_positions_df[i][2] + 3
    while a_index <= segStop:
        a_posits.append(a_index)
        a_index += 7
    while d_index <= segStop:
        d_posits.append(d_index)
        d_index += 7
    coil_i_beads = a_posits + d_posits
    coil_i_beads.sort()
    Coil1_A_beads.append(coil_i_beads)

print(f"Enumerated coil beads for the FIRST coil segment (0-indexed): {Coil1_A_beads}")
print()

# Now enumerate all the a-bead indices for each model present in the simulation. The variable MODELS_A_BEADS will have
# the 0-BASED INDICES (!!!!) for every coil model in the system. THIS IS FOR ALL COILS
MODELS_A_BEADS = []
MODEL_RANGES = []
for M in range(0, NMODELS):
    # This information stores the 0-INDEXED (!!!) ranges for each of the models in the system!
    model_additive = M * PSIZE
    model_start = 0 + (M * PSIZE)
    model_end = PSIZE + (M * PSIZE)
    MODEL_RANGES.append(range(model_start, model_end))
    # then I store the information about the coil a-bead indices (0-INDEXED!!)
    for c1i in Coil1_A_beads:       # select out each coil from Coil1 in Coil1_A_beads
        temp_abead_posits = []
        for c1i_num in c1i:         # loop through every a-bead position in coil1i selected above
            temp_abead_posits.append(c1i_num + model_additive)          # add the model addifier to the number and put
                    # all those numbers in to the temporary list to add to the final container of every coil
        MODELS_A_BEADS.append(temp_abead_posits)

if args.verbose:
    print(f"All model coils a beads (0-indexed): {MODELS_A_BEADS}")
    print(f"Model ranges (0-indexed): {MODEL_RANGES}")
    print()
# now we have all the index information we need to start doing the analysis


"""
Step 2 - do the trajectory analysis now! I will do the analysis BY FRAME
"""
"""
Here is a general overview of the analysis algorithm
set analysis variables
Iterate through every frame {
    - set the frame variables = 0
    - calculate the COMs for all the coils based on new positions
    - calculate distances between all COM points (all to all) using np.linalg.norm
    iterate through the symmetric matrix to find neighbors {
        (each array in the symmetric matrix corresponds to distances relative to a single coil)
        * find neighbors/indices using (array < CUTOFF).nonzero()
        * figure out what kind of multimer exists
        * update frame variables
        * store data about coils analyzed so I can skip those arrays and not double count
    }
    - add frame variables to master variables
    loop through again
}
"""

if args.save:
    output_file_name = f"{args.name}_contacts.txt"
    # open and turn on the output file. BUT check to make sure the file doesn't already exist because this
    # method only works by appending.
    if os.path.isfile(f"./{output_file_name}"):
        raise FileExistsError("Desired contact output file already exists. Exiting now to prevent overwriting!")
    else:
        outputFile = open(output_file_name, 'a')

# Set ANALYSIS variables
# These are the TIME SERIES/Frame data for multimers and self-vs-other interactions. 
# Each entry will correspond to the total of that category, in order of frames
FREE_COILS = []         #
DIMERS = []; TRIMERS = []; TETRAMERS = []
PENTAMERS = []; HEXAMERS = []; HIGHER = []
SELF = []; OTHER = []

# Begin the frame loop, analyze trajectory frame by frame
for frame in range(0, traj.n_frames, args.f):
    free_coil_counter = 0
    frame_dimer = 0; frame_trimer = 0; frame_tetramer = 0; 
    frame_pentamer = 0; frame_hexamer = 0; frame_higherorder = 0
    frame_self = 0; frame_other = 0
    
    # Calculate the COMs for each of the coils for the given frame
    coms = calc_COM(traj[frame], MODELS_A_BEADS, MASS)

    # now calculate the distance between ALL COMs
    # this came from https://stackoverflow.com/questions/46700326/calculate-distances-between-one-point-in-matrix-from-all-other-points
    # the first answer on the page
    dists = np.linalg.norm(coms - coms[:, None], axis=-1)

    # Loop through every coil and filter based on distances to find neighbors, the type of
    # multimer that is formed, and the indices of neighbors
    # I will also tally all the coils that are included in a multimer so that I don't 
    # double count anything
    # find the multimers, add them to the frame counters, then move onto the next frame
    excluded_coils = []
    for coili, coilD in enumerate(dists):
        if coili in excluded_coils:
            continue
        neighbors = (coilD <= CUTOFF).nonzero()[0]    # this contains the indices of other coils
                # that are neighbors, out of the total number of coils in the simulation
        multimer_size = len(neighbors)

        if args.save:
            if multimer_size > 1:
                np.savetxt(outputFile, neighbors, delimiter=" ", fmt="%i", newline=" ")
                outputFile.write(" ")

        if multimer_size == 1:
            free_coil_counter += 1
        elif multimer_size == 2:
            frame_dimer += 1
        elif multimer_size == 3:
            frame_trimer += 1
        elif multimer_size == 4:
            frame_tetramer += 1
        elif multimer_size == 5:
            frame_pentamer += 1
        elif multimer_size == 6:
            frame_hexamer += 1
        else:
            frame_higherorder += 1
        
        for indices in neighbors:
            excluded_coils.append(indices)

    if args.save:
        outputFile.write("\n")
    # Now add the frame values to the ANALYSIS variables and move onto the next frame
    # This CANNOT be a running total and instead there needs to be an entry for every 
    # frame analyzed
    FREE_COILS.append(free_coil_counter)
    DIMERS.append(frame_dimer)
    TRIMERS.append(frame_trimer)
    TETRAMERS.append(frame_tetramer)
    PENTAMERS.append(frame_pentamer)
    HEXAMERS.append(frame_hexamer)
    HIGHER.append(frame_higherorder)
    SELF.append(frame_self)
    OTHER.append(frame_other)

if args.save:
    outputFile.close()
    
""" End the timing work. This is because the majority of the time will happen above with the 
actual calculation of the haystack """
time_stop = time.perf_counter()
elapsed_time = ((time_stop - time_start) / 60) / 60 
# this reports the elapsed time in terms of hours
print("-----------------------------")
print(f"Elapsed time for the analysis fo neighbor searching and multimer analysis: {elapsed_time:.4f} hours")
print("-----------------------------")

# Calculate normalized counts (normalized to number of coil segments in a box, averaged over
#   all frames)
free_average = np.average(FREE_COILS); free_std = np.std(FREE_COILS)
norm_free = (free_average * 1.) / float(args.N)
norm_free_std = (free_std * 1.) / float(args.N)
dimer_average = np.average(DIMERS); dimer_std = np.std(DIMERS)
norm_dimer = (dimer_average * 2.) / float(args.N)
norm_dimer_std = (dimer_std * 2.) / float(args.N)
trimer_average = np.average(TRIMERS); trimer_std = np.std(TRIMERS)
norm_trimer = (trimer_average * 3.) / float(args.N)
norm_trimer_std = (trimer_std * 3.) / float(args.N)
tet_average = np.average(TETRAMERS); tet_std = np.std(TETRAMERS)
norm_tet = (tet_average * 4.) / float(args.N)
norm_tet_std = (tet_std * 4.) / float(args.N)
pent_average = np.average(PENTAMERS); pent_std = np.std(PENTAMERS)
norm_pent = (pent_average * 5.) / float(args.N)
norm_pent_std = (pent_std * 5.) / float(args.N)
hex_average = np.average(HEXAMERS); hex_std = np.std(HEXAMERS)
norm_hex = (hex_average * 6.) / float(args.N)
norm_hex_std = (hex_std * 6.) / float(args.N)


# Handle the plotting of important information
if args.distro:
    fig, ax = plt.subplots()
    ax.errorbar([1, 2, 3, 4, 5, 6], 
        [norm_free, norm_dimer, norm_trimer, norm_tet, norm_pent, norm_hex],
        yerr=[norm_free_std, norm_dimer_std, norm_trimer_std, norm_tet_std, norm_pent_std, norm_hex_std], 
        color="black", linestyle="",
        markersize=10, marker=".", ecolor="black", elinewidth=1, capsize=5)
    plt.xticks(np.arange(0, 7))
    plt.xlim(0, 6)
    plt.ylim(0, 1)
    plt.xlabel("N-mer")
    plt.ylabel("Counts")
    plt.grid(color="black", linestyle=":", alpha=0.5)
    plt.savefig(f"{args.name}_multimeric_analysis_distribution.png", dpi=600)

if args.verbose:
    print(f"Sum of normalized coil counts across multimers (should equal 1):"
        f" {np.sum(np.array([norm_free, norm_dimer, norm_trimer, norm_tet, norm_pent, norm_hex])):.4f}")

if args.single:
    arFrames = np.array([0])
else:
    arFrames = np.arange(0+args.start, (traj.n_frames*traj.timestep)+args.start, (args.f*traj.timestep))
if args.timeseries:
    fig, ax = plt.subplots(figsize=(10,5))
    ax.plot(arFrames, np.array(FREE_COILS), color="blue", label = "Free coils",
              linewidth=1)
    ax.plot(arFrames, np.array(DIMERS), color="slategrey", label = "Dimers",
             linewidth=1)
    ax.plot(arFrames, np.array(TRIMERS), color="darkorange", label = "Trimers",
             linewidth=1)
    ax.plot(arFrames, np.array(TETRAMERS), color="blueviolet", label = "Tetramers",
              linewidth=1)
    ax.plot(arFrames, np.array(PENTAMERS), color="seagreen", label = "Pentamers",
              linewidth=1)
    ax.plot(arFrames, np.array(TETRAMERS), color="goldenrod", label = "Hexamers",
              linewidth=1)
    ax.plot(arFrames, np.array(HIGHER), color="red", label = "Higher order oligos*",
              linewidth=1)
    plt.legend(loc="center left", bbox_to_anchor=(1.01, 0.5))
    plt.grid(alpha=0.5)
    plt.ylabel("Counts")
    plt.xlabel("Simulation time (ps)")
    plt.tight_layout()
    plt.savefig(f"{args.name}_multimeric_analysis_timeseries.png", dpi=600)
    plt.close()

# as of 2023 04 07 -- analysis of self to other coils is not included in the analysis
#     fig, ax = plt.subplots()
#     ax.plot(arFrames, np.array(SELF), color="black", label = "Intra-model interactions",
#               linewidth=1)
#     ax.plot(arFrames, np.array(OTHER), color="blue", label = "Inter-model interactions",
#              linewidth=1)
#     plt.legend()
#     plt.grid(alpha=0.5)
#     plt.ylabel("Counts")
#     plt.xlabel("Simulation time (ps)")
#     plt.savefig(f"{args.name}_selfvsother_analysis.png", dpi=600)
#     plt.close()


# Now save out all the data so I can use it later:
outputdata = [arFrames, np.array(FREE_COILS), np.array(DIMERS), 
    np.array(TRIMERS), np.array(TETRAMERS), np.array(PENTAMERS),
    np.array(HEXAMERS), np.array(HIGHER), 
    np.array(SELF), np.array(OTHER)]
pd.DataFrame(np.array(outputdata)).to_csv(f"{args.name}_multimerAnalysis_outData.csv", header=False)

print(f"Number of frames in trajectory after slicing from the starting point: {remaining_frames}")
print(f"The total number of analyzed frames, based on the '-f' argument: {analyzed_frames}")
if args.verbose:
    print("""
The output file is organized by rows, as follows:
Row 0 - time (ps)
Row 1 - Number of monomers
Row 2 - Number of dimers
Row 3 - Number of trimers
Row 4 - Number of tetramers
Row 5 - Number of pentamers
Row 6 - Number of hexamers
Row 7 - Number of higher order species
Row 8 - Number of self-protein interactions (not counted as of 2023 04 07)
Row 9 - Number of other-protein interactions (not counted as of 2023 04 07)
""")
if args.save:
    print("You have elected to save out the list of coil contacts in found multimers, which saves the indices of all"
          " interacting coils (0-indexed) to an output file. Each line in the file corresponds to the contacts found"
          " in a given frame, and each pair/grouping of contacts is separated by a double space. Use this fact to"
          " parse the output file for further analysis.")