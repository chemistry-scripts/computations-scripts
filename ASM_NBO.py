import argparse
import os
import sys
from scm.plams import KFReader, Settings, SingleJob, init


def main():
    """
        Main interface that reads the file and retrieve all useful geoms etc.
        Then rewrite geometries in a usable formatter_class
        Prep all NBO computations
        Run ADF on each computation
    """
    # Retrieve command line values
    args = get_input_arguments()
    inputFile = args['inputfile']

    # Init PLAMS, Setting up a plams.$PID directory in the working dir
    init()

    # Setup blank Settings, to be updated upon reading the IRC
    settings = Settings()

    geometries = IRCCoordinatesFromt21(inputFile)
    print("Writing File")
    with open("file.xyz", mode='w+') as XYZ:
        print(XYZ.name)
        for i in range(0, len(geometries)):
            for j in range(0, len(geometries[i])):
                for k in range(0, len(geometries[i][j])):
                    XYZ.write(str(geometries[i][j][k]) + "       ")
                XYZ.write("\n")
            XYZ.write("\n\n")

    # Prep a bunch of NBO computations
    NBO_jobs = [prepareNBOComputation(geom, settings) for geom in geometries]
    # Run each job in parallel as managed by plams
    NBO_results = [job.run() for job in NBO_jobs]
    # NBO_values is a list of list of charges
    NBO_values = [extractNBOCharges(output, natoms) for output in NBO_results]
    with open("NBOCharges.data", mode='w+') as output_file:
        for i in range(0, len(NBO_values)):
            output_file.write("     ".join(NBO_values[i])) + "\n"


def IRCCoordinatesFromt21(input_file):
    """
        Extract all data from a TAPE21 file.
    """
    # Read TAPE21 and get all useful data
    t21 = KFReader(input_file)

    # Number of atoms: 7
    natoms = t21.read("Geometry", "nr of atoms")
    # atom types as indexes: [1, 2, 2, 3, 3, 3, 4]
    aatoms = t21.read('Geometry', 'fragment and atomtype index')[natoms:]
    # Atom symbols as list: ['C', 'O', 'H', 'B']
    xatoms = str(t21.read('Geometry', 'atomtype')).split()
    # Actual list of atoms as used in geometry: ['C', 'O', 'O', 'H', 'H', 'B', 'H']
    satoms = [xatoms[aatoms[order - 1] - 1] for order in t21.read('Geometry', 'atom order index')[natoms:]]

    # return [[[s, x, y, z] for s, x, y, z in zip(sAtoms, xyzBlock[0::3], xyzBlock[1::3], xyzBlock[2::3])] for xyzBlock in fwdIRC + cenIRC + bwdIRC]

    nstep_fw = t21.read("IRC_Forward", "nset")
    nstep_bw = t21.read("IRC_Backward", "nset")

    geometries_fw = t21.read("IRC_Forward", "xyz")
    geometries_bw = t21.read("IRC_Backward", "xyz")
    geometries_init = t21.read("Geometry", "xyz")

    # Reformat geometries into a long list of geometries for each step
    geometries_bw = coordinatesFromList(geometries_bw, natoms, nstep_bw)
    geometries_fw = coordinatesFromList(geometries_fw, natoms, nstep_fw)
    geometries_init = coordinatesFromList(geometries_init, natoms, 1)
    geometries_fw.reverse()

    geometries = geometries_fw + geometries_init + geometries_bw

    return [[[s, x, y, z] for s, x, y, z in zip(satoms, xyzBlock[0::3], xyzBlock[1::3], xyzBlock[2::3])] for xyzBlock in geometries]


def coordinatesFromList(coordinates_list, NAtoms, NMolecules):
    """
        Rewrites a list into a list of lists of lists:
        [x1,y1,z1,x2,y2,z2,x'1,y'1,z'1,x'2,y'2,z'2] becomes:
        [ [ [x1,y1,z1], [x2,y2,z2]], [ [x'1,y'1,z'1], [x'2,y'2,z'2] ] ]
        which can be acessed as:
        list[moleculeNumber][atomNumber][coordinateNumber]
    """
    # list = [x1,y1,z1,x2,y2,z2,x'1,y'1,z'1,x'2,y'2,z'2]
    coordinates_list = list(listChunks(coordinates_list, NAtoms))
    # [[x1,y1,z1],[x2,y2,z2],[x'1,y'1,z'1],[x'2,y'2,z'2]]
    coordinates_list = list(listChunks(coordinates_list, NMolecules))
    # [ [[x1,y1,z1],[x2,y2,z2]], [[x'1,y'1,z'1],[x'2,y'2,z'2]]]
    return coordinates_list


def listChunks(list, n):
    """
        Generator that splits a list in chunks of size n
        /!\ Does not care about the end of the list if len(list) is not a multiple of n
    """
    for i in range(0, len(list), n):
        yield list[i:i + n]


def prepareNBOComputation(geometry, runparameters):
    """
        Taking a geometry and a set of ADF parameters (Basis set, Functional, ZORA, etc)
        Add the NBO Necessary keywords and put everything in a file called filename
    """
    job = NBOJob(name="NBO Computation", settings=runparameters)
    print(job.get_runscript())


def extractNBOCharges(output, natoms):
    """
        <extract NBO Charges parsing thr outFile
    """
    # Initialize charges list
    charges = []

    with open(output, mode='r') as outFile:
        line = "Foobar line"
        while line:
            line = outFile.readline()
            if line.contains("Summary of Natural Population Analysis:"):
                # We have the table we want for the charges
                # Read five lines to remove the header:
                # Summary of Natural Population Analysis:
                #
                # Natural Population
                # Natural    ---------------------------------------------
                # Atom No    Charge        Core      Valence    Rydberg      Total
                # ----------------------------------------------------------------
                for i in range(0, 5):
                    outFile.readline()
                # Then we read the actual table:
                for i in range(0, natoms):
                    # Each line follow the header with the form:
                    # C  1    0.92349      1.99948     3.03282    0.04422     5.07651
                    line = outFile.readline()
                    line = line.split(" ")
                    charges.append(line[2])
                # We have reached the end of the table, we can break the while loop
                break
        return charges


class NBOJob(SingleJob):
    def get_input(self):
        return "Input File"

    def get_runscript(self):
        return "Full run script"


def get_input_arguments():
    """
        Check command line options and accordingly set computation parameters
    """
    parser = argparse.ArgumentParser(description=help_description(),
                                     epilog=help_epilog())
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument('inputfile', type=str, nargs=1,
                        help='TAPE21 Containg IRC path')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as error:
        print(str(error))  # Print something like "option -a not recognized"
        sys.exit(2)

    # Get values from parser
    args = dict.fromkeys(['inputfile'])
    args['inputfile'] = os.path.basename(args.inputfile)

    return args


def help_description():
    """
        Returns description of program for help message
    """
    return "Help Description // To fill"


def help_epilog():
    """
        Returns additionnal help message
    """
    return "Help epilog // To Fill"


if __name__ == '__main__':
    main()
