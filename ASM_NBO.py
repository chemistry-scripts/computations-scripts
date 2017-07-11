from scm.plams import KFReader, Settings, SingleJob, init


def main():
    """
        Main interface that reads the file and retrieve all useful geoms etc.
        Then rewrite geometries in a usable formatter_class
        Prep all NBO computations
        Run ADF on each computation
    """
    # Init PLAMS, Setting up a plams.$PID directory in the working dir
    init()

    # Setup blank Settings, to be updated upon reading the IRC
    settings = Settings()

    # Read TAPE21 and get all useful data
    irc_tape21 = KFReader("CO2-HBH2/CO2-HBH2-IRC.t21")

    # Retrieve Number of atoms
    natoms = irc_tape21.read("Geometry", "nr of atoms")

    nstep_fw = irc_tape21.read("IRC_Forward", "nset")
    nstep_bw = irc_tape21.read("IRC_Backward", "nset")

    geometries_fw = irc_tape21.read("IRC_Forward", "xyz")
    geometries_bw = irc_tape21.read("IRC_Backward", "xyz")
    geometries_init = irc_tape21.read("Geometry", "xyz")

    # Reformat geometries into a long list of geometries for each step
    geometries_bw = coordinatesFromList(geometries_bw, natoms, nstep_bw)
    geometries_fw = coordinatesFromList(geometries_fw, natoms, nstep_fw)
    geometries_init = coordinatesFromList(geometries_init, natoms, 1)
    geometries_fw.reverse()

    geometries = geometries_fw + geometries_init + geometries_bw

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
        Split a list in chunks of size n
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


if __name__ == '__main__':
    main()
