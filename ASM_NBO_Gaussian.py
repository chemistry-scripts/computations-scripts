"""
Perform an NBO analysis through Gaussian on every step of an IRC path.

Takes an IRC computed from Gaussian, extract all intermediate geometries,
then create input files, run them, and analyse the output.
"""

import argparse
import os
import sys
import logging
from cclib.io import ccread
from cclib.parser.utils import PeriodicTable
from concurrent.futures import ProcessPoolExecutor


def main():
    """
    Do all the job.

    Interface that reads the file and retrieve all useful geoms etc.
    Then rewrite geometries in a usable formatter_class
    Prep all NBO computations
    Run Gaussian on each computation
    """
    # Setup logging
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    # Retrieve command line values
    logger.debug("Retrieving input arguments")
    args = get_input_arguments()
    input_files = args['input_file']
    output_file = args['output_file']
    element_list = atom_types(args['input_file'][0])
    basedir = os.path.abspath(os.curdir)
    logger.debug("Input files: " + " ".join([str(path) for path in input_files]))
    logger.debug("Output file: " + str(output_file))
    logger.debug("Current directory: " + str(basedir))

    logger.debug("Geometry extraction")
    geometries = []
    for input_file in input_files:
        geometries.extend(IRC_coordinates_from_input(input_file))
    logger.debug("Extracted geometries: " + str(len(geometries)))

    logger.debug("Number of atoms extraction")
    natoms = number_of_atoms(input_files[0])
    logger.debug("Retrieved number of atoms")
    logger.debug("natoms: " + str(natoms))

    # Get settings as a tuple
    logger.debug("Getting Gaussian input parameters")
    settings_head, settings_tail = gaussian_input_parameters(args)

    gaussian_jobs = []
    # Prep a bunch of NBO computations
    for i, geom in enumerate(geometries):
        gaussian_jobs.append(prepare_NBO_computation(basedir=basedir,
                                                     name="ASM_NBO",
                                                     geometry=geom,
                                                     job_id=i,
                                                     header=settings_head,
                                                     footer=settings_tail,
                                                     number_of_atoms=natoms,
                                                     element_list=element_list))
    # Prepare all jobs (setup directories etc.)
    for job in gaussian_jobs:
        job.setup_computation()
    # Run each job in parallel using multiple processes
    with ProcessPoolExecutor() as executor:
        for job in gaussian_jobs:
            executor.submit(job.run)

    # NBO_values is a list of list of charges
    NBO_values = [job.extract_NBO_charges() for job in gaussian_jobs]
    job_ids = [job.id for job in gaussian_jobs]
    coordinates = [job.get_coordinates() for job in gaussian_jobs]

    # Compute distances, angles and dihedrals when necessary

    # Write NBO data
    print_NBO_charges_to_file(NBO_values, output_file)


def IRC_coordinates_to_xyz_file(filename, geometries):
    """Export coordinates in geometries table to a file straight in the working directory."""
    # Open file
    with open(filename, mode='w+') as XYZ:
        # Iterate over molecules
        for i in range(0, len(geometries)):
            # For each molecule, write "New molecule", the put the geometry as C 0.00 1.00 2.00
            XYZ.write('New molecule\n')
            for j in range(0, len(geometries[i])):
                for k in range(0, len(geometries[i][j])):
                    XYZ.write(str(geometries[i][j][k]) + '       ')
                XYZ.write('\n')
            XYZ.write('\n\n')
    return


def IRC_coordinates_from_input(input_file):
    """Return a table of coordinates with all last geometries (converged or not) from an IRC."""
    file = ccread(input_file, optdone_as_list=True)
    new_indexes = [x for x, y in enumerate(file.optstatus) if y & file.OPT_NEW > 0]
    # new_indexes finishes with 0, so has to finish with -1 for the last index.
    last_indexes = [x - 1 for x in new_indexes[1:len(new_indexes)] + [new_indexes[0]]]
    # file.atomcoords is an ndarray, so can be accessed with a list!
    coordinates = file.atomcoords[last_indexes]
    return coordinates.tolist()


def distance_from_coordinates(coord1, coord2):
    """Compute distance between two points."""
    return


def angle_from_coordinates(coord1, coord2, coord3):
    """Compute angle between three points."""
    return


def dihedral_from_coordinates(coord1, coord2, coord3, coord4):
    """Compute dihedral between four points."""
    return


def prepare_NBO_computation(basedir, name, geometry, job_id, header, footer, number_of_atoms, element_list):
    """
    From geometry, header, footer, create the input file.

    Return the input file as a list of lines
    """
    input_file = []

    # Put header
    input_file.extend(header)

    # Add geometry + blank line
    input_file.extend([' '.join([element_list[i].ljust(5)] +
                                ["{:.6f}".format(s).rjust(25) for s in atom])
                       for i, atom in enumerate(geometry)])
    input_file.append('')

    # Add footer
    input_file.extend(footer)

    # Add two blank lines for the sake of Gaussian's weird behavior
    input_file.append("")
    input_file.append("")

    return Gaussian_Job(basedir, name, input_file, job_id, number_of_atoms)


def print_NBO_charges_to_file(charges_list, file):
    """Export NBO charges to a file that one can import in a spreadsheet or gnuplot."""
    with open(file, mode='w+') as output_file:
        for i in range(0, len(charges_list)):
            output_file.write('     '.join(charges_list[i]) + '\n')


def number_of_atoms(input_file):
    """Extract natoms from file."""
    file = ccread(input_file)
    return file.natom


def atom_types(input_file):
    """
    Return a list of all atom types in the right order.

    The list will look like:
    ['C', 'H', 'H', 'H', 'C', 'N', 'P']
    """
    file = ccread(input_file)
    atoms = file.atomnos.tolist()
    periodic_table = PeriodicTable()
    atom_list = [periodic_table.element[i] for i in atoms]
    return atom_list


def list_elements(input_file):
    """
    Return a list of all unique elements used in the computation.

    The list will look like:
    ['C', 'H', 'N', 'P']
    """
    file = ccread(input_file)
    atoms = dict.fromkeys(file.atomnos.tolist())
    periodic_table = PeriodicTable()
    atom_list = [periodic_table.element[i] for i in atoms]
    return atom_list


def get_input_arguments():
    """Check command line options and accordingly set computation parameters."""
    parser = argparse.ArgumentParser(description=help_description(),
                                     epilog=help_epilog())
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument('-i', '--input_file', type=str, nargs='+',
                        help='file(s) Containg IRC path')
    parser.add_argument('-o', '--output_file', type=str, nargs=1,
                        help='Output file in which to print the NBO charges')
    parser.add_argument('-f', '--functional', type=str, nargs='?', default='B3LYP-D3',
                        help='Functional used for the computation, as B3LYP-D3 or M062X\n'
                             'Hyphen will split into functional/dispersion parts when applicable')
    parser.add_argument('-b', '--basisset', type=str, nargs='?', default='6-31G*',
                        help='The basis set to use for all atoms')
    parser.add_argument('-m', '--memory', type=str, nargs='?', default='3GB',
                        help='Memory required per Gaussian calculation, i.e. per core')
    parser.add_argument('-d', '--data', type=str, nargs='?',
                        help='Useful data to extract, such as bonds or angles\n'
                             'Write it as B 1 2 (bond between atoms 1 and 2), A 3 5 4 (Angle 3-5-4)\n'
                             'or D 3 8 9 1 (Dihedral 3-8-9-1)')
    try:
        args = parser.parse_args()
    except argparse.ArgumentError as error:
        print(str(error))  # Print something like "option -a not recognized"
        sys.exit(2)

    # Get values from parser
    values = dict.fromkeys(['input_file', 'output_file', 'functional', 'dispersion',
                            'basisset', 'memory', 'data'])
    values['input_file'] = [os.path.abspath(i) for i in args.input_file]
    values['output_file'] = os.path.abspath(args.output_file[0])
    functional = args.functional.split('-')
    values['functional'] = functional[0]
    if len(functional) > 1:
        if functional[1] == 'GD3' or functional[1] == 'D3':
            values['dispersion'] = 'GD3'
    else:
        values['dispersion'] = None
    values['basisset'] = args.basisset
    values['memory'] = args.memory
    if len(args.data) > 1:
        bonds = []
        angles = []
        dihedrals = []
        iterator = iter(args.data)
        for i in iterator:
            if i == 'B':
                x = next(iterator)
                y = next(iterator)
                bonds.append([x, y])
            if i == 'A':
                m = next(iterator)
                n = next(iterator)
                o = next(iterator)
                angles.append([m, n, o])
            if i == 'D':
                a = next(iterator)
                b = next(iterator)
                c = next(iterator)
                d = next(iterator)
                dihedrals.append([a, b, c, d])
        values['data'] = dict.fromkeys('bonds', 'angles', 'dihedrals')
        values['data']['bonds'] = bonds
        values['data']['angles'] = angles
        values['data']['dihedrals'] = dihedrals

    return values


def gaussian_input_parameters(args):
    """
    Return a tuple containing the top and bottom part used for the Gaussian calculation.

    Both parts are lists of strings.
    args is the dictionary coming from parsing the command line
    """
    logger = logging.getLogger()
    header = []
    footer = []
    header.append('%NProcShared=1')
    # header.append('%Mem=' + args['memory'])
    route = '# ' + args['functional'] + " "
    if args['dispersion'] is not None:
        route += "EmpiricalDispersion=" + args['dispersion'] + " "
    # route += "gen pop=(nbo6read)"
    route += "gen pop=(npa)"
    header.append(route)
    header.append('')
    # To update probably
    header.append('Title of computation')
    header.append('')
    # This is a singlet. Careful for other systems!
    header.append('0 1')

    # Basis set is the same for all elements. No ECP either.
    elements = list_elements(args['input_file'][0])
    elements = ' '.join(elements)
    basisset = args['basisset']
    footer.append(elements + ' 0')
    footer.append(basisset)
    footer.append("****")
    footer.append('')

    # footer.append("$NBO")
    # # NBO_FILES should be updated to something more useful
    # footer.append("FILE=NBO_FILES")
    # footer.append("PLOT")
    # footer.append("$END")

    logger.debug("Header: \n" + '\n'.join(header))
    logger.debug("Footer: \n" + '\n'.join(footer))
    return header, footer


def help_description():
    """Return description of program for help message."""
    return 'Help Description // To fill'


def help_epilog():
    """Return additionnal help message."""
    return 'Help epilog // To Fill'


class Gaussian_Job():
    """
    Class that can be used as a container for Gaussian jobs.

    Attributes:
        - input (input file,  list of strings)
        - name (name of computation, string)
        - id (unique identifier, int)
        - natoms (number of atoms, int)
        - basedir (base directory, os.path object)
        - path (path in which to run current computation, os.path object)
        - input_filename (file_name.com, str)
        - output_filename (file_name.log, str)

    """

    def __init__(self, basedir, name, input, id, natoms):
        """Build  the Gaussian_job class."""
        self.name = name
        self.input = input
        self.id = id
        self.natoms = natoms
        # base directory from which all computations are started
        self.basedir = basedir
        # Set path as: /base/directory/my_name.000xx/
        self.path = os.path.join(self.basedir, self.name.replace(" ", "_") + "." + str(self.id).zfill(4))
        self.input_filename = self.name.replace(" ", "_") + ".com"
        self.output_filename = self.name.replace(" ", "_") + ".log"

    def run(self):
        """Start the job."""
        # Log computation start
        logger = logging.getLogger()
        logger.info("Starting computation " + str(self.id))
        # Get into workdir, start gaussian, then back to basedir
        os.chdir(self.path)
        os.system('g09 < ' + self.input_filename + ' > ' + self.output_filename)
        os.chdir(self.basedir)
        # Log end of computation
        logger.info("Finished computation " + str(self.id))

    def extract_NBO_charges(self):
        """Extract NBO Charges parsing the outFile."""
        # Log start
        logger = logging.getLogger()
        logger.info("Parsing results from computation " + str(self.id))

        # Get into working directory
        os.chdir(self.path)

        # Initialize charges list
        charges = []

        with open(self.output_filename, mode='r') as outFile:
            line = 'Foobar line'
            while line:
                line = outFile.readline()
                if 'Summary of Natural Population Analysis:' in line:
                    logger.debug("ID " + str(self.id) + ": Found NPA table.")
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
                    for i in range(0, self.natoms):
                        # Each line follow the header with the form:
                        # C  1    0.92349      1.99948     3.03282    0.04422     5.07651
                        line = outFile.readline()
                        line = line.split()
                        charges.append(line[2])
                    logger.debug("ID " + str(self.id) + ": " +
                                 "Charges = " + " ".join([str(i) for i in charges]))
                    # We have reached the end of the table, we can break the while loop
                    break
                # End of if 'Summary of Natural Population Analysis:'
        # Get back to the base directory
        os.chdir(self.basedir)
        return charges

    def get_coordinates(self):
        """Extract coordinates from output file."""
        # Log start
        logger = logging.getLogger()
        logger.info("Extracting coordinates " + str(self.id))

        # Get into working directory
        os.chdir(self.path)

        # Parse file with cclib
        data = ccread(self.output_filename)

        #  Return the coordinates
        return data.atomcoords

    def setup_computation(self):
        """
        Set computation up before running it.

        Create working directory, write input file
        """
        # Create working directory
        os.makedirs(self.path, mode=0o777, exist_ok=False)
        logging.info("Created directory %s", self.path)
        # Go into working directory
        os.chdir(self.path)
        # Write input file
        with open(self.input_filename, mode='w') as input_file:
            input_file.write("\n".join(self.input))
        # Get back to base directory
        os.chdir(self.basedir)


if __name__ == '__main__':
    main()
