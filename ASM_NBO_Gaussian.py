"""
Perform an NBO analysis through Gaussian on every step of an IRC path.

Takes an IRC computed from Gaussian, extract all intermediate geometries,
then create input files, run them, and analyse the output.
"""
# pylint: disable=invalid-name

import argparse
import os
import sys
import logging
from concurrent.futures import ProcessPoolExecutor
from cclib.io import ccread
from cclib.parser.utils import PeriodicTable
import numpy as np


def main():
    """
    Do all the job.

    Interface that reads the file and retrieve all useful geoms etc.
    Then rewrite geometries in a usable formatter_class
    Prep all NBO computations
    Run Gaussian on each computation
    """
    # pylint: disable=too-many-locals
    # Since it is the main, it is too annoying to remove such a number of variables.

    # Setup logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)

    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s')
    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)

    # Retrieve command line values
    logger.debug("Retrieving input arguments")
    args = get_input_arguments()
    input_files = args['input_file']
    output_file = args['output_file']
    element_list = atom_types(args['input_file'][0])
    basedir = os.path.abspath(os.curdir)
    logger.debug("Input files: %s", " ".join([str(path) for path in input_files]))
    logger.debug("Output file: %s", str(output_file))
    logger.debug("Current directory: %s", str(basedir))

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
                                                     natoms=natoms,
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
    job_ids = [job.job_id for job in gaussian_jobs]

    # Compute distances, angles and dihedrals when necessary
    measured_data = []
    logger.debug("Data to extract: %s", args['data'])
    if args['data']:
        coordinates = [job.get_coordinates() for job in gaussian_jobs]
        measured_data = [compute_measurements(coord, args['data']) for coord in coordinates]
    # Write NBO data
    print_NBO_charges_to_file(charges_list=NBO_values,
                              out_file=output_file,
                              measures=measured_data,
                              job_ids=job_ids)


def compute_measurements(coordinates, required_data):
    """Compute required measurements from coordinates."""
    # Initialize list
    extracted_data = []

    logger = logging.getLogger()
    logger.debug("Coordinates: %s", coordinates)
    logger.debug("Required data: %s", required_data)
    # Extract all bonds
    extracted_data.extend([distance_from_coordinates(coordinates[bond[0]],
                                                     coordinates[bond[1]])
                           for bond in required_data['bonds']])

    # Extract all bonds
    extracted_data.extend([angle_from_coordinates(coordinates[angle[0]],
                                                  coordinates[angle[1]],
                                                  coordinates[angle[2]])
                           for angle in required_data['angles']])

    # Extract all bonds
    extracted_data.extend([dihedral_from_coordinates(coordinates[dihedral[0]],
                                                     coordinates[dihedral[1]],
                                                     coordinates[dihedral[2]],
                                                     coordinates[dihedral[3]])
                           for dihedral in required_data['dihedrals']])

    return extracted_data


def IRC_coordinates_to_xyz_file(filename, geometries):
    """Export coordinates in geometries table to a file straight in the working directory."""
    # Open file
    with open(filename, mode='w+') as xyz_file:
        # Iterate over molecules
        for i in range(0, len(geometries)):
            # For each molecule, write "New molecule", the put the geometry as C 0.00 1.00 2.00
            xyz_file.write('New molecule\n')
            for j in range(0, len(geometries[i])):
                for k in range(0, len(geometries[i][j])):
                    xyz_file.write(str(geometries[i][j][k]) + '       ')
                xyz_file.write('\n')
            xyz_file.write('\n\n')
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
    return np.linalg.norm(coord2 - coord1)


def angle_from_coordinates(coord1, coord2, coord3):
    """Compute angle between three points."""
    bond21 = coord1 - coord2
    bond23 = coord3 - coord2

    cosine_angle = np.dot(bond21, bond23) / (np.linalg.norm(bond21) * np.linalg.norm(bond23))
    angle = np.arccos(cosine_angle)

    return np.degrees(angle)


def dihedral_from_coordinates(coord1, coord2, coord3, coord4):
    """Compute dihedral between four points."""
    # Strategy:
    # p1 <--vector1-- p2 --vector0--> p3 --vector2--> p4
    # Dihedral corresponds to angle between planes (vector0,vector1) and (vector0,vector2)

    vector0 = coord3 - coord2
    vector1 = coord1 - coord2
    vector2 = coord4 - coord3

    # normalize vector 0 in order to project properly
    vector0 /= np.linalg.norm(vector0)

    # Decomposition of vectors 1 and 2 into projection on vector 0 and other component,
    # that is kept as proj1 and proj2
    proj1 = vector1 - np.dot(vector0, vector1) * vector0
    proj2 = vector2 - np.dot(vector0, vector2) * vector0

    # Angle between proj1 and proj2 is the dihedral
    # We use here a trick with the cross product in dot2, since it allows to keep proper
    # sign with sin and cos combination, and avoids any renormalization (that can be slow)
    dot1 = np.dot(proj1, proj2)
    dot2 = np.dot(np.cross(vector0, proj1), proj2)
    return np.degrees(np.arctan2(dot2, dot1))


def prepare_NBO_computation(basedir, name, geometry, job_id, header, footer, natoms, element_list):
    """
    From geometry, header, footer, create the input file.

    Return the input file as a list of lines
    """
    # pylint: disable=too-many-arguments
    # We need them all
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

    return GaussianJob(basedir, name, input_file, job_id, natoms)


def print_NBO_charges_to_file(charges_list, out_file, measures, job_ids):
    """Export NBO charges to a file that one can import in a spreadsheet or gnuplot."""
    logger = logging.getLogger()
    with open(out_file, mode='w+') as output_file:
        for (i, job_id) in enumerate(job_ids):
            logger.debug("Printing job %s", job_id)
            logger.debug("Job_id: %s", job_id)
            logger.debug("Measures: %s", measures[i])
            logger.debug("Type: %s", type(measures[i][0]))
            logger.debug("Charges: %s", charges_list[i])
            logger.debug("Type: %s", type(charges_list[i][0]))

            line = str(job_id).ljust(5)
            # Separate job_id and measures
            line += ' '
            if measures:
                # Each measure has a maximum of three digits before decimal separation, plus a sign.
                # It is thus maximum (with rounding to three digits after decimal) 8 in length.
                line += ' '.join("{0:.3f}".format(value).rjust(8) for value in measures[i])
            # Separate measures and charges
            line += ' '
            # Each charge is rounded to three digits after the decimal separator
            # It is the same length as measures.
            line += ' '.join(["{0:.3f}".format(float(charges)).rjust(8)
                              for charges in charges_list[i]])
            # End of line
            line += '\n'
            output_file.write(line)


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
    logger = logging.getLogger()
    parser = argparse.ArgumentParser(description=help_description(),
                                     epilog=help_epilog())
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument('-i', '--input_file', type=str, nargs='+',
                        help='file(s) Containing IRC path')
    parser.add_argument('-o', '--output_file', type=str, nargs=1,
                        help='Output file in which to print the NBO charges')
    parser.add_argument('-f', '--functional', type=str, nargs='?', default='B3LYP-D3',
                        help='Functional used for the computation, as B3LYP-D3 or M062X\n'
                             'Hyphen will split into functional/dispersion parts when applicable')
    parser.add_argument('-b', '--basisset', type=str, nargs='?', default='6-31G*',
                        help='The basis set to use for all atoms')
    parser.add_argument('-m', '--memory', type=str, nargs='?', default='3GB',
                        help='Memory required per Gaussian calculation, i.e. per core')
    parser.add_argument('-d', '--data', type=str, nargs='*',
                        help='Useful data to extract, such as bonds or angles\n'
                             'Write as B 1 2 (bond between atoms 1 and 2), A 3 5 4 (Angle 3-5-4)\n'
                             'or D 3 8 9 1 (Dihedral 3-8-9-1)')
    parser.add_argument('-r', '--fragment', type=str, nargs='*',
                        help='List of atoms in one of the fragments to consider.\n'
                             'If absent, fragmentation is not considered.\n'
                             'If present, all listed atoms (as numbers in geometry) are used in\n'
                             'Frag 1, while others are added to Frag 0.\n')
    try:
        args = parser.parse_args()
    except argparse.ArgumentError as error:
        print(str(error))  # Print something like "option -a not recognized"
        sys.exit(2)

    # Get values from parser
    values = dict.fromkeys(['input_file', 'output_file', 'functional', 'dispersion',
                            'basisset', 'memory', 'data', 'frag1', 'frag2'])
    values['input_file'] = [os.path.abspath(i) for i in args.input_file]
    logger.debug("Input files: %s", values['input_file'])
    values['output_file'] = os.path.abspath(args.output_file[0])
    logger.debug("Output file: %s", values['output_file'])
    functional = args.functional.split('-')
    values['functional'] = functional[0]
    if len(functional) > 1:
        if functional[1] == 'GD3' or functional[1] == 'D3':
            values['dispersion'] = 'GD3'
    else:
        values['dispersion'] = None
    logger.debug("Functional: %s", values['functional'])
    logger.debug("Dispersion: %s", values['dispersion'])
    values['basisset'] = args.basisset
    logger.debug("Basis set: %s", values['basisset'])
    values['memory'] = args.memory
    logger.debug("Memory: %s", values['memory'])
    if len(args.data) > 1:
        bonds = []
        angles = []
        dihedrals = []
        iterator = iter(args.data)
        for i in iterator:
            if i == 'B':
                bonds.append([int(next(iterator)) - 1,
                              int(next(iterator)) - 1])
            if i == 'A':
                angles.append([int(next(iterator)) - 1,
                               int(next(iterator)) - 1,
                               int(next(iterator)) - 1])
            if i == 'D':
                dihedrals.append([int(next(iterator)) - 1,
                                  int(next(iterator)) - 1,
                                  int(next(iterator)) - 1,
                                  int(next(iterator)) - 1])
        values['data'] = dict.fromkeys(['bonds', 'angles', 'dihedrals'])
        values['data']['bonds'] = bonds
        values['data']['angles'] = angles
        values['data']['dihedrals'] = dihedrals
    logger.debug("Data to extract: %s", values['data'])
    if len(args.frag) > 1:
        values['frag1'] = args.frag
        natoms = number_of_atoms(values['input_file'])
        atoms = range(1, natoms + 1)
        values['frag0'] = [atom for atom in atoms if atom not in values['frag1']]
    logger.debug("Fragment 0: %s", values['frag0'])
    logger.debug("Fragment 1: %s", values['frag1'])
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

    logger.debug("Header: \n %s", '\n'.join(header))
    logger.debug("Footer: \n %s", '\n'.join(footer))
    return header, footer


def help_description():
    """Return description of program for help message."""
    return 'Help Description // To fill'


def help_epilog():
    """Return additionnal help message."""
    return 'Help epilog // To Fill'


class GaussianJob:
    """
    Class that can be used as a container for Gaussian jobs.

    Attributes:
        - input (input file, list of strings)
        - name (name of computation, string)
        - id (unique identifier, int)
        - natoms (number of atoms, int)
        - basedir (base directory, os.path object)
        - path (path in which to run current computation, os.path object)
        - input_filename (file_name.com, str)
        - output_filename (file_name.log, str)

    """

    # pylint: disable=too-many-instance-attributes

    def __init__(self, basedir, name, input_script, job_id, natoms):
        """Build  the GaussianJob class."""
        # pylint: disable=too-many-arguments
        # We need them all
        self.name = name
        self.input_script = input_script
        self.job_id = job_id
        self.natoms = natoms
        # base directory from which all computations are started
        self.basedir = basedir
        # Set path as: /base/directory/my_name.000xx/
        self.path = os.path.join(self.basedir,
                                 self.name.replace(" ", "_") + "." + str(self.job_id).zfill(4))
        self.input_filename = self.name.replace(" ", "_") + ".com"
        self.output_filename = self.name.replace(" ", "_") + ".log"

    def run(self):
        """Start the job."""
        # Log computation start
        logger = logging.getLogger()
        logger.info("Starting computation %s", str(self.job_id))
        # Get into workdir, start gaussian, then back to basedir
        os.chdir(self.path)
        os.system('g09 < ' + self.input_filename + ' > ' + self.output_filename)
        os.chdir(self.basedir)
        # Log end of computation
        logger.info("Finished computation %s", str(self.job_id))

    def extract_NBO_charges(self):
        """Extract NBO Charges parsing the output file."""
        # Log start
        logger = logging.getLogger()
        logger.info("Parsing results from computation %s", str(self.job_id))

        # Get into working directory
        os.chdir(self.path)

        # Initialize charges list
        charges = []

        with open(self.output_filename, mode='r') as out_file:
            line = 'Foobar line'
            while line:
                line = out_file.readline()
                if 'Summary of Natural Population Analysis:' in line:
                    logger.debug("ID %s: Found NPA table.", str(self.job_id))
                    # We have the table we want for the charges
                    # Read five lines to remove the header:
                    # Summary of Natural Population Analysis:
                    #
                    # Natural Population
                    # Natural    ---------------------------------------------
                    # Atom No    Charge        Core      Valence    Rydberg      Total
                    # ----------------------------------------------------------------
                    for _ in range(0, 5):
                        out_file.readline()
                    # Then we read the actual table:
                    for _ in range(0, self.natoms):
                        # Each line follow the header with the form:
                        # C  1    0.92349      1.99948     3.03282    0.04422     5.07651
                        line = out_file.readline()
                        line = line.split()
                        charges.append(line[2])
                    logger.debug("ID %s: Charges = %s", str(self.job_id),
                                 " ".join([str(i) for i in charges]))
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
        logger.info("Extracting coordinates for job %s", str(self.job_id))

        # Get into working directory
        os.chdir(self.path)

        # Parse file with cclib
        data = ccread(self.output_filename)

        #  Return the first coordinates, since it is a single point
        return data.atomcoords[0]

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
            input_file.write("\n".join(self.input_script))
        # Get back to base directory
        os.chdir(self.basedir)


if __name__ == '__main__':
    main()
