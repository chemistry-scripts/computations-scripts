"""
Perform an NBO analysis through Gaussian on every step of an IRC path.

Takes an IRC computed from ADF (at the moment), extract all intermediate geometries,
then create input files, run them, and analyse the output.
"""

import argparse
import os
import sys
import logging
from scm.plams import KFReader, Atom, Molecule


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
    args = get_input_arguments()
    input_file = args['input_file']
    output_file = args['output_file']
    basedir = os.path.abspath(os.curdir)

    geometries = IRC_coordinates_from_t21(input_file)
    natoms = number_of_atoms(input_file)

    # Get settings as a
    settings_head, settings_tail = gaussian_input_parameters(args)

    Gaussian_jobs = []
    # Prep a bunch of NBO computations
    for i, geom in enumerate(geometries):
        Gaussian_jobs.append(prepare_NBO_computation(basedir=basedir,
                                                     name="ASM_NBO",
                                                     geometry=geom,
                                                     id=i,
                                                     header=settings_head,
                                                     footer=settings_tail,
                                                     number_of_atoms=natoms))
    # Prepare all jobs (setup directories etc.)
    for job in Gaussian_jobs:
        job.setup_computation()
    # Run each job in parallel using threads
    # http://sametmax.com/en-python-les-threads-et-lasyncio-sutilisent-ensemble/
    # https://stackoverflow.com/users/2745756/emmanuel?tab=favorites
    for job in Gaussian_jobs:
        job.run()
    # NBO_values is a list of list of charges
    NBO_values = [job.extract_NBO_charges() for job in Gaussian_jobs]
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


def geometry_to_molecule(geometry):
    """Convert a list of XYZ coordinates to a Molecule object."""
    mol = Molecule()

    for i in range(0, len(geometry)):
        mol.add_atom(Atom(symbol=geometry[i][0], coords=(geometry[i][1], geometry[i][2], geometry[i][3])))

    return mol


def IRC_coordinates_from_t21(input_file):
    """Extract all data from a TAPE21 file."""
    # Read TAPE21 and get all useful data
    t21 = KFReader(input_file)

    # Number of atoms: 7
    natoms = t21.read('Geometry', 'nr of atoms')
    # atom types as indexes: [1, 2, 2, 3, 3, 4, 3]
    aatoms = t21.read('Geometry', 'fragment and atomtype index')[natoms:]
    # Atom symbols as list: ['C', 'O', 'H', 'B']
    xatoms = list_elements(input_file)
    # Actual list of atoms as used in geometry: ['C', 'O', 'O', 'H', 'H', 'B', 'H']
    satoms = [xatoms[aatoms[order - 1] - 1] for order in t21.read('Geometry', 'atom order index')[:natoms]]

    nstep_fw = t21.read('IRC_Forward', 'CurrentPoint')
    nstep_bw = t21.read('IRC_Backward', 'CurrentPoint')
    geometries_fw = t21.read('IRC_Forward', 'xyz')[0:nstep_fw * natoms * 3]
    geometries_bw = t21.read('IRC_Backward', 'xyz')[0:nstep_bw * natoms * 3]
    geometries_init = t21.read('IRC', 'xyz')

    # Reformat geometries into a long list of geometries for each step
    geometries_bw = coordinates_from_list(geometries_bw, natoms)
    geometries_fw = coordinates_from_list(geometries_fw, natoms)
    geometries_init = coordinates_from_list(geometries_init, natoms)
    geometries_fw.reverse()

    geometries = geometries_fw + geometries_init + geometries_bw

    return [[[s, mol[0], mol[1], mol[2]] for i, (s, mol) in enumerate(zip(satoms, molecule))] for molecule in geometries]


def coordinates_from_list(coordinates_list, natoms):
    """
    Rewrites a list into a list of lists of lists.

    Example:
    [x1,y1,z1,x2,y2,z2,x'1,y'1,z'1,x'2,y'2,z'2] becomes:
    [ [ [x1,y1,z1], [x2,y2,z2]], [ [x'1,y'1,z'1], [x'2,y'2,z'2] ] ]
    which can be acessed as:
    list[moleculeNumber][atomNumber][coordinateNumber]
    """
    # list = [x1,y1,z1,x2,y2,z2,x'1,y'1,z'1,x'2,y'2,z'2]
    coordinates_list = list(list_chunks(coordinates_list, 3))
    # [[x1,y1,z1],[x2,y2,z2],[x'1,y'1,z'1],[x'2,y'2,z'2]]
    coordinates_list = list(list_chunks(coordinates_list, natoms))
    # [ [[x1,y1,z1],[x2,y2,z2]], [[x'1,y'1,z'1],[x'2,y'2,z'2]]]
    return coordinates_list


def list_chunks(list, n):
    """
    Split a list in chunks of size n.

    Careful: does not care about the end of the list if len(list) is not a multiple of n
    It is actually a generator, so mind the usage.
    """
    for i in range(0, len(list), n):
        yield list[i:i + n]


def prepare_NBO_computation(basedir, name, geometry, id, header, footer, number_of_atoms):
    """
    From geometry, header, footer, create the input file.

    Return the input file as a list of lines
    """
    input = []

    # Put header
    input.extend(header)

    # Add geometry + blank line
    input.extend([' '.join([atom[0].ljust(5)] + [str(s).rjust(25) for s in atom[1:4]]) for atom in geometry])
    input.append('')

    # Add footer
    input.extend(footer)

    # Add two blank lines for the sake of Gaussian's weird behavior
    input.append("")
    input.append("")

    return Gaussian_Job(basedir, name, input, id, number_of_atoms)


def print_NBO_charges_to_file(charges_list, file):
    """Export NBO charges to a file that one can import in a spreadsheet or gnuplot."""
    with open(file, mode='w+') as output_file:
        for i in range(0, len(charges_list)):
            output_file.write('     '.join(charges_list[i]) + '\n')


def number_of_atoms(input_file):
    """Extract natoms from TAPE21."""
    t21 = KFReader(input_file)
    return t21.read('Geometry', 'nr of atoms')


def list_elements(input_file):
    """
    Return a list of all elements used in the computation.

    The list will look like:
    ['C', 'H', 'N', 'P']
    """
    t21 = KFReader(input_file)
    return str(t21.read('Geometry', 'atomtype')).split()


def get_input_arguments():
    """Check command line options and accordingly set computation parameters."""
    parser = argparse.ArgumentParser(description=help_description(),
                                     epilog=help_epilog())
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument('input_file', type=str, nargs=1,
                        help='TAPE21 Containg IRC path')
    parser.add_argument('output_file', type=str, nargs=1,
                        help='Output file in which to print the NBO charges')
    parser.add_argument('-f', '--functional', type=str, nargs='?', default='B3LYP-D3',
                        help='Functional used for the computation, as B3LYP-D3 or M062X\n'
                             'Hyphen will split into functional/dispersion parts when applicable')
    parser.add_argument('-b', '--basisset', type=str, nargs='?', default='6-31G*',
                        help='The basis set to use for all atoms')
    parser.add_argument('-m', '--memory', type=str, nargs='?', default='3GB',
                        help='Memory required per Gaussian calculation, i.e. per core')
    try:
        args = parser.parse_args()
    except argparse.ArgumentError as error:
        print(str(error))  # Print something like "option -a not recognized"
        sys.exit(2)

    # Get values from parser
    values = dict.fromkeys(['input_file', 'output_file', 'functional', 'dispersion',
                            'basisset', 'memory'])
    values['input_file'] = os.path.basename(args.input_file[0])
    values['output_file'] = os.path.basename(args.output_file[0])
    functional = args.functional.split('-')
    values['functional'] = functional[0]
    if len(functional) > 1:
        if functional[1] == 'GD3' or functional[1] == 'D3':
            values['dispersion'] = 'GD3'
    else:
        values['dispersion'] = None
    values['basisset'] = args.basisset
    values['memory'] = args.memory
    return values


def gaussian_input_parameters(args):
    """
    Return a tuple containing the top and bottom part used for the Gaussian calculation.

    Both parts are lists of strings.
    args is the dictionary coming from parsing the command line
    """
    header = []
    footer = []
    header.append('%NProcShared=1')
    header.append('%Mem=' + args['memory'])
    route = '# ' + args['functional'] + " "
    if args['dispersion'] is not None:
        route += "EmpiricalDispersion=" + args['dispersion'] + " "
    route += "gen pop=(nbo6read)"
    header.append(route)
    header.append('')
    # To update probably
    header.append('Title of computation')
    header.append('')
    # This is a singlet. Careful for other systems!
    header.append('0 1')

    # Basis set is the same for all elements. No ECP either.
    elements = list_elements(args['input_file'])
    elements = ' '.join(elements)
    basisset = args['basisset']
    footer.append(elements + ' 0')
    footer.append(basisset)
    footer.append("****")
    footer.append('')

    footer.append("$NBO")
    # NBO_FILES should be updated to something more useful
    footer.append("FILE=NBO_FILES")
    footer.append("PLOT")
    footer.append("$END")

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
        # Get into working directory
        os.chdir(self.path)

        # Initialize charges list
        charges = []

        with open(self.output_filename, mode='r') as outFile:
            line = 'Foobar line'
            while line:
                line = outFile.readline()
                if 'Summary of Natural Population Analysis:' in line:
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
                        # We have reached the end of the table, we can break the while loop
                break
        # Get back to the base directory
        os.chdir(self.basedir)
        return charges

    def setup_computation(self):
        """
        Set computation up before running it.

        Create working directory, write input file
        """
        # Create working directory
        os.makedirs(self.path, mode=0o777, exist_ok=False)
        logging.info("Created directory {dir}".format(dir=self.path))
        # Go into working directory
        os.chdir(self.path)
        # Write input file
        with open(self.input_filename, mode='w', buffering=-1, encoding=None, errors=None, newline=None, closefd=True) as input:
            input.write("\n".join(self.input))
        # Get back to base directory
        os.chdir(self.basedir)


if __name__ == '__main__':
    main()
