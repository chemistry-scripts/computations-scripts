import argparse
import os
import sys
from scm.plams import KFReader, Atom, Molecule


def main():
    """
        Main interface that reads the file and retrieve all useful geoms etc.
        Then rewrite geometries in a usable formatter_class
        Prep all NBO computations
        Run Gaussian on each computation
    """
    # Retrieve command line values
    args = get_input_arguments()
    input_file = args['input_file']
    output_file = args['output_file']

    geometries = IRC_coordinates_from_t21(input_file)
    natoms = number_of_atoms(input_file)

    # Get settings as a
    settings_head, settings_tail = gaussian_input_parameters(args)

    # Prep a bunch of NBO computations
    NBO_jobs = [prepare_NBO_computation(geom, settings_head, settings_tail) for geom in geometries]
    # Run each job in parallel as managed by plams
    # http://sametmax.com/en-python-les-threads-et-lasyncio-sutilisent-ensemble/
    # https://stackoverflow.com/users/2745756/emmanuel?tab=favorites
    for job in NBO_jobs:
        job.run()
    # NBO_values is a list of list of charges
    NBO_values = [extract_NBO_charges(job._filenames.get('out'), natoms) for job in NBO_jobs]
    # Write NBO data
    print_NBO_charges_to_file(NBO_values, output_file)


def IRC_coordinates_to_xyz_file(filename, geometries):
    """
        Export coordinates in geometries table to a file 'filename'
        straight in the working directory
    """
    with open(filename, mode='w+') as XYZ:
        for i in range(0, len(geometries)):
            XYZ.write('New molecule\n')
            for j in range(0, len(geometries[i])):
                for k in range(0, len(geometries[i][j])):
                    XYZ.write(str(geometries[i][j][k]) + '       ')
                XYZ.write('\n')
            XYZ.write('\n\n')
    return


def geometry_to_molecule(geometry):
    """
    Convert a list of XYZ coordinates to a Molecule object
    """
    mol = Molecule()

    for i in range(0, len(geometry)):
        mol.add_atom(Atom(symbol=geometry[i][0], coords=(geometry[i][1], geometry[i][2], geometry[i][3])))

    return mol


def IRC_coordinates_from_t21(input_file):
    """
        Extract all data from a TAPE21 file.
    """
    # Read TAPE21 and get all useful data
    t21 = KFReader(input_file)

    # Number of atoms: 7
    natoms = t21.read('Geometry', 'nr of atoms')
    # atom types as indexes: [1, 2, 2, 3, 3, 4, 3]
    aatoms = t21.read('Geometry', 'fragment and atomtype index')[natoms:]
    # Atom symbols as list: ['C', 'O', 'H', 'B']
    xatoms = str(t21.read('Geometry', 'atomtype')).split()
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
        Rewrites a list into a list of lists of lists:
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
        Generator that splits a list in chunks of size n
        /!\ Does not care about the end of the list if len(list) is not a multiple of n
    """
    for i in range(0, len(list), n):
        yield list[i:i + n]


def prepare_NBO_computation(geometry, runparameters):
    """
        Taking a geometry and a set of ADF parameters (Basis set, Functional, ZORA, etc)
        Add the NBO Necessary keywords and returns the plams job
    """
    nbo_script = '\n'.join(['"$ADFBIN/adfnbo" << eor',
                            'write',
                            'spherical',
                            'nbo-analysis',
                            'eor',
                            '',
                            '"$ADFBIN/gennbo6" FILE47',
                            '',
                            '"$ADFBIN/adfnbo" <<eor',
                            'spherical',
                            'fock',
                            'read',
                            'end input',
                            'eor'])
    runparameters.runscript.post = nbo_script

    return job


def extract_NBO_charges(output, natoms):
    """
        extract NBO Charges parsing the outFile
    """
    # Initialize charges list
    charges = []

    with open(output, mode='r') as outFile:
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
                for i in range(0, natoms):
                    # Each line follow the header with the form:
                    # C  1    0.92349      1.99948     3.03282    0.04422     5.07651
                    line = outFile.readline()
                    line = line.split()
                    charges.append(line[2])
                # We have reached the end of the table, we can break the while loop
                break
        return charges


def print_NBO_charges_to_file(charges_list, file):
    """
        Export NBO charges to a file that one can import in a spreadsheet or gnuplot
    """
    with open(file, mode='w+') as output_file:
        for i in range(0, len(charges_list)):
            output_file.write('     '.join(charges_list[i]) + '\n')


def number_of_atoms(input_file):
    """
        Extracts natoms from TAPE21
    """
    t21 = KFReader(input_file)
    return t21.read('Geometry', 'nr of atoms')


def get_input_arguments():
    """
        Check command line options and accordingly set computation parameters
    """
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
    parser.add_argument('-r', '--relativistic', type=str, nargs='?', default='None',
                        help='Relativistic effects: Scalar, Spin-Orbit or None (default)')
    parser.add_argument('-b', '--basisset', type=str, nargs='?', default='DZP',
                        help='The basis set to use for all atoms')
    parser.add_argument('-c', '--frozencore', type=str, nargs='?', default='None',
                        help='Frozen core to use: None (Default), Small, Large')
    parser.add_argument('-i', '--integrationquality', type=str, nargs='?', default='Good',
                        help='Numerical Integration Quality. Default: Good')
    try:
        args = parser.parse_args()
    except argparse.ArgumentError as error:
        print(str(error))  # Print something like "option -a not recognized"
        sys.exit(2)

    # Get values from parser
    values = dict.fromkeys(['input_file', 'output_file', 'functional', 'dispersion', 'relativistic',
                            'basisset', 'frozencore', 'integrationquality'])
    values['input_file'] = os.path.basename(args.input_file[0])
    values['output_file'] = os.path.basename(args.output_file[0])
    functional = args.functional.split('-')
    values['functional'] = functional[0]
    if len(functional) > 1:
        if functional[1] == 'GD3' or functional[1] == 'D3':
            values['dispersion'] = 'GD3'
    else:
        values['dispersion'] = None
    values['relativistic'] = args.relativistic
    values['basisset'] = args.basisset
    values['frozencore'] = args.frozencore
    values['integrationquality'] = args.integrationquality
    return values


def help_description():
    """
        Returns description of program for help message
    """
    return 'Help Description // To fill'


def help_epilog():
    """
        Returns additionnal help message
    """
    return 'Help epilog // To Fill'


if __name__ == '__main__':
    main()
