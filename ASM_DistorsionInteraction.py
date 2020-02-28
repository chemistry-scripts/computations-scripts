import os
import sys
import argparse
import logging
import scm.plams as plams
import numpy as np


def main():
    """
        The main function.
        Reads input, sends it to a parser, then sets up the computation
    """
    # Setup logging
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)

    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.DEBUG)

    formatter = logging.Formatter("%(asctime)s :: %(levelname)s :: %(message)s")
    stream_handler.setFormatter(formatter)

    logger.addHandler(stream_handler)

    # Parse command line arguments
    args = parse_args()

    geometries = get_IRC_geometries(args)
    print(np.shape(geometries))


def parse_args():
    """
        The argument parser
    """
    logger = logging.getLogger()
    parser = argparse.ArgumentParser(
        description=help_description(), epilog=help_epilog()
    )
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument(
        "-i", "--input_file", type=str, nargs=1, help="The input file to submit"
    )
    parser.add_argument(
        "-o",
        "--output_file",
        type=str,
        nargs=1,
        help="Output file in which to print the NBO charges",
    )
    parser.add_argument(
        "-f",
        "--functional",
        type=str,
        nargs="?",
        default="B3LYP-D3",
        help="Functional used for the computation, as B3LYP-D3 or M062X\n"
        "Hyphen will split into functional/dispersion parts when applicable",
    )
    parser.add_argument(
        "-b",
        "--basis_set",
        type=str,
        nargs="?",
        default="6-31G*",
        help="The basis set to use for all atoms",
    )
    parser.add_argument(
        "-d",
        "--data",
        type=str,
        nargs="*",
        help="Useful data to extract, such as bonds or angles\n"
        "Write as B 1 2 (bond between atoms 1 and 2), A 3 5 4 (Angle 3-5-4)\n"
        "or D 3 8 9 1 (Dihedral 3-8-9-1)",
    )
    parser.add_argument(
        "-r",
        "--fragment",
        type=int,
        nargs="*",
        help="List of atoms in one of the fragments to consider.\n"
        "If absent, fragmentation is not considered.\n"
        "If present, all listed atoms (as numbers in geometry) are used in\n"
        "Frag 1, while others are added to Frag 0.\n",
    )
    parser.add_argument(
        "-n",
        "--nbo",
        type=bool,
        default=False,
        nargs="?",
        help="Request NBO calculations at every point",
    )

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as error:
        print(str(error))  # Print something like "option -a not recognized"
        sys.exit(2)

    run_values = dict.fromkeys(
        [
            "input_file",
            "output_file",
            "functional",
            "basis_set",
            "data",
            "fragment",
            "NBO",
        ]
    )
    # Get values from parser
    run_values["input_file"] = [args.input_file[0]]  # os.path.basename ?
    run_values["output_file"] = args.output_file[0]
    run_values["functional"] = args.functional
    run_values["basis_set"] = args.basis_set
    # Parse data to extract
    if args.data:
        bonds = []
        angles = []
        dihedrals = []
        iterator = iter(args.data)
        for i in iterator:
            if i == "B":
                bonds.append([int(next(iterator)) - 1, int(next(iterator)) - 1])
            if i == "A":
                angles.append(
                    [
                        int(next(iterator)) - 1,
                        int(next(iterator)) - 1,
                        int(next(iterator)) - 1,
                    ]
                )
            if i == "D":
                dihedrals.append(
                    [
                        int(next(iterator)) - 1,
                        int(next(iterator)) - 1,
                        int(next(iterator)) - 1,
                        int(next(iterator)) - 1,
                    ]
                )
        run_values["data"] = dict.fromkeys(["bonds", "angles", "dihedrals"])
        run_values["data"]["bonds"] = bonds
        run_values["data"]["angles"] = angles
        run_values["data"]["dihedrals"] = dihedrals
    logger.debug("Data to extract: %s", run_values["data"])

    # Parse fragmentation
    if args.fragment:
        run_values["frag1"] = args.fragment
        natoms = number_of_atoms(run_values["input_file"][0])
        atoms = range(1, natoms + 1)
        run_values["frag0"] = [
            atom for atom in atoms if atom not in run_values["frag1"]
        ]
    else:
        run_values["frag1"] = []
        run_values["frag0"] = [
            atom for atom in range(1, number_of_atoms(run_values["input_file"][0] + 1))
        ]

    logger.debug("Fragment 0: %s", run_values["frag0"])
    logger.debug("Fragment 1: %s", run_values["frag1"])

    return run_values


def number_of_atoms(file):
    """Return number of atoms from TAPE21 file"""
    reader = plams.KFReader(file)
    return reader.read("Geometry", "nnuc")


def get_IRC_geometries(args):
    file = plams.KFFile(args["input_file"][0])
    natoms = number_of_atoms(args["input_file"][0])
    if "IRC_Forward" in file:
        geom_crude = file.read("IRC_Forward", "xyz")
        logging.info("Forward IRC size: " + str(np.shape(geom_crude)))
        geom_forward = np.reshape(geom_crude, (-1, 3, natoms))
        logging.info("Forward IRC size: " + str(np.shape(geom_forward)))

    if "IRC_Backward" in file:
        geom_crude = file.read("IRC_Backward", "xyz")
        logging.info("Backward IRC size: " + str(np.shape(geom_crude)))
        geom_backward = np.reshape(geom_crude, (-1, 3, natoms))
        logging.info("Backward IRC shape: " + str(np.shape(geom_backward)))

    if geom_backward in dir() and geom_forward in dir():
        return np.concatenate(geom_forward, geom_backward)
    elif geom_backward in dir():
        return geom_backward
    elif geom_forward in dir():
        return geom_forward
    else:
        raise ValueError("No IRC data in file")


def help_description():
    """Return description of program for help message."""
    return "Help Description // To fill"


def help_epilog():
    """Return additionnal help message."""
    return "Help epilog // To Fill"


if __name__ == "__main__":
    main()
