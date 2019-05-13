import argparse
import os
import sys


def main():
    """
        Main interface that reads the file and retrieve all useful geoms etc.
        Then rewrite geometries in a usable formatter_class
        Prep all NBO computations
        Run ADF on each computation
    """
    # Retrieve command line values
    args = get_input_arguments()
    input_file = args["input_file"]
    output_file = args["output_file"]
    nsteps = args["nsteps"]
    natoms = args["natoms"]

    # Get list of output file names
    NBO_files = [get_output_filename(input_file, i) for i in range(1, nsteps + 1)]
    # NBO_values is a list of list of charges
    print(NBO_files)
    NBO_values = [extract_NBO_charges(out_file, natoms) for out_file in NBO_files]
    # Write NBO data
    print_NBO_charges_to_file(NBO_values, output_file)


def get_output_filename(input_file, i):
    """
        Return output file name associated with input_file and step i
    """
    out_dir = os.path.basename(input_file + "." + str(i).zfill(3))
    for file in os.listdir(out_dir):
        if file.endswith(".out"):
            return out_dir + "/" + file


def extract_NBO_charges(output, natoms):
    """
        extract NBO Charges parsing the outFile
    """
    # Initialize charges list
    charges = []

    with open(output, mode="r") as outFile:
        line = "Foobar line"
        while line:
            line = outFile.readline()
            if "Summary of Natural Population Analysis:" in line:
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
    with open(file, mode="w+") as output_file:
        for i in range(0, len(charges_list)):
            output_file.write("     ".join(charges_list[i]) + "\n")


def get_input_arguments():
    """
        Check command line options and accordingly set computation parameters
    """
    parser = argparse.ArgumentParser(
        description=help_description(), epilog=help_epilog()
    )
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument(
        "-i",
        "--inputfile",
        type=str,
        nargs=1,
        help="basename for computations that have already been run",
    )
    parser.add_argument(
        "-o",
        "--outputfile",
        type=str,
        nargs=1,
        help="Output file in which to print the NBO charges",
    )
    parser.add_argument(
        "-n",
        "--numberofsteps",
        type=int,
        nargs=1,
        help="Number of computations that have been run",
    )
    parser.add_argument(
        "-a", "--numberofatoms", type=int, nargs=1, help="Number of atoms involved"
    )
    try:
        args = parser.parse_args()
    except argparse.ArgumentError as error:
        print(str(error))  # Print something like "option -a not recognized"
        sys.exit(2)

    # Get values from parser
    values = dict.fromkeys(["input_file", "output_file", "nsteps", "natoms"])
    values["input_file"] = os.path.basename(args.inputfile[0])
    values["output_file"] = os.path.basename(args.outputfile[0])
    values["nsteps"] = args.numberofsteps[0]
    values["natoms"] = args.numberofatoms[0]

    print(values)
    return values


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


if __name__ == "__main__":
    main()
