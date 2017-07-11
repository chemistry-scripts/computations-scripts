import argparse
import scm.plams as plams

reader = plams.KFReader("CO2-HBH2-NBO.t21")

reader.read('NBOs', 'Label_1')

reader.read('General', 'termination status')


def main():
    """
        The main function.
        Reads input, sends it to a parser, then sets up the computation
    """
    parse_args(args)


def parse_args(args):
    """
        The argument parser
    """
    parser = argparse.ArgumentParser(description=help_description(),
                                     epilog=help_epilog())
    parser.formatter_class = argparse.RawDescriptionHelpFormatter
    parser.add_argument('-n', '--nproc', default=24, type=int,
                        help="Number of cores used for the computation")
    parser.add_argument('-t', '--walltime', default="24:00:00", type=str,
                        help="Maximum time allowed for the computation")
    parser.add_argument('-m', '--memory', type=int,
                        help="Amount of memory allowed for the computation")
    parser.add_argument('inputfile', type=str, nargs='+',
                        help='The input file to submit')

    try:
        args = parser.parse_args()
    except argparse.ArgumentError as error:
        print(str(error))  # Print something like "option -a not recognized"
        sys.exit(2)

    runvalues = dict.fromkeys(['inputfile', 'cores', 'walltime', 'memory',
                               'chk', 'oldchk', 'rwf', 'nproc_in_input',
                               'memory_in_input'])
    # Get values from parser
    runvalues['inputfile'] = os.path.basename(args.inputfile[0])
    runvalues['cores'] = args.nproc
    runvalues['walltime'] = args.walltime
    if args.memory:
        runvalues['memory'] = args.memory

    # Initialize empty values
    runvalues['chk'] = set()
    runvalues['oldchk'] = set()
    runvalues['rwf'] = set()
    runvalues['nbo'] = False
    runvalues['nbo_basefilename'] = ''

    return runvalues
