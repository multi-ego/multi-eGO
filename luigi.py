from getopt import getopt
import os
import sys, getopt



def main(argv):
    try:
        flags, args = getopt.getopt(argv, '', ['protein=', 'path_gromacs=', 'path_pdb='])
    except getopt.GetoptError:
        sys.exit(2)

    print(flags)
        

    


    try:
        os.mkdir()

    except OSError as error:
        print(f'Folder already created {error}')


if __name__ == '__main__':
    main(sys.argv[1:])

