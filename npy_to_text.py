""" Simplest script to read npy and write to comma-separated or space-separated
"""
import sys
import numpy as np


def npy2csv(filename_npy, filename_csv):
    """ Function to convert from npy to csv
    """
    # Read the file
    x = np.load(filename_npy)
    print("Columns of the input array: {0}".format(x.dtype.names))
    # Write out the file
    # %f : float
    # %i : integer
    # for formatting options: https://docs.scipy.org/doc/numpy-1.13.0/reference
    # /generated/numpy.savetxt.html
    np.savetxt(filename_csv, x, fmt="%f,%i,%f",
               header="box_size,count,lacunarity")
    print("\nFile was saved: {0}".format(filename_csv))
    return True


def npy2txt(filename_npy, filename_txt):
    """ Function to convert from npy to space separated values
    """
    x = np.load(filename_npy)
    print("Columns of the input array: {0}".format(x.dtype.names))
    np.savetxt(filename_txt, x, fmt="%15f %15i %15f",
               header="box_size,count,lacunarity")
    print("\nFile was saved: {0}".format(filename_txt))
    return True


if __name__ == "__main__":
    print("\n{0}\n USAGE: ~> python script.py input.npy out.csv".format("="*50))
    print("\t or ~> python script.py input.npy out.txt")
    print("\t depending of the desired output")
    print("{0}\n".format("="*50))
    # The sys.argv readas the arguments from the command line and transforms
    # them to a list
    args = sys.argv

    # Checking if there are only 2 arguments
    if (len(args) > 3):
        print("Error: only 2 arguments are allowed\n")
        exit(1)
    elif (len(args) == 3):
        if (".csv" in args[2]):
            npy2csv(args[1], args[2])
        elif (".txt" in args[2]):
            npy2txt(args[1], args[2])
        else:
            print("\nUse one of the 2 allowed output extensions: csv or txt")
            exit(1)
    else:
        print("Error: input arguments are erroneous: {0}".format(args))
