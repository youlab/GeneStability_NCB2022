#!/c/Users/Andrea/Anaconda3/python

import sys
import csv
from libs.utils import read_genome, read_barcodes

"""
By: Andrea Weiss
Last edited: 15 April 2020

INPUTS
filename: sequencing reads to search through to match barcodes to 
barcode_file: file containing all known target barcodes (i.e. pattern to match sequencing reads) 
mismatch: number of mismatches to allow during matching process
"""

# parrse input arguments
def parse_args(argv):
    filename = argv[1]
    barcode_file = argv[2]
    mismatch = int(argv[3]) # number of max mismatches allowed
    return filename, barcode_file, mismatch


# naive matching algorithm - finds each occurences of a given barcode in a file of sequences allowing a predifined 
# number of mismatched basepairs and returns the number of occurences
def naive(p, t, mismatch):
    occurrences = [] 
    for i in range(len(t) - len(p) + 1):  # loop over alignments
        match = True
        k = 0
        for j in range(len(p)):  # loop over characters
            if t[i+j] != p[j]:  # compare characters
                k += 1
                if k > mismatch:
                    match = False
                    break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences

# does not store occurences... a bit faster
def naive_v2(barcode, sequ_read, max_mismatch):
    barcode_count = 0
    barcode_len = len(barcode)
    sequ_read_len = len(sequ_read)
    for i in range(sequ_read_len - barcode_len + 1):
        misatch_count = 0
        for j in range(barcode_len):
            if sequ_read[i + j] != barcode[j]:
                misatch_count += 1
            if misatch_count > max_mismatch:
                break
        if misatch_count <= max_mismatch:
            barcode_count += 1
    return barcode_count


# counts number of barcodes for multiple different barcodes and generates a dictionary containg the barcode sequence and corresponding count number
def count_barcodes(barcodes, sequ_read, max_mismatch):
    barcodes_count = {}
    for barcode in barcodes:
        barcodes_count[barcode] = naive_v2(barcode, sequ_read, max_mismatch)
    return barcodes_count


# Loops through read file once and compares each barcode instead of vise versa
# Should be a little bit faster than V1
# !!!This only works if all barcodes have the same length!!!
def count_barcodes_v2(barcodes, sequ_read, max_mismatch):
    barcodes_count = {}
    barcode_len = len(barcodes[0])
    sequ_read_len = len(sequ_read)
    # Initialize it so that the order is kept
    for barcode in barcodes:
        barcodes_count.setdefault(barcode, 0)
    for i in range(sequ_read_len - barcode_len + 1):
        for barcode in barcodes:
            misatch_count = 0
            for j in range(barcode_len):
                if sequ_read[i + j] != barcode[j]:
                    misatch_count += 1
                if misatch_count > max_mismatch:
                    break
            if misatch_count <= max_mismatch:
                barcodes_count[barcode] += 1
    return barcodes_count


def write_to_csv(filename, results):
    with open(filename,'w', newline='') as f:
        w = csv.writer(f)
        w.writerows(results.items())


def main():
    filename, barcode_file, mismatch = parse_args(sys.argv)
    reference_file = read_genome(filename)
    barcodes = read_barcodes(barcode_file)
    barcode_count = count_barcodes(barcodes, reference_file, mismatch)
    write_to_csv("%s_result.csv" % filename, barcode_count)
    print(" Table of barcodes %s " % (barcode_count))


if __name__ == '__main__':
    main()
