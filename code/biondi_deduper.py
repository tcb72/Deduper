import re
import argparse

def get_index_seq(QNAME):
    '''
    Input: str
    Output: str

    This function takes in the first column of the SAM file record, and returns the index sequence at the end of this column.
    '''
    return(QNAME.split(':')[-1])

def isReverse(FLAG):
    '''
    Input: int
    Output: bool

    This function takes in the second column of the SAM file record (the bit flag), and returns whether it is the reverse strand (or not.)
    '''
    return True if (FLAG & 16) else False

def isMapped(FLAG):
    '''
    Input: int
    Output: bool

    This function takes in the second column of the SAM file record (the bit flag), and returns whether the read is mapped (or not.)
    '''
    return False if (FLAG & 4) else True

def dedupe_umi(filename, isPairedEnd, UMI_file):

    '''
    Input: string, bool, string
    Output: None

    This function takes in the filename of the sam file, whether its paired end or not, and a file containing a whitelist of UMIs.
    It makes sure that the current record's index sequence is a valid UMI. Then, it adjusts the leftmost position of the record based on the CIGAR string. This includes softclipping.
    Finally, it will write the record to an output file given that another record with the same adjusted position, UMI, and strand has not been seen before.
    '''

    if isPairedEnd:
        print('This program does not support paired-end reads yet! Now exiting.')
        exit()

    # instantiate a dupe counter. will contain how many duplicates there are.
    dupe_counter = 0

    prev_chr = None

    # instantiate a set. will contain each unique combination of UMI, actual position, and strand
    output_set = set()

    # store UMI whitelist in list

    if UMI_file is not None:
        with open(UMI_file, 'r') as f:
            valid_UMIs = f.read().splitlines()


    # create an output file name string
    output_file_name = filename.split('.sam')[0] + '_deduped.sam'

    # open files for reading (sam) and writing (output file)

    with open(filename, 'r') as f, open(output_file_name,'a') as fo:

        # for each record
        for line in f:
            # if line starts with @, write to output, and move on to next record
            if line.startswith('@'):
                fo.write(line)
                continue

            # put record into list
            cols = line.split('\t')

            # get index sequence of current record
            index_seq = get_index_seq(cols[0])

            # if the index_seq is a valid UMI
            if index_seq in valid_UMIs:

                # get the bit flag using second column
                bit_flag = int(cols[1])

                # get whether its forward or reverse strand via bit flag
                isReverseStrand = isReverse(bit_flag)

                # only process if read is mapped
                if isMapped(bit_flag):
                    # instantiate offset
                    offset = 0

                    # split CIGAR string by the 5 letters: M, I, D, N, S. \d+ allows for 1 or more numbers prior to the letter
                    CIGAR_split = re.findall(r'\d+[MIDNS]', cols[5])

                    # get the position listed in the record
                    sam_pos = int(cols[3])

                    # if forward strand strand and 'S' is in the leftmost element of the CIGAR string list
                    # adjust offset by subtracting
                    if isReverseStrand == False:
                        if 'S' in CIGAR_split[0]:
                            offset -= int(CIGAR_split[0][:-1])

                    # if reverse strand
                    else:
                        # loop over CIGAR string list
                        # add to the SAM position if
                        for index,curr_cigar in enumerate(CIGAR_split):
                            if 'S' in curr_cigar:
                                if index != 0:
                                    offset += int(curr_cigar[:-1])
                                else:
                                    pass
                            elif ('M' in curr_cigar) | ('D' in curr_cigar) | ('N' in curr_cigar):
                                offset += int(curr_cigar[:-1])

                    # create an "actual" (adjusted position) which is the position listed in the SAM file added by the offset
                    actual_pos = sam_pos + offset
                # if the current (UMI, REAL position, strand) is in the set, its a duplicate
                if (index_seq, actual_pos, isReverseStrand) in output_set:
                    dupe_counter += 1
                # otherwise, write out to file
                else:
                    output_set.add((index_seq, actual_pos, isReverseStrand))
                    fo.write(line)

    # print how many duplicates were removed
    print('Successfully deleted ' + str(dupe_counter) + ' PCR duplicates.')


def dedupe_randomer(filename, isPairedEnd):

    '''
    Input: string, bool
    Output: None

    This function takes in the filename of the sam file, and whether its paired-end or not.
    It makes sure that the current record's index sequence is a valid randomer (no N nucleotides.) Then, it adjusts the leftmost position of the record based on the CIGAR string. This includes softclipping.
    Finally, it will write the record to an output file given that another record with the same adjusted position, UMI, and strand has not been seen before.
    '''

    if isPairedEnd:
        print('This program does not support paired-end reads yet! Now exiting.')
        exit()

    # instantiate a dupe counter. will contain how many duplicates there are.
    dupe_counter = 0

    prev_chr = None

    # instantiate a set. will contain each unique combination of UMI, actual position, and strand
    output_set = set()

    # create an output file name string
    output_file_name = filename.split('.sam')[0] + '_deduped.sam'

    # open files for reading (sam) and writing (output file)

    with open(filename, 'r') as f, open(output_file_name,'a') as fo:

        # for each record
        for line in f:
            # if line starts with @, write to output, and move on to next record
            if line.startswith('@'):
                fo.write(line)
                continue

            # put record into list
            cols = line.split('\t')

            # get index sequence of current record
            index_seq = get_index_seq(cols[0])

            # get the bit flag using second column
            bit_flag = int(cols[1])

            # get whether its forward or reverse strand via bit flag
            isReverseStrand = isReverse(bit_flag)

            # if valid randomer
            if 'N' not in index_seq:

                # only process if read is mapped
                if isMapped(bit_flag):
                    # instantiate offset
                    offset = 0

                    # split CIGAR string by the 5 letters: M, I, D, N, S. \d+ allows for 1 or more numbers prior to the letter
                    CIGAR_split = re.findall(r'\d+[MIDNS]', cols[5])

                    # get the position listed in the record
                    sam_pos = int(cols[3])

                    # if forward strand strand and 'S' is in the leftmost element of the CIGAR string list
                    # adjust offset by subtracting
                    if isReverseStrand == False:
                        if 'S' in CIGAR_split[0]:
                            offset -= int(CIGAR_split[0][:-1])

                    # if reverse strand
                    else:
                        # loop over CIGAR string list
                        # add to the SAM position if 'S' is in rightmost element of CIGAR string or 'M', 'D', or 'N' is in any element of the CIGAR string
                        for index,curr_cigar in enumerate(CIGAR_split):
                            if 'S' in curr_cigar:
                                if index != 0:
                                    offset += int(curr_cigar[:-1])
                                else:
                                    pass
                            elif ('M' in curr_cigar) | ('D' in curr_cigar) | ('N' in curr_cigar):
                                offset += int(curr_cigar[:-1])

                    # create an "actual" (adjusted position) which is the position listed in the SAM file added by the offset
                    actual_pos = sam_pos + offset
                # if the current (UMI, REAL position, strand) is in the set, its a duplicate
                if (index_seq, actual_pos, isReverseStrand) in output_set:
                    dupe_counter += 1
                # otherwise, write out to file
                else:
                    output_set.add((index_seq, actual_pos, isReverseStrand))
                    fo.write(line)

        # print how many duplicates were removed
        print('Successfully deleted ' + str(dupe_counter) + ' PCR duplicates.')

if __name__ == "__main__":
    # add arguments for argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", required=True, help="Filepath of the sorted SAM file.")
    parser.add_argument("-p", "--paired", action='store_true', help="Designates whether file is paired-end or not.")
    parser.add_argument("-u", "--umi", help="Filepath containing list of whitelisted UMIs.")
    args = parser.parse_args()

    if args.umi is None:
        # if no umi file is given, use randomers
        dedupe_randomer(args.file, args.paired)
    else:
        # otherwise, use UMI whitelist
        dedupe_umi(args.file, args.paired, args.umi)
