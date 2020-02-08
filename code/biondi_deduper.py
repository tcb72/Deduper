import argparse
from record import SAMRecord

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

    # instantiate a set. will contain each unique combination of UMI, actual position, and strand
    output_set = set()

    # store UMI whitelist in list
    valid_UMIs = []
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

            curr_record = SAMRecord(line)
            # if the curr umi is a valid UMI (if whitelist file is given)

            if (curr_record.umi in valid_UMIs) | ((UMI_file is None) & ('N' not in curr_record.umi)):

                # only process if read is mapped
                if curr_record.isMapped():
                    # instantiate offset
                    offset = 0

                    # split CIGAR string by the 5 letters: M, I, D, N, S. \d+ allows for 1 or more numbers prior to the letter
                    CIGAR_split = curr_record.cigar

                    # if forward strand strand and 'S' is in the leftmost element of the CIGAR string list
                    # adjust offset by subtracting
                    if curr_record.isReverse() == False:
                        if 'S' in CIGAR_split[0]:
                            offset = int(CIGAR_split[0][:-1])
                            curr_record.update_pos(-1 * offset)

                    # if reverse strand
                    else:
                        # loop over CIGAR string list
                        # add to the SAM position if S and index is 0
                        # add to SAM position if M, D, or N
                        for index,curr_cigar in enumerate(CIGAR_split):
                            if (('S' in curr_cigar) & (index != 0)) | (curr_cigar[-1] in ('M','D','N')):
                                offset = int(curr_cigar[:-1])
                                curr_record.update_pos(offset)

                # if the current (UMI, REAL position, strand) is in the set, its a duplicate
                if (curr_record.umi, curr_record.pos, curr_record.isReverse()) in output_set:
                    dupe_counter += 1
                # otherwise, write out to file
                else:
                    output_set.add((curr_record.umi, curr_record.pos, curr_record.isReverse()))
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

    # otherwise, use UMI whitelist
    dedupe_umi(args.file, args.paired, args.umi)
