import re

class SAMRecord:

    def __init__(self, record):
        self.record = record

        self.record_sep = record.split('\t')
        self.umi = self.record_sep[0].split(':')[-1]
        self.pos = int(self.record_sep[3])
        self.cigar = re.findall(r'\d+[MIDNS]', self.record_sep[5])
        self.flag = int(self.record_sep[1])

    def isReverse(self):
        return True if (self.flag & 16) else False

    def isMapped(self):
        return False if (self.flag & 4) else True

    def update_pos(self, offset):
        self.pos += offset
