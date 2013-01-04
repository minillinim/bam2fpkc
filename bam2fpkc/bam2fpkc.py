#!/usr/bin/env python
###############################################################################
#                                                                             #
#    bam2fpkc.py                                                              #
#                                                                             #
#    Wraps coarse workflows                                                   #
#                                                                             #
#    Copyright (C) Michael Imelfort                                           #
#                                                                             #
###############################################################################
#                                                                             #
# 888                              .d8888b.  .d888          888               #
# 888                             d88P  Y88 d88P"           888               #
# 888                                    88 888             888               #
# 88888b.   8888b.  88888b.d88b.       .d88 888888 88888b.  888  888  .d8888b # 
# 888 "88b     "88b 888 "888 "88b  .od888P" 888    888 "88b 888 .88P d88P"    #
# 888  888 .d888888 888  888  888 d88P"     888    888  888 888888K  888      #
# 888 d88P 888  888 888  888  888 888"      888    888 d88P 888 "88b Y88b.    #
# 88888P"  "Y888888 888  888  888 88888888  888    88888P"  888  888  "Y8888P #
#                                                  888                        #
#                                                  888                        #
#                                                  888                        #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

__author__ = "Michael Imelfort"
__copyright__ = "Copyright 2012"
__credits__ = ["Michael Imelfort"]
__license__ = "GPL3"
__version__ = "0.2.0"
__maintainer__ = "Michael Imelfort"
__email__ = "mike@mikeimelfort.com"
__status__ = "Development"

###############################################################################

import argparse
import sys
import pysam
import numpy as np
np.seterr(all='raise')  

# bam2fpkc imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class CouldNotOpenMappingException(BaseException): pass
class MappingNotOpenException(BaseException): pass
class BinLengthError(BaseException): pass
 
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Bam2fpkcOptionsParser():
    def __init__(self): pass
    
    def parseOptions(self, options):
        
        # get a list of contig headers
        CP = ContigParser()
        try:
            with open(options.contigs) as con_fh:
                headers = CP.getHeaders(con_fh)
        except:
             print "Could not parse contig file: %s %s" % (options.contigs,sys.exc_info()[0])
             raise
         
        # we'll need a BamParser no matter what! 
        BP = BamParser()
        BP.openBam(options.bam)
        
        if(options.subparser_name == 'contig'):
            # Make fpkc for a single fasta file
            BP.parseBam()

        elif(options.subparser_name == 'bin'):
            # Make mfpkc for a set of bins
            # we should have a tab sep'd bins file which looks like:
            #
            # cid -> bid
            #
            try:
                with open(options.bins) as bin_fh:
                    bin_assignments = self.parseBinsFile(bin_fh)
            except:
                 print "Could not parse bins file: %s %s" % (options.bins,sys.exc_info()[0])
                 raise

            # make sure that all contigs are assigned to a bin
            if len(bin_assignments.keys()) != len(headers):
                raise BinLengthError("bin file does not match headers",
                                     "B: %d != H: %d " % (len(bin_assignments.keys()),
                                                          len(headers)
                                                          )
                                     )
            BP.parseBam()
        
        # clean up
        BP.closeBam()
        
        return 0

    def parseBinsFile(self, bfh):
        """Extract contig Id to bins assignments
        
        bfh should be a tab sep'd bins file which looks like:
        
        # header
        cid -> bid
        cid -> bid
        ...
        
        """
        bin_assignments = {}
        for line in bfh:
            line = line.rstrip()
            if len(line) > 0:
                if line[0] != '#':
                    parts = line.split('\t')
                    bin_assignments[parts[0]] = parts[1]
                    
        return bin_assignments
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BamParser:
    def __init__(self):
        self.bamFile = None

    def openBam(self, bf):
        """Open the bamfile and set the file handle"""
        try:
            self.bamFile = pysam.Samfile(bf, 'rb')
        except:
            raise CouldNotOpenMappingException("Unable to open mapping file file: %s -- did you supply a SAM file?" % bf)        
    
    def closeBam(self):
        """close the open bamfile"""
        if self.bamFile is not None:
            self.bamFile.close()
        
    def parseBam(self):
        """Parse a bam file (handle) """
        if self.bamFile is None:
            raise MappingNotOpenException("No file handle to work on!")
        for reference, length in zip(self.bamFile.references, self.bamFile.lengths):
            fc = FragCounter()
            self.bamFile.fetch(reference, 0, length, callback = fc )
            print reference, fc.count

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class FragCounter:
    """AUX: Call back for counting aligned reads
    Used in conjunction with pysam.fetch
    """
    def __init__(self):
        self.count = 0
        
    def __call__(self, alignedRead):
        if not alignedRead.is_unmapped:
            self.count += 1

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ContigParser:
    """Main class for reading in and parsing contigs"""
    def __init__(self): pass

    def readfq(self, fp): # this is a generator function
        """https://github.com/lh3"""
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].split()[0], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
                    break
    def getHeaders(self, contigFile):
        """All we want id the headers"""
        headers = []
        for cid,seq,qual in self.readfq(contigFile):
            headers.append(cid)
        return headers
       
###############################################################################
###############################################################################
###############################################################################
###############################################################################
