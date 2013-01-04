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
import os

# bam2fpkc imports
 
###############################################################################
###############################################################################
###############################################################################
###############################################################################

class Sam2fpkcOptionsParser():
    def __init__(self): pass
    
    def parseOptions(self, options ):

        if(options.subparser_name == 'single'):
            # Make fpkc for a single fasta file
            pass
        elif(options.subparser_name == 'bin'):
            # Make mfpkc for a set of bins
            pass
        return 0

###############################################################################
###############################################################################
###############################################################################
###############################################################################
