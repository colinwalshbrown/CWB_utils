#!/usr/bin/env python
#
# Generate simple output for Gaussian/09 from an un-optimized pdb file
#

import sys
import re

def _main(args):
    chk_file = ""
    theory_level = ""
    gauss_options = ""
    conv_criterion = ""
    pdb = ""
    if len (args) < 7:
        print "USAGE:"
        print "pdb2gaussian.py <check_file> <n_processors> <theory_level> <gauss_options> <conv_criterion> <charge> <.pdb file>"
        sys.exit(0)
    else:
        chk_file,n_processors,theory_level, gauss_options, conv_criterion, charge, pdb = args

    print "%Chk=" + chk_file
    print "%NProcShared=" + n_processors
    print "# " + theory_level + " " + gauss_options + " SCF=" + conv_criterion
    print "\nscan rotamers\n\n",

    # Parse the charge input from standard notation (e.g. "+1" "-2")
    charge_srch = re.search("(\+|-)?(\S+)",charge)
    charge_polarity = 0

    if charge_srch.group(1) == "+":
        charge_polarity = 1
    elif charge_srch.group(1) == "-":
        charge_polarity = -1
        
    charge_multiplicity = charge_srch.group(2)
    if charge_multiplicity == "0":
        charge_multiplicity = "1"
        
    print "%s %s" % (charge_polarity,charge_multiplicity)

    # Parse and output Cartesian atom coords, maintaining PDB order
    for ln in open(pdb):
        pdb_ln_srch = re.match("ATOM\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+(\S+)\s+\S+\s+\S+\s+(\S+)",ln)
        if pdb_ln_srch:
            print "%s\t%s\t%s\t%s" % (pdb_ln_srch.group(4),pdb_ln_srch.group(1),pdb_ln_srch.group(2),pdb_ln_srch.group(3))

    # Fix Backbone Torsion Angles - I assume this will be the same for all AAs?
    print ""
    print "1 13 14 15  -150.00 F"
    print "13 14 15 7  150.00 F"
    print ""
                

if __name__ == "__main__":
    _main(sys.argv[1:])

