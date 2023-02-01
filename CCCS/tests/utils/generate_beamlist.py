"""generate_beamlet.py

Create new beamlist files with specified spacing


"""
import os.path
import datetime
import argparse

def frange(x, y, jump=1.0):
    """Range for floats."""
    i = 0.0
    x = float(x)  # Prevent yielding integers.
    x0 = x
    epsilon = jump / 2.0
    yield x  # yield always first value
    while x + epsilon < y:
        i += 1.0
        x = x0 + i * jump
        yield x

def get_beam_spec(gantry, couch):
    return '0.0 {:>7.3f} {:>7.3f}\n'.format(gantry, couch)

def get_documentation(cstart, cend, cspace, gstart, gend, gspace):
    return \
"""

##############################################################################################################
# Mandatory Format:      gantry(azimuth,theta [deg]) couch(zenith,phi:deg)           as "%f %f %f"
# Optional Modules (order insignificant):
#     <MODULE>           <DESCRIPTION OF SPECIFIED VALUES>                               <FORMAT SPECIFIER>
#   - Manual Isocenter:  isocenter.x isocenter.y isocenter.z [units: cm]             as "[...] iso: %f %f %f"
#   - Manual SAD:        SAD [units: cm]                                             as "[...] sad: %f"
#-------------------------------------------------------------------------------------------------------------
# note: spaces delimit fields, and any number of spaces may be used between fields
#
# Uses VarianIEC motion scale:
#   azimuth - Gantry/Z-axis rotation      - (0 defined as entering patient anterior, gantry above)
#   zenith  - Couch/non-coplanar rotation - (0 defined as coplanar beam, couch perp. to linac body)
#
# EXAMPLES FOR PATIENT ORIENTATION HFS:
#   azimuth = 0           : entering from patient anterior
#   azimuth = 90          : entering from patient left
#   azimuth = 180         : entering from patient posterior
#   azimuth = 270         : entering from patient right
#
#   0    < zenith  < 180  : CW couch-kick (top view)
#   -180 < zenith  <   0  : CCW couch-kick (top view)
##############################################################################################################

# Created using generate_beamlist.py on {date!s} with the following parameters:
#   couch  start:   {cstart:f}
#   couch  end:     {cend:f}
#   couch  spacing: {cspace:f}
#   ------------------------
#   gantry start:   {gstart:f}
#   gantry end:     {gend:f}
#   gantry spacing: {gspace:f}
""".format(
    date=datetime.date.today().strftime("%d%b%Y"),
    cstart=cstart,
    cend=cend,
    cspace=cspace,
    gstart=gstart,
    gend=gend,
    gspace=gspace
)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Generates beamlist files compatable with dosecalc_gpu at a discrete angular spacing")
    parser.add_argument("-v", "--verbose", action='store_true', help="print verbosely", default=False)
    parser.add_argument("--out", type=str, help="outfile name", default="beamlist")
    parser.add_argument("--spacing-couch", type=float, help="spacing of couch angles [units: degrees]", default=6)
    parser.add_argument("--spacing-gantry", type=float, help="spacing of gantry angles [units: degrees]", default=10)

    # parse args
    args = parser.parse_args()
    verbose = args.verbose
    outfile = os.path.splitext(args.out)[0] + '.txt'
    cspace = args.spacing_couch
    gspace = args.spacing_gantry

    # constants
    cstart = -90
    cend   =  90
    gstart =   0
    gend   = 360

    # open file
    N = 0
    with open(outfile, mode='w') as f:
        for c in frange(cstart, cend, cspace):
            for g in frange(gstart, gend, gspace):
                f.write(get_beam_spec(g, c))
                N = N+1
        f.write(get_documentation(cstart, cend, cspace, gstart, gend, gspace))

    if (verbose): print("wrote {:d} beams to file: \"{!s}\"".format(N, outfile))
