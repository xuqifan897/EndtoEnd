#!/usr/bin/env python3
"""Script that converts a dosecalc beamlist.txt file (or directory of files) to the latest format"""

import sys, os
import argparse

if __name__ == '__main__':
    # get input(s)
    if len(sys.argv) < 2:
        print('usage: {} <file-path>|<dir-path>'.format(sys.argv[0]))
        sys.exit(1)

    parser = argparse.ArgumentParser(description='Convert legacy dosecalc beamlists to latest format')
    parser.add_argument('input', nargs='+', type=str, help='input file or directory')
    parser.add_argument('--out', '-o', type=str, default='converted', help='output directory')
    args = parser.parse_args()

    inputs = []
    for target in args.input:
        if os.path.isfile(target):
            inputs.append(target)
        elif os.path.isdir(target):
            for f in [x for x in os.listdir(target) if os.path.isfile(os.path.join(target,x))]:
                inputs.append(os.path.join(target, f))

    # convert inputs
    os.makedirs(args.out, exist_ok=True)
    for target in inputs:
        output = os.path.join(args.out, os.path.basename(target))
        with open(output, 'w') as fout:
            with open(target, 'r') as fin:
                for line in fin:
                    tokens = list(line.rstrip('\n').split(' '))
                    if len(tokens) == 6:
                        rad, gantry, couch, isox, isoy, isoz = tokens
                        fout.write("{:<8g} {:<8g} iso: {:<8g} {:<8g} {:<8g}\n".format(
                            float(gantry), float(couch), float(isox), float(isoy), float(isoz) ))
                    elif (len(tokens) == 3):
                        rad, gantry, couch, isox, isoy, isoz = tokens
                        fout.write("{:<8g} {:<8g}\n".format( float(gantry), float(couch) ))

        print('"{!s}" --> "{!s}"'.format(target, output))
