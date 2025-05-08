#!/usr/bin/env python3
import sys
from pathlib import Path

USAGE = """\
Usage: convert_U_to_DT.py input.pdb

Produces input_T.pdb where:
 - any ALTLOC character (col 17) is set to blank
 - any residue with 3-letter code ' U '→'DT ' (T in DNA)
"""

def convert_pdb(infile, outfile):
    with open(infile, 'r') as fin, open(outfile, 'w') as fout:
        for line in fin:
            # Only process ATOM/HETATM records
            if line.startswith(('ATOM  ', 'HETATM')):
                # 1) clear altLoc (column 17, index 16)
                line = line[:16] + ' ' + line[17:]
                # 2) check resName (columns 18–20, indices 17:20)
                res3 = line[17:20]
                if res3.strip() == 'U':
                    # write 'DT ' so columns 18–20 are D, T, space
                    line = line[:17] + 'DT ' + line[20:]
            fout.write(line)

def main():
    if len(sys.argv) != 2:
        print(USAGE)
        sys.exit(1)
    inp = Path(sys.argv[1])
    if not inp.is_file():
        print(f"Error: file '{inp}' not found.")
        sys.exit(1)
    out = inp.with_name(f"{inp.stem}_T.pdb")
    convert_pdb(str(inp), str(out))
    print(f"Wrote: {out}")

if __name__ == "__main__":
    main()

