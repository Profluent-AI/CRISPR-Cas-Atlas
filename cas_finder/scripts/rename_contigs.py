#!/usr/bin/env python3

import sys

sample_id = sys.argv[1]
num = 0
for line in sys.stdin:
    if line[0] == ">":
        num += 1
        sys.stdout.write(">" + sample_id + "_" + str(num) + " " + line.split()[0][1:] + "\n")    
    else:
        sys.stdout.write(line)
