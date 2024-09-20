#
# convert a TriplClust output CSV file to gnuplot format
#
# Author: Christoph Dalitz
#         07 Nov 2022
#
# modified by Juliane Arning
#         02. Sep 2024

import sys
import os.path

usage = "Usage: " + os.path.basename(sys.argv[0]) + " [<inflie>] [-o <outfile>]\n"

opt_infile = None
opt_outfile = None
i = 1
while i < len(sys.argv):
    if sys.argv[i] == "-o":
        i += 1;
        opt_inflie = sys.argv[i]
    elif sys.argv[i][0] == "-":
        sys.stderr.write(usage)
        sys.exit(1)
    else:
        opt_infile = sys.argv[i]
    i += 1
if opt_infile and not os.path.exists(opt_infile):
    sys.stderr.write("Cannot find file '" + opt_infile + "'\n")
    sys.stderr.write(usage)
    sys.exit(1)

# data structure for labelled data
class Point:
    def __init__(self, x, y, z, label):
        self.x = x
        self.y = y
        self.z = z
        self.label = label
    def __lt__(self, other):
        return self.label < other.label
    def __gt__(self, other):
        return self.label > other.label
    def __eq__(self, other):
        return self.label == other.label

# load CSV file into memory
if (opt_infile):
    f = open(opt_infile)
else:
    f = sys.stdin
points = []
labels = set()
for line in f:
    if line.startswith("#"):
        continue
    line = line.rstrip()
    x = line.split(",")
    labels.add(x[3])
    points.append(Point(float(x[0]), float(x[1]), float(x[2]), x[3]))
if opt_infile:
    f. close()

#
# write output
#
if opt_outfile:
    f = open(opt_outfile, "w")
else:
    f = sys.stdout


f.write("unset colorbox\nunset border\nunset xtics\nunset ytics\nunset ztics\n")


# range of data
f.write("set xrange [%f:%f]\n" % (min([p.x for p in points]), max([p.x for p in points])))
f.write("set yrange [%f:%f]\n" % (min([p.y for p in points]), max([p.y for p in points])))
f.write("set zrange [%f:%f]\n" % (min([p.z for p in points]), max([p.z for p in points])))
# loop over labels with colors
first = True
if "-1" in labels:
    f.write("splot  '-' with points lc 'red' title 'noise'")
    first = False
# sort should move '-1' to beginning
labelscurves = sorted([lb for lb in labels if ";" not in lb])
labelsoverlap = sorted([lb for lb in labels if ";" in lb])
i = 0
for lb in labelscurves:
    if lb == "-1": # noise
        continue
    r = float((i * 53) % 19) / 18.0;
    g = float((i * 53) % 7) / 6.0;
    b = float((i * 53) % 5) / 4.0;
    r = int(r * 255);
    g = int(g * 255);
    b = int(b * 255);
    if not first:
        f.write(", ")
    else:
        f.write("splot ")
        first = False
    f.write("'-' with points lc '#%02x%02x%02x' title 'curve %s'" % (r,g,b,lb))
    i = i + 1
for lb in labelsoverlap:
    r = float((i * 23) % 19) / 18.0;
    g = float((i * 23) % 7) / 6.0;
    b = float((i * 23) % 5) / 4.0;
    r = int(r * 255);
    g = int(g * 255);
    b = int(b * 255);
    f.write(", '-' with points lc '#%02x%02x%02x' title 'overlap %s'" % (r,g,b,lb))
    i = i + 1
f.write(",\n")
# group points by labels
for lb in labelscurves:
    for p in points:
        if p.label == lb:
            f.write("%f %f %f\n" % (p.x, p.y, p.z))
    f.write("e\n")
for lb in labelsoverlap:
    for p in points:
        if p.label == lb:
            f.write("%f %f %f\n" % (p.x, p.y, p.z))
    f.write("e\n")

f.write("pause mouse keypress\n")
if opt_outfile:
    f.close()
