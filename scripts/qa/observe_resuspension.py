#!/usr/bin/env python3
"""Resuspension behavior observation: particle height statistics over time.

Usage: python3 scripts/qa/observe_resuspension.py <run_out_dir> [...]
Reads particles_<step>.csv snapshots written by the runner and prints
mean height, settled fraction (y<4), suspended fraction (y>16), upper-half
fraction (y>48) and max height per snapshot.
"""
import csv, glob, re, sys

for case_dir in sys.argv[1:]:
    files = sorted(glob.glob(f'{case_dir}/particles_*.csv'),
                   key=lambda p: int(re.search(r'_(\d+)\.csv', p).group(1)))
    print(f'== {case_dir}  (step | mean_y | settled<4 | susp>16 | upper>48 | max_y)')
    for f in files:
        step = int(re.search(r'_(\d+)\.csv', f).group(1))
        ys = [float(r['y']) for r in csv.DictReader(open(f))]
        n = len(ys)
        print('  %6d | %6.2f | %5.1f%% | %5.1f%% | %5.1f%% | %6.1f' % (
            step, sum(ys)/n, 100*sum(1 for y in ys if y < 4)/n,
            100*sum(1 for y in ys if y > 16)/n,
            100*sum(1 for y in ys if y > 48)/n, max(ys)))
