#!/usr/bin/env python3
"""Independent integer seed/normal inputs for base-R CCA reference fixtures.

This is a unit fixture namespace, never a calibration or confirmation stream.
Seed paths follow resample4s seed-path/v1; Python owns all integer arithmetic,
base R owns QR/SVD and Wilks statistics. No production Scala helper is called.
"""
import csv
import math
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'modules/inference/fixtures/canonical-rank-v2'
MASK = (1 << 64) - 1
ROOT_SEED = 871254

def mix(x):
    x &= MASK
    x = ((x ^ (x >> 30)) * 0xbf58476d1ce4e5b9) & MASK
    x = ((x ^ (x >> 27)) * 0x94d049bb133111eb) & MASK
    return x ^ (x >> 31)

def child(seed, domain, ordinal):
    x = mix(seed ^ 0x736565642d706174)
    x = mix(x + 0x632be59bd9b4e019 + 1)
    x = mix(x + 0x8cb92ba72f3d8dd7 + domain)
    return mix(x + 0x9e3779b185ebca87 + ordinal)

def signed(x):
    return x if x < 1 << 63 else x - (1 << 64)

def normal(seed, count):
    state = seed
    values = []
    while len(values) < count:
        uniforms = []
        for _ in range(2):
            state = (state + 0x9e3779b97f4a7c15) & MASK
            uniforms.append(((mix(state) >> 12) + 0.5) / (1 << 52))
        radius = math.sqrt(-2 * math.log(uniforms[0]))
        angle = 2 * math.pi * uniforms[1]
        values.extend((radius * math.cos(angle), radius * math.sin(angle)))
    return values[:count]

OUT.mkdir(parents=True, exist_ok=True)
with (OUT / 'gaussian-input.tsv').open('w') as stream:
    writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
    writer.writerow(('step','replicate','block','rows','columns','hypothesis_seed','value'))
    for step in range(2):
        hypothesis = child(ROOT_SEED, 1011, step)
        for replicate in range(39):
            local = child(hypothesis, 1012, replicate)
            for block, columns in enumerate((2-step,3)):
                for value in normal(child(local,1013,block), 6*columns):
                    writer.writerow((step,replicate,block,6,columns,signed(hypothesis),repr(value)))
subprocess.run(['Rscript', str(Path(__file__).with_suffix('.R')), str(ROOT)], check=True)
