"""
Gir Lion Genomic Diversity Analysis
Authors: Akshay Krishnan Pushparaj, Malarmathi Muthukumar
Institution: Dr. N.G.P. Arts and Science College, Coimbatore
Date: June 2026

Data sources:
- Mitra et al. 2019 (bioRxiv doi:10.1101/549790)
- Cho et al. 2013 (Nature Communications)
- Additional felid genomes (see manuscript Table 1)

Raw Asiatic lion data: NCBI BioProject PRJNA379375
"""

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

# Published data (Mitra et al. 2019, Cho et al. 2013, etc.)
felid_diversity = {
    'Asiatic_lion_Gir': 0.000276,
    'African_lion': 0.000480,
    'White_African_lion': 0.000350,
    'Snow_leopard': 0.000230,
    'Eurasian_lynx': 0.000290,
    'Amur_leopard': 0.000480,
    'Amur_tiger': 0.000490,
    'Bengal_tiger': 0.000730,
    'Cheetah': 0.000150,
    'Domestic_cat': 0.000900,
}

# Key statistics (Mitra et al. 2019)
atul_snvs_total = 745184
atul_homozygous = 221397
atul_heterozygous = 22790

# Inbreeding coefficient proxy
F_ROH = atul_homozygous / (atul_homozygous + atul_heterozygous)
print(f"F_ROH proxy (Asiatic lion): {F_ROH:.3f}")

# Ne estimation
mu = 1e-8  # mutation rate per site per generation
He_asiatic = felid_diversity['Asiatic_lion_Gir']
He_african = felid_diversity['African_lion']
Ne_asiatic = He_asiatic / (4 * mu * (1 - He_asiatic))
Ne_african = He_african / (4 * mu * (1 - He_african))
print(f"Ne (Asiatic): {Ne_asiatic:,.0f}")
print(f"Ne (African): {Ne_african:,.0f}")

# Z-test for diversity difference
genome_size = 2.7e9
se = np.sqrt(He_asiatic*(1-He_asiatic)/genome_size + 
             He_african*(1-He_african)/genome_size)
z = (He_african - He_asiatic) / se
p = 2 * (1 - stats.norm.cdf(abs(z)))
print(f"Z = {z:.2f}, p = {p:.2e}")

# Kuno diversity projections
t = np.arange(0, 21)
H0 = 0.000276
H_no_mgmt = H0 * (1 - 1/(2*10))**t
H_managed = H0 * (1 - 1/(2*50))**t
H_rescue = H0 * (1 - 1/(2*150))**t
print(f"\nKuno retention (10 gen):")
print(f"No management: {H_no_mgmt[10]/H0*100:.1f}%")
print(f"Managed: {H_managed[10]/H0*100:.1f}%")
print(f"Rescue: {H_rescue[10]/H0*100:.1f}%")
