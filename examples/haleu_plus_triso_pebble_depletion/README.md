# TL;DR
To run:

```bash
python kp-inspired-haleu-plus-triso-pebble.py
```

Now, sometimes, you may encounter situations where you don't want to use 
all the threads in the computer to prevent overheating. So we may want 
to use less threads. To do so, set the OMP\_NUM\_THREADS variable to 
a suitable number and run the script.

For example, if you want 10 threads:

```bash
OMP_NUM_THREADS=10 python kp-inspired-pebble-fresh-triso.py
```

For post processing, it generates a jpg file and csv file with the keff 
over time for this pebble.

```bash
python postprocessing.py
```

# Description

Basically, this is the most basic depletion for a triso pebble using OpenMC.
OpenMC currently uses surface tracking as opposed to delta tracking so this 
is extremely inefficient. However, syntax should be okay.

The fuel here is HALEU like, 19.9 atom% enriched U235. We want to test the 
extent to which U236 contamination affects fuel performance over time. 
For this, we only simulate a single pebble since simulation of triso 
pebbles in OpenMC is so hard. We then perform depletion upon it, adapting the 
OpenMC pin cell depletion for this. 

Do note that the PWR depletion chain is used as OpenMC does not have 
prebuilt depletion chains for molten salt reactors (MSR)s or 
fluoride salt cooled high temperature reactors (FHR)s. Hence, this depletion 
analysis is quite approximate.

# Future work 

In future, we want to test the effect of various U236 to U235 ratios on 
the depletion over time. From the following reference,

```bibtex
@article{zabunouglu2008determination,
  title={Determination of enrichment of recycle uranium fuels for different burnup values},
  author={Zabuno{\u{g}}lu, Okan H},
  journal={Annals of Nuclear Energy},
  volume={35},
  number={2},
  pages={285--290},
  year={2008},
  publisher={Elsevier}
}
```
We may expect about 0.72 weight percent of U236 and 0.77 weight percent of 
U235 at a burnup of 50,000 MW-day-thermal/ton U. At about 33,000 MW-day-thermal 
/ton U, the U236 was about 0.47 weight percent while U235 was about 
0.855 weight percent. The correlation between U235 and U236 composition 
over time with respect to burnup was approximately linear. 

This means we may have at most a 1:1 ratio of U235 to U236 within spent 
fuel. It is good to see the effect of various ratios of U235 to U236 
on fuel performance for the triso pebble.
