= Preliminary Testing 

This contains preliminary test results for HALEU and HALEU plus with 
50% U236 contamination. For HALEU, we have 19.9 atom% enrichment of U235.
However, U236 is difficult to separate from U235. Therefore, I class 
both of them as "light uranium" because these isotopes are lighter compared 
to U238.

For both tests, a 19.9 atom% of light uranium was set. Then, the contamination 
of U236 was set at 0 for fresh HALEU, and 50% (fraction of 0.5) for 
the HALEU plus experiment.

For convenience, the pressurised water reactor (PWR) depletion chain 
was used as it was given by the OPENMC website at the time of writing.
The material temperatures were at 600K also because of convenience as 
pregenerated cross sections at these temperatures in hdf5 format were 
available in OpenMC rather than at FHR operating temperatures especially 
with regards to S alpha beta scattering.

Also, the pebble is a non-annular TRISO pebble from my dissertation. In 
contrast, FHR pebbles have an annular design to keep them buoyant in salt.
This can be adjusted to fit FHR designs more closely in future work.

= Future Tests

Future tests can be run with regard to:

1. adjusting U236 contamination to test its effect on keff over time while 
light uranium enrichment is constant
2. using proper FHR depletion chains or MSR depletion chains
3. using a proper FHR temperature
4. keeping U235 enrichment constant, adjust the ratio of U236 to U238. 
5. using proper geometry for FHR pebbles 
6. HTGR pebble depletion studies 
7. LWR depletion studies 
8. Fast reactor depletion studies



