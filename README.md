Manipulating an Instrumental Variable in an Observational Study of
Premature Babies: Design, Bounds, and Inference
================

Code to reproduce simulation results, data analysis results and figures
from the paper “Manipulating an Instrumental Variable in an
Observational Study of Premature Babies: Design, Bounds, and Inference.”

## Organization

### code_simulations

`sim1.R` contains auxiliary functions and code used to produce results
from Simulation Study B.1 in Supplemental Material B. The expected
run-time for one simulation study using a standard desktop machine is
about 10 hours.

`sim2.R` contains auxiliary functions and code used to produce results
from Simulation Study B.2 in Supplemental Material B. The expected
run-time for one simulation study using a standard desktop machine is
about 8 hours.

`sim3.R` contains auxiliary functions and code used to produce results
from Simulation Study B.3 in Supplemental Material B. The expected
run-time for one simulation study using a standard desktop machine is
about 10 hours.

### code_matching

`util_matching.R` contains auxiliary functions used to create matched
samples M0 in `m0.R`, M1 in `m1.R` and M2 in `M2.R` from Section 2.2.

`m0.R` contains code used to create matched sample M0 from Section 2.2.
The expected run-time to construct one matched sample using a standard
desktop machine is about 8 hours.

`M1.R` contains code used to create matched sample M1 from Section 2.2.
The expected run-time to construct one matched sample using a standard
desktop machine is about 10 hours.

`M2.R` contains code used to create matched sample M2 from Section 2.2
and Figure S1 from Supplemental Material C.1. The expected run-time to
construct one matched sample using a standard desktop machine is about
10 hours.

`cpt.R` contains auxiliary functions and code used to implement the
Classification Permutation Test (CPT) from Section 5.1 and produce
Figure S2(a) in Supplemental Material C.1. The expected run-time to run
the CPT using a standard desktop machine is about 2 hours.

### code_dataanalysis

`helper_functions.R` contains auxiliary functions used for biased
randomization inference in
`M2_primary_analysis_randomization&biased_randomization_inference.R`
from Sections 5.2 and 5.3.

`M0_primary_analysis_randomization_Inference.R` contains code used to
produce results for randomization inference in Section 5.2 and Figure S4
in Supplemental Material C.2.

`M2_primary_analysis_randomization&biased_randomization_inference.R`
contains code used to produce results for biased randomization in
Sections 5.2 and 5.3. It also contains code used to generate Figure 1 in
Section 5.2; Figure 2 in Section 5.3; Figure S2(b) in Supplemental
Material C.1; and Figure S3 in Supplemental Material C.2.
