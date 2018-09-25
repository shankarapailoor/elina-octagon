# Overview
The program ```test.cpp``` demonstrates an issue we have found when using the octagon domain of ELINA. ```test.cpp``` applies a sequence of affine transformation followed by meets with (>= 0) to an initial abstract domain (essentially fully connected layer + relu). We find ```test.cpp``` will have different behaviors if ```elina_lincons0_array_clear```` is commented in line 353 of ```test.cpp```. If the line is commented, then the second linear transformation yields bottom, but if it isn't then we get nonempty constraints. Strangely we find that the constraints just before this transformation are the same regardless of whether we clear the array.

We show this by printing the linear constraints from the transformed domain, which we get from ```elina_abstract0_to_lincons_array```. Each transformation should produce two sets of constraints: 1 for the transformation and 1 for the meet with constraint (>= 0). We have provided a script to verify and usage instructions below. 
  
## Build test.cpp
g++ test.cpp -lgmp -lmpfr -lelinaux -loptoct

## Run test.cpp
```bash
$ ./a.out 1 > cleared.txt // clears constraints after meet with > 0
```
```bash
$ ./a.out 0 > notCleared.txt // doesn't clear constraints after meet
```

## Run Comparison Script
```bash
$ python compareConstraints.py cleared.txt notCleared.txt
```

This script should output nothing if the constraints are the same. 


## Expected output

### Without calling clear
```
19990: x94 - x99 + 0.067299343195221406067 >= 0
19991: -x95 - x99 + 0.067299343195221406067 >= 0
19992: x95 - x99 + 0.067299343195221406067 >= 0
19993: -x96 - x99 + 0.36762870847178780176 >= 0
19994: x96 - x99 + 0.067299343195221406067 >= 0
19995: -x97 - x99 + 0.22211343417206391715 >= 0
19996: x97 - x99 + 0.067299343195221406067 >= 0
19997: -x98 - x99 + 0.067299343195221406067 >= 0
19998: x98 - x99 + 0.067299343195221406067 >= 0
19999: -x99 + 0.067299343195221406067 >= 0
PRINTED AFTER RELU in Layer 0
opt_oct_assign_linexpr_array size: 100o->dim = 100, size = 100
opt_hmat_strong_closure returned positive
Layer: 1 before ReLU
0 disjuncts
Layer: 1 after ReLU
0 disjuncts
Layer: 2 before ReLU
0 disjuncts
Layer: 2 after ReLU
0 disjuncts
```

### Clearing array
180: x9 + 9.3590256561560938309 >= 0
181: -x0 - x9 + 16.704760200984697605 >= 0
182: x0 - x9 + 16.459121417491466844 >= 0
183: -x1 - x9 + 18.609792410238117811 >= 0
184: x1 - x9 + 16.362381296435579259 >= 0
185: -x2 - x9 + 17.848708115739377434 >= 0
186: x2 - x9 + 18.298190135296426462 >= 0
187: -x3 - x9 + 18.050018503093863132 >= 0
188: x3 - x9 + 17.5431961445769673 >= 0
189: -x4 - x9 + 16.736410361856258078 >= 0
190: x4 - x9 + 20.169857810657251918 >= 0
191: -x5 - x9 + 17.99690804300194813 >= 0
192: x5 - x9 + 21.657778485156537585 >= 0
193: -x6 - x9 + 17.485219493301762839 >= 0
194: x6 - x9 + 19.676477711972690799 >= 0
195: -x7 - x9 + 15.160692565676692212 >= 0
196: x7 - x9 + 19.339093799103338256 >= 0
197: -x8 - x9 + 24.404344589898173013 >= 0
198: x8 - x9 + 10.142582519017700804 >= 0
199: -x9 + 8.2052457123764419578 >= 0
PRINTED AFTER RELU in Layer 2
Clearing

