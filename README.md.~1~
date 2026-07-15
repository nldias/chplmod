# staeble-wrr


This repository contains 3 Chapel programs which implement the STAEBLE
Lake Evaporation Model.
You will need to install the Chapel programming language: see https://chapel-lang.org/.
The model is described [here](https://doi.org/10.1002/essoar.10511612.1)
A description of each program (a "help") can be printed with

./staeble-n --describe

(and so forth). The programs are

- staeble-n.chpl

  Implements STAEBLE-A, STAEBLE-B and STAEBLE-AB

- staeble-c.chpl

  Implements STAEBLE-C

- staeble-ch.chpl

  Implements STAEBLE-CH

The input file is lakemead-mod.dat. The outputs are staeble-a.out,
staeble-b.out, staeble-ab.out, staeble-c.out and staeble-ch.out, and
will be overwritten when you compile and run the programs.

The programs use 5 modules that are also included in this repository, as follows:

- dgrow.chpl
- sunearth.chpl
- atmgas.chpl
- angles.chpl
- evap.chpl
- ssr.chpl

All source code here comes straight from our research, and is not
directly usable at places other than Lake Mead: for example, altitude,
latitude, and various parameters for radiation models are hard-wired
into the code. For use at other lakes, the corresponding lines (at
least) will have to be modified.

All source code in this repository is released under the GPL license.
I cannot garantee that the code is correct under any circumstances of
use; use it at your own risk (see the GPL license in file GPL.txt); and I cannot
provide support for the code. However, I hope this is useful for
research and/or educational purposes.

Nelson Luís Dias
nldias@ufpr.br

