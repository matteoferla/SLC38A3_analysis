# SLC38A3 Analysis

> Please note, this repo does not contain mutant information

SLC38A3 is Sodium-coupled neutral amino acid transporter 3 ([Uniprot: S38A3_HUMAN](https://www.uniprot.org/uniprot/Q99624)).

There is no crystal structure available, but there is one for SLC38A9, PDB: 6C08
— described well in [Lei et al. 2018](https://www.nature.com/articles/s41594-018-0072-2).

## Reference

With the exception of the Extracellular part, Phyre and ITasser produced nearly identical structures.
There is no support, or need, for the extracellular fluff, so it was removed.
Swissmodel has no useful model, possibly due to this.
The Phyre model was chosen because it deviates less from the reference PDB: 6C08.
With the hindsight knowledge that the extracellular fluff was causing low similarity, 
threading EM-fluff stripped SLC38A3 against PDB: 6C08 
with `pyrosetta.rosetta.protocols.comparative_modeling.ThreadingMover` would have been better for the minimation step.
This is because membrane protein can thread weirdly in I-Tasser and Phyre.

Being a membrane protein the RosettaMP framework was used, but in conjuction with an initial constraint of the 6C08 electron density.
I.e. the SLC38A3 model was aligned against the OPM model of 6C08, membrane-ified, and then moved to the 6C08 structure as located in 
the 2Fc-Fo map.

It should be noted that `ref2015` scorefunction was used for initial steps, but subsequently `franklin2019` was.
The `ref2015` can be used with the membrane framework (cf. [Alford et al. 2020](https://www.sciencedirect.com/science/article/pii/S000634952030237X)).
Although RosettaMP was released with `mpframework_smooth_fa_2012` scorefunction, the latest scorefunction `franklin2019` is better.

...

## Ligands

...