## Instrument files

This folder contains instrument files for ISIS instruments in both NeXus and McStas formats,
and data tables needed to generate some of these files.

The NeXus files were generated by the [horace_nxs_inst.py](../scripts/horace_nxs_inst.py) in this repository.
Instrument files for LET, MAPS and MERLIN are readable by Horace and can be used by it for
resolution calculations.

The McStas `.instr` file here was used to define the guide and detector components in the NeXus files,
via [McStasScript](https://github.com/PaNOSC-ViNYL/McStasScript).
To run the McStas instruments you also need the ISIS
[moderator files](https://github.com/ISISNeutronMuon/mcstas/blob/master/docs/ISIS-Moderator-files.md).