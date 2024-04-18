# eniius

<img src="https://github.com/mducle/eniius/blob/main/ennius.png" width=150>

A program for processing neutron data in the tradition of
[Homer](http://www.libisis.org/Using_Homer_and_Getting_Started) and
[Horace](http://horace.isis.rl.ac.uk/), `eniius` is a utility for
embedding neutron instrument information using (nx)spe files.


## Aims

Proper analysis of inelastic neutron scattering (INS) data requires taking the
instrument resolution function into account
[[1]](https://ethos.bl.uk/OrderDetails.do?uin=uk.bl.ethos.239595), and this in
turn requires particular information about the neutron spectrometer
(instrument) recording the data. Hitherto this information has been obtained
post-facto from the program(s) performing the resolution convolution.

In the interest of [FAIR data](https://en.wikipedia.org/wiki/FAIR_data) and
to improve usability, this information should be embedded in the raw datafiles
produced by the instrument, and the aim of `eniius` is to provide a utility for
just this purpose.


## Features

* Instrument data is embedded in a [NeXus](https://www.nexusformat.org/) data
  file as a [NXinstrument](https://manual.nexusformat.org/classes/base_classes/NXinstrument.html)
  group.
* The embedded instrument is at a minimum enough to be imported into
  [Horace](http://horace.isis.rl.ac.uk/) for resolution calculation.
* `eniius` also has the facility to export this information as
  [McStas](http://mcstas.org/) instrument files for ray tracing simulations of
  the resolution function, and can also import certain McStas instrument files
  to embed in raw data files.
* Is compatible with the [IBEX](https://github.com/ISISComputingGroup/IBEX/)
  data acquisition program as well as the [Mantid](https://mantidproject.org/)
  data reduction suite.
* Has a GUI for visualising the instrument data.

