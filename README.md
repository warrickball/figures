Figures
=======

Code used to generate a series of plots I sometimes use for
teaching/introductory material on stellar modelling and
helio/asteroseismology.  Each figure is designed to be self-contained
with minimal dependencies.  Most should give useful information if you
run them with the `-h` flag. e.g. `python3 modelS_rays.py -h`.

In addition, each figure tries to download
the relevant datasets, some of which are large (up to about 150MB),
and some save intermediate data products to speed up frequent use.
Some datasets require are more complicated to produce: input files for
the relevant programs are included in `gen/`.  All in all, a current
run of all plots can contain up to about 1GB of data.
