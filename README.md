# README FILE FOR GRAVITY MODELING PROGRAM SUITE

## Contact information
Please contact [Dave Blair](mailto:dblair@purdue.edu) with questions.

This work (including all of the programs in this directory) is (c) David Blair,
2015, and is licensed under a [Creative Commons Attribution-NonCommercial-ShareAlike
Unported License](http://creativecommons.org/licenses/by-nc-sa/3.0/deed.en_US). 

Last updated 2014-05-01

---

## Introduction

Welcome! This is a suite of programs used to calculate the gravity anomaly over
axisymmetric computer models of geologic features. This process involves roughly
5 steps:

1. Creating a model of interest
2. Retreiving information about the model in an easily parseable format
3. Calculating the gravitational acceleration over that retrieved information
4. Subtracting the accleration of some sort of baseline or "geoid" model
5. Recording and displaying the output

I have written several programs which, between them, cover every part of this
workflow except the first. They are

- `grav_from_Abaqus.py` and `grav_from_iSALE.py` _(step 2)_
- `grav_anomaly.py` _(steps 3, 4, and 5)_

The usage of these programs is detailed below. All of these are written in
[Python](www.python.org), which should already be installed on Mac OS X or Linux
machines, and is easy to install (and broadly useful) on Windows. The final
program, `grav_anomaly.py`, depends on [NumPy](www.numpy.org) and
[SciPy](www.scipy.org) for some of its numeric and plotting routines. Finding
work-arounds for these two to reduce the number of dependencies is an end goal
of mine, but I haven't found time for that yet.

A quick word of warning, too: **none of my programs ask before overwriting
files**. Yes, I know it's bad behavior, but it's hard to implement in a
cross-platform way. You have been warned.

Finally, please note that this code was originally written for the specific
research needs of one person---myself---and thus may have unexpected bugs or may
be unable to cover cases sufficiently different from those I developed it for.
If you find that to be the case, please don't hesitate to [contact
me](mailto:dblair@purdue.edu) with questions, bug reports, and feature requests.
This README is also a work in progress, and I'll happily take input as to how I
can make it more useful.

---

## Obtaining Parseable Information from the Model

Currently, it possible to get the required data either from _Abaqus_ finite
element models or _iSALE_ hydrocode models. There are dinstinct programs for
handling each source:

- `grav_from_Abaqus.py` for _Abaqus_ models, and
- `grav_from_iSALE.py` for _iSALE_ models

Here is how to use them:

- For _Abaqus_ models:
    - Open the model's output database file in CAE, go to "File > Run
      Script...", and run `grav_from_Abaqus.py`.
    - This will generate a file called something like "foo.grav", where "foo" is
      the name of the output database minus the ".odb" extension.
    - Note that there are also other functions possile in `grav_from_Abaqus.py`,
      such as creating files meant for examining the model's density structure
      directly instead of (or in addition to) the normal ".grav" output. These
      are discussed in the **OPTIONS** section near the top of the file
      `grav_from_Abaqus.py`. You can also choose here to grab data from the
      first frame of the model as well as the last, or to grab data from a
      specific 'frame number' of an Abaqus output file.
- For _iSALE_ models:
    - Use the program `grav_from_iSALE.py`. For info on some of the options
      available here, try

            grav_from_iSALE.py --help

    - You need access to both the iSALE input file and a file containing XY data
      from the cells. The XY data will be reported for the _center_ of each
      cell, so `grav_from_iSALE.py` uses the input file to re-generate the iSALE
      mesh. It can then examine the corner positions of the cell (the boundaries
      of where that XY data applies), which will become useful in the gravity
      calculation later on.
    - To just run `grav_from_iSALE.py` on an iSALE input file called "foo.inp"
      and a data file called "bar.txt" (where the data file is X, Y, Density),
      run

            grav_from_iSALE.py foo.inp bar.txt

    - The program automatically replaces the ".txt" at the end of the data
      file's filename with ".grav", resulting in an output file called (in this
      example) "bar.grav". If you use something other than ".txt" at the end
      of your data filenames, simply change the variable _datafile_suffix_ in
      the **OPTIONS** section at the top of `grav_from_iSALE.py` from ".txt" to
      something else in order to get this same automatic-renaming behavior.

Either of these methods will create one or more plaintext ".grav" files which we
can use in the next step.


## Grabbing Parseable Information about the Geoid

To calculate a free-air or Bouguer anomaly, we need a geoid against which to
compare our model data. The procedure for doing this is much the same as for
grabbing non-geoid model data.

- For _Abaqus_ models:
    - Create a new model with the exact same exterior/overall dimensions as your
      model of interest but with a flat-layered crust equivalent in thickness to
      the far-field crustal thickness in your model
    - Make sure this model has a different name than your model of interest
      (e.g., "foo_geoid").
    - Run this model with no loads present (so that there's no absolutely no
      deformation) and repeat the procedure above
    - The result of all of this should be a set of ".grav" files for your geoid
      model, exactly as you generated for your model of interest (e.g.
      "foo_geoid.grav")
- For _iSALE_ models:
    - Grab output from time 0 in your model; this should already be a model with
      a nice flat crust-mantle boundary, although it also has the impactor
      sitting in the middle of it.
    - Run `grav_from_iSALE.py` with the option turned on to discard all cells
      above the surface (technically, cells with y > 0) so that you ignore the
      impactor; this is done using the `--discardabovezero` flag (`-0` for
      short). If your input and start-of-model data files are called "foo.inp"
      and "bar-from-time-zero.grav", then the command would look something like

            grav_from_iSALE.py --discardabovezero foo.inp bar-from-time-zero.txt

    - After doing this, you should have a new ".grav" file based on data from
      your first timestep but with the impactor removed. This is your geoid
      file. In the above example, it'd be called "bar-from-time-zero.grav".


## Doing the Gravity Calculation

Now you can move on to the actual gravity calculation. This uses the program
`grav_anomaly.py`. Let's start with a quick introduction to that program. To get
an overview, run

    grav_anomaly.py --help

This points out all of the variables in the program that can be changed at
runtime. There are quite a few, and it would be a pain to have to set something
like `--geoidname foo-very-long-annoying-name-geoid.grav` before every run, so
there's an alternative: near the top of the source file for `grav_anomaly.py`
(which, since it's in Python, is the same as the 'executable') is a section
called **OPTIONS**. In here you can set all of the same things that you can at
run time, but they'll be persistent between runs. You can even set default
behavior, like always generating text instead of visual output, or turning on
_verbose_mode_ so that you see extra visual output about the calculation as it's
happening. In the examples below, however, I'm going to assume that you have
everything in **OPTIONS** set to _False_, and generally err on the side of
over-specifying command-line options (including using long option names instead
of short, like `--verbose` insted of `-v`), for clarity. If you set some
defaults in **OPTIONS**, though, you can probably save yourself some typing. You
can check the current **OPTIONS** values by running

    grav_anomaly.py --defaults

and see some examples of usage with

    grav_anomaly.py --examples

On to actually using the program. `grav_anomaly.py` accepts input in two forms:
gravity data files (".grav"), as generated above, and XY files consisting of an
X value in kilometers in one column and a gravitational acceleration value in a
second column. These XY files are very generic, and suitable for use in external
plotting programs; in `grav_anomaly.py`, I refer to output of these files as
"GMT mode", but they could be used by anything that can make scatter plots, be
it Matlab, Mathematica, Excel, crayons and paper, whatever. `grav_anomaly.py`
is also capable of _creating_ these files, too, which makes them very useful for
recording the results of a gravity calculation that you don't want to have to
repeat. One obvious use case is for the geoid model, which will likely be shared
by a number of models that you'll want to look at. To process a geoid file (here
called "foo_geoid.grav"), use `grav_anomaly.py` in "raw acceleration" mode
(the `--acceleration` and `-a` runtime flags) so that it doesn't compare the
result to anything, like so:

    grav_anomaly.py --verbose --GMT --acceleration foo_geoid.grav

This generates a file ending in "_acc.xy" (e.g. "foo_geoid_acc.xy") chosen to
indicate that it's total raw acceleration, not an anomaly, that's stored in
the file; this is the geoid acceleration file against which you'll be comparing
your real model. The above command is also set up to generate text/XY-file
output (`--GMT`) and to provide feedback on the command line while running
(`--verbose`). You can also change the filename extensions and that "_acc"
suffix by changing "GMT_suffixes" in **OPTIONS** if you'd prefer to use
something else as the default.

It is important to note that this file, along with every other XY file generated
by `grav_anomaly.py`, is inextricably tied to the parameters with which the
calculation was done, namely spacecraft elevation, lateral range, and lateral
resolution (`-e` or `--elevation`, `-r` or `--range`, and `-n` or `--numsteps`
on the command line). If you're differing from some usual setup, you may want to
rename the files to remind you of strange parameters (e.g.,
"foo_out-to-2e6-km_acc.xy") or insert a line or two of notes at the top of the
file; the latter can be done automatically with the `--GMTheaders` command-line
flag or by setting "GMT_print_headers" to _True_ in **OPTIONS**.

Running this can take some time, depending on two things: the size of your model
(number of elements/cells), and any **subdivision** that you're doing during the
calculation. Subdivision is a technique included in this program to generate
accurate output. Essentially, inaccuracy crops up in low-resolution models
because the gravity calculation code considers each element to be a ring with
real mass and radius but with an infinitesimal cross-section. If these
"point-rings" are far enough apart, you wind up with errors across the model
that are especially prominent near the axis of symmetry, because there's
effectively a hole at at the axis (even the centermost element/cell is still
actually an annulus). What subdivision does is to take each element and break it
into two smaller elements with half of the mass of the original that are evenly
spaced within the original element/cell's boundaries. The resulting model looks
exactly the same in terms of density distribution, but the acceleration
calculated for it is much closer to analytic solutions. Generally, for the kind
of resolution I have in _Abaqus_ models I've run, one subdivision is enough, but
this is something that warrants a sensitivity check to make sure that you're not
getting bogus results from having element centroids spaced too far apart. To use
`grav_anomaly.py` with, e.g., 4 subdivisions (multiplying your element count by
16, and causing a 16-fold increase in calculation time), run

    grav_anomaly.py --verbose --GMT --subdivide 4 --acceleration foo_geoid.grav

Once you find a number of subdivisions past which your result stops changing, I
recommend setting that as the default in **OPTIONS** so you don't have to
remember to do this in future runs made with similar mesh sizes.

With that done, we can now move on to computing the free-air or Bouguer
anomalies over a model. Sensibly enough, this is done with the "freeair" and
"bouguer" modes of `grav_anomaly.py`. If our geoid file (from the procedure
above) is called "foo_geoid_acc.xy", the command to calculate the free-air
anomaly will look something like this:

    grav_anomaly.py ---verbose --freeair --geoidname foo_geoid_acc.xy bar.grav

Note that you could have also used the original ".grav" file for the geoid, like
so,

    grav_anomaly.py ---verbose --freeair --geoidname foo_geoid.grav bar.grav

but that it would've run much more slowly, as it'd have to re-process the geoid
file. Either way, you'll end up with two files: one for the acceleration over
the model (ending in "_fa_acc.xy" to distinguish it from the Bouguer one, which
we'll get to in a minute) and one for the actual free-air anomaly ("_fa.xy"). If
you decided to run your data file through `grav_anomaly.py` in raw acceleration
mode before doing this step, you can also generate a free-air anomaly even if
both the geoid and data files are in XY format; this runs in practically no time
at all, although the preceding step of turning your data file from ".grav" into
".xy" takes just as long. Basically, I recommend setting up your workflow such
that you only have to crunch each file from ".grav" to ".xy" once, as it saves a
lot of calculation time.

Calculating a Bouguer gravity anomaly is outwardly very similar to doing a
free-air calculation:

    grav_anomaly.py --verbose --bouguer --geoidname foo_geoid_acc.xy bar.grav

The ".grav" file used here is the same one as in the free-air calculation, so
there's no difference there. The files generated are also similar, with one
being the _Bouguer-corrected_ acceleration over the model ("_boug_acc.xy") and
the other being the anomaly XY file ("_boug.xy"). However, you can't use an XY
data file to generate a bouguer anomaly file, because the Bouguer
calculation relies on replacing the densities of certain elements/cells with
other assumed densities; you therefore need access to element/cell data, which
is only stored in the ".grav" file. The acceleration file coming out of this
calculation is recorded _after_ these Bouguer corrections are done, which is why
it is given a different name than the "free-air acceleration" (which I suppose
would more correctly be called the "acceleration suitable for free-air
calculation"). The assumptions used in the Bouguer correction are stored in the
"bouguer_density_cutoff" (density below which something counts as "space") and
"bouguer_assumed_crustal_density" (the assumed density of the crust) variables
in the **OPTIONS** section of `grav_anomaly.py`.


### Curved models

Doing calculations with axisymmetric models of a body with a curved surface is
largely the same as the process described above. There are only a few new
variables, which can as usual be set either in **OPTIONS** or via the command
line. First, there is a switch to turn on curved-surface mode, enabled using
the `-c` or `--curved` command line flags. Next, there is a setting for the
"theta" range (`-t` or `--thetarange`), which is the angular extent of your
calculation, measured in degrees from the axis of symmetry (the "pole" of the
model). This is equivalent to the "minimum latitude" of the calculation if the
axis of symmetry is considered to be at 90º, and replaces the lateral range
parameter for curved models. The final consideration is the radius of the planet
(`--planetradius`), which is used in determining the spacecraft path, performing
Bouguer corrections, and a few other necessary calculations. All other settings
and flags (including the spacecraft elevation) should work identically to their
flat-world counterparts.


## Plotting the Data

Now that you have this data, you may want to plot it. The XY files generated
above are plaintext, two-column files, and at their _most_ complicated they
contain a few lines of hatchmark-commented lines at the top describing the
parameters that went into the file's creation. They should be usable by your
favorite plotting software pretty easily.

If you'd like, though, there's also a plotter built into `grav_anomaly.py`. This
plotter _can_ work directly off of ".grav" data files, but unless they're very
small you won't want to do this, as it requires re-calculating the gravity
signature of the file at runtime. Instead, you can just plot up the XY files
generated previously, remembering to use `-a` or `--acceleration` (as otherwise the
program will re-compute the difference between this _anomaly_ file and the
geoid, generating bogus results):

    grav_anomaly.py --acceleration --plot foo_fa.xy

This generates a reasonably pretty Matlab-like plot using
[MatPlotLib](matplotlib.org), which can be saved to a variety of raster and
vector formats, lets you zoom in on data, allows scrollable and resizeable
windows, etc. Using `grav_anomaly.py` for plotting also has the advantage that
you can do things like

    grav_anomaly.py --plot --geoidname foo_geoid_acc.xy bar_fa_acc.xy baz_fa_acc.xy

and have it plot up the free-air anomaly of both the `bar_fa_acc.xy` and
`baz_fa_acc.xy` models against the same geoid file. It can take up to 10 files
as input and assign different colors to each (the limit is actually only the
color list it picks from, so let me know if you need to compare more than 10
files!), and if you're in free-air mode it'll calculate the difference between
each file and the geoid. Be careful about which files you're using; any time you
calculate a geoid, you want to use _raw acceleration_ files as input, or it will
have subtracted the same thing twice. Also, for plotting multiple Bouguer files,
you need to use the "foo_boug.xy" files and plot with "--acceleration", as
trying to use "--bouguer" and an XY data file simultaneously raises an error
intended to prevent unintentionally non-Bouguer files being used in a Bouguer
gravity calculation.

---

## Epilogue

Hopefully that explanation makes it pretty clear what the various pieces of code
in this package are capable of and how to use them. I welcome input as to how to
improve the code, improve this README, or otherwise make it easier for you to
benefit from the tools I've tried to create.
