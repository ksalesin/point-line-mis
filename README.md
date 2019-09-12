Combining Point and Line Samples for Direct Illumination
========================================================

This repository contains an implementation of the 2019 EGSR paper ["Combining Point and Line Samples for Direct Illumination"](https://cs.dartmouth.edu/~wjarosz/publications/salesin19combining.html) by Katherine Salesin and Wojciech Jarosz. If you use this code, please cite the paper :)

It is built on an implementation of the 2016 EGSR paper ["Line Sampling for Direct Illumination"](http://graphics.cs.kuleuven.be/publications/BD2016LSFDI/) by Niels Billen and Philip Dutré, which in turn is built on [PBRTv3](http://pbrt.org).

See PBRT's README (pasted below) for installation and build instructions.

They key classes in this implementation are `PointLineDirectIntegrator` and `Quad`.

Usage of new classes in scene files
-----------------------------------

The `scenes/` folder includes several examples of these usages.

### Quad ###

Define a quad shape by a point origin `o` and two vectors `u` and `v` representing two edges spawning from the origin (vertices of the quad will be `o`, `o+u`, `o+v`, and `o+u+v`).

### PointLineDirectIntegrator ###

Use the PointLineDirectIntegrator with the line `Integrator "pointline" [options...]` in the scene file. Options include:

Type | Name  | Acceptable Values | Description
---- | ----- | ----------------- | ----------- 
string | **method** | "render", "drawlines" | Whether to render the image normally ("render") or draw red lines on the light source to visualize the direction the lines are sampling ("drawlines").
string | **strategy1**, **strategy2** (optional), **strategy3** (optional) | "pt&#8209;surf&#8209;area", "pt&#8209;solid&#8209;angle", "brdf", "line&#8209;surf&#8209;area", "line&#8209;solid&#8209;angle" |  Which strategy the integrator will use: point or line sampling the surface area of the light source ("pt-surf-area", "line-surf-area"), point or line sampling the solid angle of the light source ("pt-solid-angle", line-solid-angle"), or point sampling the BRDF. If more than one strategy is included, the PointLineDirectIntegrator will **multiple importance sample (MIS)** the strategies.
bool | **stratified** | true, false | Whether to use stratified (true) or random (false) samples on the light source. **Important**: the `Sampler` line in the scene must also be set to "stratified." Set the number of samples to take on the light source per scene intersection on the `AreaLightSource` line with **xsamples** and **ysamples**.

The following attributes are defined for line strategies only, where the **#** following the name should match the strategy (ex. **angle1** for **strategy1**).

Type | Name  | Acceptable Values | Description
---- | ----- | ----------------- | -----------  
bool | **importance#** | true, false | Whether to importance sample lines by length (true by default since this was found to be a universally better strategy by Billen and Dutré).
float | **angle#** | [0,1] | Angle to line sample on the light source. For "line-surf-area," [0,1] maps to lines offset by [0&deg;, 180&deg;] relative to one of the quad edges. For "line-solid-angle," only 0 and 1 are acceptable values, which correspond to line sampling parallel to either edge `u` or `v` of the quad. Use the "drawlines" method to visualize the line sample direction on the light source.


pbrt, Version 3
===============

[![Build Status](https://travis-ci.org/mmp/pbrt-v3.svg?branch=master)](https://travis-ci.org/mmp/pbrt-v3)
[![Build status](https://ci.appveyor.com/api/projects/status/mlm9g91ejxlcn67s/branch/master?svg=true)](https://ci.appveyor.com/project/mmp/pbrt-v3/branch/master)

This repository holds the source code to the new version of pbrt that is
described in the third edition of *Physically Based Rendering: From
Theory to Implementation*, by [Matt Pharr](http://pharr.org/matt), [Wenzel
Jakob](http://www.mitsuba-renderer.org/~wenzel/), and Greg Humphreys.  As
before, the code is available under the BSD license.

Please see the [User's Guide](http://pbrt.org/users-guide.html) for more
information about how to check out and build the system as well as various
additional information about working with pbrt.

Example scenes
--------------

Over 10GB of example scenes are available for download. (Many are new and
weren't available with previous versions of pbrt.) We're trying an
experiment and making them available via git. Run:
```
$ git clone git://git.pbrt.org/pbrt-v3-scenes
```
to get them. We will update this repository as more scenes become
available. (See the `README.md.html file` in the scene distribution for
more information about the scenes and preview images.)

The [pbrt website](http://pbrt.org) has general information about
both *Physically Based Rendering* as well as pbrt-v2, the previous version
of the system.

Building pbrt
-------------

To check out pbrt together with all dependencies, be sure to use the
`--recursive` flag when cloning the repository, i.e.
```bash
$ git clone --recursive https://github.com/mmp/pbrt-v3/
```
If you accidentally already cloned pbrt without this flag (or to update an
pbrt source tree after a new submodule has been added, run the following
command to also fetch the dependencies:
```bash
$ git submodule update --init --recursive
```
pbrt uses [cmake](http://www.cmake.org/) for its build system.  On Linux
and OS X, cmake is available via most package management systems.  For
Windows, or to build it from source, see the [cmake downloads
page](http://www.cmake.org/download/).

* For command-line builds on Linux and OS X, once you have cmake installed,
  create a new directory for the build, change to that directory, and run
  `cmake [path to pbrt-v3]`. A Makefile will be created in that
  current directory.  Run `make -j4`, and pbrt, the obj2pbrt and imgtool
  utilities, and an executable that runs pbrt's unit tests will be built.
* To make an XCode project file on OS X, run `cmake -G Xcode [path to pbrt-v3]`.
* Finally, on Windows, the cmake GUI will create MSVC solution files that
  you can load in MSVC.

If you plan to edit the lexer and parser for pbrt's input files
(`src/core/pbrtlex.ll` and `src/core/pbrtparase.y`), you'll also want to
have [bison](https://www.gnu.org/software/bison/) and
[flex](http://flex.sourceforge.net/) installed. On OS X, note that the
version of flex that ships with the developer tools is extremely old and is
unable to process `pbrtlex.ll`; you'll need to install a more recent
version of flex in this case.

### Debug and Release Builds ###

By default, the build files that are created that will compile an optimized
release build of pbrt. These builds give the highest performance when
rendering, but many runtime checks are disabled in these builds and
optimized builds are generally difficult to trace in a debugger.

To build a debug version of pbrt, set the `CMAKE_BUILD_TYPE` flag to
`Debug` when you run cmake to create build files to make a debug build. For
example, when running cmake from the command lne, provide it with the
argument `-DCMAKE_BUILD_TYPE=Debug`. Then build pbrt using the resulting
build files. (You may want to keep two build directories, one for release
builds and one for debug builds, so that you don't need to switch back and
forth.)

Debug versions of the system run much more slowly than release
builds. Therefore, in order to avoid surprisingly slow renders when
debugging support isn't desired, debug versions of pbrt print a banner
message indicating that they were built for debugging at startup time.

### Build Configurations ###

There are two configuration settings that must be set at compile time. The
first controls whether pbrt uses 32-bit or 64-bit values for floating-point
computation, and the second controls whether tristimulus RGB values or
sampled spectral values are used for rendering.  (Both of these aren't
amenable to being chosen at runtime, but must be determined at compile time
for efficiency).

To change them from their defaults (respectively, 32-bit
and RGB.), edit the file `src/core/pbrt.h`.

To select 64-bit floating point values, remove the comment symbol before
the line:
```
//#define PBRT_FLOAT_AS_DOUBLE
```
and recompile the system.

To select full-spectral rendering, comment out the first of these two
typedefs and remove the comment from the second one:
```
typedef RGBSpectrum Spectrum;
// typedef SampledSpectrum Spectrum;
```
Again, don't forget to recompile after making this change.

