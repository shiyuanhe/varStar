## Description

This software ``gpmodel.cpp`` (come with ``Makefile``) is a simplified (and faster) standalone version of the R package for deriving frequency spectra of sparsely sample Mira light curves.

## Installation

The software has been tested on Mac and Linux operating systems. BLAS (or openBLAS), LAPACK, and armadillo packages are needed to compile it. c++11 standard is used.

(1) If armadillo is not installed in your computer, please go to [http://arma.sourceforge.net/download.html](http://arma.sourceforge.net/download.html) and download the package. Once extract the package (`tar -xcvf armadillo-x.xxx.x.tar.xz`), find the `include` directory and save it in your preferred location.

(2) Open the Makefile with your preferred text editor, switch the part `~/Programs/armadillo-6.700.4/include` with the location of your armadillo include directory (see step 1) and save the file.

(3) Type `make` in terminal to compile. An executable binary named `gpmodel` will appear after a few seconds.

## Uninstall

This software does not go to system directories. Simply delete them to uninstall.

## Usage

Type `gpmodel` (or `./gpmodel` on some machine) in terminal to start this software. You can also start it in a different working directory by calling it with the path, i.e. `/path/to/it/gpmodel`.

There are two ways to call for calculation: 

(a) To derive the frequency spectrum for a single light curve, call `./gpmodel -f filename.dat`, where `filename.dat` is the light curve file with 3 columns: time, magnitude, uncertainty. A subdirectory named `gp_spectra` will be created automatically to store the spectrum.

(b) To obtain frequency spetra for a list of light curves, call `./gpmodel -l list_file_names.lst`, where `list_file_names.lst` is text file with a list of light curve file names separated by line break, i.e.

```
001.dat
~/Work/lcs/mylc.dat
mira.dat
/Users/mira/lcs/3.dat
004.lc
```
Please be sure the light curves exist and there are enough (say, >~ 10) measurements to avoid errors.

You can always call it in a different language (IDL, R, Python, ...) within loops. Let the ``cmd`` = ``/path/to/gpmodel -f filename.dat``

```
(IDL) spawn, cmd
(R) system(cmd)
(Python) os.system(cmd)
```


## Advanced usage
Please feel free to edit the source code. 

The main function is currently at the very bottom of gpmodel.cpp file, where you can change delimiter of the light curve file, output directory, trial frequencies.