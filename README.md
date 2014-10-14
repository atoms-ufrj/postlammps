Postlammps
==========

A tool for performing post-processing with lammps log files.


Installation
------------

git clone https://github.com/atoms-ufrj/postlammps<br>
cd postlammps<br>
make<br>
sudo make install<br>


Usage
-----

Usage: postlammps [options] action [args] property-1 [property-2 ...]<br>
  * action = acfun or block or fluct or histo or ineff or print or props or sampl or stats<br>
    * acfun args = maxtime<br>
    * block args = none<br>
    * fluct args = maxtime<br>
    * histo args = nbins<br>
    * ineff args = none<br>
    * print args = none<br>
    * props args = none<br>
    * sampl args = none<br>
    * stats args = none<br>

Actions
-------

**acfun**: Computes autocorrelation functions from zero to maxtime<br>
**block**: Performs normalization group blocking analysis<br>
**fluct**: Computes normalized fluctuation autocorrelation functions from zero to maxtime<br>
**histo**: Builds histograms with specified number of bins<br>
**ineff**: Computes statistical inefficiencies and uncertainties<br>
**print**: Prints the values of the selected properties<br>
**props**: Lists all properties available in the log file<br>
**sampl**: Samples uncorrelated points from the original data<br>
**stats**: Computes basic statistics<br>

Options
-------

**-in**: Specifies the name of the log file to be processed<br>
**-p**: Tells postlammps to read a plain data file instead of a lammps log file<br>
**-e n**: Skips n lines after every property reading<br>
**-d delimiter**: Specifies the item delimiter used for output<br>
  delimiter = space or comma or semicolon or tab<br>
**-nt**: Does not print property titles<br>
**-c X**: Consider only the last X% of data<br>

Notes
-----

* If option **-in** is not used, postlammps operates in the standard input.<br>
* A plain data file contains the property titles in the first line and the property values in the subsequent lines.<br>

