Postlammps
==========

A tool for performing post-processing with lammps log files.


Installation
------------

make<br>
sudo make install<br>


Usage
-----

Usage: postlammps [options] action [args] file-name property-1 [property-2 ...]<br>
  action = acf or acfn or block or histo or ineff or print or props or stats<br>
    acf   args = maxtime<br>
    acfn  args = maxtime<br>
    block args = none<br>
    histo args = nbins<br>
    ineff args = none<br>
    print args = none<br>
    props args = none<br>
    stats args = none<br>

Actions
-------

**acf**: Computes autocorrelation functions from zero to maxtime<br>
**acfn**: Computes normalized fluctuation autocorrelation functions from zero to maxtime<br>
**block**: Performs normalization group blocking analysis<br>
**histo**: Builds histograms with specified number of bins<br>
**ineff**: Computes statistical inefficiencies and uncertainties<br>
**print**: Prints the values of the selected properties<br>
**props**: Lists all properties available in the log file<br>
**stats**: Computes basic statistics<br>

Options
-------

**-every n**: Skips n lines after every property reading<br>
**-delim d**: Specifies the item delimiter used for output<br>
  d = space or comma or semicolon<br>

