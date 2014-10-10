post_lammps
===========

A tool for performing post-processing with lammps log files


Installation
------------

make<br>
make install<br>


Usage
-----

post_lammps action [args] [file name] [property 1] [property 2] ...<br>
  action = acf or block or eacf or extract or histogram or list or stats<br>
    acf args = window<br>
    block args = none<br>
    eacf args = window<br>
    extract args = skip<br>
    histogram args = nbins<br>
    list args = none<br>
    stats args = none<br>

