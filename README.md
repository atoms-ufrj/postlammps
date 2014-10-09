post_lammps
===========

A tool for performing post-processing of lammps log files


Installation
============

make
make install


Usage
=====

post_lammps action <args> <file name> <property 1> <property 2> ...
  action = block or extract or histogram
    block args = none
    extract args = none
    histogram args = number-of-bins

