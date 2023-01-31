#!/bin/bash
source ../config.sh
# requirements: 
# cas-offinder installation
# input.txt and input2.txt with the following specifications: https://github.com/snugel/cas-offinder


cas-offinder $basedir/analysis/input2.txt C $basedir/offinder/off_target_nrch.txt
cas-offinder $basedir/analysis/input.txt C $basedir/offinder/off_target_spry.txt
