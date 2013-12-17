#!/bin/sh
# This file is only for use when compiling this software inside the HPCLinux distrubtion
# available at http://www.hpclinux.org.

if [ -r /opt/intel/composer_xe_2013_sp1.0.080/bin/compilervars.sh ]; then 
  source /opt/intel/composer_xe_2013_sp1.0.080/bin/compilervars.sh intel64
fi

if [ -r /opt/intel/impi/4.1.1.036/intel64//bin/mpivars.sh ] ; then
  source /opt/intel/impi/4.1.1.036/intel64//bin/mpivars.sh
fi 
