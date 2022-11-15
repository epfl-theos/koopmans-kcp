#!/bin/bash

echo "Running BIN2XML ..."
echo "${ESPRESSO_ROOT}/bin/bin2xml_real_space_density.x $1 $2 $3"
ls
echo "Contents of total density"
cat $2/'charge-density.xml'
echo "Contents of orbital density"
cat $2/'orbital.occ.0.00001.xml'
