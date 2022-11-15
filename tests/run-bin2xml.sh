#!/bin/bash

echo "Running BIN2XML ..."
echo "${ESPRESSO_ROOT}/bin/bin2xml_real_space_density.x $1 $2 $3"
echo "Contents of total density"
cat 'charge-density.xml'
echo "Contents of orbital density"
cat 'orbital.occ.0.00001.xml'
