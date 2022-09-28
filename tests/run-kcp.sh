#!/bin/bash
#
# Copyright (C) 2001 Quantum ESPRESSO
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

if [ $QE_USE_MPI == 1 ]; then
  export PARA_PREFIX="mpirun -np ${TESTCODE_NPROCS}"
  export PARA_SUFFIX=" "
else
  unset PARA_PREFIX
  unset PARA_SUFFIX
fi

echo "Running KCP ..."
echo "${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcp.x ${PARA_SUFFIX} -in $1 > $2 2> $3"
${PARA_PREFIX} ${ESPRESSO_ROOT}/bin/kcp.x ${PARA_SUFFIX} -in $1 > $2 2> $3
if [[ -e CRASH ]]
then
  cat $3
fi
rm -f input_tmp.in
