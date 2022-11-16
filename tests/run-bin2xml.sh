#!/bin/bash
set -e

echo "Running BIN2XML ..."
echo "${ESPRESSO_ROOT}/bin/bin2xml.x $1 tmp.xml > $2 2> $3"
${ESPRESSO_ROOT}/bin/bin2xml.x $1 tmp.xml > $2 2> $3
if [[ -e tmp.xml ]]
then
    cat tmp.xml >> $2
    rm tmp.xml
fi
