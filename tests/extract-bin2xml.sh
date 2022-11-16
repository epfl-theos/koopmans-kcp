fname=$1

echo 'maximum_charge_density'
grep 'E' $fname | grep -v 'CHARGE-DENSITY' | sort -g | tail -1
