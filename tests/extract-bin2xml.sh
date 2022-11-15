fname=$1

echo 'maximum_charge_density'
sed -e '/Contents of orbital density/,$d' $fname | grep 'E' | grep -v 'CHARGE-DENSITY' | sort -g | tail -1

echo 'maximum_orbital_density'
sed -n -e '/Contents of orbital density/,$p' $fname | grep 'E' | grep -v 'EFFECTIVE' | sort -g | tail -1
