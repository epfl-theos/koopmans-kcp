fname=$1

# sed -e "1,/TERMINATE/w $before
# /TERMINATE/,\$w $after" file

echo 'maximum_charge_density'
sed -e '/Contents of orbital density/,$d' $fname | grep 'E' | grep -v 'CHARGE-DENSITY' | sort -g | tail -1

echo 'maximum_orbital_density'
sed -n -e '/Contents of orbital density/,$p' $fname | grep 'E' | grep -v 'EFFECTIVE' | sort -g | tail -1

# tot_energy=`grep 'total energy' $fname | tail -1 | awk '{printf "%12.6f\n", $4}'`
# odd_energy=`grep 'odd energy' $fname | tail -1 | awk '{printf "%12.6f\n", $4}'`
# homo_energy=`grep -A2 'HOMO Eigenvalue' $fname | tail -1 | awk '{printf "%12.6f\n", $1}'`
# lumo_energy=`grep -A2 'LUMO Eigenvalue' $fname | tail -1 | awk '{printf "%12.6f\n", $1}'`
# 
# if test "$tot_energy" != ""; then
# 	echo tot_energy
# 	echo $tot_energy
# fi
# 
# if test "$odd_energy" != ""; then
# 	echo odd_energy
# 	echo $odd_energy
# fi
# 
# if test "$homo_energy" != ""; then
# 	echo homo_energy
# 	echo $homo_energy
# fi
# 
# if test "$lumo_energy" != ""; then
# 	echo lumo_energy
# 	echo $lumo_energy
# fi
