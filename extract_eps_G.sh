#!/bin/sh
for i in $*
do
awk '
/\/varray/ { on = 0 }
on == 1    { print $2, $3, $4, $5, $6+0 }
/varray name="epsilon_diag"/ { on=1 ; done=1 ; print " " }' <$i | sort -n
echo
done
