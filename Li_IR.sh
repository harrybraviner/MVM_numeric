#!/bin/bash

# Should run the MVM program for various values of z

z_list="1.1
1.3
1.6
2.0
2.5
3.0
3.5
4.0
5.0"

for z in $z_list
do
	./a.out $z
	mv ./output.dat Li_IR_z_$z.dat
done
