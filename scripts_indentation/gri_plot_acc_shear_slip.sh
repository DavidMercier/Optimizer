#!/bin/bash
out="$1"
gri -output ${out}.ps -p /egr/research/CMM/Code/Gri/plot line solid color own nolabel clip frame 5 linear 'increment' %g 0 950 5 2 1 log 'log(acc_shear)' %g 0.006 0.03 5 2 1 \
SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 9 'rgb 0 0 1' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 10 'rgb 0 0 1' \
SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 11 'rgb 0 0 0.5' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 12 'rgb 0 0 0.5' \
SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 13 'rgb 0 0.5 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 14 'rgb 0 0.5 0' \
SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 15 'rgb 0.5 0 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 16 'rgb 0.5 0 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 17 'rgb 0.5 0 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 18 'rgb 0.5 0 0' \
SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 19 'rgb 0 0.5 0.5' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 20 'rgb 0 0.5 0.5' \
SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 21 'rgb 0.5 0.5 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 22 'rgb 0.5 0.5 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 23 'rgb 0.5 0.5 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 24 'rgb 0.5 0.5 0' \
SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 25 'rgb 0.5 0 0.5' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 26 'rgb 0.5 0 0.5' \
SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt 1 27 'rgb 0.5 0.5 0.5' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 28 'rgb 0.5 0.5 0.5' \
SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 29 'rgb 1 0 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 30 'rgb 1 0 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 31 'rgb 1 0 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 32 'rgb 1 0 0' \
SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 33 'rgb 0 1 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 34 'rgb 0 1 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 35 'rgb 0 1 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 36 'rgb 0 1 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 37 'rgb 0 1 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 38 'rgb 0 1 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 39 'rgb 0 1 0' SX1_2015-06-06_22h42m55s__fric0.3_R1.00_cA90.0_h0.300_postdef.txt  1 40 'rgb 0 1 0'
gri_ps2pdf ${out}.ps
