#Author: Karlton Sequeira

#declare variables and paths
export PROJ_HOME=/cs/sequek/laptop/projects
export ECLAT_HOME=$PROJ_HOME/eclat
export VC_HOME=$PROJ_HOME/VC

export VC_DATA=$VC_HOME/data
export VC_SRC=$VC_HOME/src
export VC_OBJ=$VC_HOME/obj
export ECLAT_DATA=$ECLAT_HOME/data
export ECLAT_SRC=$ECLAT_HOME/src
export ECLAT_OBJ=$ECLAT_HOME/obj

g++ -O4 -o $VC_OBJ/dataGen $VC_SRC/dataGen.cc
g++ -O4 -o $VC_OBJ/convertData $VC_SRC/convertData.cc

#generate data for testing effect of number of dimensions
$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1.ha -d 50
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1 -d 50
rm $ECLAT_DATA/s1e3k5d5e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d1e2.ha -d 100
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d1e2 -d 100
rm $ECLAT_DATA/s1e3k5d1e2.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d2e2.ha -d 200
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d2e2 -d 200
rm $ECLAT_DATA/s1e3k5d2e2.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d4e2.ha -d 400
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d4e2 -d 400
rm $ECLAT_DATA/s1e3k5d4e2.ha

#generate data for testing effect of number of points
$VC_OBJ/dataGen -o $ECLAT_DATA/s1e4k5d5e1.ha -n 10000
$VC_OBJ/convertData -i $ECLAT_DATA/s1e4k5d5e1 -n 10000
rm $ECLAT_DATA/s1e4k5d5e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s5e4k5d5e1.ha -n 50000
$VC_OBJ/convertData -i $ECLAT_DATA/s5e4k5d5e1 -n 50000
rm $ECLAT_DATA/s5e4k5d5e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e5k5d5e1.ha -n 100000
$VC_OBJ/convertData -i $ECLAT_DATA/s1e5k5d5e1 -n 100000
rm $ECLAT_DATA/s1e5k5d5e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s2e5k5d5e1.ha -n 200000
$VC_OBJ/convertData -i $ECLAT_DATA/s2e5k5d5e1 -n 200000
rm $ECLAT_DATA/s2e5k5d5e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s4e5k5d5e1.ha -n 400000
$VC_OBJ/convertData -i $ECLAT_DATA/s4e5k5d5e1 -n 400000
rm $ECLAT_DATA/s4e5k5d5e1.ha

#generate datasets with different cluster strengths ($\kappa$)
$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1a2.ha -a 2.0
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1a2
rm $ECLAT_DATA/s1e3k5d5e1a2.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1a6.ha -a 6.0
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1a6
rm $ECLAT_DATA/s1e3k5d5e1a6.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1a8.ha -a 8.0
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1a8
rm $ECLAT_DATA/s1e3k5d5e1a8.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1a12.ha -a 12.0
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1a12
rm $ECLAT_DATA/s1e3k5d5e1a12.ha

#generate data for testing effect on cluster dimensionality
$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1c1.ha -c 0.1
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1c1
rm $ECLAT_DATA/s1e3k5d5e1c1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1c3.ha -c 0.3
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1c3
rm $ECLAT_DATA/s1e3k5d5e1c3.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1c7.ha -c 0.7
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1c7 
rm $ECLAT_DATA/s1e3k5d5e1c7.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1c9.ha -c 0.9
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1c9
rm $ECLAT_DATA/s1e3k5d5e1c9.ha

#generate data for testing effect of fraction of overlapping dim in clusters
$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1o1.ha -O 0.1
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1o1
rm $ECLAT_DATA/s1e3k5d5e1o1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1o3.ha -O 0.3
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1o3
rm $ECLAT_DATA/s1e3k5d5e1o3.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1o7.ha -O 0.7
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1o7 
rm $ECLAT_DATA/s1e3k5d5e1o7.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1o9.ha -O 0.9
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1o9
rm $ECLAT_DATA/s1e3k5d5e1o9.ha

#generate data for testing effect of $w$
$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1w1e1.ha -w 10
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1w1e1
rm $ECLAT_DATA/s1e3k5d5e1w1e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1w2e1.ha -w 20 
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1w2e1 
rm $ECLAT_DATA/s1e3k5d5e1w2e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1w3e1.ha -w 30
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1w3e1
rm $ECLAT_DATA/s1e3k5d5e1w3e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1w4e1.ha -w 40 
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1w4e1
rm $ECLAT_DATA/s1e3k5d5e1w4e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1m25.ha -m 1 -w 25
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1m25 
rm $ECLAT_DATA/s1e3k5d5e1m25.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1m50.ha -m 1 -w 50 
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1m50
rm $ECLAT_DATA/s1e3k5d5e1m50.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1m75.ha -m 1 -w 75
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1m75
rm $ECLAT_DATA/s1e3k5d5e1m75.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1m100.ha -m 1 -w 100
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1m100
rm $ECLAT_DATA/s1e3k5d5e1m100.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1m110.ha -m 1 -w 110
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1m110
rm $ECLAT_DATA/s1e3k5d5e1m110.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1m125.ha -m 1 -w 125
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1m125
rm $ECLAT_DATA/s1e3k5d5e1m125.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1m140.ha -m 1 -w 140 
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1m140 
rm $ECLAT_DATA/s1e3k5d5e1m140.ha

#generate data for testing effect of $k$
$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k3d5e1.ha -k 3 
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k3d5e1 -k 3
rm $ECLAT_DATA/s1e3k3d5e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k8d5e1.ha -k 8
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k8d5e1 -k 8
rm $ECLAT_DATA/s1e3k8d5e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k1e1d5e1.ha -k 10 
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k1e1d5e1 -k 10
rm $ECLAT_DATA/s1e3k1e1d5e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k12d5e1.ha -k 12 
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k12d5e1 -k 12
rm $ECLAT_DATA/s1e3k12d5e1.ha

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k15d5e1.ha -k 15
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k15d5e1 -k 15
rm $ECLAT_DATA/s1e3k15d5e1.ha

#generate data for testing effect of $\xi$
$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1x5.ha
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1x5 -x 5

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1x8.ha
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1x8 -x 8

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1x12.ha
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1x12 -x 12

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1x15.ha
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1x15 -x 15

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1x20.ha
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1x20 -x 20

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1x25.ha
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1x25 -x 25

$VC_OBJ/dataGen -o $ECLAT_DATA/s1e3k5d5e1x30.ha
$VC_OBJ/convertData -i $ECLAT_DATA/s1e3k5d5e1x30 -x 30
