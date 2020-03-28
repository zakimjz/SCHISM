#Author: Karlton Sequeira

#declare variables and paths
export PROJ_HOME=/cs/sequek/laptop/projects

export VC_HOME=$PROJ_HOME/VC
export VC_DATA=$VC_HOME/data
export VC_SRC=$VC_HOME/src
export VC_OBJ=$VC_HOME/obj

export SFA_HOME=$PROJ_HOME/sfa
export SFA_DATA=$SFA_HOME/data
export SFA_SRC=$SFA_HOME/src
export SFA_OBJ=$SFA_HOME/obj
export NUM="2"

mkdir $SFA_HOME/data
g++ -O4 -o $VC_OBJ/dataGen $VC_SRC/dataGen.cc
g++ -O4 -o $VC_OBJ/convertData $VC_SRC/convertData.cc

#generate data for testing effect of number of dimensions
$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d2e2.ha -d 200 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d2e2 -d 200 -s $NUM
rm $SFA_DATA/s1e3k7d2e2s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d1e2.ha -d 100 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d1e2 -d 100 -s $NUM
rm $SFA_DATA/s1e3k7d1e2s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1.ha -d 50 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1 -d 50 -s $NUM
rm $SFA_DATA/s1e3k7d5e1s0?.ha

#generate data for testing effect on cluster dimensionality
$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1r1.ha -c 0.1 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1r1 -s $NUM
rm $SFA_DATA/s1e3k7d5e1r1s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1r3.ha -c 0.3 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1r3 -s $NUM
rm $SFA_DATA/s1e3k7d5e1r3s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1r7.ha -c 0.7 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1r7  -s $NUM
rm $SFA_DATA/s1e3k7d5e1r7s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1r9.ha -c 0.9 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1r9 -s $NUM
rm $SFA_DATA/s1e3k7d5e1r9s0?.ha

#generate data for testing effect of fraction of overlapping dim in clusters
$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1o1.ha -O 0.1 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1o1 -s $NUM
rm $SFA_DATA/s1e3k7d5e1o1s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1o3.ha -O 0.3 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1o3 -s $NUM
rm $SFA_DATA/s1e3k7d5e1o3s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1o7.ha -O 0.7 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1o7 -s $NUM
rm $SFA_DATA/s1e3k7d5e1o7s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1o9.ha -O 0.9 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1o9 -s $NUM
rm $SFA_DATA/s1e3k7d5e1o9s0?.ha

#generate data for testing effect of $\xi$
$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1x5.ha -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1x5 -x 5 -s $NUM

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1x8.ha -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1x8 -x 8 -s $NUM

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1x12.ha -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1x12 -x 12 -s $NUM

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1x15.ha -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1x15 -x 15 -s $NUM

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1x20.ha -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1x20 -x 20 -s $NUM

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1x25.ha -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1x25 -x 25 -s $NUM

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k7d5e1x30.ha -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k7d5e1x30 -x 30 -s $NUM

#generate data for testing effect of $k$
$VC_OBJ/dataGen -o $SFA_DATA/s1e3k5d5e1.ha -k 5 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k5d5e1 -k 5 -s $NUM
rm $SFA_DATA/s1e3k5d5e1s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k1e1d5e1.ha -k 10  -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k1e1d5e1 -k 10 -s $NUM
rm $SFA_DATA/s1e3k1e1d5e1s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k12d5e1.ha -k 12  -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k12d5e1 -k 12 -s $NUM
rm $SFA_DATA/s1e3k12d5e1s0?.ha

$VC_OBJ/dataGen -o $SFA_DATA/s1e3k15d5e1.ha -k 15 -b .03 -s $NUM
$VC_OBJ/convertData -i $SFA_DATA/s1e3k15d5e1 -k 15 -s $NUM
rm $SFA_DATA/s1e3k15d5e1s0?.ha
