#Author: Karlton Sequeira

 declare variables and paths
export PROJ_HOME=$HOME/laptop/projects
export SFA_HOME=$PROJ_HOME/sfa
export SFA_OBJ=$SFA_HOME/obj
export SFA_DATA=$SFA_HOME/data
export SFA_SRC=$SFA_HOME/src

mkdir $SFA_HOME/obj
g++ -O4 -I $HOME/work/zaki/paging/boost_1_24_0/ -o $SFA_OBJ/match $SFA_SRC/Mag.cc

##############################################################################
#test effect of number of dimensions
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm 2>$SFA_DATA/d_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d1e2_out.ibm 2>>$SFA_DATA/d_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d2e2_out.ibm 2>>$SFA_DATA/d_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d3e2_out.ibm 2>>$SFA_DATA/d_t

#test effect of number of points
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm 2>$SFA_DATA/n_t
$SFA_OBJ/match -i $SFA_DATA/s1e4k7d5e1_out.ibm 2>>$SFA_DATA/n_t
$SFA_OBJ/match -i $SFA_DATA/s1e5k7d5e1_out.ibm 2>>$SFA_DATA/n_t
$SFA_OBJ/match -i $SFA_DATA/s2e5k7d5e1_out.ibm 2>>$SFA_DATA/n_t
$SFA_OBJ/match -i $SFA_DATA/s3e5k7d5e1_out.ibm 2>>$SFA_DATA/n_t

#test effect of ratio of subspace strengths
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1p2_out.ibm 2>$SFA_DATA/cn_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm 2>>$SFA_DATA/cn_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1p6_out.ibm 2>>$SFA_DATA/cn_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1p8_out.ibm 2>>$SFA_DATA/cn_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1p12_out.ibm 2>>$SFA_DATA/cn_t

#test effect of number of subspace dimensions
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1r1_out.ibm -e .1 2>$SFA_DATA/cd_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1r3_out.ibm -e .3 2>>$SFA_DATA/cd_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm -e .5 2>>$SFA_DATA/cd_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1r7_out.ibm -e .7 2>>$SFA_DATA/cd_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1r9_out.ibm -e .9 2>>$SFA_DATA/cd_t

#test effect of $\xi$ i.e. number of divisions per dimension
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1x5_out.ibm -x 5 2>$SFA_DATA/x_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1x8_out.ibm -x 8 2>>$SFA_DATA/x_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm 2>>$SFA_DATA/x_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1x12_out.ibm -x 12 2>>$SFA_DATA/x_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1x15_out.ibm -x 15 2>>$SFA_DATA/x_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1x20_out.ibm -x 20 2>>$SFA_DATA/x_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1x25_out.ibm -x 25 2>>$SFA_DATA/x_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1x30_out.ibm -x 30 2>>$SFA_DATA/x_t

#test effect of $\tau$ i.e. interestingness
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm -r 4 2>$SFA_DATA/t_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm -r 2 2>>$SFA_DATA/t_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm 2>>$SFA_DATA/t_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm -r 0.5 2>>$SFA_DATA/t_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm -r 0.25 2>>$SFA_DATA/t_t

#test effect of $o$ i.e. overlap probability
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1o1_out.ibm 2>$SFA_DATA/o_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1o3_out.ibm 2>>$SFA_DATA/o_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm 2>>$SFA_DATA/o_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1o7_out.ibm 2>>$SFA_DATA/o_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1o9_out.ibm 2>>$SFA_DATA/o_t

#test effect of widening side of hyper-rectangular subspace
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m25_out.ibm 2>$SFA_DATA/m_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m50_out.ibm 2>>$SFA_DATA/m_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m75_out.ibm 2>>$SFA_DATA/m_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m100_out.ibm 2>>$SFA_DATA/m_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m125_out.ibm 2>>$SFA_DATA/m_t

#compare effect of widening side of hyper-rectangular subspace on thresh
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m75_out.ibm 2>$SFA_DATA/comp_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m100_out.ibm 2>>$SFA_DATA/comp_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m110_out.ibm 2>>$SFA_DATA/comp_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m125_out.ibm 2>>$SFA_DATA/comp_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m135_out.ibm 2>>$SFA_DATA/comp_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1m145_out.ibm 2>>$SFA_DATA/comp_t

#test effect of constraining centres of clusters to lie within 2 std dev of each other
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1c1_out.ibm 2>$SFA_DATA/c_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1c3_out.ibm 2>>$SFA_DATA/c_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm 2>>$SFA_DATA/c_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1c7_out.ibm 2>>$SFA_DATA/c_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1c9_out.ibm 2>>$SFA_DATA/c_t

#test effect of $s$ i.e. support
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm -s .01 2>$SFA_DATA/s_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm -s .03 2>>$SFA_DATA/s_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm -s .05 2>>$SFA_DATA/s_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm -s .07 2>>$SFA_DATA/s_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1_out.ibm -s .09 2>>$SFA_DATA/s_t

#test effect of $w$
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1w1e1_out.ibm 2>$SFA_DATA/w_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1w2e1_out.ibm 2>>$SFA_DATA/w_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1w3e1_out.ibm 2>>$SFA_DATA/w_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k7d5e1w4e1_out.ibm 2>>$SFA_DATA/w_t

#test effect of $k$
$SFA_OBJ/match -i $SFA_DATA/s1e3k3d5e1w2e1_out.ibm 2>$SFA_DATA/k_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k8d5e1w2e1_out.ibm 2>>$SFA_DATA/k_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k12d5e1w2e1_out.ibm 2>>$SFA_DATA/k_t
$SFA_OBJ/match -i $SFA_DATA/s1e3k15d5e1w2e1_out.ibm 2>>$SFA_DATA/k_t

#test CLIQUE and SCHISM on real datasets
#$SFA_OBJ/match -i $SFA_DATA/pendigits_out.ibm -s .01 -r .45 2>$SFA_DATA/schism_op
#$SFA_OBJ/dna_match -i $SFA_DATA/yeast5_out.ibm -x 5 -r .5 1>Clusters 2>>$SFA_DATA/schism_op
#$SFA_OBJ/req $SFA_DATA/t4 4 Clusters .5 2>>$SFA_DATA/schism_op
