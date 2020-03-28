#Author: Karlton Sequeira

 declare variables and paths
export PROJ_HOME=$HOME/laptop/projects
export ECLAT_HOME=$PROJ_HOME/eclat
export ECLAT_DATA=$ECLAT_HOME/data

make clean
make
make dna_eclat

##############################################################################
#test effect of number of dimensions
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm 2>$ECLAT_DATA/d_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d1e2.ibm 2>>$ECLAT_DATA/d_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d2e2.ibm 2>>$ECLAT_DATA/d_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d4e2.ibm 2>>$ECLAT_DATA/d_t

$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -p 2 2>$ECLAT_DATA/d_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d1e2.ibm -p 2 2>>$ECLAT_DATA/d_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d2e2.ibm -p 2 2>>$ECLAT_DATA/d_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d4e2.ibm -p 2 2>>$ECLAT_DATA/d_t2

#test effect of number of points
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm 2>$ECLAT_DATA/n_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e4k5d5e1.ibm 2>>$ECLAT_DATA/n_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s5e4k5d5e1.ibm 2>>$ECLAT_DATA/n_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e5k5d5e1.ibm 2>>$ECLAT_DATA/n_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s2e5k5d5e1.ibm 2>>$ECLAT_DATA/n_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s4e5k5d5e1.ibm 2>>$ECLAT_DATA/n_t

$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -p 2 2>$ECLAT_DATA/n_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e4k5d5e1.ibm -p 2 2>>$ECLAT_DATA/n_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s5e4k5d5e1.ibm -p 2 2>>$ECLAT_DATA/n_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e5k5d5e1.ibm -p 2 2>>$ECLAT_DATA/n_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s2e5k5d5e1.ibm -p 2 2>>$ECLAT_DATA/n_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s4e5k5d5e1.ibm -p 2 2>>$ECLAT_DATA/n_t2

#test effect of ratio of subspace strengths
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1a2.ibm 2>$ECLAT_DATA/cn_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm 2>>$ECLAT_DATA/cn_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1a6.ibm 2>>$ECLAT_DATA/cn_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1a8.ibm 2>>$ECLAT_DATA/cn_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1a12.ibm 2>>$ECLAT_DATA/cn_t

#test effect of number of subspace dimensions
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1c1.ibm 2>$ECLAT_DATA/cd_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1c3.ibm 2>>$ECLAT_DATA/cd_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm 2>>$ECLAT_DATA/cd_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1c7.ibm 2>>$ECLAT_DATA/cd_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1c9.ibm 2>>$ECLAT_DATA/cd_t

#test effect of $o$ i.e. overlap probability
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1o1.ibm 2>$ECLAT_DATA/o_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1o3.ibm 2>>$ECLAT_DATA/o_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm 2>>$ECLAT_DATA/o_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1o7.ibm 2>>$ECLAT_DATA/o_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1o9.ibm 2>>$ECLAT_DATA/o_t

#test effect of $w$
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1w1e1.ibm 2>$ECLAT_DATA/w_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1w2e1.ibm 2>>$ECLAT_DATA/w_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1w3e1.ibm 2>>$ECLAT_DATA/w_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1w4e1.ibm 2>>$ECLAT_DATA/w_t

$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1w1e1.ibm -p 2 2>$ECLAT_DATA/w_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1w2e1.ibm -p 2 2>>$ECLAT_DATA/w_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1w3e1.ibm -p 2 2>>$ECLAT_DATA/w_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1w4e1.ibm -p 2 2>>$ECLAT_DATA/w_t2

#test effect of widening side of hyper-rectangular subspace
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m25.ibm 2>$ECLAT_DATA/m_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m50.ibm 2>>$ECLAT_DATA/m_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m75.ibm 2>>$ECLAT_DATA/m_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m100.ibm 2>>$ECLAT_DATA/m_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m110.ibm 2>>$ECLAT_DATA/m_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m125.ibm 2>>$ECLAT_DATA/m_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m140.ibm 2>>$ECLAT_DATA/m_t

#compare effect of widening side of hyper-rectangular subspace on thresh
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m75.ibm 2>$ECLAT_DATA/comp_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m100.ibm 2>>$ECLAT_DATA/comp_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m110.ibm 2>>$ECLAT_DATA/comp_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m125.ibm 2>>$ECLAT_DATA/comp_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m135.ibm 2>>$ECLAT_DATA/comp_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m145.ibm 2>>$ECLAT_DATA/comp_t

$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m75.ibm -p 2 2>$ECLAT_DATA/comp_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m100.ibm -p 2 2>>$ECLAT_DATA/comp_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m110.ibm -p 2 2>>$ECLAT_DATA/comp_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m125.ibm -p 2 2>>$ECLAT_DATA/comp_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m135.ibm -p 2 2>>$ECLAT_DATA/comp_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1m145.ibm -p 2 2>>$ECLAT_DATA/comp_t2

#test effect of $k$
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k3d5e1.ibm 2>$ECLAT_DATA/k_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k8d5e1.ibm 2>>$ECLAT_DATA/k_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k12d5e1.ibm 2>>$ECLAT_DATA/k_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k15d5e1.ibm 2>>$ECLAT_DATA/k_t

$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k3d5e1.ibm -p 2 2>$ECLAT_DATA/k_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k8d5e1.ibm -p 2 2>>$ECLAT_DATA/k_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k12d5e1.ibm -p 2 2>>$ECLAT_DATA/k_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k15d5e1.ibm -p 2 2>>$ECLAT_DATA/k_t2

#test effect of $\tau$ i.e. interestingness
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -r 4 2>$ECLAT_DATA/t_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -r 2 2>>$ECLAT_DATA/t_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm 2>>$ECLAT_DATA/t_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -r 0.5 2>>$ECLAT_DATA/t_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -r 0.25 2>>$ECLAT_DATA/t_t

#test effect of $\xi$ i.e. number of divisions per dimension
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1x5.ibm -x 5 2>$ECLAT_DATA/x_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1x8.ibm -x 8 2>>$ECLAT_DATA/x_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm 2>>$ECLAT_DATA/x_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1x12.ibm -x 12 2>>$ECLAT_DATA/x_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1x15.ibm -x 15 2>>$ECLAT_DATA/x_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1x20.ibm -x 20 2>>$ECLAT_DATA/x_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1x25.ibm -x 25 2>>$ECLAT_DATA/x_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1x30.ibm -x 30 2>>$ECLAT_DATA/x_t

#test effect of $s$ i.e. support
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -s .01 2>$ECLAT_DATA/s_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -s .03 2>>$ECLAT_DATA/s_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -s .05 2>>$ECLAT_DATA/s_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -s .07 2>>$ECLAT_DATA/s_t
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -s .09 2>>$ECLAT_DATA/s_t

$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -s .01 -p 2 2>$ECLAT_DATA/s_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -s .03 -p 2 2>>$ECLAT_DATA/s_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -s .05 -p 2 2>>$ECLAT_DATA/s_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -s .07 -p 2 2>>$ECLAT_DATA/s_t2
$ECLAT_HOME/eclat -i $ECLAT_DATA/s1e3k5d5e1.ibm -s .09 -p 2 2>>$ECLAT_DATA/s_t2

#test CLIQUE and SCHISM on real datasets
#$ECLAT_HOME/eclat -i $ECLAT_DATA/pendigits.ibm -s .01 -r .45 2>$ECLAT_DATA/schism_op
#$ECLAT_HOME/dna_eclat -i $ECLAT_DATA/yeast5.ibm -x 5 -r .5 1>Clusters 2>>$ECLAT_DATA/schism_op
#$ECLAT_HOME/req $ECLAT_DATA/t4 4 Clusters .5 2>>$ECLAT_DATA/schism_op

#$ECLAT_HOME/eclat -i $ECLAT_DATA/pendigits.ibm -s .01 -r .5 -p 2 2>$ECLAT_DATA/clique_op
#$ECLAT_HOME/dna_eclat -i $ECLAT_DATA/yeast5.ibm -x 5 -r .5 -p 2 1>Clusters 2>>$ECLAT_DATA/clique_op
#$ECLAT_HOME/req $ECLAT_DATA/t4 4 Clusters .5 2>>$ECLAT_DATA/clique_op
