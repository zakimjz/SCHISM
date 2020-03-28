clear
cd "/cs/sequek/laptop/projects/eclat/data"
set nologscale
set size 1,0.5
set terminal postscript portrait enhanced mono dashed defaultplex "Helvetica" 11
set ylabel "Score"
#set logscale y

set xlabel "Number of dimensions(d)"
set out "/cs/sequek/laptop/docs/schism/PerfVsDim.ps"
plot "/cs/sequek/laptop/projects/eclat/data/d_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/d_t" using 1:2 title "Speed(minutes)" w l

set xlabel "Number of points(n)"
set out "/cs/sequek/laptop/docs/schism/PerfVsNum.ps"
plot "/cs/sequek/laptop/projects/eclat/data/n_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/n_t" using 1:2 title "Speed(minutes)" w l

set xlabel "Number of clusters(k)"
set out "/cs/sequek/laptop/docs/schism/PerfVsk.ps"
plot "/cs/sequek/laptop/projects/eclat/data/k_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/k_t" using 1:2 title "Speed(minutes)" w l

set xlabel "Number of intervals per division(Xi)"
set out "/cs/sequek/laptop/docs/schism/PerfVsx.ps"
plot "/cs/sequek/laptop/projects/eclat/data/x_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/x_t" using 1:2 title "Speed(minutes)" w l

set xlabel "Kappa"
set out "/cs/sequek/laptop/docs/schism/PerfVsKappa.ps"
plot  "/cs/sequek/laptop/projects/eclat/data/cn_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/cn_t" using 1:2 title "Speed(minutes)" w l

set xlabel "log_{10}(Tau)"
set out "/cs/sequek/laptop/docs/schism/PerfVst.ps"
plot "/cs/sequek/laptop/projects/eclat/data/t_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/t_t" using 1:2 title "Speed(minutes)" w l

set xlabel "P(dimension in subspace is constrained)=c"
set out "/cs/sequek/laptop/docs/schism/PerfVsc.ps"
plot  "/cs/sequek/laptop/projects/eclat/data/cd_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/cd_t" using 1:2 title "Speed(minutes)" w l

set xlabel "Standard deviation of constrained dimension"
set out "/cs/sequek/laptop/docs/schism/PerfVsw.ps"
plot "/cs/sequek/laptop/projects/eclat/data/w_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/w_t" using 1:2 title "Speed(minutes)" w l
#plot "/cs/sequek/laptop/projects/eclat/data/w_t2" using 1:5 title "CLIQUE Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/w_t" using 1:5 title "SCHISM Coverage" w l, "/cs/sequek/laptop/projects/eclat/data/w_t2" using 1:2 title "CLIQUE Running time(Minutes)" w l,"/cs/sequek/laptop/projects/eclat/data/w_t" using 1:2 title "SCHISM Running time(minutes)" w l

set xlabel "P(same dimension in adjacent subspaces are constrained)=o"
set out "/cs/sequek/laptop/docs/schism/PerfVso.ps"
plot "/cs/sequek/laptop/projects/eclat/data/o_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/o_t" using 1:2 title "Speed(minutes)" w l

set xlabel "Average width of interval of constrained dimension"
set out "/cs/sequek/laptop/docs/schism/PerfVsDenseHyperRect.ps"
plot "/cs/sequek/laptop/projects/eclat/data/m_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/m_t" using 1:2 title "Running time(minutes)" w l
#plot "/cs/sequek/laptop/projects/eclat/data/comp_t2" using 1:5 title "CLIQUE Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/comp_t" using 1:5 title "SCHISM Coverage" w l, "/cs/sequek/laptop/projects/eclat/data/comp_t2" using 1:2 title "CLIQUE Running time(Minutes)" w l,"/cs/sequek/laptop/projects/eclat/data/comp_t" using 1:2 title "SCHISM Running time(minutes)" w l

set xlabel "P(adjacent subspaces with nearby means have overlapping constrained dimension)"
set out "/cs/sequek/laptop/docs/schism/PerfVsp.ps"
plot "/cs/sequek/laptop/projects/eclat/data/c_t" using 1:5 title "Coverage" w l,"/cs/sequek/laptop/projects/eclat/data/c_t" using 1:2 title "Speed(minutes)" w l

set xlabel "Average width of interval of constrained dimension"
set out "/cs/sequek/laptop/docs/schism/PerfVsDimClique.ps"
plot "/cs/sequek/laptop/projects/eclat/data/comp_t2" using 1:2 title "CLIQUE Running time(Minutes)" w l,"/cs/sequek/laptop/projects/eclat/data/comp_t" using 1:2 title "SCHISM Running time(minutes)" w l
