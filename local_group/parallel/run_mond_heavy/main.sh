ngen=20
nmod=1000
nparent=100
./init.sh ${nmod}
for ii in $( seq 1 ${ngen} ); do
./calculate.sh ${nmod}
./analyse.sh ${nmod} ${ngen} ${ii} ${nparent}
done
./final.sh ${nparent}
exit