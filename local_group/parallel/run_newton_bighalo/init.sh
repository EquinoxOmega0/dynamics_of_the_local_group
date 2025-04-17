nmod=$1
qsub -V mkinput.sh
for i in $( seq 1 ${nmod} ); do
mkdir mr${i}
cp para.dat mr${i}/para.dat
cp halo.dat mr${i}/halo.dat
cp mkmodel.out mr${i}/mkmodel.out
cp mkmodel.sh mr${i}/mkmodel.sh
cp newhexi.out mr${i}/newhexi.out
cp newhexi.sh mr${i}/newhexi.sh
done 
qsub -V real.sh
count=$(qstat -u saulder | wc -l)
while [ $count -ne 0 ]
do
sleep 1
count=$(qstat -u saulder | wc -l)
done
