nmod=$1
qsub -V mkinput.sh
for i in $( seq 1 ${nmod} ); do
mkdir mr${i}
cp para.dat mr${i}/para.dat
cp mkmodel.out mr${i}/mkmodel.out
cp mkmodel.sh mr${i}/mkmodel.sh
cp demoni.out mr${i}/demoni.out
cp demoni.sh mr${i}/demoni.sh
done 
qsub -V real.sh
count=$(qstat -u saulder | wc -l)
while [ $count -ne 0 ]
do
sleep 1
count=$(qstat -u saulder | wc -l)
done
