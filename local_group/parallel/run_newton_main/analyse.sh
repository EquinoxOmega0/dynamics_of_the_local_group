nmod=$1
ngen=$2
ii=$3
nparent=$4
qsub -V cutter.sh
count=$(qstat -u saulder | wc -l)
while [ $count -ne 0 ]
do
sleep 5
count=$(qstat -u saulder | wc -l)
done
if [ $ii -lt $ngen ]
then qsub -V genetic.sh
else qsub -V outputbest.sh
fi
count=$(qstat -u saulder | wc -l)
while [ $count -ne 0 ]
do
sleep 3
count=$(qstat -u saulder | wc -l)
done
mkdir generation${ii}
cp protocol.dat generation${ii}
cp convergence.dat convergence${ii}.dat
cp old*.dat generation${ii}
rm protocol.dat
rm old*.dat
rm mod*.dat
rm fmod*.dat
if [ $ii -lt $ngen ]
then 
for i in $( seq 1 ${nmod} ); do
cp new${i}.dat old${i}.dat
done
else
for i in $( seq 1 ${nparent} ); do
cp new${i}.dat old${i}.dat
done
fi
rm new*.dat 