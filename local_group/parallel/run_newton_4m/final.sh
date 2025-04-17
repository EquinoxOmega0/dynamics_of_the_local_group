nparent=$1
cp c2.dat c.dat
for i in $( seq 1 ${nparent} ); do
cp para2.dat mr${i}/para.dat
cp old${i}.dat mr${i}/input.dat
cd mr${i}
qsub -V mkmodel.sh
cd ..
let rest=$i%20
if [ $rest -eq 0 ]
then sleep 5
fi
done
count=$(qstat -u saulder | wc -l)
while [ $count -ne 0 ]
do
sleep 2
count=$(qstat -u saulder | wc -l)
done
for i in $( seq 1 ${nparent} ); do
cd mr${i}
cp mod${i}.dat mod.dat
rm mod${i}.dat
qsub -V newhexi.sh
cd ..
let rest=$i%20
if [ $rest -eq 0 ]
then sleep 10
fi
done
count=$(qstat -u saulder | wc -l)
while [ $count -ne 0 ]
do
sleep 3
count=$(qstat -u saulder | wc -l)
done
for i in $( seq 1 ${nparent} ); do
cd mr${i}
cp fmod.dat ../fmod${i}.dat
cd ..
done
mkdir result
cp old*.dat result
cp fmod*.dat result
cp convergence*.dat result
rm fmod*.dat
rm old*.dat
rm convergence*.dat
rm -r mr*