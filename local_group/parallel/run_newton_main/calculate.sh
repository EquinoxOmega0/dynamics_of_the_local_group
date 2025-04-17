nmod=$1
for i in $( seq 1 ${nmod} ); do
cp old${i}.dat mr${i}/input.dat
cd mr${i}
qsub -V mkmodel.sh
cd ..
let rest=$i%25
if [ $rest -eq 0 ]
then sleep 2
fi
done
count=$(qstat -u saulder | wc -l)
while [ $count -ne 0 ]
do
sleep 2
count=$(qstat -u saulder | wc -l)
done
for i in $( seq 1 ${nmod} ); do
cd mr${i}
cp mod${i}.dat mod.dat
rm mod${i}.dat
qsub -V newhexi.sh
cd ..
let rest=$i%25
if [ $rest -eq 0 ]
then sleep 5
fi
done 
count=$(qstat -u saulder | wc -l)
while [ $count -ne 0 ]
do
sleep 5
count=$(qstat -u saulder | wc -l)
done
for i in $( seq 1 ${nmod} ); do
cd mr${i}
cp fmod.dat ../fmod${i}.dat
cd ..
done 
