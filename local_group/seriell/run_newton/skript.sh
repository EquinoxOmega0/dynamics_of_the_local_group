./real.out
./mkinput.out
let ngen=50
let nmod=200
let nparent=10
for u  in $(seq 1 $ngen)
do
for i in $(seq 1 $nmod)
do
cp old${i}.dat input.dat
./mkmodel.out
done
for i in $(seq 1 $nmod)
do
cp mod${i}.dat mod.dat
./newhexi.out
cp fmod.dat fmod${i}.dat
done
for i in $(seq 1 $nmod)
do
atos fmod${i}.dat fmod${i}.snp
done
./cutter.out
if [ $u -lt $ngen ]
then ./genetic.out
else ./outputbest.out
fi
mkdir generation${u}
cp protocol.dat generation${u}
cp old*.dat generation${u}
cp fmod*.snp generation${u}
rm protocol.dat
rm input.dat
rm old*.dat
rm fmod*.snp
rm mod*.snp
rm mod*.snp2
rm mod*.dat
rm fmod*.dat
for i in $(seq 1 $nmod)
do
cp new${i}.dat old${i}.dat
done
rm new*.dat
done
for i in $(seq 1 $nparent)
do
cp old${i}.dat input.dat
./mkmodel.out
done
for i in $(seq 1 $nparent)
do
cp mod${i}.dat mod.dat
./newhexi.out
cp fmod.dat fmod${i}.dat
done
for i in $(seq 1 $nparent)
do
atos fmod${i}.dat fmod${i}.snp
done
rm fmod*.dat
rm mod*.dat
mkdir result
cp old*.dat result
cp fmod*.snp result
rm fmod*.snp
rm old*.dat
exit