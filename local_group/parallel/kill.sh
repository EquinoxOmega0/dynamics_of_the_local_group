nstart=$1 
nstop=$2
for i in $( seq ${nstart} ${nstop} ); do
qdel ${i}
done