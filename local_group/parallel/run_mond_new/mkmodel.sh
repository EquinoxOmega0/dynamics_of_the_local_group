#!/bin/sh
#echo "in script"
### Job Name
#PBS -N Make_Model

### Submits the job to the default queue (not really needed since we have only
### one queue which takes all jobs)
#PBS -q workq 

### Request nodes:cores per node
### Note: ppn=1 will not work as expected (but you get (nodes x ppn) cores for 
### sure). If you really want only one process per node, you have to request
### all cores on the nodes but use only one of them.
### for superusers:this could possibly be solved by putting:
###                JOBNODEMATCHPOLICY  EXACTNODE into /opt/maui/maui.cfg 
###                (and maybe NODEACCESSPOLICY  SHARED)  
#PBS -l nodes=1:ppn=8

### Request the nodes for estimated time (not really needed since we have only
### one queue which takes all jobs)
#PBS -l walltime=10:00:00

#echo "before cd"
#Change to directory from which job was submitted
cd $PBS_O_WORKDIR

# set number of processors to run on (list of node names is in file $PBSNODES)
NUM_CORES=`wc -l $PBS_NODEFILE | awk '{ print $1 }'`
NUM_NODES=`sort -u $PBS_NODEFILE | wc -l`

echo $NUM_CORES #for debugging
echo $NUM_NODES #for debugging
cat $PBS_NODEFILE > nodefile.tst  #for debugging

#echo "before mdpboot"
# Start mpd processes on nodes requested nodes
mpdboot -n $NUM_NODES -f $PBS_NODEFILE
#mpdboot -n 8 -f $PBS_NODEFILE #for debugging

export LD_LIBRARY_PATH=:/opt/intel/cce/10.1.015/lib:/opt/intel/fce/10.1.015/lib

echo "Job $PBS_JOBID $PBS_JOBNAME started at `date`"
echo ""
mpiexec -np $NUM_CORES ./mkmodel.out
mpdtrace #for debugging
echo ""
echo "Job $PBS_JOBID ended at `date`"

# Stop all mpd processes on nodes
mpdallexit



