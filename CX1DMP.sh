cd $HOME/DMP

#echo ". build_CX1.sh"
#echo

#echo
#. build_CX1.sh
#echo
SIMAREA=$TMPDIR/lj1913/DMP


echo
echo "mkdir -p $SIMAREA"
echo

mkdir -p $SIMAREA

echo "cd $SIMAREA"
echo
cd $SIMAREA
OUTDIR=$TMPDIR/lj1913/DMP/Output
echo "mkdir -p $OUTDIR"
echo
mkdir -p $OUTDIR 

#JOBSCRIPT=JobScripts/runFields.pbs
#JOBSCRIPT=JobScripts/runPotentials.pbs
#JOBSCRIPT=JobScripts/runOnce.pbs
#JOBSCRIPT=JobScripts/runRepeat_singlenode.pbs
#JOBSCRIPT=JobScripts/runRepeat_throughput.pbs
#JOBSCRIPT=JobScripts/runRepeat_pqplasma.pbs
JOBSCRIPT=JobScripts/jobscript.pbs

echo "Submitting Job:"
echo
echo "qsub $HOME/$JOBSCRIPT"
echo
qsub $HOME/DMP/$JOBSCRIPT
echo
echo "Script complete!" 
cd $HOME/DMP
