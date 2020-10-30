JID1=`qsub -h run_economy.sh`
JID2=`qsub -W depend=afterok:$JID1 run_economy.sh`
JID3=`qsub -W depend=afterok:$JID2 run_economy.sh`
#JID4=`qsub -W depend=afterok:$JID3 run_economy.sh`
qrls $JID1

