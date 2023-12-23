for FILE in job-*.sh; do
echo ${FILE}
sbatch ${FILE}
sleep 1 # pause to be kind to the scheduler
done
