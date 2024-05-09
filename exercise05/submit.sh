#BSUB -J pythonSC # name
#BSUB -o pythonSC_%J.out # output file
#BSUB -q gpuh100
#BSUB -n 32 ## cores
#BSUB -R "rusage[mem=1GB]" 
#BSUB -W 500 # useable time in minutes
##BSUB -N # send mail when done
#BSUB -R "span[hosts=1]"
#BSUB -gpu "num=2:mode=exclusive_process"


module load python3/3.10.13
pip3 install torch
pip3 install tqdm
pip3 install scipy

## python3 main.py 1
python3 main.py 2
## python3 main.py 3