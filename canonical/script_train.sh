#$ -q gpu
#$ -cwd
#$ -N train_spliceai
#$ -e Logs/
#$ -o Logs/
#$ -l gpus=2
#$ -l h_vmem=500g

python2.7 -u train_model.py $1 $2 > Outputs/SpliceAI${1}_c${2}.txt
