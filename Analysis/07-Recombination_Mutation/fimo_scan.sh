motif=$1
seq=$2
dir=$3

date

set -e

source /share/home/zhanglab/user/liuanguo/anaconda3/bin/activate bio

mkdir -p ${dir}

fimo --verbosity 2 --thresh 1e-4 --oc ${dir} ${motif} ${seq}
  
status=$?

echo

if [ $status -eq 0 ]; then
    echo "completed successfully."
    echo
    echo "finish!"
    echo
    echo $(date)
else
    echo "exit with error!!"
    echo
    echo $(date)
fi


