usage="$(basename "$0") [-h] [-c s] [-o s] [-s s] [-g i] [-r f]\n
-------------------------------------------------------------------\n
PROGRAM: MATCHING TEST\n
-------------------------------------------------------------------\n
where:\n
    -h  show this help message\n
    -c  [string]    config file, format: {[vcf path]\t[chromosome name]\t[chromosome length]\t[split number]} for each line\n
    -o  [string]    output directory\n
    -s  [string]    sample list\n
    -g  [int]       gap size [default: 0]\n
    -r  [float]     snp marker minor ratio cutoff [default: 0.0]\n"

# declare -f exp
SCRIPT_PATH="/mnt/disk2_workspace/zhaozicheng/zzc/swine/script"
MATCHING="/home/zhaozicheng/donkey/updated_src/matching"

config=""
out_dir=""
sample_list=""
gap_size=0
cutoff=0

while getopts ':hc:o:s:g:r:' option
do
    case "$option" in
        h) echo -e $usage
            exit
            ;;
        c) config=$OPTARG
            ;;
        o) out_dir=$OPTARG
            ;;
        s) sample_list=$OPTARG
            ;;
        g) gap_size=$OPTARG
            ;;
        r) cutoff=$OPTARG
            ;;
        :) printf "missing argument for -$s\n" "$OPTARG" >&2
            echo -e $usage >&2
            exit 1
            ;;
        \?) printf "illegal option: -%s\n" "$OPTARG" >&2
            echo -e $usage >&2
            exit 1
            ;;
    esac
done
shift $(($OPTIND - 1))

if [ -z $config ]
then
    printf "missing config file" >&2
    echo -e $usage >&2
    exit 1
fi

if [ -z $out_dir ]
then
    printf "missing output directory" >&2
    echo -e $usage >&2
    exit 1
fi

if [ -z $sample_list ]
then
    printf "missing sample list" >&2
    echo -e $usage >&2
    exit 1
fi

CUR_DIR=$PWD

if [ ! -d $out_dir ]
then
    mkdir $out_dir
else
    rm -rf $out_dir
    mkdir $out_dir
fi

cd $out_dir
echo "entering $out_dir..."

if [ ! -d simulate ]
then
    mkdir simulate
fi

cd simulate
echo "simulating..."
while read -r path name len num
do
    echo "simulate $name, written to $PWD"
    python $SCRIPT_PATH/generate_pseudo_scaff.py $path $name $num $len $cutoff
    for i in `seq 0 $((num - 1))`
    do
        echo "$PWD/"$name"_scaffold_"$i".vcf" >> ../vcf.list
    done
done < $config

cd ../

echo "simulation finished, generate simulated scaffold vcf file list at $out_dir/vcf.list, total pseudo scaffold number is $scaf_num"

if [ ! -d preprocess ]
then
    mkdir preprocess
fi

cd preprocess
scaf_num=0
echo "data format preprocessing..."
while read -r path name len num
do
    echo "processing $path"
    python $SCRIPT_PATH/preprocess.py $path $PWD/$name $num 100 $len $cutoff $gap_size
    scaf_num=$((scaf_num + num))
done < $config

cd ../
echo "preprocessing done"

# echo "generate sample list..."
# python $SCRIPT_PATH/get_sample_list.py `head -1 vcf.list` sample.list
sample_num=`wc -l < "$sample_list"`

if [ ! -d samples ]
then
    mkdir samples
fi

cd samples
echo "merge preprocessing data..."
while read sample
do
    while read -r path name len num
    do
        for i in `seq 0 $((num - 1))`
        do
            cat ../preprocess/$name/scaffold_"$i"/"$sample".snps.1 >> $sample.snps
            cat ../preprocess/$name/scaffold_"$i"/"$sample".snps.2 >> $sample.snps
        done
    done < $config
done < $sample_list

for f in *.snps
do
    echo $PWD/$f >> ../match.list
done

cd ../
echo "merge done"

echo "matching..., write result to match.result"
echo $sample_num
echo $scaf_num
$MATCHING match.list $sample_num $scaf_num > match.result

cd $CUR_DIR
echo "experiment done"
