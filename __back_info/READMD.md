# snakepipes_Hi-C
**须知**：本仓库还在构建中，暂时只作参考！！

---
**TOC**
- [snakepipes_Hi-C](#snakepipes_Hi-C)
    - [Settings](#Settings)
        - [envs](#envs)
        - [install Hi-C Pro](#install-hi-c-pro)
        - [install juicer](#install-juicer)
        - [install hicexplorer](#install-hicexplorer)
        - [enzyme cut map](#enzyme-cut-map)
        - [genome size](#genome-size)
        - [download test files](#download-test-files)
        - [file tree](#file-tree)
    - [Running](#running)
        - [[1] run mapping](#1-run-mapping)
        - [[2] call TAD](#2-call-tad)
        - [[3] call Loops](#3-call-loops)
        - [[4] call compartments](#4call-compartments)
    - [Helps](#helps)
        - [北京大学北极星slurm系统设置](#北京大学北极星slurm系统设置)
        - [config file template](#config-file-template-step01_config-hicprotxt)
---
## Settings
### envs
```shell
# https://github.com/nservant/HiC-Pro
wget https://raw.githubusercontent.com/nservant/HiC-Pro/master/environment.yml
mamba env create -f environment.yml
conda activate HiC-Pro_v3.1.0

# 2022-07-25
# bowtie 2.4.5
# samtools 1.15.1 Using htslib 1.15.1
# assert check_cmd("bowtie2")  # Bowtie 2 version 2.4.5
# assert check_cmd("samtools")  # samtools 1.15.1 Using htslib 1.15.1
# pip install iced                           
# Collecting iced
#   Downloading iced-0.5.10.tar.gz (2.3 MB)
# manually set cmd path

```

### install Hi-C Pro
```shell
tar -zxvf HiC-Pro-master.tar.gz
cd HiC-Pro-master
## Edit config-install.txt file if necessary
## do not install the app to this git repo file (itself)
make configure
make install

# add to PATH
```
### install juicer
```shell
git clone https://github.com/theaidenlab/juicer.git /lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/juicer
cd /lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/juicer/SLURM/scripts
wget https://hicfiles.tc4ga.com/public/juicer/juicer_tools.1.9.9_jcuda.0.8.jar
ln -s juicer_tools.1.9.9_jcuda.0.8.jar  juicer_tools.jar

# add to PATH
```
### install hicexplorer
```shell
mamba install hicexplorer
```

### enzyme cut map
```shell
python /lustre1/chengqiyi_pkuhpc/zhaohn/0.apps/HiC-Pro_installed/HiC-Pro_3.1.0/bin/utils/digest_genome.py \
    -r hindiii \
    -o hg38_hindiii.bed \
    /lustre1/chengqiyi_pkuhpc/zhaohn/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa
```
### genome size
```shell
samtools faidx /lustre1/chengqiyi_pkuhpc/zhaohn/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa

cut -f1-2 /lustre1/chengqiyi_pkuhpc/zhaohn/1.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.fai > chrom_hg38.sizes
```
### download test files
```shell
wget  https://zerkalo.curie.fr/partage/HiC-Pro/HiCPro_testdata.tar.gz && tar -zxvf HiCPro_testdata.tar.gz
```

### file tree
```text
./
├── fastq
│   ├── lib1
│   │   ├── SRR400264_00_R1.fastq.gz
│   │   └── SRR400264_00_R2.fastq.gz
│   └── lib2
│       ├── SRR400264_01_R1.fastq.gz
│       └── SRR400264_01_R2.fastq.gz
└── snakepipes_Hi-C
    └── ...

```
---
## Running
### [1] run mapping
```shell
########
# local
########
HiC-Pro -i ../fastq -o ../out_dir -c step01_config-hicpro.txt


########
# slurm
########
# generate sbatch 
HiC-Pro -i ../fastq -o ../out_dir -c step01_config-hicpro.txt -p
# run step1
cd ../out_dir/ && sbatch HiCPro_step1_zHiC.sh && cd -
tail -f ../out_dir/slurm-*
# run step2
cd ../out_dir/ && sbatch HiCPro_step2_zHiC.sh && cd -
tail -f ../out_dir/slurm-*
```
### [2] call TAD
HiC Explorer / TopDom overlap
重点参考了这个攻略：https://blog.csdn.net/hzau_yang/article/details/100031590

convert hicpro to h5

```shell
use step02_convert_abs_matrix_to_h5.ipynb

this notebook will convert Hi-C Pro matrix to hicexplorer h5 matrix

hicpro.matrix + hicpro.bed
↓
hicexplorer.h5
↓
hicexplorer.norm_range.h5
↓
hicexplorer.norm_range.KRcorrected.h5
```





3.校正Hi-C交互矩阵
默认使用KR标准化方法，也可以使用ICE，通过–correctionMethod参数控制；
–chromosomes控制输出的染色体，通过空格隔开，如chr1 chr2 chr3;

```shell
# KR 校正
hicCorrectMatrix correct -m hic_matrix.h5 --filterThreshold -1.5 5 -o hic_corrected.h5
# 画热图
hicPlotMatrix -m hic_corrected.h5 -o hic_plot.png --region 1:20000000-80000000 --log1p
```
该命令可以根据输出文件后缀画出不同的格式图，也可以通过–perChr参数控制每条染色体单独画；


找TAD
找TAD的软件有很多，hicexplorer的方法与Topdom有点相似，总的来说算比较简单粗暴的，实现
hicFindTADs -m hic_corrected.h5 --outPrefix hic_corrected --numberOfProcessors 16











### [3] call Loops
overlap
```shell
# Convert ValidPairs to Juicer .hic
~/0.apps/HiC-Pro_installed/HiC-Pro_3.1.0/bin/utils/hicpro2juicebox.sh -i dixon_2M.allValidPairs -g ~/0.apps/HiC-Pro_installed/HiC-Pro_3.1.0/annotation/chrom_hg38.sizes -j ~/0.apps/juicerbox/juicer_tools.jar -r ~/0.apps/HiC-Pro_installed/HiC-Pro_3.1.0/annotation/HindIII_resfrag_hg38.bed
```


call loop [juicer]
https://github.com/aidenlab/juicer

http://www.360doc.com/content/19/1224/14/68068867_881786243.shtml

juicer采用ArrowHead算法对原始的交互矩阵进行转化，并预测TAD拓扑关联结构域，采用HiCUUPS算法识别染色质环chromatin loops。和其他Hi-C数据处理软件相比，juicer的功能更为齐全

```shell
Command Line Tools Usage
Detailed documentation about the command line tools can be found on the wiki:

Annotating features with Arrowhead, HiCCUPS, MotifFinder, APA, Eigenvector, and Pearsons
Creating .hic with Pre
Extracting data from .hic files with dump
To launch the command line tools, use the shell script “juicer_tools” on Unix/MacOS or type

java -jar juicer_tools.jar (command...) [flags...] <parameters...>
There are different flavors of juicer_tools that depend on the CUDA version. If you do not use GPUs, these versions are equivalent. Otherwise, juicer_tools.X.X.jar uses CUDA version X.X

For HiCCUPS loop calling without the shell or bat script, you will need to call: java -Xms512m -Xmx2048m -Djava.library.path=path/to/natives/ -jar juicer_tools.jar hiccups [flags...] <parameters...> where path/to/natives is the path to the native libraries used for Jcuda By default, these are located in the lib/jcuda folder.

In the command line tools, there are several analysis functions:

apa for conducting aggregate peak analysis
hiccups for annotating loops
motifs for finding CTCF motifs
arrowhead for annotating contact domains
eigenvector for calculating the eigenvector (first PC) of the Pearson's
pearsons for calculating the Pearson's
The juicer_tools (Unix/MacOS) script can be used in place of the unwieldy java -Djava.library.path=path/to/natives/ -jar juicer_tools.jar
```
### [4] call compartments

5.找compartment


hicPCA -m hic_corrected.h5 --outFileName pca1.bw pca2.bw --format bigwig --pearsonMatrix pearson.h5

通常可以根据基因密度来调整第一主成分的符号，获得最终的compartment；统计完每个bin的基因数量，得到bigwig文件，然后通过–extraTrack可以直接对PCA结果调整符号；也可以通过–pearsonMatrix和–obsexpMatrix生成计算compartment的中间处理中的矩阵，如pearson矩阵；

画图

hicPlotMatrix -m pearson.h5 --outFileName pca1.png --perChr --bigwig pca1.bw

## Helps
### 北京大学北极星slurm系统设置
软件目录script下的make_slurm_script.sh需要修改
【直接复制】

```shell
#!/bin/bash

## HiC-Pro
## Copyright (c) 2015 Institut Curie                               
## Author(s): Nicolas Servant
## Contact: nicolas.servant@curie.fr
## This software is distributed without any guarantee under the terms of the BSD-3 licence.
## See the LICENCE file for details

##
## Create SLURM Torque files
##

dir=$(dirname $0)

usage()
{
    echo "usage: $0 -c CONFIG [-s STEP]"
}

MAKE_OPTS=""

while [ $# -gt 0 ]
do
    case "$1" in
    (-c) conf_file=$2; shift;;
	(-s) MAKE_OPTS=$2; shift;;
    (--) shift; break;;
    (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    (*)  suffix=$1; break;;
    esac
    shift
done

if [ -z "$conf_file" ]; then usage; exit 1; fi

CONF=$conf_file . $dir/hic.inc.sh
unset FASTQFILE

## Define input files
if [[ $MAKE_OPTS == "" || $MAKE_OPTS == *"mapping"* ]]
then
    inputfile=inputfiles_${JOB_NAME}.txt
    ifq=$(get_hic_files $RAW_DIR .fq)
    ifastq=$(get_hic_files $RAW_DIR .fastq)
    echo -e "$ifq\n$ifastq" | grep $PAIR1_EXT | sed -e "s|$RAW_DIR||" -e "s|^/||" > $inputfile
    count=$(cat $inputfile | wc -l)
elif [[ $MAKE_OPTS == *"proc_hic"* ]]
then
    inputfile=inputfiles_${JOB_NAME}.txt
    get_hic_files $RAW_DIR .bam | grep $PAIR1_EXT | sed -e "s|$RAW_DIR||" -e "s|^/||" > $inputfile
    count=$(cat $inputfile | wc -l)
fi

## Paralelle Implementation
if [[ $MAKE_OPTS == "" || $MAKE_OPTS == *"mapping"* || $MAKE_OPTS == *"proc_hic"* ]]
then
    make_target="all_sub"
    ## Remove per sample steps
    if [[ $MAKE_OPTS != "" ]]; then 
    make_target=$(echo $MAKE_OPTS | sed -e 's/,/ /g'); 
    make_target=$(echo $make_target | sed -e 's/merge_persample//g');
    make_target=$(echo $make_target | sed -e 's/build_contact_maps//g');
    make_target=$(echo $make_target | sed -e 's/ice_norm//g');
        make_target=$(echo $make_target | sed -e 's/quality_checks//g');
    fi
 
    ## step 1 - parallel
    torque_script=HiCPro_step1_${JOB_NAME}.sh
 
    cat > ${torque_script} <<EOF
#!/bin/bash
#SBATCH -N 1
#SBATCH -c ${N_CPU}
#SBATCH -p ${JOB_QUEUE}

#SBATCH --job-name=s1_${JOB_NAME}_HiCpro
#SBATCH --export=ALL
#SBATCH --no-requeue
#SBATCH -A ${JOB_ACCOUNT}
#SBATCH --qos=${JOB_QOS}
EOF
    
    if [[ $count -gt 1 ]]; then
	echo -e "#SBATCH --array=1-$count" >> ${torque_script}
    fi
    cat >> ${torque_script} <<EOF
FASTQFILE=\$SLURM_SUBMIT_DIR/$inputfile; export FASTQFILE
make --file ${SCRIPTS}/Makefile CONFIG_FILE=${conf_file} CONFIG_SYS=${INSTALL_PATH}/config-system.txt $make_target 2>&1
EOF
    
    chmod +x ${torque_script}

    ## User message
    echo "The following command will launch the parallel workflow through $count torque jobs:"
    echo sbatch ${torque_script}
fi    


## Per sample Implementation
if [[ $MAKE_OPTS == "" || $MAKE_OPTS == *"build_contact_maps"* || $MAKE_OPTS == *"ice_norm"* || $MAKE_OPTS == *"quality_checks"* ]]
then
    make_target="all_persample"
    ## Remove parallele mode
    if [[ $MAKE_OPTS != "" ]]; 
    then 
	make_target=$(echo $MAKE_OPTS | sed -e 's/,/ /g'); 
	make_target=$(echo $make_target | sed -e 's/mapping//g');
	make_target=$(echo $make_target | sed -e 's/proc_hic//g');
    fi

    torque_script_s2=HiCPro_step2_${JOB_NAME}.sh
    cat > ${torque_script_s2} <<EOF
#!/bin/bash

#SBATCH -N 1
#SBATCH -c 1
#SBATCH -p ${JOB_QUEUE}

#SBATCH --job-name=s2_${JOB_NAME}_HiCpro
#SBATCH --export=ALL
#SBATCH --no-requeue
#SBATCH -A ${JOB_ACCOUNT}
#SBATCH --qos=${JOB_QOS}

cd \$SLURM_SUBMIT_DIR

make --file ${SCRIPTS}/Makefile CONFIG_FILE=${conf_file} CONFIG_SYS=${INSTALL_PATH}/config-system.txt $make_target 2>&1
EOF
    
    chmod +x ${torque_script_s2}

    ## User message
    echo "The following command will merge the processed data and run the remaining steps per sample:"
    echo sbatch ${torque_script_s2}
fi
```
### config file template [./step01_config-hicpro.txt]
配置Hi-C Pro的模板
```shell
# Copy and edit the configuration file ‘config-hicpro.txt’ in your local folder. 
# The ‘[]’ options are optional and can be undefined. 带有[]中括号的参数可以不定义，为可选参数
# Please change the variable settings below if necessary
#########################################################################
## Paths and Settings  - Do not edit ! 输入输出文件路径，尽量不动
#########################################################################
TMP_DIR = tmp
LOGS_DIR = logs
BOWTIE2_OUTPUT_DIR = bowtie_results
MAPC_OUTPUT = hic_results
# Link to rawdata folder. The user usually not need to change this option
# 尽量不动rawdata路径
RAW_DIR = rawdata
#######################################################################
## SYSTEM AND SCHEDULER - Start Editing Here !! 从这里开始编辑
#######################################################################
# sort RAM
SORT_RAM = 1000M
# name of the main log file
LOGFILE = hicpro.log

# name of the job on the cluster【可选参数】
JOB_NAME = zHiC 

#######################################################################
## slurm squeue
#######################################################################
JOB_ACCOUNT = chengqiyi_g1
######################
# 【cnlong】
#-------------------->
N_CPU = 20
JOB_QUEUE = cn-long
JOB_QOS = chengqiyicnl
######################
######################
# 【cn-short】
#-------------------->
# N_CPU = 20
# JOB_QUEUE = cn-short
# JOB_QOS = chengqiyicns
######################
######################
# 【cn_nl】
#-------------------->
# N_CPU = 20
# JOB_QUEUE = cn_nl
# JOB_QOS = chengqiyicnnl
######################
######################
# 【fat4way】
#-------------------->
# N_CPU = 24
# JOB_QUEUE = fat4way
# JOB_QOS = chengqiyif4w
######################
######################
# 【fat8way】
#-------------------->
# N_CPU = 64
# JOB_QUEUE = fat8way
# JOB_QOS = chengqiyif8w
######################
#########################################################################
## Data
#########################################################################
# Keyword for first mate detection. Default:_R1
PAIR1_EXT = _R1
# Keywoard for seconde mate detection. Default:_R2
PAIR2_EXT = _R2
#######################################################################
## Alignment tools
#######################################################################
# Path to bowtie2 indexes
BOWTIE2_IDX_PATH = /lustre1/chengqiyi_pkuhpc/zhaohn/1.database/db_genomes/genome_fa/genome_ucsc_hg38/

# Reference genome prefix used for genome indexes. Default: hg19
REFERENCE_GENOME = genome_ucsc_hg38.fa.bowtie2_index
#######################################################################
## Alignment options
#######################################################################
# Minimum mapping quality. Reads with lower quality are discarded. Default: 0
MIN_MAPQ = 10

# bowtie2 options for mapping step1. 
# Default: –very-sensitive -L 30 –score-min L,-0.6,-0.2 –end-to-end –reorder
BOWTIE2_GLOBAL_OPTIONS = --very-sensitive -L 30 --score-min L,-0.6,-0.2 --end-to-end --reorder

# bowtie2 options for mapping step2. 
# Default: –very-sensitive -L 20 –score-min L,-0.6,-0.2 –end-to-end –reorder
BOWTIE2_LOCAL_OPTIONS =  --very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder
#######################################################################
## Annotation files
#######################################################################
# Chromsome size file. 
# Loaded from the ANNOTATION folder in the HiC-Pro installation directory. 放到安装文件夹下
# Default: chrom_hg19.sizes
GENOME_SIZE = chrom_hg38.sizes
#######################################################################
## Allele specific analysis
## http://nservant.github.io/HiC-Pro/AS.html#as
#######################################################################
# VCF file to SNPs which can be used to distinguish parental origin. See the allele specific section for more details
# 【可选参数】
# ALLELE_SPECIFIC_SNP = 
#######################################################################
## Capture Hi-C analysis
#######################################################################
# BED file of target regions to focus on (mainly used for capture Hi-C data)
# 【可选参数】
# CAPTURE_TARGET =
REPORT_CAPTURE_REPORTER = 1
#######################################################################
## Digestion Hi-C
#######################################################################
# BED file with restriction fragments. 
# Full path or name of file available in the ANNOTATION folder. 
# Default: HindIII_resfrag_hg19.bed
# 含有限制性片段的BED文件
# 放到安装文件夹下
GENOME_FRAGMENT = HindIII_resfrag_hg38.bed

# Ligation site sequence used for reads trimming. 
# Depends on the fill in strategy. 
# Example: AAGCTAGCTT
# 酶的序列，例如：
# HindIII，为AAGCTAGCTT
# MboI，为GATCGATC
LIGATION_SITE = AAGCTAGCTT

# Minimum size of restriction fragments to consider for the Hi-C processing.
# Example: 100
# 为Hi-C处理考虑的限制片段的【最小大小】
MIN_FRAG_SIZE = 100

# Maximum size of restriction fragments to consider for the Hi-C processing.
# Example: 100000
# # 为Hi-C处理考虑的限制片段的【最大大小】
MAX_FRAG_SIZE = 100000

# Minimum sequenced insert size. 
# Shorter 3C products are discarded. 
# Example: 100
# 测得【最小】插入大小。
# 较短的3C产物会被丢弃
MIN_INSERT_SIZE = 100

# Maximum sequenced insert size. 
# Larger 3C products are discarded. 
# Example: 600
# 测得【最大】插入大小。
# 较短的3C产物会被丢弃
MAX_INSERT_SIZE = 600

#######################################################################
## Hi-C processing
#######################################################################
# Filter short range contact below the specified distance. 
# Mainly useful for DNase Hi-C. 
# Example: 1000
MIN_CIS_DIST =

# Create output files with all classes of 3C products. 
# Default: 0
GET_ALL_INTERACTION_CLASSES = 1

# Create a BAM file with all aligned reads flagged according to 
# their classifaction and mapping category. 
# 是否保留BAM文件，默认不保留
# Default: 0
GET_PROCESS_SAM = 1

# Remove singleton reads. 
# Default: 1
RM_SINGLETON = 1

# Remove multi-mapped reads. 
# Default: 1
RM_MULTI = 1

# Remove duplicated reads’ pairs. 
# Default: 1
RM_DUP = 1

#######################################################################
## Contact Maps
#######################################################################
# Resolution of contact maps to generate (space separated). 
# Default: 20000 40000 150000 500000 1000000
BIN_SIZE = 5000 10000 20000 40000 100000 150000 500000 1000000

# Binning step size in ‘n’ coverage _i.e._ window step. 
# Default: 1
# BIN_STEP

# Output matrix format. 
# Must be complete, asis, upper or lower. 
# Default: upper
MATRIX_FORMAT = upper

#######################################################################
## Normalization
#######################################################################
# Maximum number of iteration for ICE normalization. 
# Default: 100
MAX_ITER = 100

# Define which pourcentage of bins with low counts should be force to zero. 
# Default: 0.02. 
# Replace SPARSE_FILTERING
FILTER_LOW_COUNT_PERC = 0.02

# Define which pourcentage of bins with low counts should be discarded 
# before normalization. 
# Default: 0
FILTER_HIGH_COUNT_PERC = 0

# The relative increment in the results before declaring convergence. 
# Default: 0.1
EPS = 0.1
```