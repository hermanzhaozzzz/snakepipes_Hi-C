# snakepipes_Hi-C
**须知**：本仓库还在构建中，暂时只作参考！！

本流程融合和更改了[HiC Pro](https://github.com/nservant/HiC-Pro)、[Juicer](https://github.com/aidenlab/juicer/)的代码

---
## 环境
```shell
# step1
pip install -r pip_env.txt
# step2
conda env create -f conda_env.yml
# or
mamba env create -f conda_env.yml
```

## 用法

搞清楚配置文件中的参数跑哪里去了

### step 0 测序质量控制
使用 [snakepipes_fastqc-multiqc](https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc)进行质量控制

### step 1 运行Snakemake Pipeline，生成Hi-C contact matrix
- **回贴Hi-C reads以及生成RAW矩阵ICE校正矩阵**
- **validPairs convert to .hic file(Juicer)**

```shell
cd HiC
git clone https://github.com/hermanzhaozzzz/snakepipes_fastqc-multiqc.git
cd snakepipes_fastqc-multiqc

# then

# use jupyterlab or runipy to run step01_generate_samples.ipynb
# get samples.json and check it

# then

# dry run, rm -n to run pipeline
snakemake -pr -j 8 -s step02_run_mapping_and_generate_matrix.py -n

# output as below
# HiC|⇒ tree . -L 1
# .
# ├── bam
# ├── fastq
# ├── hic_file
# ├── matrix
# ├── qc
# ├── quality_checks
# ├── snakepipes_fastqc-multiqc
# ├── snakepipes_Hi-C
# ├── temp_files
# └── valid_pairs
```
### step 2 Convert ValidPairs to Juicer .hic¶