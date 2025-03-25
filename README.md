## 实用全基因组分析百科全书
全基因组whole genome sequencing

## 1.	组装与评估
```bash
nohup bash -c 'for i in $(ls *.fq.gz | cut -d'_' -f1); do spades.py --isolate -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz -o ./spade/${i} -t 8 -m 300 -k 21,33,55,77; done' > spade.log &
```
##  宏基因组组装指令
```bash
nohup bash -c 'for i in $(ls *.clean.fq.gz | cut -d'_' -f1); do spades.py --meta -1 ${i}_1.clean.fq.gz -2 ${i}_2.clean.fq.gz -o ./spade/${i} -t 8 -m 600 -k 59,79,99,109,125 --careful; done' > spade.log & 
```
## 为文件夹下的双端测序的原始文件进行批量拼接，注意修改文件后缀
```bash
nohup bash -c 'ls | while read LINE; do spades.py -1 $LINE/$LINE\_1.fq.gz -2 $LINE/$LINE\_2.fq.gz -o ./spades/$LINE; done' > spades.log &
```
## 在上级文件夹下为下级文件夹下的contigs.fasta进行批量改名，并统一转移到上级文件夹下
```bash
ls | while read LINE;do cp $LINE/contigs.fasta ./$LINE.fasta;done
```
## quast评估组装的质量
```bash
quast.py *.fasta -o quast_output
```
## checkm评估组装的质量
## 安装
#mamba create -n checkM python numpy matplotlib pysam hmmer prodigal pplacer checkm-#genome
## 运行
conda activate checkM
```bash
nohup bash -c 'checkm lineage_wf -f ./checkmresult.tsv --tab_table -x fasta -t 8 --pplacer_threads 8 ./ .' > checkm.log &
```
## 阈值为5%，污染率超过5%的菌株应当舍弃
2.	基因注释
## 自建数据库查找噬菌体，接合质粒
## 对fasta文件可以直接运行以下命令，建议编写在同一个sh文件下统一运行
conda activate checkm
mlst *.fasta --threads 8 > ./pubmlst.tab
## mlst --scheme cronobacter *.fasta > ./mlst.tab
## mlst --minid 90 --scheme senterica_achtman_2 *.fasta > ./mlst.tab
abricate --threads 8 --db plasmidfinder *.fasta > ./plasmidfinder.tab
abricate --threads 8 --db BacMet2_EXP_database *.fasta > ./BacMet.tab
abricate --threads 8 --db card *.fasta > ./card.tab
abricate --threads 8 --db ncbi *.fasta > ./ncbi.tab
abricate --threads 8 --db vfdb *.fasta > ./vfdb.tab
abricate --threads 8 --db resfinder *.fasta > ./resfinder.tab
abricate --threads 8 --db ISfinder *.fasta > ./ISfinder.tab
abricate --threads 8 --db SPI *.fasta > ./SPI.tab
abricate --threads 8 --db ICE *.fasta > ./ICE.tab
abricate --threads 8 --db ice *.fasta > ./ice.tab
abricate --threads 8 --db oriT *.fasta > ./oriT.tab
abricate --threads 8 --db relaxase_t4cp *.fasta > ./relaxase_t4cp.tab
abricate --threads 8 --db prophage *.fasta > ./prophage.tab
abricate --threads 8 --db Tn *.fasta > ./Tn.tab
abricate --threads 8 --db mobileOG *.fasta > ./mobileOG.tab
abricate --threads 8 --db SGI-1 *.fasta > ./SGI-1.tab
seqkit stats -a *.fasta -j 8 -T > genome_stats.tsv
## 
## Listeria 专属数据库
abricate --threads 8 --db lis_vir *.fasta > ./lis_vir.tab
abricate --threads 8 --db lis_stress_island > ./lis_stress_island.tab
abricate --threads 8 --db lm_resistance > ./lm_resistance.tab
abricate --threads 8 --db listeria_Metal_Disinfectants_Resistance > ./lm_MDR.tab
abricate --threads 8 --db listeria_antibiotic_resistance > ./lm_AMR.tab
## 
nohup bash mlst.sh > mlst.log & 
## Transposon database: https://tncentral.ncc.unesp.br/
## Prophage database: https://phastest.ca/databases
## IS database: https://www-is.biotoul.fr/
nohup bash -c ' abricate --threads 8 --db Tn *.fasta > ./Tn.tab' > Tn.log &
## 更新pubmlst数据库
nohup bash -c 'mlst-download_pub_mlst -j 8 -d /home/student/anaconda3/envs/checkm/db/pubmlst' > download_pubmlst.log &
## iceberg 数据库介绍: Integrative and conjugative element (ICE) is bacterial mobile genetic elements (MGEs), which is integrative to the bacterial chromosome and encodes a fully functioning conjugation machinery and is thus self-transmissible between bacterial cells, chromosome-borne integrative and mobilizable elements (IMEs), cis-mobilizable elements (CIMEs), plasmids. T4SSO(Ontology of Type IV Secretion Systems) is aimed at the systematic and logic representation of Type IV Secretion Systems.
## 主流数据库统计## 
#ncbi 7010条数据，card 4805条，resfinder 3192条，MEGARes（耐药加消毒剂）6635条，#vfdb 4392条，消毒剂重金属753条，插入序列5970条，转座子加插入序列6073条
## 使用prokka进行批量注释## 
## 直接使用nohup一句话写完
nohup bash -c 'for i in *.fasta; do prokka --outdir ./prokka/${i%.fasta} --cpus 8 --prefix ${i%.fasta} ${i} --addgenes --addgenes --centre X --compliant;done' > prokka.log &
## 
## abricate 自建数据库
awk '/^>/{split($0,a,">"); print "> Tn~~~" a[2] "~~~" a[2] "~~~" a[2]} !/^>/{print}' Tn.fa > modified_Tn.fasta
awk '/^>/{split($0,a,">"); print "> mobileOG~~~" a[2] "~~~" a[2] "~~~" a[2]} !/^>/{print}' mobileOG.fasta > modified_mobileOG.fasta
head -n 1 sequences
abricate --setupdb --threads 16
abricate –list
abricate-get_db --db megares --force
vim prokka.sh   ## 编写prokka批处理文件
#按i进入插入模式
#按esc 输入:wq并回车，保存并退出文本编辑模式
cat prokka.sh # 查看文件是否写入
bash prokka.sh ## 直接运行批处理命令
#程序放服务器的后台运行
nohup bash prokka.sh > prokka.log & 
#查看正在运行的工作，查看后台有无程序运行## 
jobs
#强制停止任务运行
kill -9 %1
## 使用rgi对耐药基因进行注释，使用card数据库
nohup bash -c 'for i in *.fasta; do rgi main -a DIAMOND --input_sequence ${i} --local --clean -o card/${i%.fasta}; done' >> card.log &
## 批量生成gff3文件
nohup bash -c 'for i in *.fasta; do dante -q ${i} -o gff3/${i%.fasta}.gff3 -c 8; done' > dante.log &
## 使用bakta进行文件注释（运行较慢，注释效果优于prokka和prodigal）
conda activate bakta_env
bakta_db list
#bakta_db download --output /home/student/metagenome/bakta --type full
#bakta_db download --output /data/liushiwei/bakta_db --type full
nohup bash -c 'for i in *.fasta; do bakta --db /home/student/metagenome/bakta/db -o ./bakta/${i%.fasta} -t 8 -p ${i%.fasta} ${i};done' > bakta.log &
nohup bash -c 'for i in *.fasta; do bakta --db /home/student/metagenome/bakta/db -o ./bakta/${i%.fasta} -t 8 -p ${i%.fasta} ${i};done' > bakta.log &
nohup bash -c 'for i in *.fasta; do bakta --db /data/liushiwei/bakta_db/db -o ./bakta/${i%.fasta} -t 8 -p ${i%.fasta} ${i};done' > bakta.log &
amrfinder_update --force_update --database /data/liushiwei/bakta_db/db/amrfinderplus-db
conda install -y -c conda-forge -c bioconda --strict-channel-priority ncbi-amrfinderplus=4.0.19
（需要指定版本号，否则无法正常安装）
bakta 支持宏基因组注释，添加--meta参数
## 使用amrfinder注释耐药基因
conda activate bakta_env
mamba update ncbi-amrfinderplus
amrfinder -u -d /home/student/metagenome/bakta/db/amrfinderplus-db
amrfinder -n *.fasta -o amrfinder.tab --mutation_all --threads 8
amrfinder -g *.gff
amrfinder -v

3.	质粒
#从gff文件里直接截取一段contigs，作为fasta文件，将对应的fasta文件用gbk或者gff文件注释，最后再在easyfig或者是clinker上面进行可视化分析
这样做可视化分析会得心应手
## 使用genomad数据库预测可移动元件
genomad end-to-end --cleanup --splits 8 GCF_009025895.1.fna.gz genomad_output /data/liushiwei/genomad_db
nohup bash -c 'for i in *.fasta; do genomad end-to-end --cleanup --splits 8 ${i} genomad_output /data/liushiwei/genomad_db; done' > genomad.log &
nohup bash -c 'for i in *.fasta; do genomad end-to-end --cleanup --splits 8 ${i} genomad_output /home/student/genomad_db/genomad_db; done' > genomad.log &

## MGEfinder
for i in *.fasta; do mefinder find -t 8 --contig ./${i} ./mefinder/${i%.fasta}; done
## #质粒染色体分离
conda activate plasflow
for i in *.fasta; do PlasFlow.py --input ${i} --output ./plasflow/${i%.fasta} --threshold 0.7; done
conda activate mob_suite
for i in *.fasta; do mob_recon -i ${i} -o ./mob_suite/${i%.fasta}; done
## 使用mobile_OG对可移动元件的类型进行预测


## 
4.	泛基因组
conda activate env_roary
nohup bash -c 'roary -e --mafft -p 64 -g 10000000 -r *.gff' > roary.log &
python3 roary_plots.py core_SNP_tree.tre gene_presence_absence.csv
## coinfinder
参考资料：https://zhuanlan.zhihu.com/p/620558098
coinfinder -i gene_presence_absence.csv -I -p core_SNP_tree.tre -o coinfinder -x 8 -a -d -m
nohup bash -c 'query_pan_genome -a intersection *.gff' > pan_genome.log &
scoary -g gene_presence_absence.csv -t traits.csv -n 472_core_SNP_tree_collapsed.nwk
## 示例
scoary -t Tetracycline_resistance.csv -g Gene_presence_absence.csv -u -c I EPW
panaroo -i *.gff -o results --clean-mode strict --remove-invalid-genes
## panaroo与roary本身并不兼容，需要单独去建立panaroo的环境
## 泛基因组背景资料介绍
http://sanger-pathogens.github.io/Roary/
## 
5.	建树、SNP比较
conda activate SNP
snippy-clean_full_aln core_gene_alignment.aln > clean.full.aln
run_gubbins.py -c 32 -p gubbins clean.full.aln
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
FastTree -gtr -nt clean.core.aln > core_SNP_tree.tre
psdm -l -t 8 -o psdm_snp_clean.tab -P -d "\t" clean.core.aln
snp-dists -j 8 clean.core.aln > snp.tab
## 
#输出文件core_SNP_tree.tre为核心基因组的树文件
##  标准建树方法
run_gubbins.py -c 32 -p gubbins clean.full.aln
## 指定建树方法加快建树的进程
run_gubbins.py -c 64 --prefix enteritidis --tree-builder iqtree --first-model JC --tree-builder raxmlng --model GTR clean.full.aln
## 利用gubbins进行祖先重建
run_gubbins.py -c 24 --best-model --recon-with-dates --date 441_date_1.csv -p enteritidis clean.full.aln

#snp.tab 为两两比较的snp文件
## 
6.	血清型鉴定以及cgmlst分型
#沙门血清型鉴定
conda activate seqsero
nohup bash -c 'for i in *.fasta; do SeqSero2_package.py -m k -t 4 -i ${i} -p 8; done' > seqsero.log &
## 
#批量提取生成文件中的沙门的血清型信息
for file in SeqSero_result*/SeqSero_result.tsv; do awk 'NR==2' "$file" >> merged_second_lines.txt; done
## cgmlst_salmonella分析，cgMLSTschema99 用于做grapetree ## 
conda activate chewie
chewBBACA.py AlleleCall -i ./ -g /home/student/anaconda3/envs/chewie/db/salmonella/Salmonella_enterica_INNUENDO_cgMLST -o ./AlleleCall --cpu 8 --mode 1

chewBBACA.py ExtractCgMLST -i ./AlleleCall/results_alleles.tsv -o  ./ExtractCgMLST

chewBBACA.py AlleleCallEvaluator -i ./AlleleCall -g /home/student/anaconda3/envs/chewie/db/salmonella/Salmonella_enterica_INNUENDO_cgMLST -o ./AlleleCallEvaluator --cpu 8
## 
chewBBACA.py AlleleCall -i ./ -g /home/student/anaconda3/envs/chewie/db/listeria/Listeria_monocytogenes_Pasteur_cgMLST -o ./AlleleCall --cpu 8 --mode 1
chewBBACA.py ExtractCgMLST -i ./AlleleCall/results_alleles.tsv -o  ./ExtractCgMLST
chewBBACA.py AlleleCallEvaluator -i ./AlleleCall -g /home/student/anaconda3/envs/chewie/db/listeria/Listeria_monocytogenes_Pasteur_cgMLST -o ./AlleleCallEvaluator --cpu 8
##  利用已知数据库文件进行cgmlst建库
chewBBACA.py PrepExternalSchema -g /home/student/anaconda3/envs/chewie/db/cronobacter_alleles -o /home/student/anaconda3/envs/chewie/db/crono --cpu 18
##  cgmlst_
chewBBACA.py AlleleCall -i ./ -g /home/student/anaconda3/envs/chewie/db/crono  -o ../filtration_AlleleCall1 --cpu 8 --mode 1
chewBBACA.py ExtractCgMLST -i ../filtration_AlleleCall1/results_alleles.tsv -o  ../ExtractCgMLST
chewBBACA.py AlleleCallEvaluator -i ../filtration_AlleleCall1 -g /home/student/anaconda3/envs/chewie/db/crono  -o ../AlleleCallEvaluator --cpu 8
## 
chewBBACA.py AlleleCall -i ./ -g /home/student/anaconda3/envs/chewie/db/crono/crono  -o ../filtration_AlleleCall1 --cpu 18 --mode 1
chewBBACA.py ExtractCgMLST -i ../filtration_AlleleCall1/results_alleles.tsv -o  ../ExtractCgMLST
## Linux 隐藏文件/文件夹
在文件名前+.即可将文件/文件夹隐藏，需要对文件进行重命名才能够重新显示出来
##  计算GC含量
seqkit fx2tab -l -g -n -i -H *.fasta -j 8 > gc.tab
7.	去重、去冗余
## 使用dereplicator去除重复的样本
## （需要提前安装mash）速度和效率远远高于cd-hit
tr -d '\r' < file2.txt | xargs -I {} mv -- "{}" /data/liushiwei/salmonella/Global_Chicken_enteritidis/ST11
## dereplicator 的distance不应过大，否则大量的序列都会被删除
## 4000株筛选至165株，阈值建议设定为0.0001以下，阈值尽量设置小一点
nohup bash -c 'dereplicator.py ./ dereplicator_0.0001 --distance 0.0001 --thread 32' > dereplicator_0.0001.log &
## 使用cd-hit去冗余样本
conda activate SNP
nohup bash -c 'for i in *.fasta; do cd-hit -i ${i} -o cd-hit/${i} -aS 0.9 -c 0.95 -G 0 -g 0 -T 0 -M 0; done' > cd-hit.log &
8.	点突变
## pointfinder点突变
conda activate checkm
vim pointfinder.sh
for i in *.fasta; do python3 /home/student/pointfinder/PointFinder.py -i \./${i} -o \./pointfinder/ -p \/home/student/pointfinder/pointfinder_db -s salmonella -m blastn -m_p \/home/student/ncbi-blast-2.15.0+/bin/blastn; done
nohup bash -c 'for i in *.fasta; do python3 /data/liushiwei/pointfinder/PointFinder.py -i ${i} -o ./pointfinder -p /data/liushiwei/pointfinder_db -s salmonella -m blastn -m_p /data/liushiwei/ncbi-blast-2.16.0+/bin/blastn; done' >>  pointfinder.log &
for i in *.fasta; do python3 /data/liushiwei/pointfinder/PointFinder.py -i ${i} -o ./pointfinder -p /data/liushiwei/pointfinder_db -s campylobacter -m blastn -m_p /data/liushiwei/ncbi-blast-2.16.0+/bin/blastn; done
nohup bash pointfinder.sh >> pointfinder.log &
## pointfinder点突变结果汇总
for file in *.tsv; do awk 'NR > 1 {print FILENAME "\t" $0}' "$file" >> pointfinder.tab; done
## 注意的问题
pointfinder 无法识别包含有_下划线以及.的文件名，运行之前_.均需要去掉
9.	GO/KEGG/COG分析
## 使用prodigal进行基因注释和翻译
conda activate eggnog
nohup bash -c ' for i in *.fasta; do prodigal -i ${i} -a prodigal/${i%.fasta}.faa -o gff/${i%.fasta}.gff -d nucleotide/${i%.fasta}.fa -f gff; done' > prodigal.log &
##  eggnog基因注释，需要运行蛋白质文件，运行emapper，18m，默认diamond 1e-3; 2M,32p,1.5h, 该步骤耗时较长
cd prodigal
mkdir eggnog
nohup bash -c 'for i in *.faa; do emapper.py --data_dir /home/student/metagenome/eggnog -i ${i} --cpu 8 -m diamond --override -o eggnog/${i%.faa}; done' > eggnog.log &
##  格式化结果并显示表头
grep -v '^## ' *.emapper.annotations | sed '1 s/^#//' > output
csvtk -t headers -v output
10.	其他常见指令
## #多序列比对，需要gbk文件,由prokka生成
conda activate SNP
clinker clusters/*.gbk
clinker clusters/*.gbk -p <optional: file name to save static HTML>
clinker files/*.gbk -p plot.html
clinker -h
clinker 3-2B.gbk 7-2B.gbk -p -i 0.5 -o alignments.tab
仅适合较小的较短的序列片段
nohup bash -c 'clinker *gbk -p -i 0.5 -o alignments.tab' > dante.log &
11.	其他指令补充
##  antiSMASH（antibiotics & Secondary Metabolite Analysis Shell）识别和分析微生物中生物合成基因簇（BGCs）的工具
nohup bash -c 'for i in *.gff; do antismash ${i} --output-dir antismash --asf --pfam2go --smcog-trees --fullhmmer --output-basename ${i%.fasta}; done' > antismash.log &
## slurm系列操作指令
## 查看当前运行情况，一般不用top来查看
squeue
## 删除当前账户下所有运行的指令
scancel -u student
## 作业提交指令
sbatch *.slurm
## 其他生信代码补充
nohup bash -c 'mlst *.fasta > ./pubmlst.tab' &
## 压缩，解压指令
#压缩
tar -czvf sichuan_salmon.tar.gz *.fasta
#解压tar.gz文件
for i in *.tar.gz; do tar -zxvf ${i}; done
tar -zxvf *.tar.gz
## 解压gz文件
gunzip *.gz
##  参考资料
https://www.ncbi.nlm.nih.gov/pathogens/docs/datasets_assemblies/
12.	Download genomes from NCBI database
## 从https://www.ncbi.nlm.nih.gov/pathogens/ 下载accession号以及对应的信息表
##  datasets version: 16.43.0
## 依据NCBI给的accession号进行下载
##  The NCBI Datasets command-line tools (CLI) （默认指令，一般弃用，遇到卡顿会退出）
nohup bash -c 'datasets download genome accession --inputfile enteritidis_accessions.txt --api-key 1ef429d37c5d6103dac9cdaec0f54728d009' > ncbi.log &
## 下载速度虽然慢，但是比较稳定## 
## 编写批处理文件，避免出现单条命令出现报错会直接闪退的现象## 
## 单线程简易版本--------------------------------------------------------------------------------------------
#!/bin/bash
input_file="enteritidis_accessions.txt"
# 检查文件是否存在
if [ ! -f "$input_file" ]; then
    echo "错误: 文件 $input_file 未找到！" >&2
    exit 1
fi
# 逐行下载
while IFS= read -r accession; do
    echo "正在处理: $accession"
    datasets download genome accession "$accession" --filename "${accession}.zip" --api-key 1ef429d37c5d6103dac9cdaec0f54728d009 || echo "下载 $accession 失败"
done < "$input_file"
echo "全部任务完成。"
## 
# 检查必要依赖
check_dependency() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "错误: 未找到 $1 命令。请安装: $2" >&2
        exit 1
    fi
}
check_dependency "datasets" "https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/"
check_dependency "parallel" "apt-get install parallel 或 brew install parallel"
# 检查输入文件
if [ ! -f "$INPUT_FILE" ]; then
    echo "错误: 输入文件 $INPUT_FILE 未找到!" >&2
    exit 1
fi
## 多线程高级版本（最终版本）-----------------------------------------------------------------------------
#!/bin/bash
# NCBI基因组下载与校验脚本（支持并行和API密钥）
# ========== 用户配置 ========== 
API_KEY="1ef429d37c5d6103dac9cdaec0f54728d009"  # NCBI API密钥
INPUT_FILE="enteritidis_accessions.txt"              # 输入文件名
DOWNLOAD_DIR="downloaded_genomes"           # 下载目录
MAX_PARALLEL=8                              # 最大并行任务数
LOG_FILE="download.log"                         # 主日志文件
ERROR_LOG="error.log"                           # 错误日志
# =============================
# 创建下载目录
mkdir -p "$DOWNLOAD_DIR"
# 初始化日志
exec > >(tee -a "$LOG_FILE") 2> >(tee -a "$ERROR_LOG" >&2)
# 文件完整性检查函数
check_integrity() {
    local file="$1"
    if ! unzip -tq "$file" >/dev/null 2>&1; then
        echo "文件损坏: $(basename $file .zip)" >&2
        rm -f "$file"
        return 1
    fi
    return 0
}
# 收集需要下载的Accession
echo "[$(date +%T)] 开始文件核对..."
declare -a missing_accessions
while IFS= read -r accession; do
    accession=$(echo "$accession" | tr -d '\r' | xargs)
    [[ -z "$accession" ]] && continue    
    target_file="${DOWNLOAD_DIR}/${accession}.zip"    
    if [ -f "$target_file" ]; then
        if ! check_integrity "$target_file"; then
            missing_accessions+=("$accession")
        fi
    else
        missing_accessions+=("$accession")
    fi
done < "$INPUT_FILE"
# 退出条件检查
if [ ${#missing_accessions[@]} -eq 0 ]; then
    echo "[$(date +%T)] 所有文件均已完整存在"
    exit 0
fi
# 并行下载函数
parallel_download() {
    local acc="$1"
    echo "[开始下载] $acc"
# 下载命令（含API密钥）
    if datasets download genome accession "$acc" \
        --api-key "$API_KEY" \
        --filename "${DOWNLOAD_DIR}/${acc}.zip" 2>> "$ERROR_LOG"
    then
# 下载后验证
        if check_integrity "${DOWNLOAD_DIR}/${acc}.zip"; then
            echo "[下载成功] $acc"
            return 0
        fi
    fi
    echo "[下载失败] $acc" >&2
    return 1
}
# 导出函数和环境变量供parallel使用
export -f parallel_download check_integrity
export API_KEY DOWNLOAD_DIR ERROR_LOG
# 执行并行下载
echo "[$(date +%T)] 开始并行下载 ${#missing_accessions[@]} 个文件..."
printf "%s\n" "${missing_accessions[@]}" | parallel -j $MAX_PARALLEL \
    --progress --bar --eta \
    --joblog "${DOWNLOAD_DIR}/parallel.log" \
    --resume-failed \
    --tagstring "ACC:{}" \
    'parallel_download {}'
# 最终状态报告
success_count=$(grep -c "下载成功" "$LOG_FILE")
fail_count=$(grep -c "下载失败" "$ERROR_LOG")
echo "==============================="
echo "[最终报告] 下载完成时间: $(date)"
echo "成功: $success_count 个"
echo "失败: $fail_count 个"
echo "日志文件: $LOG_FILE"
echo "错误日志: $ERROR_LOG"
echo "并行日志: ${DOWNLOAD_DIR}/parallel.log"

## 下载后的文件立即检查ST型，剔除离群值，做一下质控，并进一步进行去重的工作
##  Simply just give it a genome file in FASTA/GenBank/EMBL format, optionally compressed with gzip, zip or bzip2. 无需解压直接运行
## 从NCBI数据库上下载fasta文件
unzip *.zip
for i in *.zip; do unzip ${i}; done
## 保留文件的前15个字符（accession ID）
for file in *.fna; do mv "$file" "${file:0:15}.fna"; done
## for file in *.gz; do mv "$file" "${file:3:50}"; done
##  # Perl语言版本批量改名
rename 's/\.all.fna/\.fasta/' *
rename 's/\.fna/\.fasta/' *
rename 's/GCA/GCA_/' *
rename 's/GCA/GCA_/' *
rename 's/\.seq/\.fasta/' *
rename 's/21L/21L-/' *
rename 's/S21./S21_/' *
rename 's/.seq/.fasta/' *
rename 's/.gff/-T.gff/' *
##  rename (util-linux 2.23.2) C语言版本
#目录中有file1.txt、file2.txt、file3.txt文件，要将所有文件名中的"file"替换为"doc"。
rename -v 'file' 'doc' *.txt
rename -v 'md5' 'MD5' *.txt

#整理并获取biosample唯一的序列号，整合到download.txt下
conda activate ncbi_datasets
nohup bash -c 'iseq -i download.txt -p 8 -g -d sra' > download.log &（不建议会卡命令，跑循环）
nohup bash -c 'cat download.txt | while read Run; do iseq -i $Run -a -g; done' > download.log &
## conda 指令增加下载渠道
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --show channels
## 单独运行某个py的方法，将py程序复制粘贴到对应环境的bin中
## 以checkm为例
## 再添加权限后即可顺利运行
chmod +x /data/liushiwei/anaconda3/envs/checkm/bin/dereplicator.py
chmod +x /data/liushiwei/anaconda3/envs/env_roary/bin/roary_plots.py
## ## #
基本指令，比如rename，git 需要在base环境下运行
## 在服务器上使用ollama进行deepseek本地化运行
ollama run DeepSeek-R1-Distill-Qwen-32B-GGUF:latest
ollama list
## ## #
13.	三代分析
#Genome assembly 三代nanopore数据分析
#Using software spades (http://cab.spbu.ru/software/spades/)
#Using software Unicycler (https://github.com/rrwick/Unicycler)
#Using software Flye (https://github.com/fenderglass/Flye)
#Using software Pilon (https://github.com/broadinstitute/pilon/)
spades.py -k 21,33,55,77,99,127 --careful --pe1-1 short_reads_1.fastq.gz --pe1-2 short_reads_2.fastq.gz -o /output/dir/ --phred-offset 33
unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -l long_reads.fastq.gz -o /output/dir/
canu -p output -d out_dir -fast genomeSize=5m -nanopore-raw reads.fastq.gz
flye -o out_dir --genome-size 5m --threads 16 --nano-raw reads.fastq.gz
## 三代数据分析
bwa index ./contigs.fa
bwa mem -t 12 ./contigs.fa ./short_reads_1.fastq.gz short_reads_2.fastq.gz | samtools view - -Sb | samtools sort - -@14 -o ./mapping.sorted.bam
samtools index mapping.sorted.bam
java -Xmx16G -jar ~/software/pilon-1.23.jar --genome ./contigs.fa --fix all --changes --frags mapping.sorted.bam --threads 12 --output ./pilon/pilon_round1 | tee ./pilon/round1.pilon
#Serotyping
#Using software SISTR (https://github.com/peterk87/sistr_cmd/)
sistr -i contigs.fasta contigs.fasta -f csv -o /output/rusults -p CGMLST_PROFILES -n NOVEL_ALLELES
cat *.csv >> ./totalsero.csv

14.	比较有用的R代码
vim Scale_zero_length_branches.R
args <- commandArgs(trailingOnly = TRUE)
# 输入文件路径
input_file <- args[1] 
# 输出文件路径
output_file <- paste0(tools::file_path_sans_ext(input_file), "_collapsed.nwk") 
##  Scale the branch length to a minimum value of 0.0001
library(ape)
tree <- read.tree(input_file)
##  判断  tree 中是否有 0 枝长 存在, 若有零枝长，即缩放至最小枝长为0.0001
condition_str <- '0' %in% tree$edge.length 
if ( condition_str == TRUE ) {
  tree <- compute.brlen(tree, minlength = 0.0001)
}
##  collapsing zero length branches
#library(phytools)
#tree <- read.newick(input_file)
#tree <- di2multi.simmap(tree, tol = 0.0001 ) 
#tree <- di2multi(tree, tol = 0.0001 ) 
# 将折叠后的树保存为Newick文件
write.tree(tree, file = output_file)
Rscript Scale_zero_length_branches.R core_SNP_tree.tre

15.	分子钟
##  Estimation of time scaled phylogenies (核心语句，用于筛选离群值)
treetime --tree 441_core_SNP_tree.tre --dates 441_date.csv --aln 441_clean.core.aln --outdir salmonella --coalescent skyline
## 方案二，使用全长aln或者全长进行建树，利用全长去做treetime

其他输出文件包括与重建祖先序列的比对、带注释的 nexus 格式树，其中分支长度对应于年份和突变，并且节点日期（采用上面详述的数字格式）作为注释添加到每个节点。此外，绘制从根到尖与时间的回归和树并保存到文件中。
##  Ancestral sequence reconstruction （祖先序列重建）
treetime ancestral --aln 441_clean.core.aln --tree 441_core_SNP_tree.tre --outdir ancestral_results

treetime homoplasy --aln 441_clean.core.aln --tree 441_core_SNP_tree.tre
重建祖先序列并将突变映射到树中。输出由一个包含祖先序列的文件“ancestral.fasta”和一个树“annotated_tree.nexus”组成，其中添加了突变作为注释，例如 A45G、G136T、…，SNP 中的数字默认使用基于 1 的索引。推断出的 GTR 模型写入标准输出。
treetime ancestral [-h] --aln ALN [--vcf-reference VCF_REFERENCE]
                   [--tree TREE] [--rng-seed RNG_SEED] [--gtr GTR]
                   [--gtr-params GTR_PARAMS [GTR_PARAMS ...]] [--aa]
                   [--custom-gtr CUSTOM_GTR] [--marginal] [--keep-overhangs]
                   [--zero-based] [--reconstruct-tip-states]
                   [--report-ambiguous]
                   [--method-anc {parsimony,fitch,probabilistic,ml}]
                   [--verbose VERBOSE] [--outdir OUTDIR]
##  Estimation of evolutionary rates and tree rerooting
treetime clock --tree 441_core_SNP_tree.tre --dates 441_date.csv --outdir clock_results --aln 441_clean.core.aln
 --sequence-len 1400 
计算从根到尖的回归并量化树的“时钟似然性”。除非使用 –keep-root 运行，否则它将重新根植树以最大化时钟似信号并重新计算分支长度
treetime clock [-h] --tree TREE [--rng-seed RNG_SEED] [--dates DATES]
               [--name-column NAME_COLUMN] [--date-column DATE_COLUMN]
               [--sequence-length SEQUENCE_LENGTH] [--aln ALN]
               [--vcf-reference VCF_REFERENCE] [--clock-filter CLOCK_FILTER]
               [--clock-filter-method {residual,local}]
               [--reroot REROOT [REROOT ...] | --keep-root]
               [--tip-slack TIP_SLACK] [--covariation] [--allow-negative-rate]
               [--plot-rtt PLOT_RTT] [--verbose VERBOSE] [--outdir OUTDIR]

##  Inference of transition between discrete characters and ‘mugration’ models
treetime mugration --tree *.nwk --states *.csv --attribute country
离散字符与“迁移”模型之间的转换推断
此命令将生成带注释的 nexus 树，其中属性的状态作为注释添加到每个节点（例如[&country="brazil"]）。此外，推断出的不同状态之间的 GTR 模型被写入文件。

##  Analyzing homoplasies and recurrent mutations
treetime homoplasy --aln 441_clean.core.aln --tree 441_core_SNP_tree.tre
## 同源性
重建祖先序列并将突变映射到树上。然后扫描树以查找同源性。过多的同源性可能表明存在污染、重组、文化适应或类似情况。
treetime homoplasy [-h] --aln ALN [--vcf-reference VCF_REFERENCE]
                   [--tree TREE] [--rng-seed RNG_SEED] [--const CONST]
                   [--rescale RESCALE] [--detailed] [--gtr GTR]
                   [--gtr-params GTR_PARAMS [GTR_PARAMS ...]] [--aa]
                   [--custom-gtr CUSTOM_GTR] [--zero-based] [-n N]
                   [--drms DRMS] [--verbose VERBOSE] [--outdir OUTDIR]
##  arg
计算从根到尖的回归并量化树的“时钟似性”。除非使用 –keep_root 运行，否则它将重新根植树以最大化时钟似信号并重新计算分支长度。

分子钟分析总流程
核苷酸/氨基酸序列数据(ClustalW, Muscle, MAFFT)
多序列比对对齐
寻找最佳进化模型(ModelFinder)
Beauti配置参数(生成XML文件)
(I)	分子钟校准方式:节点校准(Taxa，文献报道或者化石记录)，端点校准(Tps，使用采样时间校准)
(II)	配置进化模型参数:核苷酸氨基酸替换模型位点变异速率型
(III)	选择分子钟类型:严格or宽松选择树先验模型(种群增长模型 or 物种增长模型)
(IV)	设置参数先验(本质为取值范围)
(V)	配置MCMC参数:运行总代数，样本容量等 
Beast运行XML文件（估算节点分化时间和平均替换速率）
收敛诊断直看(Tracer，ESS>200)
生成带分化时间的最大可信进化分枝树(TreeAnnotator)（FigTree查看树文件）
BEAUti术语简介
partitions 序列数据(支持数据分区)
Taxa tips 分子钟校准方法
Traits 性状数据(比如地理状态)
Sites 进化模型
Clocks 分子钟类型
Trees 树模型
Priors 先验
MCMC参数
建议总样本数量大于等于1万个T)官网推荐的处理参数不收敛4)(ESS<200)的策略:增加总代数;②或合并多个独立运行链的结果;或增大采样频率，即减少Log parameters every值。
总样本数量 = Length of chain (总代数)÷Log parameters every (样本容量)
建议勾选Create tree log file with branch length in substitutions
treetime软件做时间信号检(需要准备三个文件)
提取采样时间做时间信号检测
严格来说，要用脚本做日期随机化检验(date-randomization test)，但这个步骤十分之繁琐，故通常我们可选择用Treetime之类的去做时间信号检测，并剔除掉离群点
Sites(位点替换模型)使用phylosuite软件中的ModelFinder进行模型选择，选择生成的log文件中BIC最佳得分模型即可
clocks和Trees模型选择(对于病毒序列一般直接选择uncorrelated宽松分子钟模型和Bayesian skyline贝叶斯天际线树模型)，Bayesian Skygrid模型需要填写Time at least transition point
参数
Priors中的参数不建议修改，勾选Use classic priors/operators参数即可
nohup beast -beagle_GPU *.xml&
BEAST软件运行完成使用Tracer软件进行收敛诊断（将上一步生成的log文件直接拖入Tracer软件中）如果有参数ESS值小于200没有收敛，可考虑将多次运行结果拖入Tracer软件进行Combined,或者调整运行总代数重新运行xml文件
使用TreeAnnotator生成MCC树(上步生成的*.(time).trees.txt拖入)
生成的tree文件在FigTree中添加时间标尺
Time Scale ---Scale by facter ---Offset by (填写最近采样日期)
取消勾选Scale Bar ---勾选Scale Axis ---勾选Reverse axis
勾选Node Bars 添加置信区间
Chiplot的node ages 和length axis选项也可以进行带分歧时间树的美化
生成贝叶斯天际线种群增长模型图
上传完成过返回第一项，点击Create partition from trait…
其余选项参数配置同上，得到的文件，依次将log文件导入Tracer软件查看是否收敛
将trees.txt文件导入TreeAnnotator软件合并生成一棵树
在SpreaD3软件中加载合并后的树，并选择location 参数
MCC tree with discrete traits
加载经纬度信息（纬度在前，经度在后）
设置最近采样时间及导入json格式地图文件，保存到输出文件
使用Rendering模式渲染输出文件
输出的html文件无法正常打开，需要使用convert-spread3.sh脚本处理
将脚本复制到文件夹中，右键点击你的文件夹，在菜单中选择 “Git Bash Here”。
在Git Bash 中运行以下命令bash convert-spread3.sh
运行完成后生成new.html文件即可正常显示
生成BF值传播图
依次上传location.rates.log,（用宿主或地区分log文件，不用总的log文件）经纬度文件，地图.json文件






##  翻译核酸为对应蛋白序列, --trim去除结尾的*
mkdir protein_seq
nohup bash -c ' for i in *.fasta; do seqkit translate --trim ${i} -j 8 > ./protein_seq/${i}; done' > translate.log &
##  格式化结果并显示表头
grep -v '^## ' *.emapper.annotations | sed '1 s/^#//' > output
csvtk -t headers -v output
summarizeAbundance.py -i ../../salmon/gene.TPM -m output --dropkeycolumn -c '7,12,19' -s '*+,+,' -n raw -o eggnog
sed -i 's#^ko:## ' eggnog.KEGG_ko.raw.txt
sed -i '/^-/d' eggnog*
head -n3 eggnog*
## 添加注释生成STAMP的spf格式
cd /home/student/lilanqi/aquaproduct/spade/prodigal/nucleotide/translate/eggnog
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{print a[$1],$0}' /home/student/metagenome/EasyMicrobiome/kegg/KO_description.txt eggnog.KEGG_ko.raw.txt | sed 's/^\t/Unannotated\t/' > eggnog.KEGG_ko.TPM.spf
head -n 5 eggnog.KEGG_ko.TPM.spf
## KO to level 1/2/3
summarizeAbundance.py -i eggnog.KEGG_ko.raw.txt -m /home/student/metagenome/EasyMicrobiome/kegg/KO1-4.txt -c 2,3,4 -s ',+,+,' -n raw --dropkeycolumn -o KEGG
head -n3 KEGG*
## 
## 多个文件寻找相同的序列
seqkit common [flags]
# By ID (default,>后面，空格之前的名字)输出ID名字相同的
seqkit common test1.fa test2.fa -o common.fasta
# By full name（整个序列的名字，包含description部分）。输出序列名字相同的。
seqkit common test1.fa test2.fa -n -o common.fasta
# 输出要比较的文件中序列相同的序列
seqkit common test1.fa test2.fa -s -i -o common.fasta
# 输出要比较的文件中序列相同的序列 (for large sequences)
seqkit common *.fasta -s -i -o common.fasta --md5
https://bioinf.shenwei.me/seqkit/
## ## 
## 使用samtools，bcftools，bedtools，seqkit提取上下游基因片段
samtools faidx *.fasta
最后生成的.fai文件如下， 共5列，\t分隔；
one 66 5 30 31
two 28 98 14 15
第一列 NAME   :   序列的名称，只保留“>”后，第一个空白之前的内容；
第二列 LENGTH:   序列的长度， 单位为bp；
第三列 OFFSET :   第一个碱基的偏移量， 从0开始计数，换行符也统计进行；
第四列 LINEBASES : 除了最后一行外， 其他代表序列的行的碱基数， 单位为bp；
第五列 LINEWIDTH : 行宽， 除了最后一行外， 其他代表序列的行的长度， 包括换行符， 在windows系统中换行符为\r\n, 要在序列长度的基础上加2；
## #提取gff文件的所有基因位置,并转换成bed格式
## #将标准注释gff3文件转换得到bed格式的gff文件
## #方法1  #调用的bedops，也可以自己
awk '{if($3~/^gene$/)print}' EVM.gff3 > genes.gff && convert2bed  --input=gff --output=bed  < gene.gff >genes.bed #调用的bedops，也可以自己用awk
awk '{if($3~/^gene$/)print}' EVM.gff3 > genes.gff && gff2bed <genes.gff> genes.bed #结果同上
## #方法2
less EVM.final.gene.gff3 |grep -w gene|awk '{print$1"\t"$4-1"\t"$5"\t"$9"\t"".""\t"$7}'  >genes.gff
2.2 计算染色体长度
samtools faidx *.fasta
cut -f 1,2 final.genome.fasta.fai >final.genome.fasta.len
2.3 创建包含promoter位置的bed文件,up or down位置
## 注意事项
## slop-根据已有特征区间向外延展，可分别指定上下游延伸长度；
## flank-根据现有区间，指定侧翼延伸长度，得到两侧翼位置的新区间，不包含现有区间；
## 一般默认启动子区域应该是上下游2kb，最大不超过5kb。
#up or down 
# 提取基因上游3k 区间
bedtools flank -i genes.bed -g final.genome.fasta.len  -l 3000 -r 0 -s > up_3k.promoters.bed
# 提取基因下游2k的区间
bedtools flank -i genes.bed -g final.genome.fasta.len  -l 0 -r 2000 -s > down_2k.promoters.bed
# 提取基因上游3k，下游2k区间
bedtools flank -i genes.bed -g final.genome.fasta.len  -l 3000 -r 2000 -s > up_down.promoters.bed
#提取up+gene+down区间
bedtools slop  -i genes.bed -g final.genome.fasta.len  -l 3000 -r 2000 -s > up_3k.promoters.slop.bed
## ##  参数说明
# -l 基因起始位置前多少bp
# -r 基因后多少bp
# -s 考虑正负链
2.4 根据bed中的位置信息，在基因组序列中提取指定序列
#分上下游提取：结果1示例cmd
bedtools getfasta -s -fi final.genome.fasta -bed up_3k.promoters.bed -fo up_3k.promoters.fa -name
bedtools getfasta -s -fi final.genome.fasta -bed down_2k.promoters.bed -fo down_2k.promoters.fa -name
#上下游一起提取，fa文件中id会一致：结果2示例cmd
bedtools getfasta -s -fi final.genome.fasta -bed up_down_2k.promoters.bed -fo down_2k.promoters.fa -name
#提取上游+gene+下游，重点：结果3示例cmd
bedtools getfasta -s -fi final.genome.fasta -bed up_3k.promoters.slop.bed -fo J493.up_genes_down.fa -name

#提取gene序列，不考虑延长基因坐标左右的边位置。
bedtools getfasta -s -fi ../01.data/03.Assembly_final/final.genome.fasta -bed genes.bed -fo genes.fa -name
## #参数说明：
-name 显示名字，即bed的第四列的名字。不写则显示坐标的范围
-s strand #考虑正负链条
-fi 基因组fa文件
-bed 准备好的bed格式文件
-fo 输出文件名 
export PATH=/share/nas1/pengzw/software/bedops-v2.4.38/bin/:$PATH
#step0:数据准备
gff=Chr_genome_final_gene.gff3
genome=Lachesis_assembly_changed.fa
fai=Lachesis_assembly_changed.fa.fai
#step1:将标准注释gff3文件转换得到bed格式的gff文件
convert2bed  --input=gff --output=bed  <$gff>all.bed && awk '{if($8~/^gene$/)print}' all.bed  >genes.bed
#awk '{if($3~/^gene$/)print}' $gff > genes.gff && gff2bed <genes.gff> genes.bed #结果同上
#step2准备基因组长度文件
#samtools faidx $genome >fai
cut -f 1,2 $fai >genome.len
len=genome.len
#step3:getfasta
#flank 分别提取上游，下游序列，因为写在仪器，id会一样，无法区分上下游
l=2000
r=2000
bedtools flank  -i genes.bed -g $len -l $l -r 0 -s > up.flank.bed
bedtools getfasta -s -fi $genome -bed up.flank.bed -fo up.gene.flank.promoter.fa -name 
bedtools flank  -i genes.bed -g $len -l 0 -r $r -s > down.flank.bed
bedtools getfasta -s -fi $genome -bed down.flank.bed -fo down.gene.flank.promoter.fa -name 
#slop 上下游延伸提取
l=2000
r=2000
bedtools slop  -i genes.bed -g $len -l $l -r $r -s > slop.bed
bedtools getfasta -s -fi $genome -bed slop.bed -fo all.gene.slop.promoter.fa -name
提取序列：
samtools faidx input.fa chr1 > chr1.fa
samtools faidx input.fa chr1:100-200 > chr1.fa
##  输出A和B有交集的区域
bedtools intersect -a cpg.bed -b exons.bed  > a_int_b.txt
##  在有重叠区域，输出文件A中的原始特征
bedtools intersect -a cpg.bed -b exons.bed  -wa |head
##  在有重叠区域，输出文件A和文件B的原始特征
bedtools intersect -a cpg.bed -b exons.bed  -wa -wb |head
##  对文件A中的每个特征输出与文件B的重叠，如果没有重叠，则为B输出为NULL
bedtools intersect -a cpg.bed -b exons.bed -loj |head
##  输出文件A和B的特征以及它们之间的碱基对重叠数量
bedtools intersect -a cpg.bed -b exons.bed -wao |head
##  只输出文件A中重叠一次的特征
bedtools intersect -a cpg.bed -b exons.bed -u |head
##  对文件A中的每个条目，输出与文件B重叠的次数
bedtools intersect -a cpg.bed -b exons.bed -c |head
##  只输出文件A中不与文件B重叠的特征
bedtools intersect -a cpg.bed -b exons.bed -v |head
##  A文件与一个或多个B文件取交集
bedtools intersect -a exons.bed -b cpg.bed gwas.bed hesc.chromHmm.bed -sorted | head
##  用标签标识，"A"文件（例如外显子）与多个"B"文件（例如CpG岛、GWAS SNPs和ChromHMM注释）相交时，每个交集来自哪个"B"文件
bedtools intersect -a exons.bed -b cpg.bed gwas.bed hesc.chromHmm.bed -sorted -wa -wb -names cpg gwas chromhmm |head -n 10000 |tail
##  对bed文件排序
sort -k1,1 -k2,2n foo.bed > foo.sort.bed
##  合并重叠区间
bedtools merge -i exons.bed |head
##  合并重叠区间，同时输出合并成新区间的原始区间数量
bedtools merge -i exons.bed -c 1 -o count|head
-c #参数用于指定输入文件中你想要总结的列，
-o #参数定义了你希望应用到`-c`参数列上的操作
##  合并那些虽然不重叠但彼此接近的区间。例如，如果你想合并所有相距不超过1000个碱基对（bp）的特征区间，你可以使用 `-d 1000` 参数
bedtools merge -i exons.bed -d 1000 -c 1 -o count|head
##  合并后，同时列出合并后区间包含的每个特征（例如外显子）的名称
bedtools merge -i exons.bed -d 90  -c 1,4  -o count,collapse|head
bedtools genomecov -i ./demo_date/exons.bed -g ./demo_date/genome.txt |head
-d #生成的输出详细列出了每个位置的覆盖次数
-bg #输出BEDGRAPH格式，它合并连续的具有相同覆盖深度的区域
-bga #参数提供了类似于BEDGRAPH的输出，但是对于未覆盖的区域也会显示出来，覆盖次数为0
-d -split #参数考虑了剪接对覆盖度的影响，剪接或分裂的读段被计算在内
-bga -split #纳入剪接读段的考虑，并输出EDGRAPH格式，同样会显示未覆盖区域
如果没有bed文件，使用gff，vcf文件也可以
bedtools intersect -a file1.bed -b file2.bed > output.bed
vcftools --vcf input.vcf --freq --out output
samtools sort -o output.bam input.bam
plotHeatmap -m matrix.gz -out output.pdf
nohup bash -c 'for i in $(ls *.fq.gz | cut -d'_' -f1); do spades.py --isolate -1 ${i}_R1.fq.gz -2 ${i}_R2.fq.gz -o ./spade/${i} -t 8 -m 300 -k 21,33,55,77; done' > spade.log &
conda activate SNP
nohup bash -c 'for i in *.fasta; do transposon_classifier_RFSB -mode classify -fastaFile ${i} -outputPredictionFile transposon/${i%.fasta}.txt; done' > transposon.log &
transposon_classifier_RFSB -mode evaluate -predLabelFile demoFiles/demo2_predLabel.txt -trueLabelFile demoFiles/demo2_trueLabel.txt -outputPickleFile True
transposon_classifier_RFSB -mode trainModel -fastaFile demoFiles/demo3_transposonDB.fasta -labelFile demoFiles/demo3_labels.txt -outputModelFile demoFiles/demo3_model.pickle -eThreshold 5.0
#MLST typing
#Using software MLST (https://github.com/peterk87/sistr_cmd/)
mlst --csv contigs.fasta > /output/rusults.csv
#AMR analysis
#ARG and plasmid replicon analysis
#Using software ABRicate (https://github.com/tseemann/abricate)
abricate --db resfinder --quiet contigs.fasta
abricate --db plasmidfinder --quiet contigs.fasta
#Point mutation analysis
#Using software AMRFinder (https://github.com/ncbi/amr)
amrfinder -n contigs.fasta --plus -O Salmonella
#Sequence alignment
#Using software Clinker (https://github.com/gamcil/clinker)
clinker files/*.gbk -p plot.html
#Phylogenetic analysis
#Using software snippy (https://github.com/tseemann/snippy)
#Using software gubbins (https://github.com/nickjcroucher/gubbins#generating-input-files)
#Using software RAxML (https://github.com/stamatak/standard-RAxML)
snippy --cpus 64 --outdir ./outdir/contigs.fa-SNP --ref ./SO4698-09.fna --ctgs contigs.fa
snippy-core --ref ./SO4698-09.fna *-SNP
snp-sites -b -c -o phylo.aln core.full.aln
snippy-clean_full_aln core.full.aln > clean.full.aln
run_gubbins.py --threads 64 -p gubbins clean.full.aln
raxmlHPC -f a -x 12345 -# 1000 -p 12345 -m GTRGAMMA -s core.aln -n ex -T 64
library(rhierbaps)
snp.matrix <- load_fasta(fasta.file.name)
hb.results <- hierBAPS(snp.matrix, max.depth = 2, n.pops = 20, quiet = TRUE)
#Genome annotation
#Using software Prokka (https://github.com/tseemann/prokka)
prokka --kingdom Bacteria --outdir mydir --prefix mygenome contigs.fa
#Pan genome and GWAS analysis
#Using software Roary (https://github.com/tseemann/prokka)
#Using software Scoary (https://github.com/tseemann/prokka)
roary -e --mafft -p 64 *.gff
scoary -g scaex -t example.csv -o mydir -s 11 -c BH -p 0.05
emapper.py --cpu 64 --data_dir /ref/eggNOG/ --tax_scope bacteria -i /protein.fasta --output_dir mydir
#Phylogeography analysis
#Transitions between states are estimated under the asymmetric model of Bayesian Stochastic Search Variable Selection -BSSVS. 
#Using software Beast (https://beast.community/index.html)
software/BEASTv1.10.4/bin/beast -threads 64 align.xml
software/BEASTv1.10.4/bin/logcombiner -burnin 10000000 -trees align-1.trees.txt align-2.trees.txt align-3.trees.txt aligncombinetree.txt
software/BEASTv1.10.4/bin/logcombiner -burnin 10000000 align-1.log.txt align-2.log.txt align-3.log.txt aligncombinelog.txt
software/BEASTv1.10.4/bin/treeannotator -limit 90 -heights mean aligncombinelog.txt output.txt
library(TipDatingBeast)
RandomDates(name="align", reps=20, writeTrees=F)
PlotDRT(name="align", reps=20, burnin=0.1)
