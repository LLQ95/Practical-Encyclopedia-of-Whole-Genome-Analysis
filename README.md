# 实用全基因组分析百科全书
全基因组whole genome sequencing
### 代码说明
将代码中的${your_path}替换成自己系统的路径即可使用
## 1.	组装与评估
```bash
nohup bash -c 'for i in $(ls *.fq.gz | cut -d'_' -f1); do spades.py --isolate -1 ${i}_1.fq.gz -2 ${i}_2.fq.gz -o ./spade/${i} -t 8 -m 300 -k 21,33,55,77; done' > spade.log &
```
宏基因组组装指令
```bash
nohup bash -c 'for i in $(ls *.clean.fq.gz | cut -d'_' -f1); do spades.py --meta -1 ${i}_1.clean.fq.gz -2 ${i}_2.clean.fq.gz -o ./spade/${i} -t 8 -m 600 -k 59,79,99,109,125 --careful; done' > spade.log & 
```
为文件夹下的双端测序的原始文件进行批量拼接，注意修改文件后缀
```bash
nohup bash -c 'ls | while read LINE; do spades.py -1 $LINE/$LINE\_1.fq.gz -2 $LINE/$LINE\_2.fq.gz -o ./spades/$LINE; done' > spades.log &
```
在上级文件夹下为下级文件夹下的contigs.fasta进行批量改名，并统一转移到上级文件夹下
```bash
ls | while read LINE;do cp $LINE/contigs.fasta ./$LINE.fasta;done
```
## quast评估组装的质量
```bash
quast.py *.fasta -o quast_output
```
## checkm评估组装的质量
```bash
nohup bash -c 'checkm lineage_wf -f ./checkmresult.tsv --tab_table -x fasta -t 8 --pplacer_threads 8 ./ .' > checkm.log &
```
阈值为5%，污染率超过5%的菌株应当舍弃
## 2.	基因注释
对fasta文件可以直接运行以下命令，建议编写在同一个sh文件下统一运行
```bash
mlst *.fasta --threads 8 > ./pubmlst.tab
abricate --threads 8 --db plasmidfinder *.fasta > ./plasmidfinder.tab
abricate --threads 8 --db BacMet2_EXP_database *.fasta > ./BacMet.tab
abricate --threads 8 --db card *.fasta > ./card.tab
abricate --threads 8 --db ncbi *.fasta > ./ncbi.tab
abricate --threads 8 --db vfdb *.fasta > ./vfdb.tab
abricate --threads 8 --db resfinder *.fasta > ./resfinder.tab
abricate --threads 8 --db ISfinder *.fasta > ./ISfinder.tab
abricate --threads 8 --db SPI *.fasta > ./SPI.tab
abricate --threads 8 --db Tn *.fasta > ./Tn.tab
abricate --threads 8 --db mobileOG *.fasta > ./mobileOG.tab
seqkit stats -a *.fasta -j 8 -T > genome_stats.tsv
```
## Transposon database
https://tncentral.ncc.unesp.br/
## Prophage database
https://phastest.ca/databases
## IS database
https://www-is.biotoul.fr/
## mobileOG
https://github.com/clb21565/mobileOG-db/tree/main
## BacMet
http://bacmet.biomedicine.gu.se/

## 更新pubmlst数据库
```bash
nohup bash -c 'mlst-download_pub_mlst -j 8 -d ${your_path}/anaconda3/envs/checkm/db/pubmlst' > download_pubmlst.log &
```
## 使用prokka进行批量注释 
```bash
nohup bash -c 'for i in *.fasta; do prokka --outdir ./prokka/${i%.fasta} --cpus 8 --prefix ${i%.fasta} ${i} --addgenes --addgenes --centre X --compliant;done' > prokka.log &
```
## abricate 自建数据库
拿到数据库的fasta文件以后，将其放入对应环境的db下，例如`anaconda3/envs/abricate/db/mobileOG/`
对fasta文件执行以下操作
```bash
awk '/^>/{split($0,a,">"); print "> Tn~~~" a[2] "~~~" a[2] "~~~" a[2]} !/^>/{print}' Tn.fa > modified_Tn.fasta
```
案例二
```bash
awk '/^>/{split($0,a,">"); print "> mobileOG~~~" a[2] "~~~" a[2] "~~~" a[2]} !/^>/{print}' mobileOG.fasta > modified_mobileOG.fasta
```
核心建库语句
```bash
abricate --setupdb --threads 8
```
查看abricate数据库
```bash
abricate –list
```
强制更新abricate自带的数据库
```bash
abricate-get_db --db megares --force
```
### 使用vim编写批处理文件
```bash
vim prokka.sh  
```
按i进入插入模式
按esc 输入:wq并回车，保存并退出文本编辑模式
查看文件是否写入
```bash
cat prokka.sh
```
直接运行批处理命令
```bash
bash prokka.sh
```
程序放服务器的后台运行
```bash
nohup bash prokka.sh > prokka.log &
```
查看正在运行的工作，查看后台有无程序运行
```bash
jobs
```
强制停止任务运行
```bash
kill -9 %1
```
### 使用rgi对耐药基因进行注释，使用card数据库
```bash
nohup bash -c 'for i in *.fasta; do rgi main -a DIAMOND --input_sequence ${i} --local --clean -o card/${i%.fasta}; done' >> card.log &## 批量生成gff3文件
```
```bash
nohup bash -c 'for i in *.fasta; do dante -q ${i} -o gff3/${i%.fasta}.gff3 -c 8; done' > dante.log &
```
### 使用bakta进行文件注释（运行较慢，注释效果优于prokka和prodigal）
```bash
nohup bash -c 'for i in *.fasta; do bakta --db ${your_path}/metagenome/bakta/db -o ./bakta/${i%.fasta} -t 8 -p ${i%.fasta} ${i};done' > bakta.log &
nohup bash -c 'for i in *.fasta; do bakta --db ${your_path}/metagenome/bakta/db -o ./bakta/${i%.fasta} -t 8 -p ${i%.fasta} ${i};done' > bakta.log &
nohup bash -c 'for i in *.fasta; do bakta --db ${your_path}/bakta_db/db -o ./bakta/${i%.fasta} -t 8 -p ${i%.fasta} ${i};done' > bakta.log &
amrfinder_update --force_update --database ${your_path}/bakta_db/db/amrfinderplus-db
conda install -y -c conda-forge -c bioconda --strict-channel-priority ncbi-amrfinderplus=4.0.19
```
（需要指定版本号，否则无法正常安装）
bakta 支持宏基因组注释，添加--meta参数
### 使用amrfinder注释耐药基因
```bash
amrfinder -u -d ${your_path}/metagenome/bakta/db/amrfinderplus-db
amrfinder -n *.fasta -o amrfinder.tab --mutation_all --threads 8
amrfinder -g *.gff
```
# 3.	质粒
从gff文件里直接截取一段contigs，作为fasta文件，将对应的fasta文件用gbk或者gff文件注释，最后再在easyfig或者是clinker上面进行可视化分析
### 使用genomad数据库预测可移动元件
```bash
genomad end-to-end --cleanup --splits 8 GCF_009025895.1.fna.gz genomad_output ${your_path}/genomad_db
nohup bash -c 'for i in *.fasta; do genomad end-to-end --cleanup --splits 8 ${i} genomad_output ${your_path}/genomad_db; done' > genomad.log &
nohup bash -c 'for i in *.fasta; do genomad end-to-end --cleanup --splits 8 ${i} genomad_output ${your_path}/genomad_db/genomad_db; done' > genomad.log &
```
### MGEfinder
```bash
for i in *.fasta; do mefinder find -t 8 --contig ./${i} ./mefinder/${i%.fasta}; done
```
### 质粒染色体分离
```bash
for i in *.fasta; do PlasFlow.py --input ${i} --output ./plasflow/${i%.fasta} --threshold 0.7; done
```
```bash
for i in *.fasta; do mob_recon -i ${i} -o ./mob_suite/${i%.fasta}; done
```
# 4.	泛基因组
```bash
nohup bash -c 'roary -e --mafft -p 8 -g 10000000 -r *.gff' > roary.log &
python3 roary_plots.py core_SNP_tree.tre gene_presence_absence.csv
```
### coinfinder
参考资料：https://zhuanlan.zhihu.com/p/620558098
```bash
coinfinder -i gene_presence_absence.csv -I -p core_SNP_tree.tre -o coinfinder -x 8 -a -d -m
nohup bash -c 'query_pan_genome -a intersection *.gff' > pan_genome.log &
scoary -g gene_presence_absence.csv -t traits.csv -n core_SNP_tree_collapsed.nwk
```
### panaroo
```bash
scoary -t Tetracycline_resistance.csv -g Gene_presence_absence.csv -u -c I EPW
panaroo -i *.gff -o results --clean-mode strict --remove-invalid-genes -t 8
```
panaroo与roary本身并不兼容，需要单独去建立panaroo的环境
泛基因组背景资料介绍
http://sanger-pathogens.github.io/Roary/

# 5.	建树、SNP比较
```bash
snippy-clean_full_aln core_gene_alignment.aln > clean.full.aln
run_gubbins.py -c 8 -p gubbins clean.full.aln
snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
FastTree -gtr -nt clean.core.aln > core_SNP_tree.tre
psdm -l -t 8 -o psdm_snp_clean.tab -P -d "\t" clean.core.aln
snp-dists -j 8 clean.core.aln > snp.tab
```
输出文件core_SNP_tree.tre为核心基因组的树文件,snp.tab 为两两比较的snp文件
## 指定建树方法
```bash
run_gubbins.py -c 8 --prefix enteritidis --tree-builder iqtree --first-model JC --tree-builder raxmlng --model GTR clean.full.aln
```
## 利用gubbins进行祖先重建
```bash
run_gubbins.py -c 8 --best-model --recon-with-dates --date 441_date_1.csv -p enteritidis clean.full.aln
```
# 6.	血清型鉴定以及cgmlst分型
```bash
nohup bash -c 'for i in *.fasta; do SeqSero2_package.py -m k -t 4 -i ${i} -p 8; done' > seqsero.log &
```
批量提取生成文件中的沙门的血清型信息
```bash
for file in SeqSero_result*/SeqSero_result.tsv; do awk 'NR==2' "$file" >> merged_second_lines.txt; done
```
## cgmlst_salmonella分析，cgMLSTschema99 用于做grapetree ## 
```bash
chewBBACA.py AlleleCall -i ./ -g ${your_path}/anaconda3/envs/chewie/db/salmonella/Salmonella_enterica_INNUENDO_cgMLST -o ./AlleleCall --cpu 8 --mode 1
chewBBACA.py ExtractCgMLST -i ./AlleleCall/results_alleles.tsv -o  ./ExtractCgMLST
chewBBACA.py AlleleCallEvaluator -i ./AlleleCall -g ${your_path}/anaconda3/envs/chewie/db/salmonella/Salmonella_enterica_INNUENDO_cgMLST -o ./AlleleCallEvaluator --cpu 8
chewBBACA.py AlleleCall -i ./ -g ${your_path}/anaconda3/envs/chewie/db/listeria/Listeria_monocytogenes_Pasteur_cgMLST -o ./AlleleCall --cpu 8 --mode 1
chewBBACA.py ExtractCgMLST -i ./AlleleCall/results_alleles.tsv -o  ./ExtractCgMLST
chewBBACA.py AlleleCallEvaluator -i ./AlleleCall -g ${your_path}/anaconda3/envs/chewie/db/listeria/Listeria_monocytogenes_Pasteur_cgMLST -o ./AlleleCallEvaluator --cpu 8
```
利用已知数据库文件进行cgmlst建库
```bash
chewBBACA.py PrepExternalSchema -g ${your_path}/anaconda3/envs/chewie/db/cronobacter_alleles -o /home/student/anaconda3/envs/chewie/db/crono --cpu 18
```
cgmlst
```bash
chewBBACA.py AlleleCall -i ./ -g ${your_path}/anaconda3/envs/chewie/db/crono  -o ../filtration_AlleleCall1 --cpu 8 --mode 1
chewBBACA.py ExtractCgMLST -i ../filtration_AlleleCall1/results_alleles.tsv -o  ../ExtractCgMLST
chewBBACA.py AlleleCallEvaluator -i ../filtration_AlleleCall1 -g ${your_path}/anaconda3/envs/chewie/db/crono  -o ../AlleleCallEvaluator --cpu 8
chewBBACA.py AlleleCall -i ./ -g ${your_path}/anaconda3/envs/chewie/db/crono/crono  -o ../filtration_AlleleCall1 --cpu 18 --mode 1
chewBBACA.py ExtractCgMLST -i ../filtration_AlleleCall1/results_alleles.tsv -o  ../ExtractCgMLST
```
## Linux 隐藏文件/文件夹
在文件名前+.即可将文件/文件夹隐藏，需要对文件进行重命名才能够重新显示出来
##  计算GC含量
```bash
seqkit fx2tab -l -g -n -i -H *.fasta -j 8 > gc.tab
```
# 7.	去重、去冗余
## 使用dereplicator去除重复的样本（需要提前安装mash）
```bash
nohup bash -c 'dereplicator.py ./ dereplicator_0.0001 --distance 0.0001 --thread 32' > dereplicator_0.0001.log &
```
## 使用cd-hit去冗余样本
```bash
nohup bash -c 'for i in *.fasta; do cd-hit -i ${i} -o cd-hit/${i} -aS 0.9 -c 0.95 -G 0 -g 0 -T 0 -M 0; done' > cd-hit.log &
```
# 8.	点突变
## pointfinder点突变
```bash
vim pointfinder.sh
for i in *.fasta; do python3 ${your_path}/pointfinder/PointFinder.py -i \./${i} -o \./pointfinder/ -p \${your_path}/pointfinder/pointfinder_db -s salmonella -m blastn -m_p \${your_path}/ncbi-blast-2.15.0+/bin/blastn; done
```
```bash
nohup bash -c 'for i in *.fasta; do python3 ${your_path}/pointfinder/PointFinder.py -i ${i} -o ./pointfinder -p ${your_path}/pointfinder_db -s salmonella -m blastn -m_p ${your_path}/ncbi-blast-2.16.0+/bin/blastn; done' >>  pointfinder.log &
```
## pointfinder点突变结果汇总
```bash
for file in *.tsv; do awk 'NR > 1 {print FILENAME "\t" $0}' "$file" >> pointfinder.tab; done
```
注意的问题:pointfinder 无法识别包含有_下划线以及.的文件名，运行之前_.均需要去掉
# 9.	GO/KEGG/COG分析
## 使用prodigal进行基因注释和翻译
```bash
nohup bash -c ' for i in *.fasta; do prodigal -i ${i} -a prodigal/${i%.fasta}.faa -o gff/${i%.fasta}.gff -d nucleotide/${i%.fasta}.fa -f gff; done' > prodigal.log &
```
##  eggnog基因注释，需要运行蛋白质文件
```bash
mkdir eggnog
nohup bash -c 'for i in *.faa; do emapper.py --data_dir ${your_path}/metagenome/eggnog -i ${i} --cpu 8 -m diamond --override -o eggnog/${i%.faa}; done' > eggnog.log &
```
格式化结果并显示表头
```bash
grep -v '^## ' *.emapper.annotations | sed '1 s/^#//' > output
csvtk -t headers -v output
```
# 10.	其他常见指令
多序列比对，需要gbk文件
```bash
clinker files/*.gbk -p plot.html
```
```bash
nohup bash -c 'clinker *gbk -p -i 0.5 -o alignments.tab' > dante.log &
```
antiSMASH（antibiotics & Secondary Metabolite Analysis Shell）识别和分析微生物中生物合成基因簇（BGCs）的工具
```bash
nohup bash -c 'for i in *.gff; do antismash ${i} --output-dir antismash --asf --pfam2go --smcog-trees --fullhmmer --output-basename ${i%.fasta}; done' > antismash.log &
```
slurm系列操作指令
查看当前运行情况，一般不用top来查看
```bash
squeue
```
删除当前账户下所有运行的指令
```bash
scancel -u student
```
作业提交指令
```bash
sbatch *.slurm
```
其他生信代码补充
```bash
nohup bash -c 'mlst *.fasta > ./pubmlst.tab' &
```
压缩，解压指令
压缩
```bash
tar -czvf sichuan_salmon.tar.gz *.fasta
```
解压tar.gz文件
```bash
for i in *.tar.gz; do tar -zxvf ${i}; done
tar -zxvf *.tar.gz
```
解压gz文件
```bash
gunzip *.gz
```
参考资料
https://www.ncbi.nlm.nih.gov/pathogens/docs/datasets_assemblies/
# 11.	Download genomes from NCBI database
从https://www.ncbi.nlm.nih.gov/pathogens/ 下载accession号以及对应的信息表
依据NCBI给的accession号进行下载
he NCBI Datasets command-line tools (CLI) （默认指令，一般弃用，遇到卡顿会退出）
```bash
nohup bash -c 'datasets download genome accession --inputfile enteritidis_accessions.txt --api-key 1ef429d37c5d6103dac9cdaec0f54728d009' > ncbi.log &
``` 
从NCBI数据库上下载fasta文件
```bash
unzip *.zip
for i in *.zip; do unzip ${i}; done
``` 
## 保留文件的前15个字符（accession ID）
```bash
for file in *.fna; do mv "$file" "${file:0:15}.fna"; done
``` 
##  Perl语言版本批量改名
```bash
rename 's/\.all.fna/\.fasta/' *
rename 's/\.fna/\.fasta/' *
rename 's/GCA/GCA_/' *
rename 's/GCA/GCA_/' *
rename 's/\.seq/\.fasta/' *
rename 's/21L/21L-/' *
rename 's/S21./S21_/' *
rename 's/.seq/.fasta/' *
rename 's/.gff/-T.gff/' *
``` 
##  rename (util-linux 2.23.2) C语言版本
```bash
rename -v 'file' 'doc' *.txt
rename -v 'md5' 'MD5' *.txt
``` 
##  整理并获取biosample唯一的序列号，整合到download.txt下
```bash
nohup bash -c 'iseq -i download.txt -p 8 -g -d sra' > download.log &
nohup bash -c 'cat download.txt | while read Run; do iseq -i $Run -a -g; done' > download.log &
``` 
## conda 指令增加下载渠道
```bash
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --show channels
``` 
单独运行某个py的方法，将py程序复制粘贴到对应环境的bin中,再添加权限后即可顺利运行
```bash
chmod +x ${your_path}/anaconda3/envs/checkm/bin/dereplicator.py
chmod +x ${your_path}/anaconda3/envs/env_roary/bin/roary_plots.py
``` 
## 在服务器上使用ollama进行deepseek本地化运行
```bash
ollama run DeepSeek-R1-Distill-Qwen-32B-GGUF:latest
ollama list
``` 
# 12.	三代分析
Genome assembly 三代nanopore数据分析
## Using software spades 
(http://cab.spbu.ru/software/spades/)
## Using software Unicycler 
(https://github.com/rrwick/Unicycler)
## Using software Flye 
(https://github.com/fenderglass/Flye)
## Using software Pilon 
(https://github.com/broadinstitute/pilon/)
```bash
spades.py -k 21,33,55,77,99,127 --careful --pe1-1 short_reads_1.fastq.gz --pe1-2 short_reads_2.fastq.gz -o /output/dir/ --phred-offset 33
unicycler -1 short_reads_1.fastq.gz -2 short_reads_2.fastq.gz -l long_reads.fastq.gz -o /output/dir/
canu -p output -d out_dir -fast genomeSize=5m -nanopore-raw reads.fastq.gz
flye -o out_dir --genome-size 5m --threads 16 --nano-raw reads.fastq.gz
``` 
## 三代数据分析
```bash
bwa index ./contigs.fa
bwa mem -t 12 ./contigs.fa ./short_reads_1.fastq.gz short_reads_2.fastq.gz | samtools view - -Sb | samtools sort - -@14 -o ./mapping.sorted.bam
samtools index mapping.sorted.bam
java -Xmx16G -jar ~/software/pilon-1.23.jar --genome ./contigs.fa --fix all --changes --frags mapping.sorted.bam --threads 12 --output ./pilon/pilon_round1 | tee ./pilon/round1.pilon
``` 
## Serotyping Using software SISTR 
(https://github.com/peterk87/sistr_cmd/)
```bash
sistr -i contigs.fasta contigs.fasta -f csv -o /output/rusults -p CGMLST_PROFILES -n NOVEL_ALLELES
cat *.csv >> ./totalsero.csv
``` 
# 13.	分子钟
##  Estimation of time scaled phylogenies (核心语句，用于筛选离群值)
```bash
treetime --tree 441_core_SNP_tree.tre --dates 441_date.csv --aln clean.core.aln --outdir salmonella --coalescent skyline
``` 
其他输出文件包括与重建祖先序列的比对、带注释的 nexus 格式树，其中分支长度对应于年份和突变，并且节点日期（采用上面详述的数字格式）作为注释添加到每个节点。此外，绘制从根到尖与时间的回归和树并保存到文件中。
##  Ancestral sequence reconstruction （祖先序列重建）
```bash
treetime ancestral --aln lean.core.aln --tree core_SNP_tree.tre --outdir ancestral_results
```
```bash
treetime homoplasy --aln clean.core.aln --tree core_SNP_tree.tre
```
重建祖先序列并将突变映射到树中。输出由一个包含祖先序列的文件“ancestral.fasta”和一个树“annotated_tree.nexus”组成，其中添加了突变作为注释，例如 A45G、G136T、…，SNP 中的数字默认使用基于 1 的索引。推断出的 GTR 模型写入标准输出。
```bash
treetime ancestral [-h] --aln ALN [--vcf-reference VCF_REFERENCE]
                   [--tree TREE] [--rng-seed RNG_SEED] [--gtr GTR]
                   [--gtr-params GTR_PARAMS [GTR_PARAMS ...]] [--aa]
                   [--custom-gtr CUSTOM_GTR] [--marginal] [--keep-overhangs]
                   [--zero-based] [--reconstruct-tip-states]
                   [--report-ambiguous]
                   [--method-anc {parsimony,fitch,probabilistic,ml}]
                   [--verbose VERBOSE] [--outdir OUTDIR]
```                  
##  Estimation of evolutionary rates and tree rerooting
```bash
treetime clock --tree core_SNP_tree.tre --dates date.csv --outdir clock_results --aln clean.core.aln
```   
计算从根到尖的回归并量化树的“时钟似然性”。除非使用 –keep-root 运行，否则它将重新根植树以最大化时钟似信号并重新计算分支长度
```bash
treetime clock [-h] --tree TREE [--rng-seed RNG_SEED] [--dates DATES]
               [--name-column NAME_COLUMN] [--date-column DATE_COLUMN]
               [--sequence-length SEQUENCE_LENGTH] [--aln ALN]
               [--vcf-reference VCF_REFERENCE] [--clock-filter CLOCK_FILTER]
               [--clock-filter-method {residual,local}]
               [--reroot REROOT [REROOT ...] | --keep-root]
               [--tip-slack TIP_SLACK] [--covariation] [--allow-negative-rate]
               [--plot-rtt PLOT_RTT] [--verbose VERBOSE] [--outdir OUTDIR]
```  
##  Inference of transition between discrete characters and ‘mugration’ models
```bash
treetime mugration --tree *.nwk --states *.csv --attribute country
```  
离散字符与“迁移”模型之间的转换推断
此命令将生成带注释的 nexus 树，其中属性的状态作为注释添加到每个节点（例如[&country="brazil"]）。此外，推断出的不同状态之间的 GTR 模型被写入文件。

##  Analyzing homoplasies and recurrent mutations
```bash
treetime homoplasy --aln 441_clean.core.aln --tree 441_core_SNP_tree.tre
```  
## 同源性
重建祖先序列并将突变映射到树上。然后扫描树以查找同源性。过多的同源性可能表明存在污染、重组、文化适应或类似情况。
```bash
treetime homoplasy [-h] --aln ALN [--vcf-reference VCF_REFERENCE]
                   [--tree TREE] [--rng-seed RNG_SEED] [--const CONST]
                   [--rescale RESCALE] [--detailed] [--gtr GTR]
                   [--gtr-params GTR_PARAMS [GTR_PARAMS ...]] [--aa]
                   [--custom-gtr CUSTOM_GTR] [--zero-based] [-n N]
                   [--drms DRMS] [--verbose VERBOSE] [--outdir OUTDIR]
```  
