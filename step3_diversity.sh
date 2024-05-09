conda activate py312

samples_input=/scratch365/zhuang8/R16S_data/workspace/input_data
declare -a prefixes
for file in $samples_input/*.final.clean.fq; do
    prefix=$(echo "${file}" | sed -r 's/(_1|_2)\.final\.clean\.fq$//')
    if [[ ! " ${prefixes[@]} " =~ " ${prefix} " ]]; then
        prefixes+=("$prefix")
    fi
done

# bowtie2 (download script)
# cd /scratch365/zhuang8/R16S_data/database/b_index && wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip && unzip GRCh38_noalt_as.zip
# bowtie2-build GCA_000001405.15_GRCh38_genomic.fna GRCH38
# cd /scratch365/zhuang8/R16S_data/database/b_index && wget https://genome-idx.s3.amazonaws.com/kraken/protocol/k2protocol_bowtie2indices.tgz &&  tar -xzvf k2protocol_bowtie2indices.tgz

GRCH38=/scratch365/zhuang8/R16S_data/database/b_index/GRCh38_noalt_as/GRCh38_noalt_as
bowtied_input=/scratch365/zhuang8/R16S_data/workspace/bowtied_input

cd $bowtied_input && rm -rf * && cd -
# bowtie2 for human reads
for prefix in "${prefixes[@]}"; do
    f_fq="${prefix}_1.final.clean.fq"
    r_fq="${prefix}_2.final.clean.fq"
    base=$(basename $prefix)
    echo "Processing $base"
    bowtie2 \
        -x $GRCH38 \
        -p 8 \
        -1 $f_fq \
        -2 $r_fq \
        --un-conc $bowtied_input/${base}_%.bowtied.fq \
        -S human_reads.sam
done


# kraken2 database (download script)
# cd /scratch365/zhuang8/R16S_data/database/kraken/k2protocol_db && wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20240112.tar.gz && tar -xvf k2_standard_16gb_20240112.tar.gz
# cd /scratch365/zhuang8/R16S_data/database/kraken/16S_SILVA138_k2db && wget https://genome-idx.s3.amazonaws.com/kraken/16S_Silva138_20200326.tgz && tar -xvf 16S_Silva138_20200326.tgz

prefixes=()
for file in $bowtied_input/*.bowtied.fq; do
    prefix=$(echo "${file}" | sed -r 's/(_1|_2)\.bowtied\.fq$//')
    if [[ ! " ${prefixes[@]} " =~ " ${prefix} " ]]; then
        prefixes+=("$prefix")
    fi
done

k2_protocol_db=/scratch365/zhuang8/R16S_data/database/kraken/k2protocol_db

# kraken2
kraken_output=/scratch365/zhuang8/R16S_data/workspace/kraken_output
cd $kraken_output && rm -rf * && cd -
mkdir -p $kraken_output/kreports

for prefix in "${prefixes[@]}"; do
    f_fq="${prefix}_1.bowtied.fq"
    r_fq="${prefix}_2.bowtied.fq"
    base=$(basename $prefix)
    kraken2 \
        --db $k2_protocol_db \
        --threads 8 \
        --report-minimizer-data \
        --minimum-hit-groups 3 \
        $f_fq $r_fq > $kraken_output/$base.kraken2 \
        --report $kraken_output/kreports/$base.k2report
done


# bracken
bracken_output=/scratch365/zhuang8/R16S_data/workspace/bracken_output
# cd $bracken_output && rm -rf * && cd -
mkdir -p $bracken_output/breports

for prefix in "${prefixes[@]}"; do
    f_fq="${prefix}_1.bowtied.fq"
    r_fq="${prefix}_2.bowtied.fq"
    base=$(basename $prefix)
    echo "running bracken for $f_fq and $r_fq with base $base"
    bracken -d $k2_protocol_db \
        -i $kraken_output/kreports/$base.k2report \
        -r 100 \
        -l S \
        -t 10 \
        -o $bracken_output/$base.bracken \
        -w $bracken_output/breports/$base.breport

    python $KrakenTools/filter_bracken.out.py \
        -i $bracken_output/$base.bracken \
        --exclude [Pseudomonas,Pseudomonadaceae,Pseudomonadales] \
        -o $bracken_output/${base}_filtered.bracken

done


# Alpha diversity using KrakenTools

# Alpha diversity using KrakenTools
KrakenTools=/scratch365/zhuang8/R16S_data/software/KrakenTools 
alpha_diversity_output=/scratch365/zhuang8/R16S_data/workspace/alpha_diversity_output
rm -rf $alpha_diversity_output
mkdir -p $alpha_diversity_output

for prefix in "${prefixes[@]}"; do
    base=$(basename $prefix)
    declare -a metrics=("BP" "Sh" "F" "Si" "ISi")
    for metric in "${metrics[@]}"; do
        python $KrakenTools/DiversityTools/alpha_diversity.py \
        -f $bracken_output/${base}.bracken \
        -a $metric > $alpha_diversity_output/${base}_${metric}.txt
    done
done

# Beta diversity using KrakenTools
beta_diversity_output=/scratch365/zhuang8/R16S_data/workspace/beta_diversity_output
mkdir -p $beta_diversity_output
bracken_files=()
for prefix in "${prefixes[@]}"; do
    base=$(basename $prefix)
    bracken_files+=("$bracken_output/${base}.bracken")
done

# Join the array into a string of filenames separated by spaces
bracken_files_string=$(IFS=" "; echo "${bracken_files[*]}")
# Run beta diversity calculation
python $KrakenTools/DiversityTools/beta_diversity.py \
    -i $bracken_files_string \
    --type bracken > $beta_diversity_output/output.txt


# Pavian plots
# 1. Open https://fbreitwieser.shinyapps.io/pavian/
# 2. Click ‘Browse…’, and upload the four separate *.breport files
# 3. Once files are uploaded, click ‘Sample’ to see the hierarchical classification visualization results.
# 4. Click on the drop-down list to select the sample that you are viewing.
# 5. (Optional) Choose ‘Configure Sankey’ to change what taxonomical ranks to display, number of taxa \
#   at each level, etc. There are ten different settings that can be changed to customize the plots
# 6. Save the network by clicking ‘Save Network’.


# Compare samples to controls using Pavian
# 1. Open the Pavian Shiny app via https://fbreitwieser.shinyapps.io/pavian/
# 2. Click ‘Browse…’, and upload the ten separate *.k2report files. Once files are uploaded, click
#   ‘Comparison’ to see a table that summarizes the total read counts per organism in each sample.
# 3. Choose ‘Species’ and ‘Z-score (reads)’ to focus on species-level reads and calculate z-scores.
# 4. Sort by maximum z-scores
# 