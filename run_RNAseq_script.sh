#!/bin/bash
#SBATCH -p fast                    # Partition with a Nvidia Tesla V100S GPU 
#SBATCH -n 12
#SBATCH -c 8
#SBATCH --job-name="RNA_seq"      
#SBATCH -o log.RNAseq.%j.out   # File to which STDOUT will be written
#SBATCH -e log.RNAseq.%j.err   # File to which STDERR will be written

module load snakemake
module load fastqc
module load hisat2
module load samtools
# activate conda env with subread

ref_idx="/home/jpeter/DATA/ReferenceGenomes/Athaliana/TAIR10/hisat2_indexes/TAIR10_chr_all"
gff="/home/jpeter/DATA/ReferenceGenomes/Athaliana/Araport11/Araport11_GTF_genes_transposons.Jul2023.gtf"

time_id=$(date +"%FT%H%M%S")
outdir="/home/jpeter/DATA/RNAseq/RNAseq072023/"$time_id"_Analysis"
mkdir $outdir

edger_script="/home/jpeter/Scripts/RNA_seq/RNAseq_script/edgeR_analysis.R"

echo "###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@ hisat2"

#for read1 in $outdir/*/trimming/*.1.fastq
for read1 in ~/DATA/RNAseq/RNAseq072023/230707/*/*.R1.fastq.gz
do
    echo $read1
    bn_dn=$(basename $(dirname $read1))
    read2=${read1%.R1.fastq.gz}.R2.fastq.gz
    run_bname="$(basename ${read1%%.*})"
    tissue_dir=$outdir/$bn_dn/$run_bname/hisat2
    mkdir -p $tissue_dir
    output_bam=$tissue_dir/$run_bname.sorted.bam
    summary="$tissue_dir/$run_bname.hisat_summary.txt"

    srun --exclusive -N1 -n1 hisat2 -p 4 -k 20 --max-intronlen 2000 --rna-strandness RF -x $ref_idx --summary-file $summary --new-summary -1 $read1 -2 $read2 | samtools sort -@ 4 -o $output_bam &

done
wait

echo "###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@ index"
#for read1 in $outdir/*/trimming/*.1.fastq
for read1 in ~/DATA/RNAseq/RNAseq072023/230707/*/*.R1.fastq.gz
do
    bn_dn=$(basename $(dirname $read1))
    run_bname="$(basename ${read1%%.*})"
    tissue_dir=$outdir/$bn_dn/$run_bname/hisat2
    output_bam=$tissue_dir/$run_bname.sorted.bam
    srun --exclusive -N1 -n1 samtools index $output_bam &

done
wait

echo "###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@ featurecounts"

featureCounts -O -s 2 -T 4 -t "gene" -g "gene_id" -p -a $gff -o $outdir/seeds/seeds.gene_models.fCounts.txt $outdir/seeds/*/hisat2/*.sorted.bam
featureCounts -O -s 2 -T 4 -t "gene" -g "gene_id" -p -a $gff -o $outdir/flowers/flowers.gene_models.fCounts.txt $outdir/flowers/*/hisat2/*.sorted.bam

echo "###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@###@@@ edgeR script"

echo Rscript $edger_script $outdir/seeds "seeds" ".gene_models.fCounts.txt"
echo Rscript $edger_script $outdir/flowers "flowers" ".gene_models.fCounts.txt"

