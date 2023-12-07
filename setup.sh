#!/bin/bash
echo "Welcome to muSTEC workflow setup script!"
ref=$1 #Passing reference file as argument
#Check if reference file exists
if [ ! -f "$ref" ]; then
    echo "Reference file does not exist!"
    exit 1
fi
# less $ref
#Add a flag for single and paired end reads
echo "Are your reads single or paired end? (s/p)"
read flag
if [ "$flag" == "s" ]; then
    echo "Single end reads selected"
#Here each sample for single end file reads and continues for the pipeline
# ------------------------------------------------------------------------------------------------------------------------
#Analysis starts heregit remote add origin https://github.com/ajaykumarmizzou/muSTEC_workflow.git
#Count total number of files in the data directory
    total_files=$(ls -l data/tura_data/*.fsa_nt | wc -l)
    echo "Total files: $total_files"
#read files in i variable
    for ((i=1;i<=total_files;i++))
    do
        # echo "Sample $i:"
        file=$(ls data/tura_data/*.fsa_nt | sed -n "$i"p)
        # less $file
        # echo "Read: $file"
        # echo "Is this correct? (y/n)"
        # read ans
        # if [ "$ans" == "y" ]; then
        #     echo "Correct!"
        # else
        #     echo "Incorrect!"
        #     exit 1
        # fi #Files successfully readed in i variable so far
    
#Alignment of sequences samples to reference genome - samtools alignment
#run bwa mem for each sample
    f=$(echo "./data/tura_data/"${file%%.*}.1.fsa_nt"" | sed 's|/data/tura_data||')
    s=$(echo "./results/aligned/sam/"${file%%.*}.sam"" | sed 's|/data/tura_data||')
    bwa mem $ref $f > $s
    echo "Alignment done!"
#Converting sam file into bam
    bam=$(echo "./results/aligned/bam/"${file%%.*}.bam"" | sed 's|/data/tura_data||')
    samtools view -bS $s > $bam
    echo "Converting sam into bam, done!"
#Sorting bam file
    sorted=$(echo "./results/aligned/sorted_bam/"${file%%.*}.sorted.bam"" | sed 's|/data/tura_data||')
    samtools sort $bam -o $sorted
    echo "Sorting bam, done!"
#samtools flagstats
    flagstats=$(echo "./results/aligned/flagstats/"${file%%.*}.flagstats.txt"" | sed 's|/data/tura_data||')
    samtools flagstat $sorted > $flagstats
    echo "Samtools flagstats, done!"
#bcf mpileup
    mpileup=$(echo "./results/variant_calling/bcf/"${file%%.*}.bcf"" | sed 's|/data/tura_data||')
    bcftools mpileup -O b -o $mpileup -f $ref $sorted
    echo "bcf mpileup, done!"
#Identify SNVs using call ploidy
    vcf=$(echo "./results/variant_calling/vcf/"${file%%.*}.vcf"" | sed 's|/data/tura_data||')
    bcftools call --ploidy 1 -m -v -o $vcf $mpileup
    echo "Identifying SNVs, done!"
#Filter using varFilter
    filtered=$(echo "./results/variant_calling/filtered_vcf/"${file%%.*}.vcf"" | sed 's|/data/tura_data||')
    vcfutils.pl varFilter $vcf > $filtered
    echo "Filtering SNVs, done!"
#bcftools stats to generate vcf stats
    stats=$(echo "./results/variant_calling/vcf_stats/"${file%%.*}.vcf.txt"" | sed 's|/data/tura_data||')
    bcftools stats -F $ref -s - $filtered > $stats
    echo "Generating vcf stats, done!"
#plot-vcfstats
    # plot=$(echo "./results/variant_calling/vcf_stats/"${file%%.*}.vcf.plot.pdf"" | sed 's|/data/tura_data||')
    # plot-vcfstats -p $plot $stats
    # echo "Plotting vcf stats, done!"
#sort variants
    sorted_vcf=$(echo "./results/variant_calling/sorted_vcf/"${file%%.*}.sorted.vcf"" | sed 's|/data/tura_data||')
    bcftools sort -o $sorted_vcf $filtered
    echo "Sorting vcf, done!"
#bgzip and tabix
    bgzip $sorted_vcf
    tabix -p vcf $sorted_vcf.gz
    echo "bgzip and tabix, done!"
    done
#Merge into single vcf file
bcftools merge ./results/variants/sorted_vcf/*.gz > ./results/variant_calling/merged_variants/merged.vcf
echo "Merging vcf, done!"



#For loop to merge all sorted_vcf.gz files into merged.vcf file using bcftools
#Analysis ends here
# ------------------------------------------------------------------------------------------------------------------------


elif [ "$flag" == "p" ]; then
    echo "Paired end reads selected"
#Here each sample for two paired end files reads and continues for the pipeline
    #count total files in the data directory
    total_files=$(ls -l data/abdalhamid_data/*.fastq | wc -l)
    echo "Total files: $total_files"
    #check if total files are even
    if [ $((total_files % 2)) -eq 0 ]; then
        echo "Total files are even"
    else
        echo "Total files are odd"
        exit 1
    fi
    #read files in i and j variables
    for ((i=1,j=2;i<=total_files;i+=2,j+=2))
    do
        echo "Sample $((i/2+1)):"
        file_1=$(ls data/abdalhamid_data/*.fastq | sed -n "$i"p)
        # less $file_1
        echo "Read 1: $file_1"
        file_2=$(ls  data/abdalhamid_data/*.fastq | sed -n "$j"p)
        # less $file_2
        echo "Read 2: $file_2"
        echo "Is this correct? (y/n)"
        read ans
        if [ "$ans" == "y" ]; then
            echo "Correct!"
        else
            echo "Incorrect!"
            exit 1
        fi #Files successfully readed in i and j variables so far
# ------------------------------------------------------------------------------------------------------------------------
#Analysis starts here
#Trimming and quality control
        perl ./tools/trim-galore/trim_galore.pl --paired --fastqc --illumina --output_dir results/abdalhamid_results_2/trimmed/ $file_1 $file_2
#indexing reference genome
        # bwa index $ref
#Alignment of sequences samples to reference genome - samtools alignment
#run bwa mem for each sample
        f1=$(echo "./results/abdalhamid_results_2/trimmed//"${file_1%%.*}_val_1.fq"" | sed 's|/data/abdalhamid_data/||')
        f2=$(echo "./results/abdalhamid_results_2/trimmed//"${file_2%%.*}_val_2.fq"" | sed 's|/data/abdalhamid_data/||')
        s=$(echo "./results/abdalhamid_results_2/aligned/sam//"${file_1%%.*}.sam"" | sed 's|/data/abdalhamid_data/||')
        # echo $s
        bwa mem $ref $f1 $f2 > $s
        echo "Alignment done!"
#Converting sam file into bam
        bam=$(echo "./results/abdalhamid_results_2/aligned/bam//"${file_1%%.*}.bam"" | sed 's|/data/abdalhamid_data/||')
        samtools view -bS $s > $bam
        echo "Converting sam into bam, done!"
#Sorting bam file
        sorted=$(echo "./results/abdalhamid_results_2/aligned/sorted_bam//"${file_1%%.*}.sorted.bam"" | sed 's|/data/abdalhamid_data/||')
        samtools sort $bam -o $sorted
        echo "Sorting bam, done!"
#samtools flagstats
        flagstats=$(echo "./results/abdalhamid_results_2/aligned/flagstats//"${file_1%%.*}.flagstats.txt"" | sed 's|/data/abdalhamid_data/||')
        samtools flagstat $sorted > $flagstats
        echo "Samtools flagstats, done!"
#bcf mpileup
        mpileup=$(echo "./results/abdalhamid_results_2/variant_calling/bcf//"${file_1%%.*}.bcf"" | sed 's|/data/abdalhamid_data/||')
        bcftools mpileup -O b -o $mpileup -f $ref $sorted
        echo "bcf mpileup, done!"
#Identify SNVs using call ploidy
        vcf=$(echo "./results/abdalhamid_results_2/variant_calling/vcf//"${file_1%%.*}.vcf"" | sed 's|/data/abdalhamid_data/||')
        bcftools call --ploidy 1 -m -v -o $vcf $mpileup
#         echo "Identifying SNVs, done!"
#Filter using varFilter
        filtered=$(echo "./results/abdalhamid_results_2/variant_calling/filtered_vcf//"${file_1%%.*}.vcf"" | sed 's|/data/abdalhamid_data/||')
#         vcfutils.pl varFilter $vcf > $filtered
#         echo "Filtering SNVs, done!"
#bcftools stats to generate vcf stats
        stats=$(echo "./results/abdalhamid_results_2/variant_calling/vcf_stats//"${file_1%%.*}.vcf.txt"" | sed 's|/data/abdalhamid_data/||')
        bcftools stats -F $ref -s - $filtered > $stats
        echo "Generating vcf stats, done!"
#plot-vcfstats
        # plot=$(echo "./results/variant_calling/vcf_stats/"${file_1%%.*}.vcf.plot.pdf"" | sed 's|/data/demo||')
        # plot-vcfstats -p $plot $stats
        # echo "Plotting vcf stats, done!"
#sort variants
        sorted_vcf=$(echo "./results/abdalhamid_results_2/variant_calling/sorted_vcf//"${file_1%%.*}.sorted.vcf"" | sed 's|/data/abdalhamid_data/||')
        bcftools sort -o $sorted_vcf $filtered
        echo "Sorting vcf, done!"
#bgzip and tabix
        bgzip $sorted_vcf
        tabix -p vcf $sorted_vcf.gz
        echo "bgzip and tabix, done!"
    done
#Merge into single vcf file
    bcftools merge ./results/abdalhamid_results_2/variant_calling/sorted_vcf/*.gz > ./results/abdalhamid_results_2/variant_calling/merged_variants/merged.vcf
    echo "Merging vcf, done!"
#Annotating multiVCF file using snpEff
    merged=$("./results/abdalhamid_results_2/variant_calling/merged_variants/merged.vcf" | sed 's|/data/tura_data||')
    java -Xmx8g -jar /data/akt5b/.Github_repositories/muSTEC_workflow/muSTEC_workflow/tools/snpEff/snpEff/snpEff.jar  -v e.coli.o157.h7.sakai  $merged > test.chr22.ann.vcf
#Phylogeny creation
/data/akt5b/.Github_repositories/muSTEC_workflow/tools/vcf2phylip/vcf2phylip.py -i /data/akt5b/.Github_repositories/muSTEC_workflow/results/abdalhamid_results_2/annotated_vcf/merged_annotated.vcf -o /data/akt5b/.Github_repositories/muSTEC_workflow/results/abdalhamid_results_2/phylogeny/merged_annotated.phy
#.phy file generated

#tree command 

#multiqc
# multiqc ./results/abdalhamid_results_2/trimmed/ -o ./results/abdalhamid_results_2/multiqc_trimmed/
#merging all vcf into one using bcftools merge from sorted_vcf
    bcftools merge ./results/variants/merged_variants/*.gz > ./results/variant_calling/merged_variants/merged.vcf.gz
    echo "Merging vcf, done!"

#Analysis ends here
# ------------------------------------------------------------------------------------------------------------------------

else
    echo "Invalid input!"
    exit 1
fi

