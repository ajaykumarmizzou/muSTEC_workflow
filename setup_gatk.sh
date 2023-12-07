#!/bin/bash
filename="/data/akt5b/.Github_repositories/muSTEC_workflow/data/tura_SRR_samples/accesion_list.txt"
chdir="/data/akt5b/.Github_repositories/muSTEC_workflow/results/tura_gatk_results/running_on_all_tura_samples/filtered_snp_vcf/"
# Check if the file exists
if [ -e "$filename" ]; then
    # Open the file for reading
    exec 3< "$filename"

    # Read the file line by line using a while loop
    while IFS= read -r line <&3; do
        # Process each line (replace this with your actual processing)
        l=$(echo "${line}_1.fastq")
        r=$(echo "${line}_2.fastq")
        sam=$(echo "${line}.sam")
        fastqc-dump --split-files $line 
        # echo $l $r 
        # # Run the GATK pipeline
        # #Aligning reads to reference genome
        bwa mem ref_seq.fa $l $r > $sam
        # #Sort the mapped reads
        gatk SortSam -I $(echo "${line}.sam") -O $(echo "${line}_sorted.bam") -SO coordinate --CREATE_INDEX true --MAX_RECORDS_IN_RAM 5000000 --VALIDATION_STRINGENCY LENIENT
        # #Deduplicate the sorted bam file
        gatk MarkDuplicates -I $(echo "${line}_sorted.bam") -O $(echo "${line}_deduped_reads.bam") -M $(echo "${line}_deduped.metrics")  --CREATE_INDEX true --MAX_RECORDS_IN_RAM 5000000 --VALIDATION_STRINGENCY LENIENT
        # #Assign the reads to a read group
        gatk AddOrReplaceReadGroups -I $(echo "${line}_deduped_reads.bam") -O $(echo "${line}_addrepl.bam") -PL ILLUMINA -LB $line -PU $line -SM $line -SO coordinate --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT
        # #Call SNPs and INDELs via local re-assembly of haplotypes
        gatk HaplotypeCaller -I $(echo "${line}_addrepl.bam") -O $(echo "${line}.g.vcf") -R ref_seq.fa -ERC GVCF
        # #Perform joint genotyping on one or more samples pre-called with HaplotypeCaller
        gatk GenotypeGVCFs -V $(echo "${line}.g.vcf") -O $(echo "${line}.vcf") -R ref_seq.fa
        # #Select only a certain type of variants from the input file : INDEL, SNP, MIXED, SYMBOLIC, NO_VARIATION
        gatk SelectVariants -V $(echo "${line}.vcf") -R ref_seq.fa -O $(echo "${line}_snp.vcf") --select-type-to-include SNP
        gatk SelectVariants -V $(echo "${line}.vcf") -R ref_seq.fa -O $(echo "${line}_indel.vcf") --select-type-to-include INDEL
        # #Filter variant calls based on INFO and/or FORMAT field in the vcf file
        gatk VariantFiltration -V $(echo "${line}_snp.vcf") -R ref_seq.fa -O $(echo "${line}_snp.filter.vcf") -filter-expression 'QD < 2.0 || FS > 60.0 || MQ <40.0' --filter-name my_snp_filter
        gatk VariantFiltration -V $(echo "${line}_indel.vcf")  -R ref_seq.fa -O $(echo "${line}_indel.filter.vcf") -filter-expression 'QD < 2.0 || FS > 200.0 || MQ < 40.0' --filter-name my_indel_filter       
        #Run a grep for loop to extract the filtered variants from the vcf file
        in=$(echo "$chdir${line}_snp.filter.vcf")
        outdir='/data/akt5b/.Github_repositories/muSTEC_workflow/results/tura_gatk_results/running_on_all_tura_samples/filtered_snp_vcf/filtered_only'
        out=$(echo "$outdir/${line}_snp.filter.vcf")
        grep -v "my_snp_filter" $in > $out
        
    done
    # Close the file descriptor
    exec 3<&-
# #bgzip to compress and tabix to index then merge into single file
for i in /data/akt5b/.Github_repositories/muSTEC_workflow/results/tura_gatk_results/running_on_all_tura_samples/merged_snp_vcf/*;do bgzip $i;tabix -p vcf $i.gz;done
bcftools merge *.gz > merged_snp.vcf
# #Annotate the merged vcf file
# merged=$("/data/akt5b/.Github_repositories/muSTEC_workflow/results/tura_gatk_results/running_on_all_tura_samples/merged_snp_vcf/merged_snp.vcf")
java -Xmx8g -jar /data/akt5b/.Github_repositories/muSTEC_workflow/muSTEC_workflow/tools/snpEff/snpEff/snpEff.jar  -v e.coli.o157.h7.sakai  $merged > merged_annotated.vcf
# #Phylogeny creation
/data/akt5b/.Github_repositories/muSTEC_workflow/tools/vcf2phylip/vcf2phylip.py -i /data/akt5b/.Github_repositories/muSTEC_workflow/results/abdalhamid_results_2/annotated_vcf/merged_annotated.vcf -o /data/akt5b/.Github_repositories/muSTEC_workflow/results/abdalhamid_results_2/phylogeny/merged_annotated.phy
#Fasttree
fasttree ./phylogeny/merged_annotated.min4.phy > ./tree/merged_ann_vcf_phy.tree
#Trim # snpann file
grep -v "##" merged_annotated.vcf > merged_annotated_trim.vcf

#snpsift tool
#filter missense variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'missense_variant')" merged_annotated.vcf > missense_.vcf
#filter synonymous variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'synonymous_variant')" merged_annotated.vcf > synonymous_.vcf
#filter high impact variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'HIGH')" merged_annotated.vcf > high_impact_.vcf
#filter moderate impact variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'MODERATE')" merged_annotated.vcf > moderate_impact_.vcf
#filter low impact variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'LOW')" merged_annotated.vcf > low_impact_.vcf
#filter modifier impact variants

java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'MODIFIER')" merged_annotated.vcf > modifier_impact_.vcf
#filter stop gained variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'stop_gained')" merged_annotated.vcf > stop_gained_.vcf
#filter stop lost variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'stop_lost')" merged_annotated.vcf > stop_lost_.vcf
#filter start lost variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'start_lost')" merged_annotated.vcf > start_lost_.vcf
#filter splice acceptor variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'splice_acceptor_variant')" merged_annotated.vcf > splice_acceptor_.vcf
#filter splice donor variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] has 'splice_donor_variant')" merged_annotated.vcf > splice_donor_.vcf
#filter frameshift variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] contains 'frameshift')" merged_annotated.vcf > frameshift_.vcf
#filter nonframeshift variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] contains 'nonframeshift')" merged_annotated.vcf > nonframeshift_.vcf
#filter nonframeshift deletion variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] contains 'nonframeshift_deletion')" merged_annotated.vcf > nonframeshift_deletion_.vcf
#filter nonframeshift insertion variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] contains 'nonframeshift_insertion')" merged_annotated.vcf > nonframeshift_insertion_.vcf
#filter frameshift deletion variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] contains 'frameshift_deletion')" merged_annotated.vcf > frameshift_deletion_.vcf
#filter frameshift insertion variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] contains 'frameshift_insertion')" merged_annotated.vcf > frameshift_insertion_.vcf
#filter splice region variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] contains 'splice_region')" merged_annotated.vcf > splice_region_.vcf
#filter intron variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] contains 'intron')" merged_annotated.vcf > intron_.vcf
#filter 5 prime UTR variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] contains '5_prime_UTR')" merged_annotated.vcf > 5_prime_UTR_.vcf
#filter 3 prime UTR variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(ANN[*] contains '3_prime_UTR')" merged_annotated.vcf > 3_prime_UTR_.vcf
#filter intergenic region variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(isIntergenic)" merged_annotated.vcf > intergenic_.vcf
#filter upstream gene variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(isUpstream)" merged_annotated.vcf > upstream_.vcf
#filter downstream gene variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(isDownstream)" merged_annotated.vcf > downstream_.vcf
#filter intragenic variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(isIntragenic)" merged_annotated.vcf > intragenic_.vcf
#filter intergenic variants
java -jar /data/akt5b/.Github_repositories/muSTEC_workflow/tools/snpEff/SnpSift.jar filter "(isIntergenic)" merged_annotated.vcf > intergenic_.vcf




else
    echo "File not found: $filename"
fi
