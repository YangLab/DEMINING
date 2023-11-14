
#!/usr/bin/env bash
#touch $0\.arg 2>/dev/null && args_path=`dirname $0` ||args_path=`pwd`
set -eo pipefail

echo -e `date` dir:$PWD  command:$0 "$*" #>>$args_path/`basename $0`\.arg  ###### Sun Dec 8 10:27:47 CST 2019 fzc introduce $args_path
######!!!!!!!!!Warining!!!!!!!!!!############
##USE bash instead of sh to run this shell script!
######!!!!!!!!!Warining!!!!!!!!!!############
###############Usage############




################################



DEMINING_main_flow(){
    ###### Tue Dec 3 14:35:27 CST 2019
    ####################init#################################################
    ##Usage:
    ###paired end:
    ###single end:

    #ref_genome="hg19" or "hg38"

    ####################init END#################################################
    while getopts :1:2:s:o:n:g:t:p: ARGS  
    do  
    case $ARGS in   
        1)  
            fq1=$OPTARG
            ;;  
        2)  
            fq2=$OPTARG  
            ;;
        s)  
            fq0=$OPTARG
            ;;  
        t)  
            threads=$OPTARG
            ;;  
        o)  
            work_path=$OPTARG
            ;;  
        n)  
            sample_name=$OPTARG
            ;; 
        g)  
            ref_genome=$OPTARG
            ;;  
        p) ###### Sat Feb 22 15:48:50 CST 2020 fzc; specify run step when fix 
            step=$OPTARG
            ;;  
        *)  
            echo "Unknown option: $ARGS"
            ;;
        \?)
        echo "Invalid option: -$OPTARG" 
        ;;
    esac
    done
    tmp_work_path_0=$work_path/${sample_name}_DEMINING_tmp
    tmp_a="__${fq0}__"
    test "$tmp_a" == "____" && {
        layout=paired
        fq_path=`dirname ${fq1}`
        }|| {
            layout=single
            fq_path=`dirname ${fq0}`
        }
    threads_tmp_a="__${threads}__"
    test "$threads_tmp_a" == "____" && threads=10
    step_tmp_a="__${step}__"
    test "$step_tmp_a" == "____" && step=01234

    rRNAdeplete_path="${tmp_work_path_0}/0-0-0remove_rRNA"
    HISAT2_mapping_path="${tmp_work_path_0}/1-1-0HISAT_map"
    BWA_mapping_path="${tmp_work_path_0}/2-0-0BWA_map"
    fine_tune_path="${tmp_work_path_0}/3-0-0Combine_bam"
    MutationCalling_wp=${tmp_work_path_0}/4-0-0Editing_sites
    CMCconst_wp="${tmp_work_path_0}/5-0-0CMC_construction"
    
    echo "#####[`date`]${sample_name} BEGIN..."

    Y_or_N=$(echo "$step"|grep -q "0" && echo "Y"||echo "N")  ###### Sat Feb 22 16:12:59 CST 2020 fzc; specify step
    if [ "$Y_or_N" == "Y" ];then
    ###0. remove rRNA, remove_rRNA_mapped_reads
    echo "#####[`date`]BEGIN remove_rRNA_mapped_reads"
    if [ "$layout" == "paired" ];then
    out_fq1=${rRNAdeplete_path}/${sample_name}_R1.fastq.gz
	out_fq2=${rRNAdeplete_path}/${sample_name}_R2.fastq.gz
    STEP0_remove_rRNA_mapped_reads $fq1 $fq2 $out_fq1 $out_fq2  ###### Wed Dec 23 17:07:21 CST 2020  
    fq1=$out_fq1
    fq2=$out_fq2
    elif [ "$layout" == "single" ];then
    out_fq0=${rRNAdeplete_path}/${sample_name}.fastq.gz
    STEP0_remove_rRNA_mapped_reads $fq0 $out_fq0
    fq0=$out_fq0
    fi
    fi
    Y_or_N=$(echo "$step"|grep -q "1" && echo "Y"||echo "N")  ###### Sat Feb 22 16:12:59 CST 2020 fzc; specify step
    if [ "$Y_or_N" == "Y" ];then
    ###1. mapping, HISAT2_2mismatch_following_BWA_6mismatch_mapping
    if [ "$layout" == "paired" ];then
        echo "#####[`date`]BEGIN STEP1_HISAT2_2mismatch_following_BWA_6mismatch_mapping"
        STEP1_HISAT2_2mismatch_following_BWA_6mismatch_mapping $sample_name $layout $HISAT2_mapping_path $BWA_mapping_path $fine_tune_path $fq1 $fq2 ###### Tue Dec 3 14:35:36 CST 2019 fzc  support set threads number when running mapping
    elif [ "$layout" == "single" ];then
        echo "#####[`date`]BEGIN STEP1_HISAT2_2mismatch_following_BWA_6mismatch_mapping"
        STEP1_HISAT2_2mismatch_following_BWA_6mismatch_mapping $sample_name $layout $HISAT2_mapping_path $BWA_mapping_path $fine_tune_path $fq0
    fi
    fi

    Y_or_N=$(echo "$step"|grep -q "2" && echo "Y"||echo "N")  ###### Sat Feb 22 16:12:59 CST 2020 fzc; specify step
    if [ "$Y_or_N" == "Y" ];then
    ###2. sam_fine_tune, Markduplicate, BQSR
    echo "######[`date`]BEGIN STEP2_sam_fine_tune"
    bam_file=${fine_tune_path}/${sample_name}_combine.bam  ###### Sat Dec 21 06:25:53 CST 2019 fzc; change ${combine_bam} to ${fine_tune_path}
    STEP2_sam_fine_tune $bam_file $fine_tune_path 
    fi

    Y_or_N=$(echo "$step"|grep -q "3" && echo "Y"||echo "N")  ###### Sat Feb 22 16:12:59 CST 2020 fzc; specify step
    if [ "$Y_or_N" == "Y" ];then
    ###3. Samtools_mpileup Calling
    echo "#####[`date`]BEGIN STEP3_Variant_calling"
    bam_file_prefix=${sample_name}_combine
    bam_file=$fine_tune_path/${bam_file_prefix}_rgadd_dedupped_split_recal.bam
    STEP3_Variant_calling $bam_file "_recal" $MutationCalling_wp ${sample_name} $ref_genome
    fi

    Y_or_N=$(echo "$step"|grep -q "4" && echo "Y"||echo "N")  ###### Sat Feb 22 16:12:59 CST 2020 fzc; specify step
    if [ "$Y_or_N" == "Y" ];then
    ###4. CMC_construction 
    echo "#####[`date`]BEGIN STEP4_CMC_construction"
    bam_file_prefix=${sample_name}_combine
    bam_file=$fine_tune_path/${bam_file_prefix}_rgadd_dedupped_split_recal.bam
    RADAR_out=${MutationCalling_wp}/${sample_name}.BQ20o6ES95v2.allvariants_recal
    predict_out_file="$work_path/${sample_name}_DEMINING.tsv"
    STEP4_DeepDDR_predict $ref_genome $sample_name $RADAR_out $bam_file $predict_out_file $CMCconst_wp $threads
    fi
    echo "#####[`date`]${sample_name} END"
}
STEP0_remove_rRNA_mapped_reads(){
    thread=$threads
    test -d $rRNAdeplete_path ||mkdir -p $rRNAdeplete_path
    echo "[`date`] remove rRNA mapped reads" 
    if [ "$layout" == "paired" ];then
    in_fq1=$1
    in_fq2=$2
    out_fq1=$3
    out_fq2=$4
    bwa mem -t ${thread} ${rDNA_index_bwa_mem} $in_fq1  $in_fq2|samtools view -bh  -f 4 |samtools sort -n  -o ${rRNAdeplete_path}/${sample_name}-rRNA_unmapped_sort.bam 
    samtools fastq -c 6 -1 ${out_fq1} -2 ${out_fq2} -s /dev/null ${rRNAdeplete_path}/${sample_name}-rRNA_unmapped_sort.bam
    # fastqc_out=`dirname ${out_fq1}`/fastqc_out
    # test -d $fastqc_out ||mkdir  $fastqc_out
    # $fastqc -o ${fastqc_out} -d ~/tmp -q ${out_fq1} ${out_fq2}  

    elif [ "$layout" == "single" ];then
    in_fq0=$1
    out_fq0=$2
    bwa mem -t ${thread} ${rDNA_index_bwa_mem} $in_fq0|samtools view  -bh  -f 4 |samtools sort -n  -o ${rRNAdeplete_path}/${sample_name}-rRNA_unmapped_sort.bam 
    samtools fastq  -s ${rRNAdeplete_path}/${sample_name}_singleton.fq ${rRNAdeplete_path}/${sample_name}-rRNA_unmapped_sort.bam |gzip > $out_fq0
	fi
}
STEP1_HISAT2_2mismatch_following_BWA_6mismatch_mapping(){
    #HISAT2_2mismatch_following_BWA_6mismatch_mapping $sample_name $layout $fq_path $HISAT2_mapping_path $BWA_mapping_path $combine_bam
    ###### Wed Sep 25 19:04:01 CST 2019
    #used for HISAT2 and BWA two pass mapping.
    #First, HISAT2 2 mismatches
    #HISAT2 unmapped reads are extracted and send to BWA 6 mismatches mapping.
    #for paired ends reads, mismatches are all count at read level not fragment level, for example, HISAT2 2 mismatches means read1's mismatches maximum count is 2,read2's mismatches maximum count is 2, read1 mismatches plus read2 mismatches can up to 4.

    sample_name=$1
    layout=$2
    HISAT2_mapping_path=$3
    BWA_mapping_path=$4
    combine_bam=$5

    test -d $HISAT2_mapping_path || mkdir -p $HISAT2_mapping_path
    test -d $BWA_mapping_path || mkdir -p $BWA_mapping_path
    test -d ${combine_bam} ||mkdir -p ${combine_bam}
    HISAT_map=$HISAT2_mapping_path
    bwa_map=${BWA_mapping_path}

    if [ "$layout" == "paired" ];then
        fq1=$6
        fq2=$7
        ### 1. HISAT2 2 mismatches mapping
        echo "hisat2  --rna-strandness RF --no-mixed --secondary --no-temp-splicesite --known-splicesite-infile ${annotation_splice_sites} --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p ${threads} -x ${genome_index_hisat2} -1  $fq1  -2 $fq2  --un-conc-gz ${HISAT_map}/${sample_name}_un_conc_%.fastq.gz -S ${HISAT_map}/${sample_name}_HISAT2_mapped.sam"
        hisat2  --rna-strandness RF --no-mixed --secondary --no-temp-splicesite --known-splicesite-infile ${annotation_splice_sites} --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p ${threads} -x ${genome_index_hisat2} -1  $fq1  -2 $fq2  --un-conc-gz ${HISAT_map}/${sample_name}_un_conc_%.fastq.gz -S ${HISAT_map}/${sample_name}_HISAT2_mapped.sam  
        
        samtools view -h -F 4 ${HISAT_map}/${sample_name}_HISAT2_mapped.sam|awk 'BEGIN{FS="XM:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "XM"){split($2,a,"\t");if ( a[1] <= 2 ) print $0 } else print $0 " not have XM tag"}}'|awk 'BEGIN{FS="NH:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "NH"){split($2,a,"\t");if ( a[1] == 1 ) print $0 } else print $0 " not have NH tag"  }}' >${HISAT_map}/${sample_name}_unique_mismatch2.sam &
        
        samtools view -bS -f 4 -o ${HISAT_map}/${sample_name}_HISAT2_unmapped.bam ${HISAT_map}/${sample_name}_HISAT2_mapped.sam

        samtools view -f 4 -S ${HISAT_map}/${sample_name}_HISAT2_mapped.sam |awk 'BEGIN{FS="\t"}{print $1}'|sort|uniq >${HISAT_map}/${sample_name}_hisat2_unmap.readid 

        seqtk subseq ${HISAT_map}/${sample_name}_un_conc_1.fastq.gz ${HISAT_map}/${sample_name}_hisat2_unmap.readid |gzip > ${HISAT_map}/${sample_name}_unmapped_1.fastq.gz  &
        seqtk subseq ${HISAT_map}/${sample_name}_un_conc_2.fastq.gz ${HISAT_map}/${sample_name}_hisat2_unmap.readid |gzip > ${HISAT_map}/${sample_name}_unmapped_2.fastq.gz  
        wait

        ### 2. BWA 6 mismatches mapping
        
        bwa mem -t ${threads}  -A 1 -B 4  ${genome_index_bwa_mem}  ${HISAT_map}/${sample_name}_unmapped_1.fastq.gz  ${HISAT_map}/${sample_name}_unmapped_2.fastq.gz > ${bwa_map}/${sample_name}_bwa_mapped.sam |tee 2>${bwa_map}/log_BWA_6mismatch_`date +%Y_%m_%d`.log
    elif [ "$layout" == "single" ];then

        fq0=$6
        hisat2 --secondary --no-temp-splicesite --known-splicesite-infile ${annotation_splice_sites} --no-softclip --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 --dta -t -p ${threads} -x ${genome_index_hisat2} -U $fq0 -S ${HISAT_map}/${sample_name}_HISAT2_mapped.sam |tee 2>${HISAT_map}/log_HISAT2_2mismatch_${sample_name}_`date +%Y_%m_%d`.log ###### Tue Dec 17 19:58:15 CST 2019 fzc; change log path

        samtools view -h -F 4 ${HISAT_map}/${sample_name}_HISAT2_mapped.sam|awk 'BEGIN{FS="XM:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "XM"){split($2,a,"\t");if ( a[1] <= 2 ) print $0 } else print $0 " not have XM tag" }}'|awk 'BEGIN{FS="NH:i:"}{if($0 ~/^@/){print $0}else{if ($0 ~ "NH"){split($2,a,"\t");if ( a[1] == 1 ) print $0 } else print $0 " not have NH tag" }}' >${HISAT_map}/${sample_name}_unique_mismatch2.sam &
        samtools view -bS -f 4 -o ${HISAT_map}/${sample_name}_HISAT2_unmapped.bam ${HISAT_map}/${sample_name}_HISAT2_mapped.sam 
        bedtools bamtofastq -i  ${HISAT_map}/${sample_name}_HISAT2_unmapped.bam -fq /dev/stdout | gzip >${HISAT_map}/${sample_name}_HISAT2_unmapped.fastq.gz
        ### 2. BWA 6 mismatches mapping
        bwa mem -t ${threads} ${genome_index_bwa_mem}  ${HISAT_map}/${sample_name}_HISAT2_unmapped.fastq.gz > ${bwa_map}/${sample_name}_bwa_mapped.sam
    else
        echo "Unknown layout:$layout, please specify single or paired." 

    fi
    python ${DEMINING_path}/src/bwa_unique_mismatch6.py ${bwa_map}/${sample_name}_bwa_mapped.sam ${bwa_map}/${sample_name}_bwa_unique_mis6_mapq0.sam 
    samtools_log_file=${bwa_map}/samtools_${sample_name}.log 
    samtools view -bT ${ref_genome_path} -o ${bwa_map}/${sample_name}_unmapped.nbam ${bwa_map}/${sample_name}_bwa_unique_mis6_mapq0.sam |tee 2>>$samtools_log_file
    samtools sort ${bwa_map}/${sample_name}_unmapped.nbam -o ${bwa_map}/${sample_name}_unmapped.sort.bam |tee 2>>$samtools_log_file
    samtools view -H ${bwa_map}/${sample_name}_unmapped.sort.bam > ${bwa_map}/${sample_name}_bwa.header |tee 2>>$samtools_log_file 
    
    wait
    cat ${bwa_map}/${sample_name}_bwa.header ${HISAT_map}/${sample_name}_unique_mismatch2.sam > ${combine_bam}/${sample_name}_accepted_hits.nsam |tee 2>>$samtools_log_file
    samtools view -bT ${ref_genome_path} -o  ${combine_bam}/${sample_name}_accepted_hits.nbam ${combine_bam}/${sample_name}_accepted_hits.nsam |tee 2>>$samtools_log_file
    samtools sort ${combine_bam}/${sample_name}_accepted_hits.nbam -o ${combine_bam}/${sample_name}_accepted_hits.sort.bam |tee 2>>$samtools_log_file

    samtools merge -f ${combine_bam}/${sample_name}_combine.bam ${combine_bam}/${sample_name}_accepted_hits.sort.bam ${bwa_map}/${sample_name}_unmapped.sort.bam 2>>$samtools_log_file
    samtools flagstat  ${combine_bam}/${sample_name}_combine.bam |tee >${combine_bam}/${sample_name}_combine_flagstat.log 
    #samtools index ${combine_bam}/${sample_name}_combine.bam  2>>$samtools_log_file
    samtools index ${combine_bam}/${sample_name}_combine.bam || samtools index -c ${combine_bam}/${sample_name}_combine.bam |tee 2>>$samtools_log_file ###### Sat Oct 8 13:54:25 CST 2022
    need_rm_internal_file_list=(${combine_bam}/${sample_name}_accepted_hits.nsam ${combine_bam}/${sample_name}_accepted_hits.nbam ${combine_bam}/${sample_name}_accepted_hits.sort.bam ${bwa_map}/${sample_name}_bwa.header ${bwa_map}/${sample_name}_unmapped.nbam ${bwa_map}/${sample_name}_unmapped.nbam ${bwa_map}/${sample_name}_unmapped.sort.bam ${bwa_map}/${sample_name}_bwa_unique_mis6_mapq0.sam ${bwa_map}/${sample_name}_bwa_mapped.sam ${HISAT_map}/${sample_name}_HISAT2_mapped.sam ${HISAT_map}/${sample_name}_unique_mismatch2.sam ${HISAT_map}/${sample_name}_unmapped_1.fastq.gz ${HISAT_map}/${sample_name}_unmapped_2.fastq.gz ${HISAT_map}/${sample_name}_HISAT2_unmapped.fastq.gz  ${HISAT_map}/${sample_name}_hisat2_unmap.readid ${HISAT_map}/${sample_name}_un_conc_2.fastq.gz ${HISAT_map}/${sample_name}_un_conc_1.fastq.gz) #${HISAT_map}/${sample_name}_HISAT2_unmapped.bam

    echo "#####[`date`]${sample_name} satrt rm internal files"
    for need_rm_internal_file in ${need_rm_internal_file_list[@]}
    do
        test -e $need_rm_internal_file && rm $need_rm_internal_file
    done
    
    echo "#####[`date`]${sample_name} end rm internal files"
}
STEP2_sam_fine_tune(){
    #sam_fine_tune $bam_file $fine_tune_path $bam_file_basename_prefix
    bam_file=$1
    fine_tune_path=$2
    bam_file_basename_prefix=$3  ###### Tue Oct 8 19:27:03 CST 2019 fzc added
    test -d $fine_tune_path || mkdir -p $fine_tune_path
    if [ "$bam_file_basename_prefix" == "" ];then
        bam_file_basename_prefix=`basename $bam_file|awk -F".bam" '{print $1}'`
    fi
    bam_file_prefix=${fine_tune_path}/${bam_file_basename_prefix}

    
    knownSNP_for_BQSR=$dbSNP_all 
    
    picard AddOrReplaceReadGroups I=${bam_file} O=${bam_file_prefix}_rgadd.bam SO=coordinate RGID=1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20 

    picard MarkDuplicates I=${bam_file_prefix}_rgadd.bam O=${bam_file_prefix}_rgadd_dedupped.bam CREATE_INDEX=false VALIDATION_STRINGENCY=SILENT M=${bam_file_prefix}_rgadd_MarkDuplicates_output.metrics 
    samtools index ${bam_file_prefix}_rgadd_dedupped.bam ||{
    test -s ${bam_file_prefix}_rgadd_dedupped.bai && rm ${bam_file_prefix}_rgadd_dedupped.bai
    samtools index -c ${bam_file_prefix}_rgadd_dedupped.bam ###### Sun Oct 9 13:35:07 CST 2022
    }


    gatk SplitNCigarReads -R ${ref_genome_path} -I ${bam_file_prefix}_rgadd_dedupped.bam -O ${bam_file_prefix}_rgadd_dedupped_split.bam --create-output-bam-index false #  change create-output-bam-index from true (default) to false, because index bam failed when contig longer than 512M. ###### Sun Oct 9 13:36:57 CST 2022
    samtools index ${bam_file_prefix}_rgadd_dedupped_split.bam ||{
    test -s ${bam_file_prefix}_rgadd_dedupped_split.bai && rm ${bam_file_prefix}_rgadd_dedupped_split.bai
    samtools index -c ${bam_file_prefix}_rgadd_dedupped_split.bam ###### Sun Oct 9 13:35:07 CST 2022
    }

    gatk BaseRecalibrator -R ${ref_genome_path}  -I ${bam_file_prefix}_rgadd_dedupped_split.bam --known-sites ${knownSNP_for_BQSR}  -O ${bam_file_prefix}_rgadd_dedupped_split_recal_gatk.grv 
    gatk ApplyBQSR  -R ${ref_genome_path}  -I ${bam_file_prefix}_rgadd_dedupped_split.bam  --bqsr-recal-file ${bam_file_prefix}_rgadd_dedupped_split_recal_gatk.grv -O ${bam_file_prefix}_rgadd_dedupped_split_recal.bam --create-output-bam-index false # change CREATE_INDEX from true to false, because index bam failed when contig longer than 512M. ###### Sun Oct 9 13:36:57 CST 2022
    #gatk ApplyBQSR  -R ${ref_genome_path} -I ${bam_file_prefix}_rgadd_dedupped_split.bam  --bqsr-recal-file ${bam_file_prefix}_rgadd_dedupped_split_recal_gatk.grv -O ${bam_file_prefix}_rgadd_dedupped_split_recal.bam
    samtools index ${bam_file_prefix}_rgadd_dedupped_split_recal.bam ||{
    test -s ${bam_file_prefix}_rgadd_dedupped_split_recal.bai && rm ${bam_file_prefix}_rgadd_dedupped_split_recal.bai
    samtools index -c ${bam_file_prefix}_rgadd_dedupped_split_recal.bam ###### Sun Oct 9 13:35:07 CST 2022
    }

    samtools flagstat ${bam_file_prefix}_rgadd_dedupped_split_recal.bam >${bam_file_prefix}_rgadd_dedupped_split_recal_flagstat.log & 

    tmp_remove_file_list=(${bam_file_prefix}_rgadd.bam ${bam_file_prefix}_rgadd_dedupped.bam ${bam_file_prefix}_rgadd_dedupped.bam.bai ${bam_file_prefix}_rgadd_dedupped.bai ${bam_file_prefix}_rgadd_dedupped_split.bam ${bam_file_prefix}_rgadd_dedupped_split.bai ${bam_file_prefix}_rgadd_dedupped_split.bam.bai ${bam_file_prefix}_rgadd_dedupped_split_recal_gatk.grv ${bam_file_prefix}_rgadd_MarkDuplicates_output.metrics)

    for file1 in ${tmp_remove_file_list[@]}
    do
    test -e $file1 && rm $file1
    done
}
STEP3_Variant_calling(){
    #bmc_Variant_filter $bam_file $tag $work_path $sample_name
    #Used for: Calling and filtering, Calling mutation use samtools mpileup, filtering same as gatk_Variant_filter.
    #tag="_recal"
    local from_bam=$1
    local tag=$2
    local MutationCalling_wp=$3
    local sample_name=$4
    local ref_genome=$5
    test -d $MutationCalling_wp || mkdir -p $MutationCalling_wp
    #####All editing sites (BQ20; overhang6; ES95)
    
    perl ${DEMINING_path}/src/npileupBam_sszhu.pl.backup -i $from_bam -s ${ref_genome_path} -depth 10000000 -minBQ 20 -o 6 -HPB 0 -eSignal 0.95 -v 0 --cRatio 0 >${MutationCalling_wp}/${sample_name}.BQ20o6ES95v0${tag}
    perl -ane 'print "$F[0]:$F[1]\t",join("\t",@F[2..$#F]),"\n"' ${MutationCalling_wp}/${sample_name}.BQ20o6ES95v0${tag} >${MutationCalling_wp}/${sample_name}.BQ20o6ES95v0.new${tag}
    
    python ${DEMINING_path}/src/npileup_to_ES_variant_allvariants.py ${MutationCalling_wp}/${sample_name}.BQ20o6ES95v0.new${tag} 0.95 2 >${MutationCalling_wp}/${sample_name}.BQ20o6ES95v2.allvariants${tag}
    test -e  ${MutationCalling_wp}/${sample_name}.BQ20o6ES95v0${tag} && rm  ${MutationCalling_wp}/${sample_name}.BQ20o6ES95v0${tag}
    rm ${MutationCalling_wp}/${sample_name}.BQ20o6ES95v0.new${tag} &
}
STEP4_DeepDDR_predict(){
    # RADAR_out="${MutationCalling_wp}/${sample_name}.BQ20o6ES95v2.allvariants_recal"
    # BamPath="${project_path}/hg38/3-0-0Combine_bam/${sample_name}_combine_rgadd_dedupped_split_recal.bam"
    ref_genome=$1
    sample_name=$2
    RADAR_out=$3
    BamPath=$4
    predict_out_file=$5
    CMCconst_wp=$6 #/5-0-0CMC_construction
    Treads=$7
    memkdir $CMCconst_wp
    Patch=False
    length0=601
   
    if [ "$ref_genome" == "hg38" ];then
    DeepDDR_path="${DEMINING_path}/DeepDDR"
    else 
    DeepDDR_path="${DEMINING_path}/DeepDDR_Transfer"
    fi

    CMC_file="${CMCconst_wp}/2-1-1-CMC_ER5HPB3MR2_Length-${length0}_Sample-${sample_name}.txt.gz"
    echo "python3 ${DEMINING_path}/src/CMC_construction_func.py $sample_name $RADAR_out $BamPath $CMCconst_wp $Treads $Patch $ref_genome_path $DeepDDR_path $CMC_file $predict_out_file"
    python3 ${DEMINING_path}/src/CMC_construction_func.py $sample_name $RADAR_out $BamPath $CMCconst_wp $Treads $Patch $ref_genome_path $DeepDDR_path $CMC_file $predict_out_file


    


}

memkdir(){  
    for path in $@  
    do  
    test -d $path ||mkdir -p $path  
    done  
}  
DEMINING_src_path=$(cd "$(dirname "$0")";pwd)
DEMINING_path=$(dirname $DEMINING_src_path)
ref_genome_path=$genome_fasta
ref_genome=$genome_build_version


"$@"
