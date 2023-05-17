#!/bin/bash
#!/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/software/miniconda3/envs/py37/bin/python3.7


batch_no=$1
project_url=$2

dataset_url="${project_url}/data/${batch_no}"
results_url="${project_url}/results"
check_csv="${results_url}/${organism}_${batch_no}.csv"
jobs_stat_txt="${batch_no}_jobs_stat.txt"

echoAll(){
    echo -n ${SRP},${SampleID},${size_of_fq} >> ${check_csv}
    for txt in $(sed -e '2,3d;/>1/d;s/[ ]*[0-9]*[ ]*(//;s/reads;/ /;s/overall/ /;s/%//;s/)//' ${alignment_results_file} | awk '{print $1}')
    do
        echo -n ,$txt >> ${check_csv}
    done
    echo -n ,${uniquely_reads},${index},${n_known},${n_unknown},${ADAR},${ADARB1},${ADARB2},${AG},${AI} >> ${check_csv}
    echo >> ${check_csv}
}
human_getAdarExp(){
    if test -z "$(grep -w ADAR ${ADAR_exp})"; 
    then
        ADAR='NA'
    else
        ADAR=$(grep -w ADAR ${ADAR_exp} | awk '{print $6}')
    fi
    if test -z "$(grep -w ADARB1 ${ADAR_exp})"; 
    then
        ADARB1='NA'
    else
        ADARB1=$(grep -w ADARB1 ${ADAR_exp} | awk '{print $6}')
    fi
    if test -z "$(grep -w ADARB2 ${ADAR_exp})"; 
    then
        ADARB2='NA'
    else
        ADARB2=$(grep -w ADARB2 ${ADAR_exp} | awk '{print $6}')
    fi
}
mouse_getAdarExp(){
    if test -z "$(grep -w Adar ${ADAR_exp})"; 
    then
        ADAR='NA'
    else
        ADAR=$(grep -w Adar ${ADAR_exp} | awk '{print $6}')
    fi
    if test -z "$(grep -w Adarb1 ${ADAR_exp})"; 
    then
        ADARB1='NA'
    else
        ADARB1=$(grep -w Adarb1 ${ADAR_exp} | awk '{print $6}')
    fi
    if test -z "$(grep -w Adarb2 ${ADAR_exp})"; 
    then
        ADARB2='NA'
    else
        ADARB2=$(grep -w Adarb2 ${ADAR_exp} | awk '{print $6}')
    fi
}
getCheckCsv(){
    if test -e "${results_url}/*${batch_no}.csv"
    then
        $(rm ${results_url}/*${batch_no}.csv)
    else 
        echo no_necessary_to_remove
    fi                         
    for organism in `ls ${dataset_url}`
    do
        echo SRP,SampleID,size_of_fq,reads,1_time_rate,overall_alignment_rate,uniquely_reads,index,n_known,n_unknown,ADAR,ADARB1,ADARB2,AG,AI >> ${check_csv}
        for SRP in $(ls "${dataset_url}/${organism}")
        do
            # SraRunTable.csv url samples:
            # /hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/data/batch1_cell_data/human/SRP362070/SraRunTable.csv
            SraRunTable_csv=${dataset_url}/${organism}/${SRP}/SraRunTable.csv
            #get SampleID
            for SampleID in `awk -F, '{if (NR>1) {print $1}}' ${SraRunTable_csv}|cut -d '_' -f 1|sort|uniq`
            do
                ## 1.size of fq file
                ## /hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/human/SRP334277/${SampleID}/fq_file
                fq_file="${results_url}/${organism}/${SRP}/${SampleID}/fq_file"
                size_of_fq=$(du -sh ${fq_file} | awk '{print $1}')
                ## 2.results of bam file
                ## /hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/human/SRP334277/${SampleID}/bam_file/${SampleID}4.Map2GenomeStat.xls
                bam_file="${results_url}/${organism}/${SRP}/${SampleID}/bam_file"
                alignment_results_file="${bam_file}/${SampleID}.Map2GenomeStat.xls"
                ## 3.uniquely_reads
                uniquely_reads=$(sed -e '1,2d;/0[ ]times/d;/>1/d;/overall/d' ${alignment_results_file} | awk '{print $1}')
                ## 4.index ; n_known ; n_unknown
                ## /hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/human/SRP334277/${SampleID}/REDI_table/second/DnaRna_REDI/overall_redilevel.txt
                overall_redilevel="${results_url}/${organism}/${SRP}/${SampleID}/REDI_table/second/DnaRna_REDI/1_overall_redilevel.txt"
                index=$(awk '{print $1}' ${overall_redilevel}) && n_known=$(awk '{print $2}' ${overall_redilevel}) && n_unknown=$(awk '{print $3}' ${overall_redilevel})
                ## 5.ADAR,ADARB1,ADARB2
                ##  /hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/human/SRP334277/${SampleID}/ADAR_exp/${SampleID}.genes.results
                ADAR_exp="${results_url}/${organism}/${SRP}/${SampleID}/ADAR_exp/${SampleID}.genes.results"
                ${organism}_getAdarExp && echoAll 
            done
        done    
    done
}


$(qstat > ${results_url}/${jobs_stat_txt})
$(sed -i 1,2d ${jobs_stat_txt})
