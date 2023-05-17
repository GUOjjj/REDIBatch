#!/bin/bash
##########################################################################################
##########################################################################################
### ATTENTION(1/2):
### The arguments need to be modified:
# 1.absolute url and project name
#   e.g. projtect_url=="/user/project/RNAediting"
projtect_url="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/A-to-I"
# 2.Enter the name of the batch data, which corresponds to the name of the [batch_no] under the path: [your_url]/[your_project_name]/data/data_info
#   e.g. batch_no="batch1_cell_data"
batch_no="batch1_cell_data"
# 3.Submit job setting properties
#   e.g. q="qsub -cwd -l vf=16.0g -l num_proc=4  -binding linear:4 -q stt.q -P aaa"
qsub="qsub -cwd -l vf=16.0g -l num_proc=4  -binding linear:4 -q st.q -P P18H19700N0356"
# 4.Set the version number of the reference sequence
#   e.g. oganismVersion=( ["human"]="hg19" ["mouse"]="mm10" )
declare -A oganismVersion
oganismVersion=( ["human"]="hg19" ["mouse"]="mm10" )
# 5.pipeline
#   e.g. echo_preprocess_sh="${projtect_url}/tool/RE/pipeline_fq2redML_test.sh"
echo_preprocess_sh="${projtect_url}/tool/RE/pipeline_fq2redML_test.sh"
##########################################################################################
##########################################################################################

dataset_url="${projtect_url}/data/data_info/${batch_no}"
src_url="${projtect_url}/src"
results_url="${projtect_url}/results"

for organism in `ls ${dataset_url}`
do 
    ref_no=${oganismVersion[${organism}]}
    for SRP in `ls "${dataset_url}/${organism}"`
    do
        # SraRunTable.csv url samples:
        # /user/project/RNAediting/data/data_info/batch4_cell_data/human/SRP362070/SraRunTable.csv
        # content:
        # Sample_Run	Organism	Cell_type	PHENOTYPE
        # GSM3617760_SRR8606234 Homo sapiens	induced neuron	wild type
        SraRunTable_csv=${dataset_url}/${organism}/${SRP}/SraRunTable.csv
        #get SampleID
        for SampleID in `awk -F, '{if (NR>1) {print $1}}' ${SraRunTable_csv}|cut -d '_' -f 1|sort|uniq`
            do
            $(mkdir -p ${src_url}/${organism}/${SRP}/${SampleID})
            $(mkdir -p ${results_url}/${organism}/${SRP}/${SampleID}/fq_file/)
            SampleID_1_fq="${results_url}/${organism}/${SRP}/${SampleID}/fq_file/${SampleID}_1.fastq"
            SampleID_2_fq="${results_url}/${organism}/${SRP}/${SampleID}/fq_file/${SampleID}_2.fastq"  
            if test -e "${SampleID_1_fq}.gz"
            then
                rm ${SampleID_1_fq}.gz && echo old_${SampleID_1_fq}.gz_removed
            fi
            if test -e "${SampleID_2_fq}.gz"
            then
                rm ${SampleID_2_fq}.gz && old_echo ${SampleID_2_fq}.gz_removed
            fi
            if test -e "${SampleID_1_fq}"
            then
                rm ${SampleID_1_fq} && echo old_${SampleID_1_fq}_removed
            fi
            if test -e "${SampleID_2_fq}"
            then
                rm ${SampleID_2_fq} && echo old_${SampleID_2_fq}_removed
            fi
            # get SampleID corresponds to SRR(s)
            for SRR in $(awk -F, '{if ($1 ~ /'"$SampleID"'/) {print $1}}' ${SraRunTable_csv}|cut -d '_' -f 2)
            do  
                ## 0. srr.fastq.gz>gsm.fq.gz>bam>redi
                ## 2. generate bam2reditable script
                ## 3. submit job
                preprocess_sh="${src_url}/${organism}/${SRP}/${SampleID}/${SampleID}_${SRP}_${preprocess}.sh"
                if test -e "${preprocess_sh}"
                then
                    rm ${preprocess_sh} && echo old_${preprocess_sh}_removed
                fi
                if test -e "${src_url}/${organism}/${SRP}/${SampleID}/${preprocess}.q*"
                then
                    rm ${src_url}/${organism}/${SRP}/${SampleID}/${preprocess}.q* && echo ${preprocess}_qOut_qError_removed
                fi
                q="${qsub} -e ${src_url}/${organism}/${SRP}/${SampleID}/preprocess.qError.txt -o ${src_url}/${organism}/${SRP}/${SampleID}/preprocess.qOut.txt"
                `sh ${echo_preprocess_sh} ${organism} ${SRP} ${SampleID} ${SRR} ${ref_no} ${projtect_url} ${preprocess_sh}` && ${q} ${preprocess_sh}    
            done  
        done
    done
done