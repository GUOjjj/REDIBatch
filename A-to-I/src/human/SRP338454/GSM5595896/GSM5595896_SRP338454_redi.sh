#!/bin/bash

### parameter 
organism=human
srp=SRP338454
gsm=GSM5595896
srr=SRR16027460
ref_no=hg19
project_url=/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/A-to-I


# perl5
#!/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/software/miniconda3/envs/pl5/bin/perl
# sratoolkit3.0.0
export PATH="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/software/sratoolkit/sratoolkit.3.0.0-centos_linux64/bin:$PATH"
# hisat2-2.1.0
export PATH="/hwfssz5/ST_EARTH/P18H19700N0356/HuangJR/localshare/hisat2-2.1.0:$PATH"
# samtools-0.1.19
export PATH="/hwfssz5/ST_EARTH/P18H19700N0356/HuangJR/localshare/samtools-0.1.19:$PATH"
# blat
blat="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/software/miniconda3/pkgs/blat-36-0/bin/blat"
# py27
py27="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/software/miniconda3/envs/py27/bin/python2.7"
# py37
py37="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/software/miniconda3/envs/py37/bin/python3.7"
#bowtie2
bowtie2="/hwfssz5/ST_EARTH/P18H19700N0356/HuangJR/share/GeneExp/bin/bowtie2-2.2.5/bowtie2"
#samtools2
samtools2="/hwfssz5/ST_EARTH/P18H19700N0356/HuangJR/share/GeneExp/bin/rsem-1.2.12/sam/samtools"
#rsem
rsem="/hwfssz5/ST_EARTH/P18H19700N0356/HuangJR/share/GeneExp/bin/rsem-1.2.12/rsem-calculate-expression"
REDItoolDnaRna="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/software/REDItools/main/RE/REDItoolDnaRna.py"
readPsl="/hwfssz5/ST_EARTH/P18H19700N0356/HuangJR/localshare/REDItools/REDItools-1.3/accessory/readPsl.py"


# REDI_method="REDItools"
# raw data url
sra_url="/hwfssz1/pub/database/ftp.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/"
# reference url
reference_url="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/database/reference"
# reference genome
ref_prefix="${reference_url}/reference_genome/${ref_no}/${ref_no}"
# Indexed transcript reference file
rnaref_file="${reference_url}/reference_rna/${ref_no}/bowtie_index/${ref_no}.fa"


HISAT_uniquely="${project_url}/tool/RE/HISAT_uniquelyMapped.pl"
REDIcalculator="${project_url}/tool/RE/REDIcalculator.py"
#sra
sra_file_url="${sra_url}/SRP/${srp:0:6}/${srp}/${srr}"
sra_file_name=$(ls $sra_file_url)
sra_file="${sra_file_url}/${sra_file_name}"
#results url
results="${project_url}/results/${organism}/${srp}/${gsm}"
#fq file
fq_file_url="${results}/fq_file"
#bam file
bam_file_url="${results}/bam_file"
rnabam="${bam_file_url}/${gsm}.sort.rmdup.unique.bam"
#REDI table
REDI_table_url="${results}/REDI_table"
#ADAR_exp
ADAR_exp_url="${results}/ADAR_exp"
#runtime txt
runtime_txt="${project_url}/results/${organism}/runtime.txt"


sra2fq(){
    echo sra2fq_start && fastq-dump --gzip --split-3 "${sra_file}" --outdir ${fq_file_url} && echo sra2fq_done
}

renameFq(){
    echo renameFq_start
    if [ "$(ls -l ${fq_file_url}|grep "^-"|wc -l)" == 1 ]
    then
        $(mv * ${fq_file_url}/${srr}_2.fq.gz) && echo renameSingleFq_done:${file} to ${srr}_2.fq.gz
    else
        for file in $(ls ${fq_file_url})
        do
            new_file_name=${file/\.[0-9]/}
            if [ "${file}" == "${new_file_name}" ]
            then
                echo no_necessary_to_rename
            else    
                $(mv ${fq_file_url}/"${file}" ${fq_file_url}/"${new_file_name}") && echo renamePairedFq_done:${file} to ${new_file_name}
            fi
        done  
    fi
}

cpFq(){
    $(cp ${sra_file_url}/* ${fq_file_url})
}


# cat SRR.fq.gz file by SampleID------->SampleID.fq.gz file
catFq(){                         
    $(cat ${fq_file_url}/${srr}_1.fastq.gz >> ${fq_file_url}/${gsm}_1.fastq.gz) && echo srr1_to_gsm1_done
    $(cat ${fq_file_url}/${srr}_2.fastq.gz >> ${fq_file_url}/${gsm}_2.fastq.gz) && echo srr2_to_gsm2_done
}
                        

#single edge -U read2
bamProcess_SE(){
    echo hisat2_samtools_start && hisat2 -x ${ref_prefix} --no-mixed --no-discordant -U ${fq_file_url}/${gsm}_2.fastq.gz 2>${bam_file_url}/${gsm}.Map2GenomeStat.xls | samtools view -b -S -o ${bam_file_url}/${gsm}.bam - && echo view_done && samtools sort ${bam_file_url}/${gsm}.bam ${bam_file_url}/${gsm}.sort && echo sort_done && samtools rmdup ${bam_file_url}/${gsm}.sort.bam ${bam_file_url}/${gsm}.sort.rmdup.bam && echo rmdup_done && samtools index ${bam_file_url}/${gsm}.sort.rmdup.bam echo index_done && echo hisat2_samtools_over && echo use_perl && perl ${HISAT_uniquely} -samtools /hwfssz5/ST_EARTH/P18H19700N0356/HuangJR/localshare/samtools-0.1.19/samtools -bam ${bam_file_url}/${gsm}.sort.rmdup.bam -outdir ${bam_file_url} -suffix sort.rmdup.bam && echo perl_hisat_uniquely_mapped_over
}
#paired
bamProcess_PE(){
    echo hisat2_samtools_start && hisat2 -x ${ref_prefix} --no-mixed --no-discordant -1 ${fq_file_url}/${gsm}_1.fastq.gz -2 ${fq_file_url}/${gsm}_2.fastq.gz 2>${bam_file_url}/${gsm}.Map2GenomeStat.xls | samtools view -b -S -o ${bam_file_url}/${gsm}.bam - && echo view_done && samtools sort ${bam_file_url}/${gsm}.bam ${bam_file_url}/${gsm}.sort && echo sort_done && samtools rmdup ${bam_file_url}/${gsm}.sort.bam ${bam_file_url}/${gsm}.sort.rmdup.bam && echo rmdup_done && samtools index ${bam_file_url}/${gsm}.sort.rmdup.bam echo index_done && echo hisat2_samtools_over && echo use_perl && perl ${HISAT_uniquely} -samtools /hwfssz5/ST_EARTH/P18H19700N0356/HuangJR/localshare/samtools-0.1.19/samtools -bam ${bam_file_url}/${gsm}.sort.rmdup.bam -outdir ${bam_file_url} -suffix sort.rmdup.bam && echo perl_hisat_uniquely_mapped_over
}


REDI(){
    echo reditools_start && ${py27} ${REDItoolDnaRna} -o ${REDI_table_url} -i ${rnabam} -f ${ref_prefix}.fa -c 10,10 -q 30,30 -a 6-6 -v 3 -n 0.01 -R -z --gzip && echo reditools_end && echo ${srp}/${gsm}_precrocess_over
}

REDI_blat_REDI2(){
    echo reditools_start && ${py27} ${REDItoolDnaRna} -o ${REDI_table_url}/first -i ${rnabam} -f ${ref_prefix}.fa -c 10,10 -q 30,30 -a 6-6 -v 3 -n 0.01 -R -z --reads --addP --gzip && echo REDI2_1_over && echo touch_file && $(touch ${REDI_table_url}/first/DnaRna_REDI/reads.psl) && echo touch_file_over && echo blat_start && blat -t=dna -q=dna -minIdentity=95 ${ref_prefix}.fa ${REDI_table_url}/first/DnaRna_REDI/outReads_REDI ${REDI_table_url}/first/DnaRna_REDI/reads.psl && echo blat_done && echo read_psl_start && ${py37} ${readPsl} ${REDI_table_url}/first/DnaRna_REDI/reads.psl ${REDI_table_url}/first/DnaRna_REDI/badreads.txt && echo read_psl_over && echo REDI2_2_start && $py27 $REDItoolDnaRna -o ${REDI_table_url}/second -i ${rnabam} -f ${ref_prefix}.fa -b ${REDI_table_url}/first/DnaRna_REDI/badreads.txt -c 10,10 -q 30,30 -a 6-6 -v 3 -n 0.01 -R -z --gzip && echo REDI2_2_over
}

calculate(){
    echo claculate_start && echo overall_known_geneAA_type_func && ${py37} ${REDIcalculator} ${organism} ${srp} ${gsm} && echo calculate_end
}

#single
adarExp_SE(){
    echo adarExp_SE_start && echo bowtie2start && ${bowtie2} -q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -p 3 -k 200 -x ${rnaref_file} -U ${fq_file_url}/${gsm}_2.fastq.gz 2>${ADAR_exp_url}/${gsm}.adarMap2GenomeStat.xls | ${samtools2} view -S -b -o ${ADAR_exp_url}/${gsm}.bam - && echo bowtie2end && echo rsemstart && ${rsem} --forward-prob 0.5 -p 3 --bam --no-bam-output ${ADAR_exp_url}/${gsm}.bam ${rnaref_file} ${ADAR_exp_url}/${gsm} && echo rsemend && echo adarExp_SE_end
}
#paired
adarExp_PE(){
    echo adarExp_PE_start && echo bowtie2start && ${bowtie2} -q --phred33 --sensitive --dpad 0 --gbar 99999999 --mp 1,1 --np 1 --score-min L,0,-0.1 -I 1 -X 1000 --no-mixed --no-discordant -p 3 -k 200 -x ${rnaref_file} -1 ${fq_file_url}/${gsm}_1.fastq.gz -2 ${fq_file_url}/${gsm}_2.fastq.gz 2>${ADAR_exp_url}/${gsm}.adarMap2GenomeStat.xls | ${samtools2} view -S -b -o ${ADAR_exp_url}/${gsm}.bam - && echo bowtie2end && echo rsemstart && ${rsem} --forward-prob 0.5 -p 3 --paired-end --bam --no-bam-output ${ADAR_exp_url}/${gsm}.bam ${rnaref_file} ${ADAR_exp_url}/${gsm} && echo rsemend && echo adarExp_PE_end
}


echo start
starttime=$(date +"%s")
starttime_human=$(date +"%Y-%m-%d %H:%M:%S")
echo "${starttime_human}"

mkdir -p ${fq_file_url}
mkdir -p ${bam_file_url}
mkdir -p ${REDI_table_url}
if [[ ${sra_file} =~ ".sra" ]]
then
    REDI
    
else
    REDI
fi

endtime=$(date +"%s")
endtime_human=$(date +"%Y-%m-%d %H:%M:%S")
echo "${endtime_human}"
echo end

elapsed_time=$(($endtime-$starttime))
hour=$(( ${elapsed_time}/3600 ))
min=$(( (${elapsed_time}-${hour}*3600)/60 ))
sec=$(( ${elapsed_time}-${hour}*3600-${min}*60 ))
elapsed_time_human=$($echo ${hour}:${min}:${sec})
echo "${starttime_human}_${srp}_${gsm}_${elapsed_time_human}" >> ${runtime_txt}


