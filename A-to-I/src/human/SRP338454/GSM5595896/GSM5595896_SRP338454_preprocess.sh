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
# redML
redML="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/software/RED-ML/bin/red_ML.pl"

# REDI_method="redML"
# raw data url
sra_url="/hwfssz1/pub/database/ftp.ncbi.nih.gov/sra/sra-instant/reads/ByStudy/sra/"
# reference url
reference_url="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/database/reference"
# reference genome
ref_prefix="${reference_url}/reference_genome/${ref_no}/${ref_no}"
# dbsnp
dbsnp="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/database/snpdb/${ref_no}/00-All.vcf"
# simpleRepeat
simpleRepeat="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/database/simpleRepeat/${ref_no}/simpleRepeat.merge.bed"
# Alu
alu_bed="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/database/Alu/${ref_no}/${ref_no}.alu.bed"
# REDML threhold
threhold=0.65


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
REDI_table_url="${results}/REDI_table/${REDI_method}"
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

#redML
redML(){
    perl ${redML} --rnabam ${rnabam} --reference ${ref_prefix}.fa --dbsnp ${dbsnp} --simpleRepeat ${simpleRepeat} --alu ${alu_bed} --outdir ${REDI_table_url} --p ${threhold}
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
    redML
    
else
    redML
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

