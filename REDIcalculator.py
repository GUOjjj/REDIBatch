#!/bin/bash
#!/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/software/miniconda3/envs/py37/bin/python3.7


import gzip
from sys import argv
import re
import json

## input content : organism, srp, srr
organism, srp, srr = argv[1],argv[2],argv[3]
## output content :   0_redi_results_dic.json
                    # 1_overall_redilevel.txt
                    # 2_known_redi_table.txt
                    # 3_geneAA_redi_level.txt
                    # 4_contribution_type.txt
                    # 5_contribution_func.txt
# INPUT
## 1 outTable_REDI_url examples : /hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/human/SRP334277/SRR15629334/REDI_table/second/DnaRna_REDI/outTable_REDI.gz
outTable_REDI_url= "/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/{organism}/{srp}/{srr}/REDI_table/second/DnaRna_REDI/outTable_REDI.gz".format(organism=organism,srp=srp,srr=srr)
## 2 DB
if organism=='human':
    ref_no = 'hg19'
    redi_db_url="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/database/REDIportal/REDIportalV2.0/AnnotatedPositions/TABLE1_{ref_no}.txt.gz".format(ref_no=ref_no)
elif organism=='mouse':
    ref_no = 'mm10'
    redi_db_url="/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/database/REDIportal/REDIportalV2.0/AnnotatedPositions/TABLE1_{ref_no}.txt.gz".format(ref_no=ref_no)
# output
## 0_json_dictionary
redi_results_dic_json = "/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/{organism}/{srp}/{srr}/REDI_table/second/DnaRna_REDI/0_redi_results_dic.json".format(organism=organism,srp=srp,srr=srr)
## 1_overall_redilevel.txt
overall_redilevel_url = "/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/{organism}/{srp}/{srr}/REDI_table/second/DnaRna_REDI/1_overall_redilevel.txt".format(organism=organism,srp=srp,srr=srr)
## 2_known_redi_table.txt
known_redi_table_url = "/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/{organism}/{srp}/{srr}/REDI_table/second/DnaRna_REDI/2_known_redi_table.txt".format(organism=organism,srp=srp,srr=srr)
## 3_geneAA_redi_level.txt
geneAA_index_table_url = "/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/{organism}/{srp}/{srr}/REDI_table/second/DnaRna_REDI/3_geneAA_index_table.txt".format(organism=organism,srp=srp,srr=srr)
## 4_contribution_type.txt
contribution_type_url = "/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/{organism}/{srp}/{srr}/REDI_table/second/DnaRna_REDI/4_contribution_type.txt".format(organism=organism,srp=srp,srr=srr)
## 5_contribution_func.txt
contribution_func_url = "/hwfssz5/ST_EARTH/P18H19700N0356/guosijia/project/RNAediting/results/{organism}/{srp}/{srr}/REDI_table/second/DnaRna_REDI/5_contribution_func.txt".format(organism=organism,srp=srp,srr=srr)




def getGenename(AAchange_content):
    b=AAchange_content.split(',')
    b=list(set(b))
    if len(b)>1:
        #Q/R
        res=re.search(r'.*:p.(.)\d*(.)',b[0])
        c='({Q}/{R})'.format(Q=res.group(1),R=res.group(2))
    else:
        #Q67R
        res=re.search(r'.*:p.(.*)',b[0])
        c='({AAchange})'.format(AAchange=res.group(1))
    return c

def saveAsJson(jf,dic):
    with open(jf,'wt',encoding='utf-8')as f:
        json.dump(dic,f)

def main():
    redi_table_lines=[]
    try:
        with gzip.open(outTable_REDI_url, "rb") as redi_table:
            redi_table_content = redi_table.read()
    except FileNotFoundError:
        print("Sorry, file : {outTable_REDI_url} does not exist.".format(outTable_REDI_url=outTable_REDI_url))
    else:
        redi_table_content = str(redi_table_content, encoding="utf-8")
        redi_table_lines = redi_table_content.split('\n')
        with gzip.open(redi_db_url,"rb") as redi_db:
            db_content = redi_db.read()
        db_content = str(db_content, encoding="utf-8")
        db_lines = db_content.split('\n')
        ### count outTable
        # N_known_sides='Numeber of known editing sites:{N}'.format(N=N)
        N_known = 0
        redi_results_dic = {}
        # redi outTable
        ### Region    Position    Reference    Strand    Coverage-q30    MeanQ    BaseCount[A,C,G,T]    AllSubs Frequency    gCoverage-q30    gMeanQ    gBaseCount[A,C,G,T]    gAllSubs    gFrequency
        ### chr16_KE145796_random    8550    T    2    22    38.14    [0, 8, 0, 14]    TC    0.36    -    -    -    -    -
        for i in range(1,len(redi_table_lines)-1):
        # for i in range(1,3):
            columns = redi_table_lines[i].split('\t')
            # current dictionary: redi_results_dic
            # {
            #   'chr1_18329':'['AG',[4,0,8,0]]'
            # }
            redi_region = columns[0].split('_')[0]
            redi_pos = columns[1]
            redi_results_key = redi_region + "_" + redi_pos
            redi_variation = columns[7]
            redi_basecounts = columns[6]
            redi_results_value = []
            redi_results_value.extend([redi_variation,redi_basecounts])
            redi_results_dic[redi_results_key] = redi_results_value
        N_total_sites = len(redi_results_dic)
        # redi DB
        ### Region  Position        Ref     Ed      Strand  db      type    dbsnp   repeat  Func.wgEncodeGencodeBasicV34lift37      Gene.wgEncodeGencodeBasicV34lift37      GeneDetail.wgEncodeGencodeBasicV34lift37        ExonicFunc.wgEncodeGencodeBasicV34lift37        AAChange.wgEncodeGencodeBasicV34lift37  Func.refGene    Gene.refGene    GeneDetail.refGene      ExonicFunc.refGene      AAChange.refGene        Func.knownGene  Gene.knownGene  GeneDetail.knownGene    ExonicFunc.knownGene    AAChange.knownGene      phastConsElements100way
        ### chr1    87158   T       C       -       A       ALU     -       SINE/AluJo      intergenic      OR4F5;AL627309.3        -       -       intergenic      OR4F5;LOC729737 -       -       intergenic      OR4F5;LOC729737 -       -       -
        for i in range(1, len(db_lines) - 1):
        # for i in range(1, 3):
            db_columns = db_lines[i].split('\t')
            # dictionary
            # {
            #   'chr1_18329':'['AG','[4,0,8,0]','ALU','SINE/AluJo','intergenic','OR4F5;LOC729737',12,8,0.67]',
            # nonsynonymous RNA editing sites (AAchange) in dic[key][5] have brackets
            # e.g. ZNF189(Q67R)
            #   'chr9_104169775':'['AG','[5,0,6,0]','ALU','SINE/AluJb','exonic','ZNF189(Q67R)',11,6,0.55]'
            # }
            db_region = db_columns[0].split('_')[0]
            db_pos = db_columns[1]
            db_side = db_region + "_" + db_pos
            db_variation = db_columns[2] + db_columns[3]
            db_type = db_columns[6]
            db_repeat = db_columns[8]
            db_func = db_columns[13]
            if db_side in redi_results_dic.keys():
                new_results = redi_results_dic[db_side]
                if db_variation in new_results:
                    # db_columns[14]: Gene.refGene
                    # content: OR4F5;LOC729737
                    # db_columns[15]: ExonicFunc.refGene
                    # content: nonsynonymous SNV
                    # db_columns[16]: AAChange.refGene
                    # content: CAVIN1:NM_012232:exon1:c.A200G:p.Q67R (0 or more, comma separated)
                    if db_columns[15] == 'nonsynonymous SNV':
                        # ZNF189(Q67R)
                        db_gene = db_columns[14] + getGenename(db_columns[16])
                    else:
                        # OR4F5;LOC729737
                        db_gene = db_columns[14]
                    new_results.extend([db_type, db_repeat, db_func, db_gene])
                    redi_results_dic[db_side] = new_results
                    N_known += 1
        n_total_reads = 0
        n_edited_reads = 0
        for known_side_info in redi_results_dic.values():
            if len(known_side_info) > 2:
                redi_basecounts = known_side_info[1]
                A = int(re.split('[\[,\]]', redi_basecounts)[1])
                G = int(re.split('[\[,\]]', redi_basecounts)[3])
                T = int(re.split('[\[,\]]', redi_basecounts)[4])
                C = int(re.split('[\[,\]]', redi_basecounts)[2])
                if known_side_info[0].startswith('A'):
                    edited_reads = G
                    total_reads = A + G
                    redi_level_per_side = edited_reads / total_reads
                    known_side_info.extend([edited_reads, total_reads, redi_level_per_side])
                elif known_side_info[0].startswith('T'):
                    edited_reads = C
                    total_reads = T + C
                    redi_level_per_side = edited_reads / total_reads
                    known_side_info.extend([edited_reads, total_reads, redi_level_per_side])
                n_edited_reads += int(known_side_info[6])
                n_total_reads += int(known_side_info[7])
        redi_level = n_edited_reads / n_total_reads
        # 1. overall_redilevel.txt
        # format： 0.25 5452 1213
        # redi_level：Overall RNA editing level；
        # N_known：number of known RNA editing sites
        # N_unknown：number of unknown RNA editing sites
        N_unknown = N_total_sites - N_known
        with open(overall_redilevel_url, 'w') as known_redi_level:
            known_redi_level.write(str(redi_level) + '\t' + str(N_known) + '\t' + str(N_unknown))
        # 2.known_redi_table.txt
        # 将已知位点单独输出成table
        # redi_results_dic的value长度>2：注释好的已知位点
        # redi_resutts-dic的value长度=2：未知的位点，没有注释
        # known_redi_table.txt
        with open(known_redi_table_url, 'w') as known_txt:
            for key in redi_results_dic.keys():
                if len(redi_results_dic[key]) > 2:
                    known_txt.write(key + '\t')
                    for value in redi_results_dic[key]:
                        known_txt.write(str(value) + '\t')
                    known_txt.write('\n')
        # 3.geneAA_redilevel.txt
        # 将geneAAchange的位点单独输出
        # format： sites  geneAA  redilevel
        # {
        #   'chr9_104169775':['ZNF189(Q67R)',0.55,[5,0,6,0]]
        # }
        # 原内容
        # {
        #   'chr1_18329':'['AG',[4,0,8,0],'ALU','SINE/AluJo','intergenic','OR4F5;LOC729737',12,8,0.67]',
        #   # 是AAchange的dic[key][5]gene列有括号
        #   'chr9_104169775':'['AG',[5,0,6,0],'ALU','SINE/AluJb','exonic','ZNF189(Q67R)',11,6,0.55]'
        # }
        with open(geneAA_index_table_url, 'w') as known_txt:
            for key, value in redi_results_dic.items():
                if len(value) > 2:
                    res = re.search(r'.*[(].*[)]', value[5])
                    try:
                        gene_AA = res.group(0)
                    except:
                        pass
                    else:
                        if gene_AA == value[5]:
                            known_txt.write(key + '\t' + value[5] + '\t' + str(value[8]) + '\t' + str(value[1]) + '\n')
        # 4. contribution_type.txt
        # Count Type: ALU、L1-Line、others in known RNA editing sites
        types_count = {}
        for value in redi_results_dic.values():
            if len(value) > 2:
                type = value[2] + ':' + value[3]
                if type in types_count.keys():
                    types_count[type] += 1
                else:
                    types_count[type] = 1
        with open(contribution_type_url, 'w') as contribution_type_txt:
            for key in types_count:
                contribution_type_txt.write(key + '\t' + str(types_count[key]) + '\n')
        # 5. contribution_func.txt
        # Conut Func: intergenic, exonic, intronic, UTR3…… in known RNA editing sites
        funcs_count = {}
        for value in redi_results_dic.values():
            if len(value) > 2:
                func = value[4]
                if func in funcs_count.keys():
                    funcs_count[func] += 1
                else:
                    funcs_count[func] = 1
        with open(contribution_func_url, 'w') as contribution_func_txt:
            for key in funcs_count:
                contribution_func_txt.write(key + '\t' + str(funcs_count[key]) + '\n')

        # 0. redi_results_dic.json
        ## saveDicAsJson
        saveAsJson(redi_results_dic_json, redi_results_dic)

if __name__ == "__main__":
   main()
