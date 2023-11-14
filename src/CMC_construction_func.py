## ###### Sat Jun 25 10:11:24 CST 2022 add mm10,fzc

import pandas as pd
import numpy as np
import multiprocessing as mp
import os,sys,time,gzip,pysam
from tensorflow.keras import models

print(" [pid:%s; %s ], Start!"%(os.getpid(),time.ctime()))
def memkdir(Path):
    if not os.path.exists(Path):os.makedirs(Path)

# pysam MSA

def Extract_Aln_per_Site_601bp_v6(BamPath,Site,outpath=None):
    # print(time.ctime(),"Extract_Aln_per_Site_601bp_v4 Starting1")
    Chrn=Site.split(":")[0]
    Start_0base=int(Site.split(":")[1])-1-300
    End_0base=int(Site.split(":")[1])-1+300+1
    Start_1base=int(Site.split(":")[1])-300
    End_1base=int(Site.split(":")[1])+300
    if Start_0base<=0:return("") #Extract_Aln_per_Site_601bp_v6 ###### Tue Jun 7 14:17:55 CST 2022 
    t1={} #n-read-base
    samfile = pysam.AlignmentFile(BamPath, "rb")
    for pileupcolumn in samfile.pileup(Chrn, Start_0base, End_0base,flag_filter=1024,truncate=True):
        t1.setdefault(pileupcolumn.pos,{})
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                readName=pileupread.alignment.query_name
                readSeq=pileupread.alignment.query_sequence
                t1[pileupcolumn.pos][readName]=readSeq[pileupread.query_position]
    # print(time.ctime(),"Formating")
    df1=pd.DataFrame.from_dict(t1).fillna('-')
    Pos_list=range(Start_1base,End_1base+1)
    df_empty=pd.DataFrame([['-']*len(Pos_list)]*df1.shape[0],columns=Pos_list,index=df1.index.values)
    remain_columns=[x for x in Pos_list if x not in df1.columns.values]
    df1_compeleted=pd.concat([df1,df_empty.loc[:,remain_columns]],axis=1).loc[:,Pos_list]
    array1=df1_compeleted.to_numpy()
    readSeq_list=[''.join(x) for x in array1]
    # print(time.ctime(),"ref_fa")
    if pysam_version_gt084 ==True:
        ref_seq="".join([x for x in pysam.faidx(ref_genome_path, "%s:%s-%s"%(Chrn,Start_1base,End_1base)).strip().split("\n") if ">" not in x])
    else:
        ref_seq="".join([x.strip() for x in pysam.faidx(ref_genome_path, "%s:%s-%s"%(Chrn,Start_1base,End_1base)) if ">" not in x])
    result_MSA_list=[ref_seq]+readSeq_list
    # print(time.ctime(),"writing")
    if outpath!=None:
        outfile=os.path.join(outpath,Site+"_MSA.txt")
        with open(outfile,'w') as file1:
            file1.write("\n".join(result_MSA_list))
    else:
        ti=">"+Site+'\n'+"\n".join(result_MSA_list)+"\n"
        return(ti)

def MSA_of_dataset_per_Sample_v2_601bp_BigFile_ER5HPB3MR2_v3(df1,MSA_OneSample_Summary_file,Treads):
    print(SampleName,time.ctime()," Starting")
    ok_token1=MSA_OneSample_Summary_file+"_MSA_extract_ok"
    running_token1=MSA_OneSample_Summary_file+"_MSA_extract_running"
    if Patch=="True" and (os.path.exists(ok_token1) or os.path.exists(running_token1)):
       return(None)
    with open(running_token1,'w')as file1:
        file1.write("Starting_at_%s"%time.ctime())
    with mp.Pool(Treads) as pool:
        MSA_string_list=pool.starmap(Extract_Aln_per_Site_601bp_v6,[(BamPath,Site) for Site in df1.Site.unique()])
    with gzip.open(MSA_OneSample_Summary_file,"wt") as file1:
        file1.write("\n".join(MSA_string_list))
    with open(ok_token1,'w')as file1:
        file1.write("Ending_at_%s"%time.ctime())
    print(SampleName,time.ctime()," Ending")


# Cal prob of base pair

### fine-tune length_BigFile_ER5HPB3MR2

#### reads saved MSA

### reads saved dict
def read_MSA_to_dict_FromSummary_601bp_BigFile_ER5HPB3MR2(SampleName,MSA_OneSample_Summary_file):
    t1={}
    t1.setdefault(SampleName,{})
    with gzip.open(MSA_OneSample_Summary_file,'rt') as file1:
        MSA=[]
        for i in file1:
            if i.startswith(">"):
                if len(MSA)!=0:
                    t1[SampleName][Site]=MSA
                Site=i[1:-1]
                MSA=[]
            elif i.strip().upper()!="":
                MSA.append(i.strip().upper())
        t1[SampleName][Site]=MSA
    return(t1)
def MSA_prob_Cal_Save_to_df_v2_BigFile_ER5HPB3MR2(SampleName,Site,MSA_string,MSA_length=501):
#     Site,SampleName="chr1:16458508,Human_Heart_8Week_PostConception_Female_rep1".split(",")
#     global t1_601bp_BigFile_ER5HPB3MR2_OneSample
    MSA_prob_length=MSA_length-1
    df_empty=pd.DataFrame([[0]*MSA_prob_length]*len(index_all_v2),index=index_all_v2)
    half_length=int(MSA_prob_length/2)
    MSA_list=[list(x.strip().upper()) for x in MSA_string]
    array1=np.array(MSA_list)
    start_token=300-half_length ### 300 was determined by 601bp MSA
    try:
        nuls_list=[[array1[row,col]+array1[row,300] for row in range(array1.shape[0]) if "-" not in array1[row,col]+array1[row,300]  and "N" not in array1[row,col]+array1[row,300]]for col in range(array1.shape[1]) if col!=300 and col >=300-half_length and col <=300+half_length]   ### 300 was determined by 601bp MSA
    except:
        if int(t_chrom_size[Site.split(':')[0]])-int(Site.split(':')[1])<300:
            return(None)
        print(SampleName,Site)
    df_nul_occurrence=pd.DataFrame(nuls_list).T.apply(pd.value_counts,normalize=True).fillna(0)
    ### flat index
    remain_index=[x for x in index_all_v2 if x not in df_nul_occurrence.index.values]
    df_nul_occurrence_flat=pd.concat([df_nul_occurrence,df_empty.loc[remain_index]]).loc[index_all_v2]
    df_nul_occurrence_flat.loc[:,"SampleName Site".split()]=SampleName,Site
    return(df_nul_occurrence_flat)
def Read_RADAR_out_ER5HPB3MR2(RADAR_out):
    df_RADAR_out=pd.read_csv(RADAR_out,sep='\s+',names="Site Ref_nul Total_hits #a_hits %a_hits #c_hits %c_hits #g_hits %g_hits #n_hits %n_hits #t_hits %t_hits #mismatch_hits Total_hits(BQ_not_filtered) HPB(BQ_not_filtered) Max_mismatch_nul %Max_mismatch_nul_hits".split())
    if "hg38" in ref_genome_path or "rheMac10" in ref_genome_path or  "mm10" in ref_genome_path:
        df_RADAR_out_NoChrM_contig=df_RADAR_out[~df_RADAR_out.Site.str.contains("_|M")]
    else:
        df_RADAR_out_NoChrM_contig=df_RADAR_out
    df_RADAR_out_NoChrM_contig_ER5_HPB3_MutReads2=df_RADAR_out_NoChrM_contig.query("`HPB(BQ_not_filtered)`>=3&`%Max_mismatch_nul_hits`>=0.05&`Total_hits`*`%Max_mismatch_nul_hits`>=2") # 20230720,fzc, `Total_hits(BQ_not_filtered)`*`%Max_mismatch_nul_hits` to `Total_hits`*`%Max_mismatch_nul_hits`
    df_RADAR_out_NoChrM_contig_ER5_HPB3_MutReads2.loc[:,"Sample_name"]=SampleName
    return(df_RADAR_out_NoChrM_contig_ER5_HPB3_MutReads2)
def Gen_CMC_main(RADAR_out):
    global ref_genome_path,SampleName,t_chrom_size,index_all_v2
    ref_genome_fai_path=ref_genome_path+".fai"
    if not os.path.exists(ref_genome_fai_path):raise ValueError("Not Exists: "+ref_genome_fai_path)
    if ".BQ20o6ES95v2.allvariants_recal" in RADAR_out:
        df1=Read_RADAR_out_ER5HPB3MR2(RADAR_out)
    else:
        df1=pd.read_csv(RADAR_out,sep='\t',names=["Site"],usecols=[0])

    MSA_OneSample_Summary_file=os.path.join(CMCconst_wp,"1-1-1-MSA_out_pysam_601bp_ER5HPB3MR2_"+SampleName+".txt.gz")
    MSA_of_dataset_per_Sample_v2_601bp_BigFile_ER5HPB3MR2_v3(df1,MSA_OneSample_Summary_file,Treads)
    #### cal
    index_all_v2=[x+y for x in "ACTG" for y in "ACTG" ]
    with open(ref_genome_fai_path) as file1:
        t_chrom_size={x.split("\t")[0]:x.split("\t")[1]for x in file1}
    print(" [pid:%s; %s ], read_MSA_to_dict_FromSummary_601bp_BigFile_ER5HPB3MR2 Start!"%(os.getpid(),time.ctime()))
    t1_601bp_BigFile_ER5HPB3MR2_OneSample=\
    read_MSA_to_dict_FromSummary_601bp_BigFile_ER5HPB3MR2(SampleName,MSA_OneSample_Summary_file)
    print(" [pid:%s; %s ], read_MSA_to_dict_FromSummary_601bp_BigFile_ER5HPB3MR2 End!"%(os.getpid(),time.ctime()))
    t1=t1_601bp_BigFile_ER5HPB3MR2_OneSample
    with mp.Pool(Treads) as pool:
        df_CMC=pd.concat(pool.starmap(MSA_prob_Cal_Save_to_df_v2_BigFile_ER5HPB3MR2,[(row[0], row[1],t1[row[0]][row[1]],length0) for row in df1.query("Sample_name==\"%s\""%SampleName)["Sample_name Site".split()].to_numpy()if row[1] in t1[row[0]]]))
        print(" [pid:%s; %s ], MSA_prob_Cal_Save_to_df_v2_BigFile_ER5HPB3MR2 End!"%(os.getpid(),time.ctime()))
        df_CMC.to_csv(CMC_file,sep='\t')
        print(" [pid:%s; %s ], write csv End!"%(os.getpid(),time.ctime()))
        print(SampleName,time.ctime()," Ending")
    return(df_CMC)
def DeepDDR_predict(df_CMC,DeepDDR_path,predict_out_file,MSALen1=600):   
    CMC_array_Len600=df_CMC.iloc[:,:MSALen1].to_numpy().reshape(int(df_CMC.shape[0]/16),16,MSALen1) # 600*16
    MSA_len1=50
    half_cut=round((600-MSA_len1)/2)
    CMC_array_Len50=CMC_array_Len600[...,half_cut:-1*half_cut] # 50*16
    DeepDDR_model=models.load_model(DeepDDR_path)
    y_predict_2D = DeepDDR_model.predict(CMC_array_Len50)
    y_predict_1D=1-y_predict_2D[:,0]
    df_Predict=df_CMC.drop_duplicates(["Site"])[["Site"]]
    df_Predict["DeepDDR_rawScore"]=y_predict_1D
    df_Predict["DeepDDR_DM(DNA_mutation)_or_RM(RNA_mutation)"]=np.where(y_predict_1D>=0.5,"DM","RM")
    df_Predict.to_csv(predict_out_file,sep='\t',index=False)
    return(predict_out_file)

if __name__ == "__main__":
    SampleName=sys.argv[1]
    RADAR_out=sys.argv[2]
    BamPath=sys.argv[3]
    CMCconst_wp=sys.argv[4]
    Treads=int(sys.argv[5])
    Patch=sys.argv[6]
    ref_genome_path=sys.argv[7]
    DeepDDR_path=sys.argv[8]
    CMC_file=sys.argv[9]
    predict_out_file=sys.argv[10]
    length0=601
    length1=length0-1
    try:
        _=pysam.faidx(ref_genome_path, "chr22:1-601").strip()
        pysam_version_gt084=True
    except:
        pysam_version_gt084=False
    df_CMC=Gen_CMC_main(RADAR_out)
    DeepDDR_predict(df_CMC,DeepDDR_path,predict_out_file,MSALen1=length1)

        
    




