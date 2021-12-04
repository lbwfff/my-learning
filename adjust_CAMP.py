import pandas as pd
import numpy as np
import pickle

def aa_ss_concat(aa,ss):
    if len(aa)!= len(ss):
        return 'string length error!'
    else:
        new_str = ''
        for i in range(len(aa)):
            concat_str = aa[i]+ss[i]+','
            new_str = new_str+concat_str
    final_str = new_str[:-1]
    return final_str


df_org = pd.read_csv('./ss/test.out.ss',sep='#',header = None) #SCRATCH1D预测得到的结果，将氨基酸序列做了二级结构的预测
# print(df_org)
df_org=df_org[0]
# print(df_org)
df_org.columns = ['col_1']
# print(df_org)

ss_idx=[]
seq_idx=[]
seq_idx = [2*x for x in list(range(int(df_org.shape[0]/2)))]

ss_idx = [x+1 for x in seq_idx]


df_seq = df_org.iloc[seq_idx]
print(df_seq)
df_seq.columns = ['seq_id']
df_ss = df_org.iloc[ss_idx]
print(df_ss)
df_ss.columns = ['seq_ss']

df_seq = df_seq.reset_index(drop=True)
df_ss = df_ss.reset_index(drop=True)

# join sequence & sse together
df_seq_ss = pd.merge(df_seq, df_ss,left_index=True, right_index=True)

# load id mapping file
df_id = pd.read_csv('./ss/test.fasta',sep='#',header = None) #the input asta file used for SCRATCH1D SSPro
df_id=df_id[0]
df_id.columns = ['col_1']

ss_idx=[]
seq_idx=[]
seq_idx = [2*x for x in list(range(int(df_id.shape[0]/2)))]
ss_idx = [x+1 for x in seq_idx]

# subset sequence dataframe and sse dataframe
df_seq = df_id.iloc[seq_idx]
df_seq.columns = ['seq_id']
df_ss = df_id.iloc[ss_idx]
df_ss.columns = ['seq']

df_seq = df_seq.reset_index(drop=True)
df_ss = df_ss.reset_index(drop=True)


# join sequence &  sse together
df_idx = pd.merge(df_seq, df_ss,left_index=True, right_index=True)
print(df_idx)
print(df_seq_ss)
df_output_ss = pd.merge(df_idx, df_seq_ss, left_on=['0_x'], right_on=['0_x'])#把含序列信息的fasta文件和含二级结构信息的预测文件进行了合并
df_output_ss.columns = ['seq_id','seq','seq_ss']#我不大明白作者想要做的文件到底是什么样子，这里是我猜测的栏名字
print(df_output_ss)
df_output_ss['concat_seq'] = df_output_ss.apply(lambda x: aa_ss_concat(x['seq'],x['seq_ss']),axis=1)
df_output_ss.to_csv('./test.tsv', encoding = 'utf-8', index = False, sep = '\t') # 'output_ss_filename' is the name of the output tsv you like
#居然顺滑的跑出来了

### Load Protein PSSM Files (first change the value of protein_number)
# prot_pssm_dict : key is protein sequence, value is protein PSSM Matrix
prot_pssm_dict_all={}
prot_pssm_dict={}
protein_num = 4 ### NEED TO BE CHANGED TO the total number of protein sequences#在test文件中我只包括了4个序列，所以这里用的4
for i in range(protein_num):
    filename_pssm = 'new_prot_'+str(i)+'.pssm' # need to name each individual fasta and pssm file with the same prefix
    filename_fasta = 'new_prot_'+str(i)+'.fasta'#需要把fasta拆分，然后分别得到pssm矩阵，这里fasta拆分使用R语言做的，python实在不熟悉
    prot_key = 'new_prot_'+str(i)
    pssm_line_list= []
    print('./pssm/prot_file/'+filename_fasta)
    with open('./pssm/prot_file/'+filename_fasta,'r') as f: # directory to store fasta files (single file of each protein)
        for line in f.readlines():
            prot_seq = line.strip()
    
    with open('./pssm/pssm_result/'+filename_pssm,'r') as f:  # directory to store pssm files (single file of each protein)
        for line in f.readlines()[3:-6]:
            line_list = line.strip().split(' ')
            line_list = [x for x in line_list if x!=''][2:22]
            line_list = [int(x) for x in line_list]
            if len(line_list)!=20:
                print('Error line:')
                print(line_list)
            pssm_line_list.append(line_list)
        pssm_array = np.array(pssm_line_list)
        if pssm_array.shape[1]!=20:
            print('Error!')
            print(filename_pssm)
        else:
            prot_pssm_dict_all[prot_key] = (prot_seq,pssm_array)
            prot_pssm_dict[prot_seq]=pssm_array

with open('./output_pssm_dict','wb') as f:  # 'output_pssm_dict' is the name of the output dict you like
    pickle.dump(prot_pssm_dict,f)#之后就得到了dict文件
            
