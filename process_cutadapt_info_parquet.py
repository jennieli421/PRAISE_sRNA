import csv
import pandas as pd
import sys
import os
# 解除 csv 字段长度限制
csv.field_size_limit(sys.maxsize)

# under env sports 

"""
    解析 Cutadapt --info-file 输出
    逻辑：
    1. 如果第2列是 -1，则视为无 Adapter。
    2. 如果包含 3AD 或 5AD，则视为有 Adapter。
"""
def parse_cutadapt_info(input_file):
    headers = [ "read_name", "errors", "start_coord", "end_coord", 
        "seq_left", "seq_adapter_match", "seq_right", "adapter_name", 
        "qual_left", "qual_adapter", "qual_right", "is_revcomp" ]
    results = []
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            # 使用 csv.reader 处理 tab 分隔，quotechar 设为 None 防止引号干扰
            reader = csv.reader(f, delimiter='\t', quotechar=None)
            for row in reader:
                if not row: continue # 跳过空行
                # --- 情况 1: 无 Adapter 的行 ---
                # 规则：第二列 (index 1) 是 '-1'
                if len(row) == 4 and row[1] == '-1':
                    # 格式: Read name(0), -1(1), Sequence(2), Quality(3)
                    entry = {
                        "read_name": row[0],          
                        "errors": row[1],      
                        "start_coord": None,    
                        "end_coord": None,      
                        "seq_left": "",         
                        "seq_adapter_match": "",
                        "seq_right": row[2], 
                        "adapter_name": "None", 
                        "qual_left": "",        
                        "qual_adapter": "",     
                        "qual_right": row[3], 
                        "is_revcomp": ""      
                    }
                # --- 情况 2: 有 Adapter 的行 ---
                # 规则：包含 '3AD' 或 '5AD'
                elif any(x in row for x in ['3AD', '5AD']):
                    if len(row) < 12:
                        print(f"Warning: row has less than 12 fields {row}...")
                    entry = {
                        "read_name": row[0],        
                        "errors": row[1],           
                        "start_coord": row[2],      
                        "end_coord": row[3],        
                        "seq_left": row[4],         
                        "seq_adapter_match": row[5],
                        "seq_right": row[6],        
                        "adapter_name": row[7],       
                        "qual_left": row[8],          
                        "qual_adapter": row[9],       
                        "qual_right": row[10],         
                        "is_revcomp": row[11]          
                    }
                else:
                    # 兜底：不符合上述规则的异常行，跳过
                    continue
                results.append(entry)
    except FileNotFoundError:
        print(f"Error: file not found {input_file}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: an unknown error occurred while processing the file {input_file}: {e}")
        sys.exit(1)
    df = pd.DataFrame(results, columns=headers)
    return df


# --------- 处理3AD ---------   
if os.path.exists("info_3AD_parsed.parquet"):
    df3 = pd.read_parquet("info_3AD_parsed.parquet")
else:
    df3 = parse_cutadapt_info('info_3AD.tsv')
    df3['adapter_name'].unique() # sanity check

# df3.to_parquet("info_3AD_parsed.parquet", compression="zstd", index=False)

# --------- 处理5AD ---------
if os.path.exists("info_5AD_parsed.parquet"):
    df5 = pd.read_parquet("info_5AD_parsed.parquet")
else:
    df5 = parse_cutadapt_info('info_5AD.tsv')
    df5['adapter_name'].unique() # sanity check

# df5.to_parquet("info_5AD_parsed.parquet", compression="zstd", index=False) # 保存


# ---------- 合并 3AD 和 5AD 信息 ----------
# 提取 read_name 中的基础标识（空格前的部分）
df3['read_base'] = df3['read_name'].str.split().str[0]
df5['read_base'] = df5['read_name'].str.split().str[0]


df3_tmp = df3.drop(columns=['read_name','errors', 'start_coord', 'end_coord', 'seq_adapter_match', 'qual_left','qual_adapter','qual_right','is_revcomp'])
df5_tmp = df5.drop(columns=['read_name','errors', 'start_coord', 'end_coord', 'seq_adapter_match', 'qual_left','qual_adapter','qual_right','is_revcomp'])

# 只保留同时包含两种 adapter 的 read（inner join）
merged_df = pd.merge(df3_tmp, df5_tmp, on='read_base', how='inner', suffixes=('_R1', '_R2'))

# 将 read_base 移到第一列
merged_df.insert(0, 'read_base', merged_df.pop('read_base'))

merged_df["adapter_name_R1"] = merged_df["adapter_name_R1"].astype("category")
merged_df["adapter_name_R2"] = merged_df["adapter_name_R2"].astype("category")
merged_df.to_parquet("info_bothAD_parsed.parquet", compression="zstd", index=False) # 保存


# ----------- 读取保存的文件 -------------
df = pd.read_parquet("info_bothAD_parsed.parquet")

# # 按第一列（read_base）排序（稳定排序，保留原有行顺序）
# df = df.sort_values(by='read_base', kind='stable')


# val = df['seq_left_R1'][4]
# print(f"值: '{val}'")
# print(f"类型: {type(val)}")
# print(f"是否为空字符串: {val == ''}")
# print(f"是否为 None: {val is None}")
# print(f"是否为 NaN: {pd.isna(val)}")

df[['read_base', 'adapter_name_R1', 'seq_left_R1', 'seq_right_R1']]
df[['read_base', 'adapter_name_R2', 'seq_left_R2', 'seq_right_R2']]

# view examples of seq start from inner adapter 
tmp = df[(df['seq_left_R1']=='') & (df['seq_right_R1']!='')]
tmp[['read_base', 'seq_left_R1', 'seq_right_R1']]

# view examples of seq start from inner adapter 
tmp = df[(df['seq_left_R2']=='') & (df['seq_right_R2']!='')]
tmp[['read_base', 'seq_left_R2', 'seq_right_R2']]


combination_counts = df.groupby(['adapter_name_R1', 'adapter_name_R2'], dropna=False).size().reset_index(name='count')
combination_counts

# control_1
#   adapter_name_R1 adapter_name_R2     count
# 0             3AD             5AD  12,495,406
# 1             3AD            None     57,372
# 2            None             5AD    738,099
# 3            None            None   1,311,321

# praise_1
#   adapter_name_R1 adapter_name_R2    count
# 0             3AD             5AD  8,280,567
# 1             3AD            None    58,321
# 2            None             5AD   815,525
# 3            None            None   553,196

# ---------- 分类 ----------

# DataFrame 拆分：“both”（两端都有接头）和 “others”（其他所有情况）
both_mask = (df['adapter_name_R1'] == '3AD') & (df['adapter_name_R2'] == '5AD')
both_df = df[both_mask].copy()
other_df = df[~both_mask].copy()

both_df[['read_base', 'seq_left_R2', 'seq_left_R1']]


# 对 other_df 进一步细分
na_mask = (other_df['adapter_name_R1'] == "None") & (other_df['adapter_name_R2'] == "None")
na_df = other_df[na_mask].copy()
one_df = other_df[~na_mask].copy()

# ---------- 保存 ----------
na_df.to_csv("valid_reads_noAD.tsv", sep="\t", index=False) 
na_df[['read_base']].to_csv('valid_ids_noAD.txt', index=False, header=False)

one_df.to_csv("valid_reads_singleAD.tsv", sep="\t", index=False) 
one_df[['read_base']].to_csv('valid_ids_singleAD.txt', index=False, header=False)

both_df.to_csv("valid_reads_bothAD.tsv", sep="\t", index=False) 
both_df[['read_base']].to_csv('valid_ids_bothAD.txt', index=False, header=False)






# ---------- 检查 RNA 区间 ----------
# 定义反向互补函数（处理 N）
def revcomp(seq):
    comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join(comp.get(base.upper(), 'N') for base in seq[::-1])

final_df['seq_R2_rc'] = final_df['seq_R2'].apply(revcomp)




