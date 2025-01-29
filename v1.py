### import libraries

import pandas as pd
import os
import xlsxwriter

### import data

comm_abnorm_reg = pd.read_csv("Common_Abnormal_Regions.csv")
calls_dir = "/facility/nfdata-omics/projects/GEDI_LPS/sample_prioritization/controlFREEC_calls/"

sample_cnv_call_path = "/facility/nfdata-omics/projects/GEDI_LPS/sample_prioritization/controlFREEC_calls/KOLF2_1J_TET2_het_C14.markdup.bam_CNVs" 
sample_cnv_call = pd.read_csv(sample_cnv_call_path, sep="\t", names=["chr", "start", "end", "copy_number", "category"], header=None)

# region overlap

def get_overlap(a, b):
   return max(0, min(a[1], b[1]) - max(a[0], b[0]))

# test one file

matching_cnv = []

for chromosome in sample_cnv_call["chr"].unique():
    if comm_abnorm_reg["Chromosome"].isin([chromosome]).any() == True:
        comm_abnorm_reg_chr = comm_abnorm_reg[comm_abnorm_reg["Chromosome"] == str(chromosome)]
        sample_cnv_call_chr = sample_cnv_call[sample_cnv_call["chr"] == str(chromosome)]
        for _, row in sample_cnv_call_chr.iterrows():
            sample_range = [row["start"], row["end"]]
            for _, row in comm_abnorm_reg_chr.iterrows():
                abnorm_reg_range = [row["start_clean"], row["end_clean"]]
                overlap = get_overlap(sample_range, abnorm_reg_range)
                if overlap > 0:
                    matching_cnv.append({
                        "comm_abnorm_reg_ID": row["comm_abnorm_reg_ID"],
                        "chromosome": chromosome,
                        "sample_start": sample_range[0],
                        "sample_end": sample_range[1],
                        "overlap_start": max(sample_range[0], abnorm_reg_range[0]),
                        "overlap_end": min(sample_range[1], abnorm_reg_range[1]),
                        "overlap_length": overlap, 
                        "percent_overlap": overlap/(sample_range[1]-sample_range[0])*100
                    })
                else:
                    continue
    else:
       continue

# final function (for all files in directory)

matched_calls = []

for file in os.listdir(calls_dir):
    if file=="A549.markdup.bam_CNVs": ### exclude cancer sample A549.markdup.bam_CNVs
        continue
    else:
        sample_cnv_call = pd.read_csv(calls_dir+str(file), sep="\t", names=["chr", "start", "end", "copy_number", "category"], header=None)
        for chromosome in sample_cnv_call["chr"].unique():
            if comm_abnorm_reg["Chromosome"].isin([chromosome]).any() == True:
                comm_abnorm_reg_chr = comm_abnorm_reg[comm_abnorm_reg["Chromosome"] == str(chromosome)]
                sample_cnv_call_chr = sample_cnv_call[sample_cnv_call["chr"] == str(chromosome)]
                for _, row in sample_cnv_call_chr.iterrows():
                    sample_range = [row["start"]+1, row["end"]]
                    for _, row in comm_abnorm_reg_chr.iterrows():
                        abnorm_reg_range = [row["start_clean"], row["end_clean"]]
                        overlap = get_overlap(sample_range, abnorm_reg_range)
                        if overlap > 0:
                            matched_calls.append({
                                "sample": file,
                                "comm_abnorm_reg_ID": row["comm_abnorm_reg_ID"],
                                "chromosome": chromosome,
                                "sample_start": sample_range[0],
                                "sample_end": sample_range[1],
                                "overlap_start": max(sample_range[0], abnorm_reg_range[0]),
                                "overlap_end": min(sample_range[1], abnorm_reg_range[1]),
                                "overlap_length": overlap, 
                                "percent_overlap": overlap/(sample_range[1]-sample_range[0])*100,
                                # add parental or clone col, quality, 
                            })
                        else:
                            continue
        else:
            continue

# save excel

matched_calls_df = pd.DataFrame(matched_calls)
writer = pd.ExcelWriter('hotspot_in_sample.xlsx', engine='xlsxwriter')
matched_calls_df.to_excel(writer, sheet_name='Sheet1')
writer.save()

