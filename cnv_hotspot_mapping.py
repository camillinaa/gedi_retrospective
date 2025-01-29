# import libraries

import pandas as pd
import os
import re
import xlsxwriter

# import data

Common_Abnormal_Regions_df = pd.read_csv("Common_Abnormal_Regions.csv")
cnv_calls_dir = "/facility/nfdata-omics/projects/GEDI_LPS/sample_prioritization/controlFREEC_calls/"
testfile_dir = "/facility/nfdata-omics/projects/GEDI_LPS/sample_prioritization/testfile/" 
blacklist_bed = "/facility/nfdata-omics/reference/human/encode/hg38_v2/hg38-blacklist.v2.bed"

# overlap function

def get_overlap(a, b):
   return max(0, min(a[1], b[1]) - max(a[0], b[0])+1)

# cnv-hotspot matching function

def match_cnvs_to_hotspots(cnvs_dir, hotspots_df, skip_list):
    cnvs_hotspots_matched = []
    for file in os.listdir(cnvs_dir):
        if file in skip_list:
            continue
        sample_cnv_call = pd.read_csv(cnvs_dir+str(file), sep="\t", names=["chr", "start", "end", "copy_number", "category"], header=None)
        sample_type = "clone" if re.search(r"C\d+(?=\.markdup\.bam_CNVs)", file) else "parental"
        for chromosome in hotspots_df["Chromosome"].unique():
            if sample_cnv_call["chr"].isin([chromosome]).any() == True:
                hotspots_chr = hotspots_df[hotspots_df["Chromosome"] == str(chromosome)]
                sample_cnv_call_chr = sample_cnv_call[sample_cnv_call["chr"] == str(chromosome)]
                for _, row in hotspots_chr.iterrows():
                    hotspot_range = [row["start_clean"], row["end_clean"]]
                    for _, row in sample_cnv_call_chr.iterrows():
                        sample_range = [row["start"], row["end"]]
                        overlap = get_overlap(hotspot_range, sample_range)
                        if overlap > 0:
                            cnvs_hotspots_matched.append({
                                "sample": file,
                                "hotspot_ID": hotspots_chr["hotspot_ID"],
                                "chromosome": chromosome,
                                "sample_type": sample_type,
                                "cnv_start": sample_range[0],
                                "cnv_end": sample_range[1],
                                "hotspot_start": hotspot_range[0],
                                "hotspot_end": hotspot_range[1],
                                "overlap_start": max(sample_range[0], hotspot_range[0]),
                                "overlap_end": min(sample_range[1], hotspot_range[1]),
                                "overlap_length": overlap, 
                                "percent_overlap": overlap/(hotspot_range[1]-hotspot_range[0])*100
                                # add quality metrics
                            })
    return cnvs_hotspots_matched

def add_total_percentage_overlap(matched_matrix):
    matched_df = pd.DataFrame(matched_matrix)
    matched_df["percent_overlap_total"] = matched_df.groupby(["sample", "chromosome", "hotspot_start", "hotspot_end"])["percent_overlap"].transform("sum")
    return matched_df

# apply function to all files in directory

not_included = ["A549.markdup.bam_CNVs"] # define samples to exclude

result = match_cnvs_to_hotspots(cnv_calls_dir, Common_Abnormal_Regions_df, not_included)

result_total_overlap = add_total_percentage_overlap(result)

# save excel

writer = pd.ExcelWriter('cnvs_hotspots_mapped.xlsx', engine='xlsxwriter')
result_total_overlap.to_excel(writer, sheet_name='Sheet1')
writer.save()

