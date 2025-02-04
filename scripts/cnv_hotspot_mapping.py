# import libraries

import pandas as pd
import numpy as np
import os
import re
import xlsxwriter

# set paths

Common_Abnormal_Regions_path = "~/GEDI_retrospective/data/Common_Abnormal_Regions.csv"
cnv_calls_dir = "/facility/nfdata-omics/projects/GEDI_LPS/sample_prioritization/controlFREEC_calls/"
testfile_dir = "/facility/nfdata-omics/projects/GEDI_LPS/sample_prioritization/testfile/" 
blacklist_bed = "/facility/nfdata-omics/reference/human/encode/hg38_v2/hg38-blacklist.v2.bed"
gapregions_bed = "/facility/nfdata-omics/reference/human/gencode/v44/gap_regions.bed"

# overlap function

def get_overlap(a, b):
   return max(0, min(a[1], b[1]) - max(a[0], b[0])+1)

# cnv-hotspot matching function

def match_cnvs_to_hotspots(cnvs_dir, hotspots_path, skip_list):
    
    hotspots_df = pd.read_csv(hotspots_path)
    matched = []

    for file in os.listdir(cnvs_dir):
        
        if file in skip_list:
            continue

        sample_cnv_call = pd.read_csv(cnvs_dir+str(file), sep="\t", names=["chr", "start", "end", "copy_number", "category"], header=None)
        sample_type = "clone" if re.search(r"C\d+(?=\.markdup\.bam_CNVs)", file) else "parental"

        for _, hotspot in hotspots_df.iterrows():
            hotspot_range = [hotspot["start_clean"], hotspot["end_clean"]]

            for _, cnv in sample_cnv_call.iterrows():
                sample_range = [cnv["start"], cnv["end"]]

                if hotspot["Chromosome"] != cnv["chr"]:
                    continue

                overlap = get_overlap(hotspot_range, sample_range)

                if overlap > 0:
                    matched.append({
                        "sample": file,
                        "hotspot_ID": hotspot["hotspot_ID"],
                        "chromosome": hotspot["Chromosome"],
                        "sample_type": sample_type,
                        "cnv_start": sample_range[0],
                        "cnv_end": sample_range[1],
                        "hotspot_start": hotspot_range[0],
                        "hotspot_end": hotspot_range[1],
                        "overlap_start": max(sample_range[0], hotspot_range[0]),
                        "overlap_end": min(sample_range[1], hotspot_range[1]),
                        "overlap_length": overlap, 
                        "percent_overlap": overlap/(hotspot_range[1]-hotspot_range[0])*100
                    })

    matched_df = pd.DataFrame(matched)
    
    return matched_df

# add percentage overlap for all cnvs covering the same hotspot in each sample 

def add_total_overlap(matched_df):
    total_overlap_series = matched_df.groupby(["sample", "chromosome", "hotspot_start", "hotspot_end"])["percent_overlap"].sum().rename("total_overlap")
    df = matched_df.merge(total_overlap_series, on=["sample", "chromosome", "hotspot_start", "hotspot_end"], how="left")
    return df

# general region annotation % function 

def apply_annotation(df, annotation_df, colname):

    annotation_df = annotation_df.drop_duplicates(subset=["Chromosome", "Start", "End"])
    annotation_df["Chromosome"] = annotation_df["Chromosome"].str.replace("chr", "", regex=True)

    df_merged = df.merge(annotation_df, left_on="chromosome", right_on="Chromosome", how="left")

    print("df after merge: \n", df_merged[(df_merged["chromosome"] == "11") & (df_merged["cnv_start"] == 7800000) & (df_merged["cnv_end"] == 7900000) & (df_merged["hotspot_ID"] == "11_2800000_10700000")])

    df_merged["overlap_start"] = df_merged[["cnv_start", "Start"]].max(axis=1).where(df_merged["Start"].notna(), other=pd.NA)
    df_merged["overlap_end"] = df_merged[["cnv_end", "End"]].min(axis=1).where(df_merged["End"].notna(), other=pd.NA)
    df_merged["overlap"] = (df_merged["overlap_end"] - df_merged["overlap_start"] + 1).where(df_merged["End"].notna() & df_merged["Start"].notna(), other=pd.NA)
    df_merged[colname] = ((df_merged["overlap"] / (df_merged["cnv_end"] - df_merged["cnv_start"]+1)) * 100).where(df_merged["overlap"] > 0, other=0) # includes case if annot_df is NA
    
    print("df after overlap percentage calculation: \n", df_merged[(df_merged["chromosome"] == "11") & (df_merged["cnv_start"] == 7800000) & (df_merged["cnv_end"] == 7900000) & (df_merged["hotspot_ID"] == "11_2800000_10700000")])

    print("df after column drop: \n", df_merged[(df_merged["chromosome"] == "11") & (df_merged["cnv_start"] == 7800000) & (df_merged["cnv_end"] == 7900000) & (df_merged["hotspot_ID"] == "11_2800000_10700000")])

    #df_debug = df[(df["chromosome"] == "1") & (df["cnv_start"] == 16200000) & (df["cnv_end"] == 17500000) & (df["hotspot_ID"] == "1_16200000_17074942")]
    #print(df_debug)

    #df_debug = df[(df["chromosome"] == "11") & (df["cnv_start"] == 7800000) & (df["cnv_end"] == 7900000) & (df["hotspot_ID"] == "11_2800000_10700000")]
    #print(df_debug)

    #df_debug = df_debug.groupby(["sample", "sample_type", "chromosome", "hotspot_ID", "cnv_start", "cnv_end", "total_overlap"], as_index=False).sum()
    #df_debug[colname] = ((df_debug["overlap"] / (df_debug["cnv_end"] - df_debug["cnv_start"]+1)) * 100).where(df_debug["overlap"] > 0, other=0)
    
    return df_merged[colname]


# group by hotspots

def group_hotspots(df_test):

    df_test["cnv_length"] = df_test["cnv_end"]-df_test["cnv_start"]+1
    df_length = df_test.groupby(["sample","chromosome","hotspot_ID"], as_index=False)["cnv_length"].sum().rename(columns={"cnv_length": "total_hotspot"})
    df_test = df_test.merge(df_length, left_on=["sample","chromosome","hotspot_ID"], right_on=["sample","chromosome","hotspot_ID"], how="left")
    df_test["cnv_fraction"] = df_test["cnv_length"]/df_test["total_hotspot"]

    df_test["% High Signal Region"] = df_test["% High Signal Region"]*df_test["cnv_fraction"]
    df_test["% Low Mappability Region"] = df_test["% Low Mappability Region"]*df_test["cnv_fraction"]
    df_test["Gap Fraction"] = df_test["Gap Fraction"]*df_test["cnv_fraction"]

    df_test = df_test.groupby(["sample","chromosome","hotspot_ID","total_overlap"], as_index=False)["% High Signal Region", "% Low Mappability Region", "Gap Fraction"].sum()

    df_test["sample_type"] = df_test["sample"].apply(lambda x: "clone" if re.search(r"C\d+(?=\.markdup\.bam_CNVs)", x) else "parental")

    return df_test

# apply functions to all files in directory

not_included = ["A549.markdup.bam_CNVs"] # define samples to exclude

cnv_hotspot_match = match_cnvs_to_hotspots(cnv_calls_dir, Common_Abnormal_Regions_path, not_included)
total_overlap = add_total_overlap(cnv_hotspot_match)
df = total_overlap[["sample","sample_type", "chromosome","hotspot_ID","cnv_start","cnv_end","total_overlap"]] 

blacklist_df = pd.read_csv(blacklist_bed, sep="\t", names=["Chromosome", "Start", "End", "Motivation"], header=None)
highsignal_annot = blacklist_df[blacklist_df["Motivation"]=="High Signal Region"].drop("Motivation", axis=1)
lowmap_annot = blacklist_df[blacklist_df["Motivation"]=="Low Mappability"].drop("Motivation", axis=1)
gapregions_annot = pd.read_csv(gapregions_bed, sep="\t", names=["Chromosome", "Start", "End"], header=None)

annotations = [
    (highsignal_annot, "% High Signal Region"),
    (lowmap_annot, "% Low Mappability Region"),
    (gapregions_annot, "Gap Fraction")
]
df_annotated = df.copy()
for annot, colname in annotations:
    df_annotated[colname] = apply_annotation(df_annotated, annot, colname)

#df_hotspots = group_hotspots(df_annotated)

# save excels

writer = pd.ExcelWriter('merged_df.xlsx', engine='xlsxwriter')
df.to_excel(writer, sheet_name='Sheet1')
writer._save()

writer = pd.ExcelWriter('cnvs_hotspots_mapped.xlsx', engine='xlsxwriter')
df_annotated.to_excel(writer, sheet_name='Sheet1')
writer._save()

writer = pd.ExcelWriter('hotspots_mapped.xlsx', engine='xlsxwriter')
df_hotspots.to_excel(writer, sheet_name='Sheet1')
writer.save()