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

# compute total percentage overlap for all cnvs covering the same hotspot in each sample (output: dictionary)

def get_total_overlap(matched_df):
    total_overlap_dict = matched_df.groupby(["sample", "chromosome", "hotspot_start", "hotspot_end"])["percent_overlap"].sum().to_dict()
    return total_overlap_dict

# general region annotation % function 

def apply_annotation(df, annotation_df, colname):
    
    df = df.merge(annotation_df.assign(Chromosome=annotation_df["Chromosome"].str.replace("chr","")), left_on="chromosome", right_on="Chromosome", how="inner")
    df["overlap_start"] = df[["cnv_start", "Start"]].max(axis=1)
    df["overlap_end"] = df[["cnv_end", "End"]].min(axis=1)
    df["overlap"] = df["overlap_end"]-df["overlap_start"]

    df[colname] = (df["overlap"] / (df["cnv_end"] - df["cnv_start"])) * 100

    return df

# apply functions to all files in directory

not_included = ["A549.markdup.bam_CNVs"] # define samples to exclude

cnv_hotspot_match = match_cnvs_to_hotspots(cnv_calls_dir, Common_Abnormal_Regions_path, not_included)
total_overlap = get_total_overlap(cnv_hotspot_match)

blacklist_df = pd.read_csv(blacklist_bed, sep="\t", names=["Chromosome", "Start", "End", "Motivation"], header=None)
highsignal_annot = blacklist_df[blacklist_df["Motivation"]=="High Signal Region"].drop("Motivation", axis=1)
lowmap_annot = blacklist_df[blacklist_df["Motivation"]=="Low Mappability"].drop("Motivation", axis=1)
gapregions_annot = pd.read_csv(gapregions_bed, sep="\t", names=["Chromosome", "Start", "End"], header=None)

df_annotated = annotate(cnv_hotspot_match, highsignal_annot, "% High Signal Region")
xxx = annotate(cnv_hotspot_match, lowmap_annot, "% Low Mappability Region")
xxx = annotate(cnv_hotspot_match, gapregions_annot, "Gap Fraction")

writer = pd.ExcelWriter('cnvs_hotspots_mapped.xlsx', engine='xlsxwriter')
final_df.to_excel(writer, sheet_name='Sheet1')
writer.save()
