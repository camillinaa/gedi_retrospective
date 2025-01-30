# import libraries

import pandas as pd
import os
import re
import xlsxwriter

# set paths

Common_Abnormal_Regions_path = "Common_Abnormal_Regions.csv"
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
    matched_df = pd.DataFrame(cnvs_hotspots_matched)
    return matched_df

# compute total percentage overlap for all cnvs covering the same hotspot in each sample (output: dictionary)

def get_total_overlap(matched_df):
    total_overlap_dict = matched_df.groupby(["sample", "chromosome", "hotspot_start", "hotspot_end"])["percent_overlap"].sum().to_dict()
    return total_overlap_dict

# compute blacklist/gap/lowsignal regions (output: dictionary)

def get_genome_regions(matched_df, blacklist_bed, gapregions_bed):
    blacklist_df = pd.read_csv(blacklist_bed, sep="\t", names=["Chromosome", "Start", "End", "Motivation"], header=None)
    genome_regions_dict = {}
    #gapregions_df = pd.read_csv(gapregions_bed, sep="\t", names=["Chromosome", "Start", "End"], header=None)
    for _, call in matched_df.iterrows():
        cnv_range = [call["cnv_start"], call["cnv_end"]]
        regions = set()
        for _, region in blacklist_df.iterrows():
            blacklist_range = [region["Start"], region["End"]]
            overlap = get_overlap(cnv_range, blacklist_range)
            if overlap > 0:
                regions.add(region["Motivation"])
        genome_regions_dict[call["sample"]] = ", ".join(regions) if regions else "None"
    return genome_regions_dict

# apply function to all files in directory

def main():
    not_included = ["A549.markdup.bam_CNVs"] # define samples to exclude

    cnv_hotspot_match = match_cnvs_to_hotspots(cnv_calls_dir, Common_Abnormal_Regions_path, not_included)
    
    total_overlap = get_total_overlap(cnv_hotspot_match)
    genome_regions = get_genome_regions(cnv_hotspot_match, blacklist_bed, gapregions_bed)

    final_df = cnv_hotspot_match.assign(
        total_overlap=cnv_hotspot_match.apply(
            lambda row: total_overlap.get(
                (row["sample"], row["chromosome"], row["hotspot_start"], row["hotspot_end"]), 0
            ), axis=1
        ),
        genome_region=cnv_hotspot_match["sample"].map(genome_regions)
    )

    writer = pd.ExcelWriter('cnvs_hotspots_mapped.xlsx', engine='xlsxwriter')
    final_df.to_excel(writer, sheet_name='Sheet1')
    writer.save()

    return final_df



if __name__ == "__main__":
    final_result=main()
