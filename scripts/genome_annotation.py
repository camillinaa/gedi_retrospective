# import libraries
import pandas as pd

# import merged cnv-hotspot data 
df = pd.read_excel("merged_df.xlsx")
# df = df[(df["chromosome"] == "11") & (df["cnv_start"] == 7800000) & (df["cnv_end"] == 7900000) & (df["hotspot_ID"] == "11_2800000_10700000")]
df.shape

# set paths
blacklist_bed = "/facility/nfdata-omics/reference/human/encode/hg38_v2/hg38-blacklist.v2.bed"
gapregions_bed = "/facility/nfdata-omics/reference/human/gencode/v44/gap_regions.bed"

# import and prepare annotation dataframes
blacklist_df = pd.read_csv(blacklist_bed, sep="\t", names=["Chromosome", "Start", "End", "Motivation"], header=None)
highsignal_annot = blacklist_df[blacklist_df["Motivation"]=="High Signal Region"].drop("Motivation", axis=1).assign(Chromosome=lambda df: df["Chromosome"].str.replace("chr", ""))
lowmap_annot = blacklist_df[blacklist_df["Motivation"]=="Low Mappability"].drop("Motivation", axis=1).assign(Chromosome=lambda df: df["Chromosome"].str.replace("chr", ""))
gapregions_annot = pd.read_csv(gapregions_bed, sep="\t", names=["Chromosome", "Start", "End"], header=None).assign(Chromosome=lambda df: df["Chromosome"].str.replace("chr", ""))

# function to calculate overlap percentage
def calculate_overlap(cnv_start, cnv_end, annot_regions):
    overlap_bp = 0
    for _, row in annot_regions.iterrows():
        annot_start, annot_end = row["Start"], row["End"]
        overlap = max(0, min(cnv_end, annot_end) - max(cnv_start, annot_start))
        overlap_bp += overlap
    return (overlap_bp / (cnv_end - cnv_start)) * 100 if (cnv_end - cnv_start) > 0 else 0

# annotation datasets with corresponding column names to be added
annotations = [
    (highsignal_annot, "% High Signal Region"),
    (lowmap_annot, "% Low Mappability Region"),
    (gapregions_annot, "Gap Fraction")
]

# calculate percentage overlap for each CNV across all annotation datasets
for annot_df, column_name in annotations:
    percent_overlaps = []
    for _, row in df.iterrows():
        chrom = str(row["chromosome"])  # Ensure consistent format
        cnv_start, cnv_end = row["cnv_start"], row["cnv_end"]
        
        # filter annotation regions for matching chromosome
        annot_regions = annot_df[annot_df["Chromosome"] == chrom]
        
        # computer overlap percentage
        overlap_pct = calculate_overlap(cnv_start, cnv_end, annot_regions)
        percent_overlaps.append(overlap_pct)
    
    # Add results to dataframe
    df[column_name] = percent_overlaps


# save xlsx

df.to_excel("cnv_hotspots_annotated.xlsx", index=False)