import pandas as pd
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scikit_posthocs as sp
import os
from matplotlib.patches import Rectangle
from matplotlib import rc



class analyzePTR:

    def __init__(self, PROJECT_DIR, analysis_name, project, samples_list):
        self.project = project
        self.analysis_name = analysis_name
        self.PROJECT_DIR = PROJECT_DIR


        # Metadata for the whole PTR project
        self.NZ_TYPE_FILE="/Users/stebliankin/Desktop/FLINTproject/data/annotations/nz_type.txt"
        self.TAX_FILE="/Users/stebliankin/Desktop/FLINTproject/data/annotations/kraken_taxonomy.txt"
        self.tax_levels = ["phylum", "class", "order", "family", "genus", "species"]

        # Inferred variables:
        self.SAMPLES_LIST_FILE = PROJECT_DIR + "/metadata/" + samples_list
        self.ptr_dir = PROJECT_DIR + "/ptr_samples"
        self.figures_dir = PROJECT_DIR + "/figures"
        self.ptr_time_series_file = PROJECT_DIR + "/PTR_time_series.csv"
        self.abundance_time_series_file = PROJECT_DIR + "/abundance_time_series.csv"

        if not os.path.exists(self.figures_dir):
            os.mkdir(self.figures_dir)

        if not os.path.exists(self.ptr_time_series_file):
            self.merge_time_series()

    def merge_samples_append(self, type):
        #
        # Merge all samples into one dataframe (no time series)
        #   Assuming that ptr and abundance are in the same folder
        #
        # Output with the following columns:
        #   genome | abundance or ptr | sample
        SAMPLES_LIST_FILE = self.SAMPLES_LIST_FILE
        ptr_dir = self.ptr_dir

        if type == "PTR":
            type="bPTR"
            names = ["genome", "ori", "ter", "PTR", "type"]
            measure = "PTR"
        if type == "abundance":
            names = ["genome", "coverage"]
            measure = "coverage"

        merged_genome_list = []
        merged_ptr_list = []
        merged_samples_list = []

        with open(SAMPLES_LIST_FILE, "r") as samples:
            for i, sample in enumerate(samples.readlines()):
                sample = sample.strip("\n")
                ptr_file = ptr_dir + "/" + type + "_" + sample + ".txt"
                ptr_df = pd.read_csv(ptr_file, sep="\t", names=names)

                if type == "abundance":
                    summ = ptr_df[measure].sum()
                    ptr_df[measure] = ptr_df[measure].apply(lambda x: x / summ)
                    ptr_df["genome"] = ptr_df["genome"].apply(lambda x: x.split("_")[1] + "_" + x.split("_")[2])
                genomes_current_list = list(ptr_df["genome"])

                ptr_current_list = list(ptr_df[measure])
                current_list_len = len(ptr_current_list)

                samples_current_list = [sample] * current_list_len

                merged_ptr_list = merged_ptr_list + ptr_current_list
                merged_samples_list = merged_samples_list + samples_current_list
                merged_genome_list = merged_genome_list + genomes_current_list

        merged_df = pd.DataFrame(
            {"genome": merged_genome_list, measure: merged_ptr_list, "sample": merged_samples_list})
        return merged_df


    def plot_distributions(self, measure, groups, merged_df, figures_dir, affix):
        # Plot distribution of Groups in merged_dataframe
        # Input:
        #   measure - abundance or ptr
        #   groups - dictionary or list with groups

        # Groups dict can be a list of groups

        analysis_name = self.analysis_name

        if not os.path.exists(figures_dir):
            os.mkdir(figures_dir)

        plt.ion()

        for group in groups:
            group_df = merged_df[merged_df["group"] == group]
            hist = group_df[measure].hist(bins=120, alpha=0.7, label=group)
            #plt.savefig(figures_dir + "/" + "distribution_" + group + ".png")

            # create legend
            hist.legend(prop={'size': 10})

        hist.set_xlabel(measure)
        hist.title.set_text(measure + " distribution " +analysis_name + " gropus")
        plt.savefig(figures_dir + "/" + "distribution_"+ affix + ".png")
            # plt.clf()
        plt.clf()

    def bar_plot(self, measure, merged_df, figures_dir, affix):
        analysis_name = self.analysis_name
        if not os.path.exists(figures_dir):
            os.mkdir(figures_dir)

        sns.set(style="whitegrid")
        tips = sns.load_dataset("tips")
        ax = sns.barplot(x="group", y=measure, data=merged_df, capsize=0.12)
        # plt.boxplot([ptr_listA, ptr_listB])
        plt.ylim(1, 2.1)
        plt.savefig(figures_dir + "/" + "hist_" + affix + ".png")
        plt.clf()


    def plot_violin(self, df, x, y, hue=False, order=None):
        fig = plt.subplots(figsize=(12, 5))
        if hue!=False:
            ax = sns.violinplot(data=df, y=y, x=x, scale="width", ci=68, hue=hue, order=order)
        else:
            ax = sns.violinplot(data=df, y=y, x=x, scale="width", ci=68, order=order)
        plt.tight_layout()
        plt.show()
    #def scatter_plot(self,):

    def box_plot(self, measure, merged_df, figures_dir, affix):
        analysis_name = self.analysis_name

        sns.set(style="whitegrid")
        tips = sns.load_dataset("tips")
        ax = sns.boxplot(x="group", y=measure, data=merged_df)
        # plt.boxplot([ptr_listA, ptr_listB])
        if measure == "abundance":
            plt.ylim(0, 1)
        else:
            plt.ylim(1, 5)
        plt.savefig(figures_dir + "/" + "box_" + affix + ".png")
        plt.clf()

    def calculate_p_values_groups(self, measure, merged_df, figures_dir, affix, plot=True):
        merged_df = merged_df[merged_df[measure].notnull()]
        p_values = sp.posthoc_ttest(merged_df, val_col=measure, group_col="group")
        if plot:
            p_values.to_csv(figures_dir + "/p_values_" +affix+ ".csv")
        return p_values

    def merge_time_series(self):
        def read_as_df(sample, PROJECT_DIR):
            sample = sample.strip("\n")
            ptr_file = PROJECT_DIR + "/" + "ptr_samples/" + "bPTR_" + sample + ".txt"
            abundance_file = PROJECT_DIR + "/" + "ptr_samples/" + "abundance_" + sample + ".txt"

            current_ptr_df = pd.read_csv(ptr_file, sep="\t", names=["genome", "ori", "ter", "ptr"])
            current_abundance_df = pd.read_csv(abundance_file, sep="\t", names=["genome", "coverage"])

            return current_ptr_df, current_abundance_df

        print("Merging time series samples..")
        SAMPLES_LIST_FILE = self.SAMPLES_LIST_FILE
        PROJECT_DIR = self.PROJECT_DIR
        OUT_PTR = PROJECT_DIR + "/PTR_time_series.csv"
        OUT_ABUNDANCE = PROJECT_DIR + "/abundance_time_series.csv"

        with open(SAMPLES_LIST_FILE, "r") as samples:
            columns_ptr = ["genome"]  # Initialize array with columns
            columns_abundance = ["genome"]
            for i, sample in enumerate(samples.readlines()):
                # ----------------------------
                # Read PTR and abundance
                # ----------------------------

                current_ptr_df, current_abundance_df = read_as_df(sample, PROJECT_DIR)

                # ------------------------------
                # Extract names for bPTR df
                # ------------------------------
                # Extract NZ number:
                # From PTR
                current_ptr_df["NZ"] = current_ptr_df["genome"].apply(
                    lambda row: row.split("/")[-1].split("_")[0] + "_" + row.split("/")[-1].split("_")[1])
                current_ptr_df = current_ptr_df[["NZ", "ptr"]]

                # From Abundance
                current_abundance_df["NZ"] = current_abundance_df["genome"].apply(
                    lambda row: row.split("_")[1] + "_" + row.split("_")[2])

                # Rename columns
                current_ptr_df[sample] = current_ptr_df["ptr"]
                current_abundance_df[sample] = current_abundance_df["coverage"]

                # Normalize abundance
                sum_coverage = current_abundance_df[sample].sum()
                current_abundance_df[sample] = current_abundance_df[sample].apply(lambda row: row/sum_coverage)

                current_ptr_df = current_ptr_df[["NZ", sample]]
                current_abundance_df = current_abundance_df[["NZ", sample]]

                if i == 0:
                    previous_ptr_df = current_ptr_df
                    previous_abundance_df = current_abundance_df
                else:
                    previous_ptr_df = previous_ptr_df.merge(current_ptr_df, on="NZ", how="outer")
                    previous_abundance_df = previous_abundance_df.merge(current_abundance_df, on="NZ", how="outer")
                pass
                print(i)

        previous_abundance_df.to_csv(OUT_ABUNDANCE, index=False)
        previous_ptr_df.to_csv(OUT_PTR, index=False)

    def merge_strains(self, ptr_merged_path, taxonomy_file, abundance_path, out_path):
        ################################
        # Merge strains to species level
        ################################

        ptr_df = pd.read_csv(ptr_merged_path)
        samples_list = list(ptr_df.columns)
        samples_list.remove("NZ")
        samples_list = [x.strip("\n") for x in samples_list]

        ptr_df.rename(columns=lambda x: x.strip("\n") + "_PTR", inplace=True)
        ptr_df.rename(columns={"NZ_PTR": "NZ"}, inplace=True)
        # Step 1 - Merge with abundance file
        abundance_df = pd.read_csv(abundance_path)

        abundance_df.rename(columns=lambda x: x.strip("\n") + "_abundance", inplace=True)
        abundance_df.rename(columns={"NZ_abundance": "NZ"}, inplace=True)

        nz_list = list(ptr_df["NZ"])
        abundance_df = abundance_df[abundance_df["NZ"].isin(nz_list)]
        # ptr_df = ptr_df.fillna(0) # convenient for computation. Note: Change for 1 later
        abundance_df = abundance_df.fillna(0)
        abundance_df = abundance_df.reset_index(drop=True)

        ptr_df = ptr_df.merge(abundance_df, on="NZ", how="left")

        # Step 2 - Merge NZ with taxa species name
        taxonomy_df = pd.read_csv(taxonomy_file, sep="\t")

        taxonomy_df["NZ"] = taxonomy_df["NZ"].apply(lambda x: x.split(".")[0])
        taxonomy_cols = list(taxonomy_df.columns)
        ptr_df = taxonomy_df.merge(ptr_df, on="NZ", how="right")

        species_list = list(ptr_df["species"].unique())

        # Step 3 - Combine taxa

        # create dictionary with new values
        ptr_dict = {}
        ptr_dict["sample"] = samples_list

        previous_species = ""
        k = -1

        for specie in species_list:
            print(specie)
            if isinstance(specie, str):
                tmp_df = ptr_df[ptr_df["species"] == specie]
                ptr_dict[specie + "#abundance"] = []
                ptr_dict[specie + "#PTR"] = []
                for sample in ptr_dict["sample"]:
                    tmp_df_sample = tmp_df[["NZ", "species", sample + "_PTR", sample + "_abundance"]]
                    tmp_df_sample = tmp_df_sample.dropna()
                    total_abundance = tmp_df_sample[sample + "_abundance"].sum()

                    if len(tmp_df_sample) > 0:
                        tmp_df_sample[sample + "_cumulative"] = tmp_df.apply(
                            lambda row: row[sample + "_abundance"] * row[sample + "_PTR"] / total_abundance, axis=1)
                        ptr = tmp_df_sample[sample + "_cumulative"].sum()
                        ptr_dict[specie + "#abundance"].append(total_abundance)
                        ptr_dict[specie + "#PTR"].append(ptr)

                    else:
                        ptr_dict[specie + "#abundance"].append(0)
                        ptr_dict[specie + "#PTR"].append(1)

        new_ptr_df = pd.DataFrame.from_dict(ptr_dict)
        new_ptr_df.to_csv(out_path, index=False)

        return


