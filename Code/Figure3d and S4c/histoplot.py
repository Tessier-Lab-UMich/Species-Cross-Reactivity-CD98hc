'''
This code is used to plot the PSERM scores of sequences from the generated virtual libraries, showing the distribution of PSERM scores across the library sequences.
'''
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

file_path1 = './HCDR3_PSERMs_average_opt_S1_10000.csv'
file_path10 = './HCDR3_PSERMs_average_unopt_S1_10000.csv'

file_path2 = './HCDR3_PSERMs_average_opt_S2_10000.csv'
file_path20 = './HCDR3_PSERMs_average_unopt_S2_10000.csv'


def draw_histoplot(file_opt, file_un):
    dataset_un = pd.read_csv(file_un) 
    dataset_opt = pd.read_csv(file_opt)

    #merged_df = pd.merge(dataset_un, dataset_opt, how='outer')
    merged_df = pd.concat([dataset_un, dataset_opt], ignore_index=True)

    custom_palette = ["#E76F51", "#4C72B0"]

    sns.histplot(
        data=merged_df,
        x="PSERM average",
        bins=50,
        kde=True,
        hue="Library",
        fill=True,
        palette=custom_palette,
        alpha=0.5
    )

    plt.show()


draw_histoplot(file_path1, file_path10)
draw_histoplot(file_path2, file_path20)