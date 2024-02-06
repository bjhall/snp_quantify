#!/usr/bin/env python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_table("simulated.tsv")

sample_ids = df["sample"].unique()

g = sns.FacetGrid(df, col="total_reads")
g.map_dataframe(sns.scatterplot, x="expected_fraction", y="observed_fraction", hue="sample", hue_order=sample_ids)
g.add_legend()
plt.savefig("evaluation.pdf", bbox_inches = "tight")

