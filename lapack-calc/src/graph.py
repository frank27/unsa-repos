import pandas as pd
import matplotlib.pyplot as plt

class MatrixAlgorithmAnalysis:

    ALGORITHM_LEGEND = {
        0: "Standard Algorithm",
        1: "Divide And Conquer",
        2: "Strassen Algorithm",
    }

    CSV_FILE = "result_data.csv"

    def __processData(self, raw_df):
        df_stats = raw_df.groupby(["Algorithm", "Size"], as_index=False).agg(
                MeanTime=("TimeSeconds", "mean"),
                MedianTime=("TimeSeconds", "median")
            )
        print(df_stats)
        return df_stats

    def __export_data(self,raw_df, df_stats):
        repeats_df = (
            raw_df
            .pivot_table(
                index=["Algorithm", "AlgorithmName", "Size"],
                columns="RepeatNumber",
                values="TimeSeconds"
            )
            .rename(columns=lambda x: f"Repeat {x}")
            .reset_index()
        )
        df_stats["AlgorithmName"] = df_stats["Algorithm"].map(self.ALGORITHM_LEGEND)

        result = repeats_df.merge(df_stats,on=["Algorithm", "AlgorithmName", "Size"])

        result["MeanTime"]= result["MeanTime"] * 1_000_000
        result["MedianTime"]= result["MedianTime"] * 1_000_000
        result["Repeat 0"]= result["Repeat 0"] * 1_000_000
        result["Repeat 1"]= result["Repeat 1"] * 1_000_000
        result["Repeat 2"]= result["Repeat 2"] * 1_000_000

        result.to_csv("matrix_multiplication_comparison.csv", index=False)

    def __plot_from_csv(self, df):
        plt.figure(figsize=(12, 6))
        ax = plt.gca()

        # Color maps per algorithm
        ALGO_COLORMAPS = {
            0: plt.cm.Blues,
            1: plt.cm.Purples,
            2: plt.cm.Oranges,
            3: plt.cm.Purples
        }

        # Line styles
        line_styles = {
            0: ("o-", 1.5),
            1: ("s--", 1.5),
            2: ("d-.", 1.5),
            3: ("d-.", 1.5),
        }

        for algo, group in df.groupby("Algorithm"):
            x = group["Size"].values
            y = group["MeanTime"].values

            linestyle, lw = line_styles.get(algo, ("o-", 1.5))
            color = ALGO_COLORMAPS.get(algo, plt.cm.Greys)(0.85)

            ax.plot(
                x,
                y,
                linestyle,
                lw=lw,
                color=color,
                label=self.ALGORITHM_LEGEND.get(algo, f"Algorithm {algo}")
            )

        ax.set_xticks(sorted(df["Size"].unique()))
        ax.set_title("Matrix Multiplication Execution Time")
        ax.set_xlabel("Matrix Size (N x N)")
        ax.set_ylabel("Time (seconds)")
        ax.grid(True)
        plt.legend()
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig("matrix_multiplication_analysis_cpp.png")
        plt.show()

    def execute(self):
        raw_df = pd.read_csv(self.CSV_FILE)
        stats_df = self.__processData(raw_df)
        self.__plot_from_csv(stats_df)
        self.__export_data(raw_df, stats_df)


MatrixAlgorithmAnalysis().execute()
