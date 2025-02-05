import polars as pl


def show_data():
    csvfile = "metric.csv"
    data = pl.read_csv(source=csvfile, separator="\t").filter(
        pl.col("longIndel").str.contains("-")
    )
    print(data.select("longIndel").head(20))


if __name__ == "__main__":
    show_data()
    pass
