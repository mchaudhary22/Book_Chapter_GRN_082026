import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import dask
dask.config.set({"dataframe.query-planning": False})

import numpy as np
import pandas as pd
from distributed import Client, LocalCluster
from arboreto.algo import grnboost2


def main():
    np.random.seed(42)

    n_cells, n_genes, n_tfs = 50, 20, 5
    gene_names = [f"gene{i}" for i in range(n_genes)]
    tf_names = gene_names[:n_tfs]

    expr = pd.DataFrame(
        np.random.rand(n_cells, n_genes),
        columns=gene_names,
    )

    cluster = LocalCluster(n_workers=4, threads_per_worker=1, processes=True)
    client = Client(cluster)

    network = grnboost2(
        expression_data=expr,
        tf_names=tf_names,
        client_or_address=client,
        seed=42,
    )

    client.close()
    cluster.close()

    print("EDGE_COUNT", len(network))
    print(network.head())


if __name__ == "__main__":
    main()
