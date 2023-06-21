#!/usr/bin/env python3
"""
To run:
```
conda activate metabolike
pip install -r requirements.txt
./inspect_paths.py
```
"""

import random
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# (this is importing from reaction_distance.py in the current folder)
from reaction_distance import print_path


def main():
    SEED = None
    random.seed(SEED)

    # acetate, dGDP, etc not in blacklist here
    #data_dir = Path('2023-04-03_outputs')

    # Added a few things to compound name blacklist for this, but it crashed partway
    # through run
    #data_dir = Path('2023-04-18_partial_outputs')
    data_dir = Path('2023-04-19_outputs')
    #data_dir = Path('.')

    print('loading distances and paths...', flush=True, end='')

    min_dists = pd.read_pickle(data_dir / 'min_dists.p')
    pathway_min_dists = pd.read_pickle(data_dir / 'pathway_min_dists.p')

    min_paths = pd.read_pickle(data_dir / 'min_paths.p')
    pathway_min_paths = pd.read_pickle(data_dir / 'pathway_min_paths.p')

    print(' done', flush=True)


    def sample_path(path_df: pd.DataFrame, seed=SEED) -> dict:
        """
        Args:
            path_df: chemicals x chemicals DataFrame, where each element is a list of
                paths, and each path (each list element) is a dict with the keys
                'compounds' and 'reactions'

        Returns a single path between two chemicals ()
        """
        n = 1
        path_df_subset = path_df.sample(axis=0, n=n, random_state=seed
            ).sample(axis=1, n=n, random_state=seed)

        assert path_df_subset.shape == (n, n)
        # (assumes n=1 above)
        paths: list[dict]|float = path_df_subset.iat[0, 0]

        if type(paths) is not list:
            assert pd.isna(paths)
            # as long as we have some non-null entries, should terminate eventually...
            return sample_path(path_df, seed=seed)

        assert type(paths) is list and len(paths) > 0
        assert set([type(path) for path in paths]) == {dict}

        # (random.sample will return a list even if only sampling one)
        path: dict = random.sample(paths, 1)[0]
        return path

    min_dists = min_dists.replace(float('inf'), float('nan')
        ).dropna(how='all', axis=0).dropna(how='all', axis=1)

    min_paths = min_paths.loc[min_dists.index, min_dists.columns]

    for _ in range(5):
        path = sample_path(min_paths)
        print_path(path)
        print()
        import ipdb; ipdb.set_trace()

    vmin = 0
    # Actual max goes above this (which maybe it shouldn't?), but some of Google's
    # histograms topped out ~12 steps.
    vmax = 12

    fig, ax = plt.subplots()
    sns.heatmap(min_dists, ax=ax, cmap='crest', vmin=vmin, vmax=vmax,
        xticklabels=True, yticklabels=True
    )

    # TODO assert nothing outside of min_dists that is notna here
    pathway_min_dists = pathway_min_dists.loc[min_dists.index, min_dists.columns]
    pathway_min_paths = pathway_min_paths.loc[min_dists.index, min_dists.columns]

    fig, ax = plt.subplots()
    # This one is dramatically sparser.
    sns.heatmap(pathway_min_dists, ax=ax, cmap='crest', vmin=vmin, vmax=vmax,
        xticklabels=True, yticklabels=True
    )
    plt.show()
    import ipdb; ipdb.set_trace()


if __name__ == '__main__':
    # TODO TODO TODO also take optional args to pass in chemical names to check
    # TODO also arg to check all paths sequentially? just always do in above case?
    main()

