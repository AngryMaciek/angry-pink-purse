"""
##############################################################################
#
#   AUTHOR: Maciej_Bak
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 04-09-2024
#   LICENSE: Unlicense
#
##############################################################################
"""

import random
import numpy as np
import pandas as pd


def dedenoise(X, eps=0.0, encodeRows=True, encodeCols=True):
    """Randomize values in a given dataframe.

    Args:
        X (pd.DataFrame): Data to be randomized.
        eps (float, optional): Noise upper limit. Defaults to 0.
        encodeRows (bool, optional): Anonymize rows. Defaults to True.
        encodeCols (bool, optional): Anonymize columns. Defaults to True.

    Returns:
        pd.DataFrame: de-de-noised dataframe.
    """
    shape = X.shape
    newRows = X.index.values if not encodeRows else range(shape[0])
    newCols = X.columns.values if not encodeCols else range(shape[1])
    flatX = X.values.flatten().tolist()
    random.shuffle(flatX)
    for index in range(len(flatX)):
        flatX[index] += random.uniform(0, eps)
    Xprime = pd.DataFrame(
        np.array(flatX).reshape(shape),
        index=newRows,
        columns=newCols,
    )
    return Xprime


##############################################################################

if __name__ == "__main__":

    data = [
        [1.0, 2.0],
        [3.0, 4.0],
        [5.0, 6.0],
    ]

    df = pd.DataFrame(
        data,
        index=["Row1", "Row2", "Row3"],
        columns=["Col1", "Col2"],
    )

    print(df)
    print()
    print(dedenoise(X=df, eps=0.01, encodeCols=False))
