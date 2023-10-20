# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.15.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import numpy as np
import pandas as pd
import scipy as sp

nA=1000,
nB=1000,
mA=[0, 0, 0]
mB=[0, 0, 0]
covA=[[1.0, 0.2, 0.2], [0.2, 1.0, 0.2], [0.2, 0.2, 1.0]]
covB=[[1.0, 0.2, 0.2], [0.2, 1.0, 0.2], [0.2, 0.2, 1.0]]
px1c=[0.5, 0.5]
px2c=[0.333, 0.334, 0.333]
px3c=[0.25, 0.25, 0.25, 0.25]
p=0.6
aA=[1, 1, 1, 1, 1, 1]
aB=[1, 1, 1, 1, 1, 1]
eps=0

# +
X_glob1 = np.random.multivariate_normal(mA, covA, nA)
X_glob2 = np.random.multivariate_normal(mB, covB, nB)

X_glob1 = pd.DataFrame(X_glob1, columns=["X1", "X2", "X3"])
X_glob2 = pd.DataFrame(X_glob2, columns=["X1", "X2", "X3"])
X_glob = pd.concat([X_glob1, X_glob2])

px1cc = np.cumsum(px1c)[:-1]
px2cc = np.cumsum(px2c)[:-1]
px3cc = np.cumsum(px3c)[:-1]

px1cc, px2cc, px3cc
# -

qx1c = sp.stats.norm.ppf(px1cc, loc=0, scale=1)
qx2c = sp.stats.norm.ppf(px2cc, loc=0, scale=1)
qx3c = sp.stats.norm.ppf(px3cc, loc=0, scale=1)
qx1c, qx2c, qx3c

bins11 = [min(X_glob1["X1"]) - 100] + list(qx1c) + [max(X_glob1["X1"]) + 100]
bins12 = [min(X_glob1["X2"]) - 100] + list(qx2c) + [max(X_glob1["X2"]) + 100]
bins13 = [min(X_glob1["X3"]) - 100] + list(qx3c) + [max(X_glob1["X3"]) + 100]

# +
X = pd.DataFrame(np.random.multivariate_normal(mA, covA, nA), columns=["X1", "X2", "X3"])

X1 = np.digitize(X["X1"], bins=bins11)
X2 = np.digitize(X["X2"], bins=bins12)
X3 = np.digitize(X["X3"], bins=bins13)
# -

X1c = pd.get_dummies(X1).iloc[:,1:]
X2c = pd.get_dummies(X2).iloc[:,1:]
X3c = pd.get_dummies(X3).iloc[:,1:]
X3c

# +
X_dumm1 = np.c_[X1c, X2c, X3c]

Y1 = X_dumm1 @ aA

b11 = np.quantile(Y1, [0.25, 0.5, 0.75])
b22 = np.quantile(Y1, [1 / 3, 2 / 3])
binsY11 = [min(Y1) - 100] + list(b11) + [max(Y1) + 100]
binsY22 = [min(Y1) - 100] + list(b22) + [max(Y1) + 100]

binsY11, binsY22

# +
X_glob1 = pd.DataFrame(np.random.multivariate_normal(mA, covA, nA), columns=["X1", "X2", "X3"])
X_glob2 = pd.DataFrame(np.random.multivariate_normal(mB, covB, nB), columns=["X1", "X2", "X3"])

X11 = np.digitize(X_glob1["X1"], bins=bins11)
X21 = np.digitize(X_glob2["X1"], bins=bins11)
X12 = np.digitize(X_glob1["X2"], bins=bins12)
X22 = np.digitize(X_glob2["X2"], bins=bins12)
X13 = np.digitize(X_glob1["X3"], bins=bins13)
X23 = np.digitize(X_glob2["X3"], bins=bins13)

# +
X11c = pd.get_dummies(X11).iloc[:,1:]
X21c = pd.get_dummies(X21).iloc[:,1:]
X12c = pd.get_dummies(X12).iloc[:,1:]
X22c = pd.get_dummies(X22).iloc[:,1:]
X13c = pd.get_dummies(X13).iloc[:,1:]
X23c = pd.get_dummies(X23).iloc[:,1:]

Y1 = np.c_[X11c, X12c, X13c] @ aA
Y2 = np.c_[X21c, X22c, X23c] @ aB

Y11 = np.digitize(Y1, bins=binsY11)
Y12 = np.digitize(Y1, bins=binsY22)
binsY11eps = [x + params.eps for x in binsY11]
Y21 = np.digitize(Y2, bins=binsY11eps)
binsY22eps = [x + params.eps for x in binsY22]
Y22 = np.digitize(Y2, bins=binsY22eps)

YY = np.r_[Y11, Y21]
ZZ = np.r_[Y12, Y22]

# +
X1 = np.r_[X11, X21]
X2 = np.r_[X12, X22]
X3 = np.r_[X13, X23]

X1 = np.r_[X11, X21]
X2 = np.r_[X12, X22]
X3 = np.r_[X13, X23]

py = [
    1 - p,
    (p - p**2) / 2,
    (p - p**2) / 2,
    p**2 / 4,
    p**2 / 4,
    p**2 / 4,
    p**2 / 4,
]
pz = [1 - p, (p - p**2) / 2, (p - p**2) / 2, p**2 / 3, p**2 / 3, p**2 / 3]

U = np.random.multinomial(1, py, nA + nB)
V = np.random.multinomial(1, pz, nA + nB)
UU = U[:, 1] + U[:, 2] * 2 + U[:, 3] * 3 + U[:, 4] * 4 + U[:, 5] * 5 + U[:, 6] * 6
VV = V[:, 1] + V[:, 2] * 2 + V[:, 3] * 3 + V[:, 4] * 4 + V[:, 5] * 5

Y = np.arange(nA + nB)
Z = np.arange(nA + nB)
Y[UU == 0] = YY[UU == 0]
Y[UU == 1] = YY[UU == 1] - 1
Y[UU == 2] = YY[UU == 2] + 1
Y[UU == 3] = 1
Y[UU == 4] = 2
Y[UU == 5] = 3
Y[UU == 6] = 4
Y[Y > 4] = 4
Y[Y < 1] = 1
Z[VV == 0] = ZZ[VV == 0]
Z[VV == 1] = ZZ[VV == 1] - 1
Z[VV == 2] = ZZ[VV == 2] + 1
Z[VV == 3] = 1
Z[VV == 4] = 2
Z[VV == 5] = 3
Z[Z > 3] = 3
Z[Z < 1] = 1

database = np.array([1] * nA + [2] * nB)
data = pd.DataFrame(
    np.c_[X1, X2, X3, Y, Z, database],
    columns=["X1", "X2", "X3", "Y", "Z", "database"],
)
data[["X1", "X2", "X3"]] -= 1
# -

params = DataParameters(nA=1000, nB=1000, mB=[0, 0, 0], eps=0.00, p=0.2)
params

params.nA

data = generator_xcat_ycat(params)

df = pd.DataFrame(np.c_[data.X1, data.X2, data.X3, 
                                    data.Y, data.Z, data.database],
                              columns = ["X1", "X2", "X3", "Y", "Z", "database"])
    

df


