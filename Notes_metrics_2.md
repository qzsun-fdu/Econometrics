# Econometrics by Prof. Ju (Continue)

## Chapter 11: IV and GMM

### 11.1 IV Estimation

How to solve endogeneity?

* Proxy

* IV

* Structural Modelling

* Panel Method: differencing to erase FE
$$
y_i=X_{i_{m\times1}}'\alpha+u_i,i=1,2,...,n,\quad E(u_i|X_i)\neq 0 \\
IV:Z_{i_{l\times1}},\quad l\geq m
$$

Two requirement of IV:

* Relevance condition: $\frac{1}{n}Z'X=\frac{1}{n}\sum^n_{i=1}Z_iX_i'\to E(ZX)\neq0$

* Exogeneity condition: $\frac{1}{n}Z'u\to_a0$

**1st stage:**
$$
X=Z\gamma+\mu \quad\Rightarrow\quad X=Z\hat{\gamma}+\hat{u}
$$

**2nd stage:**
$$
y=X\alpha+u=\hat{X}\alpha+(X-\hat{X})\alpha+u,\quad (X-\hat{X})\alpha+u\equiv\varepsilon
$$

Note that
$$
\hat{X}\perp X-\hat{X},\\
\hat{X}'u=X'Z(Z'Z)^{-1}(Z'u)\to_a0\Rightarrow \hat{X}\perp u,\\
\Rightarrow \hat{X}\perp \varepsilon.
$$

Thus, we can use the **least square method** to estimate:
$$
\hat{\alpha}_{2SLS}=\hat{\alpha}_{IV}=(\hat{X}'\hat{X})^{-1}\hat{X}'y\\
\overset{\hat{X}=P_ZX}{=}(X'P_ZX)^{-1}X'P_Zy\\
=\alpha+(X'P_ZX)^{-1}X'P_Z\varepsilon\\
=\alpha+(X'P_ZX)^{-1}X'P_Zu\\
=\alpha+(\frac{X'Z}{n}\frac{(Z'Z)^{-1}}{n}\frac{Z'X}{n})^{-1}\frac{X'Z}{n}\frac{(Z'Z)^{-1}}{n}\frac{Z'u}{n}\\
\to_a \alpha
$$

The last equation uses LLN (each term converges to finite matrix & $\frac{Z'u}{n}=o_p(1)$).

Assume that $Var(u|Z,X)=\Omega$ is known. We have
$$
Var(\hat{\alpha}_{2SLS}|Z,X)=(X'P_ZX)^{-1}X'P_Z\Omega P_ZX(X'P_ZX)^{-1}.
$$

In reality, we don't estimate $\Omega$ directly alone; instead, we estimate $Z'\Omega Z$ (lower dimension, so that less parameters are needed to estimate).

**If $\;l= m$:**

$Z'X$是方阵。
$$
\hat{\alpha}=(\hat{X}'\hat{X})^{-1}\hat{X}'y\\
=(X'Z(Z'Z)Z'X)^{-1}X'(Z(Z'Z)Z')y\\
=(Z'X)^{-1}Z'y
$$

is SAE (sample analogue estimator). Why?

Proof:
$$
y_i=X_i'\alpha+u_i\\
Z_iy_i=Z_iX_i'\alpha+Z_iu_i\\
\overset{E(\cdot)}{\Rightarrow} E(Z_iy_i)=E(Z_iX_i')\alpha+E(Z_iu_i)=E(Z_iX_i')\alpha\\
\Rightarrow \hat{\alpha}=(E(Z_iX_i'))^{-1}E(Z_iy_i)
$$

$$
\Rightarrow \hat{\alpha}_{SAE}=(\frac{1}{n}\sum_{i=1}^nZ_iX_i')^{-1}(\frac{1}{n}\sum_{i=1}^nZ_iy_i')\\
=(Z'X)^{-1}Z'y
$$

**If $\;l> m$:**
$$
Z'y=Z'X\alpha+Z'u
$$

要求$\hat{\alpha}$，左乘$m\times l$阶矩阵
$$
\Rightarrow \frac{X'Z}{n}\frac{(Z'Z)^{-1}}{n}\frac{Z'y}{n}=\frac{X'Z}{n}\frac{(Z'Z)^{-1}}{n}\frac{Z'X\alpha}{n}+\frac{X'Z}{n}\frac{(Z'Z)^{-1}}{n}\frac{Z'u}{n}\\
\to \frac{X'Z}{n}\frac{(Z'Z)^{-1}}{n}\frac{Z'X\alpha}{n}
$$

$$
\Rightarrow \hat{\alpha}_{SAE}=(X'P_ZX)^{-1}X'P_Zy
$$

但是左乘$m\times l$阶矩阵不是唯一的，怎么找到最好的取法？
$$
X_{n\times l}\to X_{n\times l}T_{l\times m},
$$

$T_{l\times m}$ is called **"selection matrix"**.

In 2SLS,
$$
T=(Z'Z)^{-1}Z'X.
$$

But what is the **most efficient** one (consistent & the least variance)? Let's find it!
$$
Z'y=Z'X\alpha+Z'u
$$

引入 selection matrix $T$, 转化成 square matrix 求逆:
$$
T'Z'y=T'Z'X\alpha+T'Z'u
$$

$$
\Rightarrow \quad \hat{\alpha}=(T'Z'X)^{-1}T'Z'y=\alpha+(T'Z'X)^{-1}T'Z'u
$$

$$
\Rightarrow Var(\hat{\alpha})=(T'Z'X)^{-1}\cdot T'Z'\Omega ZT\cdot(X'ZT)^{-1}
$$

夹心估计量的优化技巧：当$T'Z'X,\;T'Z'\Omega ZT,\;(X'ZT)^{-1} $三者相等时，方差最小。
凑得
$$
T^*=(Z'\Omega Z)^{-1}Z'X
$$

is the most efficient.

$$
Var(\hat{\alpha})=(X'Z(Z'\Omega Z)^{-1}Z'X)^{-1}
$$

In order to prove it, we introduce below

#### Theorem 11.1

$A$ and $B$ are p.s.d.. Then (using eigen decomposition)
$$
A\geq B\quad\Leftrightarrow\quad A^{-1}\leq B^{-1}.
$$

Thus, we only need to prove
$$
X'Z(Z'\Omega Z)^{-1}Z'X\geq Var(\alpha)
$$

which is true. (We don't present the proof here.)

### 11.2 GMM
