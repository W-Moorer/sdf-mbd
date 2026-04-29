# 无压缩 SDF-NCP 接触动力学理论推导与实现修正指南

> 面向第一篇传统计算力学论文：**Signed-Distance-Field-Based Nonlinear Complementarity Contact Dynamics**  
> 目标：不使用 AI、不使用神经网络，建立一个可验证、可复现、可由 Codex 修正实现错误的 **无压缩 SDF-NCP 接触求解框架**。

---

## 0. 本文档用途

本文档用于指导 Codex 检查、修正和完善当前代码实现。重点是：

1. 推导无压缩 SDF-NCP 接触动力学的完整理论；
2. 明确每个变量的维度、符号和物理意义；
3. 给出时间步求解残差的正确形式；
4. 给出平滑 Fischer-Burmeister NCP 残差及其导数；
5. 给出无压缩多点接触的装配流程；
6. 给出常见实现错误和修正标准；
7. 给出 Codex 可直接执行的实现检查清单。

本文档中的“无压缩”指：进入候选接触集合的每一个接触采样点都保留独立的 gap、normal、Jacobian row 和 multiplier，不做聚类、合并、最深点筛选或代表点约化。

---

## 1. 基本设定

考虑一个离散机械系统，其广义坐标为

\[
q \in \mathbb{R}^{n_q},
\]

广义速度为

\[
v = \dot q \in \mathbb{R}^{n_q}.
\]

系统质量矩阵为

\[
M(q) \in \mathbb{R}^{n_q \times n_q}.
\]

外力为

\[
Q(q,v,t) \in \mathbb{R}^{n_q}.
\]

偏置力或内力项记为

\[
h(q,v,t) \in \mathbb{R}^{n_q}.
\]

如果当前原型只实现点质量或简单刚体，可以令

\[
h(q,v,t)=0.
\]

---

## 2. SDF 几何定义

设目标障碍物或目标体的符号距离场为

\[
\phi(x):\mathbb{R}^d \rightarrow \mathbb{R},
\]

通常二维中

\[
d=2,
\]

三维中

\[
d=3.
\]

本文采用以下符号约定：

\[
\phi(x)>0 \quad \text{表示点在目标体外部},
\]

\[
\phi(x)=0 \quad \text{表示点在目标体边界上},
\]

\[
\phi(x)<0 \quad \text{表示点在目标体内部，即发生穿透}.
\]

如果 \(\phi\) 是严格 SDF，则满足 Eikonal 条件：

\[
\|\nabla \phi(x)\|=1.
\]

实际数值实现中，需要允许 \(\|\nabla\phi\|\) 有微小误差。因此接触法向建议写为

\[
n(x)=\frac{\nabla \phi(x)}{\|\nabla \phi(x)\|+\delta},
\]

其中

\[
\delta \approx 10^{-12}.
\]

---

## 3. 接触采样点与无压缩候选集

设主动体表面有一组采样点。对于第 \(i\) 个采样点，其世界坐标为

\[
x_i(q) \in \mathbb{R}^d.
\]

一共有

\[
N_s
\]

个表面采样点。

### 3.1 无压缩原则

如果第 \(i\) 个点进入候选接触集，则必须保留以下独立变量：

\[
g_i,\quad n_i,\quad J_i,\quad \lambda_i.
\]

不得做以下操作：

- 不得将多个点聚类；
- 不得只保留最深穿透点；
- 不得将多个接触点合成一个代表点；
- 不得对接触约束降阶；
- 不得用一个乘子代表多个接触点；
- 不得把多个 gap 平均成一个 gap。

### 3.2 窄带筛选不是压缩

可以使用窄带筛选来剔除明显远离接触的点：

\[
\mathcal{C}_k = \{i\mid \phi(x_i(q_k)) < d_{\text{band}}\}.
\]

其中 \(d_{\text{band}}>0\) 是窄带半径。

进入 \(\mathcal{C}_k\) 的点全部保留。若

\[
N_c = |\mathcal{C}_k|,
\]

则该时间步中有 \(N_c\) 个接触候选点、\(N_c\) 个法向乘子和 \(N_c\) 条 NCP 方程。

建议第一版实现中：

1. 每个时间步开始时筛选一次候选集；
2. Newton 迭代过程中固定候选集；
3. 不在 Newton 迭代中改变候选点数量。

这样可以保证残差维度固定。

---

## 4. SDF gap 定义

对每个候选接触点 \(i\in\mathcal{C}_k\)，定义法向间隙：

\[
g_i(q)=\phi(x_i(q)).
\]

堆叠所有候选点：

\[
g(q)=
\begin{bmatrix}
g_1(q)\\
g_2(q)\\
\vdots\\
g_{N_c}(q)
\end{bmatrix}
\in \mathbb{R}^{N_c}.
\]

物理意义：

- \(g_i>0\)：点 \(i\) 与目标体分离；
- \(g_i=0\)：点 \(i\) 处于接触边界；
- \(g_i<0\)：点 \(i\) 穿透目标体。

---

## 5. 接触法向

对每个候选点：

\[
n_i(q)=\nabla\phi(x_i(q)).
\]

若数值 SDF 不严格满足单位梯度，则建议归一化：

\[
n_i(q)=\frac{\nabla\phi(x_i(q))}{\|\nabla\phi(x_i(q))\|+\delta}.
\]

注意：如果使用归一化后的法向计算 Jacobian，需要保持理论一致。严格地说，gap 的导数是

\[
\frac{\partial g_i}{\partial q}=\nabla\phi(x_i)^T\frac{\partial x_i}{\partial q}.
\]

因此：

- 若 `sdf.grad(x)` 返回真实 \(\nabla\phi\)，Jacobian 应使用真实梯度；
- 若 `sdf.grad(x)` 返回归一化法向，则相当于假设 \(\|\nabla\phi\|=1\)，这对严格 SDF 成立；
- 对数值 SDF，建议保留两个接口：`grad_phi(x)` 和 `normal(x)`。

推荐实现：

```python
raw_grad = sdf.grad(x)
normal = raw_grad / (np.linalg.norm(raw_grad) + 1e-12)
J_i = raw_grad.T @ dx_dq_i
```

若所有 SDF primitive 都是解析严格 SDF，则 `raw_grad` 和 `normal` 在接触区应几乎相同。

---

## 6. 接触 Jacobian 推导

### 6.1 一般形式

法向间隙为

\[
g_i(q)=\phi(x_i(q)).
\]

对 \(q\) 求导：

\[
J_i(q)=\frac{\partial g_i(q)}{\partial q}.
\]

由链式法则：

\[
J_i(q)=\nabla\phi(x_i(q))^T\frac{\partial x_i(q)}{\partial q}.
\]

其中

\[
J_i(q)\in\mathbb{R}^{1\times n_q}.
\]

堆叠所有候选点：

\[
J_c(q)=
\begin{bmatrix}
J_1(q)\\
J_2(q)\\
\vdots\\
J_{N_c}(q)
\end{bmatrix}
\in\mathbb{R}^{N_c\times n_q}.
\]

### 6.2 接触广义力

每个点的法向乘子为

\[
\lambda_i\ge 0.
\]

所有乘子堆叠：

\[
\lambda=
\begin{bmatrix}
\lambda_1\\
\lambda_2\\
\vdots\\
\lambda_{N_c}
\end{bmatrix}
\in\mathbb{R}^{N_c}.
\]

接触广义力为

\[
Q_c(q,\lambda)=J_c(q)^T\lambda.
\]

维度检查：

\[
J_c^T\in\mathbb{R}^{n_q\times N_c},
\]

\[
\lambda\in\mathbb{R}^{N_c},
\]

所以

\[
Q_c\in\mathbb{R}^{n_q}.
\]

逐点展开：

\[
Q_c=\sum_{i=1}^{N_c}J_i(q)^T\lambda_i.
\]

---

## 7. 2D 点质量的特殊情况

点质量广义坐标：

\[
q=x=\begin{bmatrix}x\\y\end{bmatrix}.
\]

此时

\[
x_i(q)=q.
\]

所以

\[
\frac{\partial x_i}{\partial q}=I_{2\times2}.
\]

因此

\[
J_i(q)=\nabla\phi(q)^T.
\]

如果只有一个点质量接触点，则

\[
J_c\in\mathbb{R}^{1\times2}.
\]

接触力：

\[
Q_c=J_c^T\lambda=\lambda\nabla\phi(q).
\]

---

## 8. 2D 刚体接触 Jacobian

二维刚体广义坐标：

\[
q=\begin{bmatrix}p_x&p_y&\theta\end{bmatrix}^T.
\]

局部表面采样点：

\[
r_i=\begin{bmatrix}r_{ix}\\r_{iy}\end{bmatrix}.
\]

世界坐标：

\[
x_i(q)=p+R(\theta)r_i,
\]

其中

\[
p=\begin{bmatrix}p_x\\p_y\end{bmatrix},
\]

\[
R(\theta)=
\begin{bmatrix}
\cos\theta & -\sin\theta\\
\sin\theta & \cos\theta
\end{bmatrix}.
\]

旋转导数：

\[
\frac{\partial}{\partial\theta}R(\theta)r_i
=
R(\theta)
\begin{bmatrix}
-r_{iy}\\
r_{ix}
\end{bmatrix}.
\]

因此

\[
\frac{\partial x_i}{\partial q}
=
\begin{bmatrix}
1 & 0 & \left[R(\theta)(-r_{iy},r_{ix})^T\right]_x\\
0 & 1 & \left[R(\theta)(-r_{iy},r_{ix})^T\right]_y
\end{bmatrix}.
\]

令

\[
n_i=\nabla\phi(x_i)=\begin{bmatrix}n_x\\n_y\end{bmatrix}.
\]

则

\[
J_i=
\begin{bmatrix}
n_x & n_y & n_i^T R(\theta)
\begin{bmatrix}
-r_{iy}\\
r_{ix}
\end{bmatrix}
\end{bmatrix}.
\]

实现中务必检查：

```python
J_i.shape == (3,)
# or as row:
J_i.reshape(1, 3)
```

---

## 9. Signorini 接触条件

无摩擦法向接触的 Signorini 条件为

\[
g_i(q)\ge 0,
\]

\[
\lambda_i\ge 0,
\]

\[
g_i(q)\lambda_i=0.
\]

写成互补形式：

\[
0\le g_i(q)\perp \lambda_i\ge0.
\]

对所有候选点：

\[
0\le g(q)\perp \lambda\ge0.
\]

含义：

1. 若 \(g_i>0\)，则 \(\lambda_i=0\)：分离点没有接触力；
2. 若 \(\lambda_i>0\)，则 \(g_i=0\)：有压力时必须处于接触边界；
3. \(g_i<0\) 不允许出现在收敛解中；
4. \(\lambda_i<0\) 不允许，因为接触不能提供拉力。

---

## 10. Fischer-Burmeister NCP 函数

标准 Fischer-Burmeister 函数：

\[
\Phi_{FB}(a,b)=\sqrt{a^2+b^2}-a-b.
\]

其性质是：

\[
\Phi_{FB}(a,b)=0
\quad\Longleftrightarrow\quad
0\le a\perp b\ge0.
\]

为便于普通 Newton、自动微分和数值实现，使用平滑形式：

\[
\Phi_\epsilon(a,b)=\sqrt{a^2+b^2+\epsilon^2}-a-b.
\]

其中

\[
\epsilon>0.
\]

对接触点 \(i\)：

\[
R_{\lambda,i}=\Phi_\epsilon(g_i(q),\lambda_i)
=\sqrt{g_i(q)^2+\lambda_i^2+\epsilon^2}-g_i(q)-\lambda_i.
\]

堆叠：

\[
R_\lambda(q,\lambda)=
\begin{bmatrix}
R_{\lambda,1}\\
R_{\lambda,2}\\
\vdots\\
R_{\lambda,N_c}
\end{bmatrix}.
\]

---

## 11. 平滑 FB 函数导数

令

\[
s_i=\sqrt{g_i^2+\lambda_i^2+\epsilon^2}.
\]

则

\[
R_{\lambda,i}=s_i-g_i-\lambda_i.
\]

对 \(g_i\) 的导数：

\[
\frac{\partial R_{\lambda,i}}{\partial g_i}
=\frac{g_i}{s_i}-1.
\]

对 \(\lambda_i\) 的导数：

\[
\frac{\partial R_{\lambda,i}}{\partial \lambda_i}
=\frac{\lambda_i}{s_i}-1.
\]

由于

\[
g_i=g_i(q),
\]

对 \(q\) 的导数为

\[
\frac{\partial R_{\lambda,i}}{\partial q}
=\left(\frac{g_i}{s_i}-1\right)J_i(q).
\]

如果时间步未知量使用 \(v_{k+1}\)，且

\[
q_{k+1}=q_k+\Delta t v_{k+1},
\]

则

\[
\frac{\partial R_{\lambda,i}}{\partial v_{k+1}}
=\left(\frac{g_i}{s_i}-1\right)J_i(q_{k+1})\Delta t.
\]

对所有点，\(R_\lambda\) 对 \(\lambda\) 的导数是对角矩阵：

\[
\frac{\partial R_\lambda}{\partial \lambda}
=\operatorname{diag}\left(\frac{\lambda_i}{s_i}-1\right).
\]

---

## 12. 动力学方程

连续时间动力学写为

\[
M(q)\dot v+h(q,v)=Q(q,v,t)+J_c(q)^T\lambda.
\]

或

\[
M(q)\ddot q+h(q,\dot q)=Q(q,\dot q,t)+J_c(q)^T\lambda.
\]

注意接触项的符号必须与 gap 定义和 Jacobian 定义一致。

如果

\[
g=\phi(x),\quad \phi>0 \text{ outside},\quad \nabla\phi \text{ outward from obstacle},
\]

则对于点质量落在平面 \(y=0\) 上方，取

\[
\phi(x,y)=y,
\]

\[
\nabla\phi=[0,1]^T.
\]

接触力应向上：

\[
Q_c=J^T\lambda=[0,1]^T\lambda.
\]

因此动力学中应写成

\[
M\dot v=Q+J^T\lambda.
\]

如果实现中出现接触力向下，通常是 SDF 符号或 Jacobian 符号错了。

---

## 13. 隐式 Euler 时间离散

给定时间步 \(k\)：

\[
q_k,\quad v_k.
\]

时间步长：

\[
\Delta t.
\]

未知量：

\[
v_{k+1},\quad \lambda_{k+1}.
\]

位置通过隐式 Euler 更新：

\[
q_{k+1}=q_k+\Delta t v_{k+1}.
\]

动力学离散为：

\[
M(q_{k+1})(v_{k+1}-v_k)
=
\Delta t\left[Q_{k+1}-h_{k+1}+J_c(q_{k+1})^T\lambda_{k+1}\right].
\]

移到左边定义残差：

\[
R_v=
M(q_{k+1})(v_{k+1}-v_k)
-
\Delta t\left[Q_{k+1}-h_{k+1}+J_c(q_{k+1})^T\lambda_{k+1}\right].
\]

如果当前原型没有 \(h\)：

\[
R_v=
M(v_{k+1}-v_k)
-
\Delta t\left[Q+J_c(q_{k+1})^T\lambda_{k+1}\right].
\]

---

## 14. 完整无压缩时间步非线性系统

定义时间步未知量：

\[
z=
\begin{bmatrix}
v_{k+1}\\
\lambda_{k+1}
\end{bmatrix}
\in\mathbb{R}^{n_q+N_c}.
\]

其中

\[
\lambda_{k+1}\in\mathbb{R}^{N_c}.
\]

完整残差：

\[
R(z)=
\begin{bmatrix}
R_v(z)\\
R_\lambda(z)
\end{bmatrix}
=
0.
\]

维度：

\[
R_v\in\mathbb{R}^{n_q},
\]

\[
R_\lambda\in\mathbb{R}^{N_c},
\]

\[
R\in\mathbb{R}^{n_q+N_c}.
\]

这是无压缩 SDF-NCP 求解的核心系统。

---

## 15. 残差装配流程

每次给定 \(z\)，执行以下操作。

### 15.1 拆分未知量

\[
v = z[0:n_q],
\]

\[
\lambda = z[n_q:n_q+N_c].
\]

### 15.2 更新位置

\[
q = q_k+\Delta t v.
\]

### 15.3 对每个候选点计算几何量

对 \(i=1,\dots,N_c\)：

\[
x_i=x_i(q),
\]

\[
g_i=\phi(x_i),
\]

\[
\nabla\phi_i=\nabla\phi(x_i),
\]

\[
D_i=\frac{\partial x_i}{\partial q},
\]

\[
J_i=\nabla\phi_i^T D_i.
\]

### 15.4 堆叠

\[
g=[g_i],
\]

\[
J_c=\begin{bmatrix}J_i\end{bmatrix}.
\]

### 15.5 动力学残差

\[
R_v=M(v-v_k)-\Delta t(Q+J_c^T\lambda)
\]

若有 \(h\)：

\[
R_v=M(v-v_k)-\Delta t(Q-h+J_c^T\lambda).
\]

### 15.6 NCP 残差

\[
R_{\lambda,i}=\sqrt{g_i^2+\lambda_i^2+\epsilon^2}-g_i-\lambda_i.
\]

### 15.7 返回联合残差

\[
R=\begin{bmatrix}R_v\\R_\lambda\end{bmatrix}.
\]

---

## 16. 解析 Jacobian 的结构

第一版可以使用有限差分 Jacobian，但为了检查实现和后续优化，需要知道解析结构。

令未知量为

\[
z=[v,\lambda].
\]

Jacobian 块结构：

\[
K=\frac{\partial R}{\partial z}
=
\begin{bmatrix}
\frac{\partial R_v}{\partial v} & \frac{\partial R_v}{\partial \lambda}\\
\frac{\partial R_\lambda}{\partial v} & \frac{\partial R_\lambda}{\partial \lambda}
\end{bmatrix}.
\]

### 16.1 动力学残差对 \(\lambda\)

若

\[
R_v=M(v-v_k)-\Delta t(Q+J^T\lambda),
\]

且暂时忽略 \(J\) 对 \(\lambda\) 的依赖，则

\[
\frac{\partial R_v}{\partial \lambda}=-\Delta t J^T.
\]

维度：

\[
J^T\in\mathbb{R}^{n_q\times N_c}.
\]

### 16.2 NCP 残差对 \(\lambda\)

\[
\frac{\partial R_\lambda}{\partial \lambda}
=\operatorname{diag}\left(\frac{\lambda_i}{s_i}-1\right).
\]

### 16.3 NCP 残差对 \(v\)

\[
\frac{\partial R_\lambda}{\partial v}
=
\operatorname{diag}\left(\frac{g_i}{s_i}-1\right)J\Delta t.
\]

因为

\[
\frac{\partial g}{\partial v}=J\frac{\partial q}{\partial v}=J\Delta t.
\]

### 16.4 动力学残差对 \(v\)

完整形式包含

- \(M(q)\) 对 \(q\) 的导数；
- \(Q(q,v)\) 对 \(q,v\) 的导数；
- \(J(q)^T\lambda\) 对 \(q\) 的导数；
- \(h(q,v)\) 对 \(q,v\) 的导数。

第一版点质量常质量矩阵、常外力、无 \(h\) 时，可近似为

\[
\frac{\partial R_v}{\partial v}\approx M.
\]

如果考虑 \(J(q)^T\lambda\) 对 \(v\) 的依赖：

\[
\frac{\partial}{\partial v}\left(J(q)^T\lambda\right)
=
\frac{\partial}{\partial q}\left(J(q)^T\lambda\right)\Delta t.
\]

该项涉及 SDF Hessian 和运动学二阶导数。第一版可用有限差分整体 Jacobian 避免手推错误。

---

## 17. Newton / LM 求解流程

给定初值 \(z^{(0)}\)。对 \(m=0,1,2,\dots\)：

1. 计算残差

\[
R^{(m)}=R(z^{(m)}).
\]

2. 若

\[
\|R^{(m)}\|<\text{tol}_R,
\]

则收敛。

3. 构造 Jacobian

\[
K^{(m)}=\frac{\partial R}{\partial z}.
\]

4. 求解 Newton 步

\[
K^{(m)}\Delta z=-R^{(m)}.
\]

如果 \(K\) 病态，使用 LM 步：

\[
(K^TK+\eta I)\Delta z=-K^TR.
\]

5. 线搜索：

\[
z^{(m+1)}=z^{(m)}+\alpha\Delta z,
\]

选择 \(\alpha\in(0,1]\) 使残差下降。

---

## 18. 初值策略

速度预测：

\[
v_{k+1}^{(0)}=v_k+\Delta t M^{-1}Q.
\]

位置预测：

\[
q_{k+1}^{(0)}=q_k+\Delta t v_{k+1}^{(0)}.
\]

乘子初值可以选择：

### 18.1 零初值

\[
\lambda_i^{(0)}=0.
\]

优点：简单。缺点：接触冲击时可能收敛慢。

### 18.2 罚函数初值

\[
\lambda_i^{(0)}=k_{init}\max(0,-g_i(q_{k+1}^{(0)})).
\]

优点：更容易收敛。注意最终方程仍是 NCP，不是 penalty。

推荐第一版使用 penalty warm start。

---

## 19. 收敛判据

建议同时检查以下量：

### 19.1 联合残差

\[
\|R\|_2<\text{tol}_R.
\]

### 19.2 最大穿透

\[
\max_i\max(0,-g_i)<\text{tol}_g.
\]

### 19.3 最大互补误差

\[
\max_i |g_i\lambda_i|<\text{tol}_c.
\]

### 19.4 乘子非负性

\[
\min_i \lambda_i > -\text{tol}_\lambda.
\]

### 19.5 NCP 残差

\[
\|R_\lambda\|_\infty<\text{tol}_{NCP}.
\]

---

## 20. 诊断量定义

每个时间步应保存以下标量：

\[
\text{max\_penetration}=\max_i\max(0,-g_i),
\]

\[
\text{mean\_penetration}=\frac{1}{N_c}\sum_i\max(0,-g_i),
\]

\[
\text{max\_complementarity}=\max_i |g_i\lambda_i|,
\]

\[
\text{ncp\_residual\_inf}=\|R_\lambda\|_\infty,
\]

\[
\text{total\_normal\_force}=\sum_i \lambda_i.
\]

对刚体或多体系统，更重要的是总广义接触力：

\[
Q_c=J^T\lambda.
\]

---

## 21. 面积权重的可选扩展

如果候选点来自表面积分点，可引入权重 \(w_i\)。

此时 \(\lambda_i\) 更接近接触压力，接触广义力为

\[
Q_c=\sum_i w_i J_i^T\lambda_i.
\]

矩阵形式：

\[
Q_c=J^T W\lambda,
\]

其中

\[
W=\operatorname{diag}(w_i).
\]

NCP 条件仍为

\[
0\le g_i\perp \lambda_i\ge0.
\]

第一版可以不启用权重，或令

\[
w_i=1.
\]

但代码结构可以预留 `weights`。

---

## 22. 无压缩多点接触的不唯一性

当多个接触点同时满足 \(g_i=0\) 时，乘子分布 \(\lambda_i\) 可能不唯一。尤其是刚体平面接触刚性平面时，多个点可共同承担同一个总反力。

这不是实现错误，而是无压缩点接触离散的自然现象。

如果需要稳定化，可以使用：

1. LM 正则化；
2. 弱乘子范数正则；
3. 面积权重；
4. 表面采样均匀化；
5. 使用上一时间步 \(\lambda\) warm start。

第一篇基础版本可以讨论该限制，但不要引入压缩。

---

## 23. Penalty baseline

罚函数基线用于对比，不是主方法。

定义：

\[
\lambda_i^{pen}=k_n\max(0,-g_i).
\]

接触广义力：

\[
Q_c^{pen}=J^T\lambda^{pen}.
\]

如果加入阻尼：

\[
\dot g_i=J_i v,
\]

\[
\lambda_i^{pen}=k_n\max(0,-g_i)+c_n\max(0,-\dot g_i)\mathbf{1}_{g_i<0}.
\]

第一篇建议先实现无阻尼 penalty，避免额外参数。

Penalty 与 NCP 的核心差别：

- penalty 允许 \(g_i<0\)，用刚度惩罚穿透；
- NCP 通过互补约束逼近 \(g_i\ge0\)。

---

## 24. 点质量-平面接触的符号验证

这是最重要的单元测试。

平面：

\[
y=0.
\]

SDF：

\[
\phi(x,y)=y.
\]

梯度：

\[
\nabla\phi=[0,1]^T.
\]

点质量：

\[
q=[x,y]^T.
\]

Jacobian：

\[
J=[0,1].
\]

重力：

\[
Q=[0,-mg]^T.
\]

动力学残差：

\[
R_v=m(v_{k+1}-v_k)-\Delta t([0,-mg]^T + [0,1]^T\lambda).
\]

如果接触发生，\(\lambda>0\)，则接触力为向上。

Codex 必须检查：

```python
assert np.allclose(sdf.phi(np.array([0.0, 1.0])), 1.0)
assert np.allclose(sdf.grad(np.array([0.0, 1.0])), np.array([0.0, 1.0]))
assert J.shape == (1, 2)
assert np.allclose(J, np.array([[0.0, 1.0]]))
```

---

## 25. 常见实现错误与修正

### 错误 1：接触力符号反了

症状：点质量落到平面后被继续向下拉，穿透越来越大。

检查：

\[
Q_c=J^T\lambda
\]

是否被写成

\[
-Q_c.
\]

正确：若 \(\phi=y\)，\(J=[0,1]\)，则 \(J^T\lambda=[0,\lambda]^T\)。

---

### 错误 2：gap 使用了负号

错误写法：

\[
g=-\phi(x).
\]

如果 SDF 约定外部为正，则正确是：

\[
g=\phi(x).
\]

除非整个代码统一采用相反约定。第一篇建议统一使用外部为正。

---

### 错误 3：Jacobian 维度错

正确：

\[
J\in\mathbb{R}^{N_c\times n_q}.
\]

错误常见为：

```python
J.shape == (nq, Nc)
```

修正：每个 \(J_i\) 是一行，最终 `np.vstack(rows)`。

---

### 错误 4：用 normal 代替 gap 导数导致不一致

如果 `normal = grad / norm(grad)`，而 SDF 不是严格单位梯度，则

\[
J_i=normal^T dx/dq
\]

不是严格的 \(\partial g/\partial q\)。

修正：使用 raw gradient 构造 \(J_i\)，normal 只用于可视化或摩擦切向投影。

---

### 错误 5：Newton 过程中候选集变化

症状：残差维度变化，Newton 崩溃。

修正：每个时间步开始确定候选集，Newton 内固定。

---

### 错误 6：所有点都压缩成一个最深点

这违反无压缩版本。

修正：所有候选点均创建独立 \(\lambda_i\) 和 NCP 残差。

---

### 错误 7：NCP 残差没有逐点对应

错误：

\[
\Phi(\sum_i g_i, \sum_i \lambda_i)=0.
\]

正确：

\[
\Phi(g_i,\lambda_i)=0,\quad i=1,\dots,N_c.
\]

---

### 错误 8：平滑 FB 函数写错符号

正确：

\[
\Phi_\epsilon(g,\lambda)=\sqrt{g^2+\lambda^2+\epsilon^2}-g-\lambda.
\]

不是

\[
g+\lambda-\sqrt{g^2+\lambda^2+\epsilon^2}.
\]

后者只是负号相反，理论上根相同，但导数和收敛实现要一致。建议统一使用前者。

---

### 错误 9：平滑参数过大

如果 \(\epsilon\) 太大，互补条件会被明显放松。

建议测试：

\[
\epsilon=10^{-2},10^{-4},10^{-6}.
\]

第一版默认：

\[
\epsilon=10^{-6}\text{ 或 }10^{-5}.
\]

---

### 错误 10：未检查 \(\lambda\ge0\)

平滑 FB 理论上逼近非负约束，但数值误差可能导致微小负值。

诊断中必须保存：

```python
lambda_min = np.min(lambda)
```

并检查：

```python
assert lambda_min > -tol_lambda
```

---

## 26. Codex 实现检查清单

Codex 应逐项检查当前代码。

### 26.1 SDF 检查

- `PlaneSDF.phi([0,1]) == 1` for plane `y=0`；
- `PlaneSDF.grad([0,1]) == [0,1]`；
- `CircleSDF.phi(center + [radius,0]) == 0`；
- `CircleSDF.grad` 在非中心点单位化；
- `BoxSDF` 内部负、边界零、外部正。

### 26.2 Gap 检查

- `g_i = sdf.phi(x_i)`；
- 不允许默认取负号；
- 所有候选点都有 gap。

### 26.3 Jacobian 检查

- 点质量：`J = grad_phi.reshape(1, d)`；
- 刚体：`J_i = grad_phi.T @ dx_dq_i`；
- `J.shape == (Nc, nq)`；
- 用有限差分验证 `J_i`。

### 26.4 NCP 检查

- `smooth_fb(g, lam, eps) = sqrt(g*g + lam*lam + eps*eps) - g - lam`；
- 支持 scalar 和 vector；
- `dPhi_dg = g/s - 1`；
- `dPhi_dlam = lam/s - 1`；
- NCP 残差逐点计算。

### 26.5 残差检查

- 未知量是 `[v_next, lambda]`；
- `q_next = q_k + dt * v_next`；
- `R_dyn = M @ (v_next - v_k) - dt * (Q - h + J.T @ lambda)`；
- 若无 `h`，使用 `Q + J.T @ lambda`；
- `R_ncp.shape == (Nc,)`；
- `R.shape == (nq + Nc,)`。

### 26.6 无压缩检查

- 不得出现 `deepest_only=True` 作为默认；
- 不得只保留 `argmin(g)`；
- 不得聚类 contact points；
- 不得把多个 gap 取平均；
- 不得用单个 lambda 代表多个点；
- `len(lambda) == len(candidates)`。

### 26.7 时间步检查

- 每个时间步开始筛选候选集；
- Newton 内候选集固定；
- 若无候选点，则执行自由运动步；
- 若有候选点，求解联合残差。

---

## 27. 推荐代码接口

### 27.1 SDF 接口

```python
class SDF:
    def phi(self, x: np.ndarray) -> float:
        ...

    def grad(self, x: np.ndarray) -> np.ndarray:
        ...

    def normal(self, x: np.ndarray) -> np.ndarray:
        g = self.grad(x)
        return g / (np.linalg.norm(g) + 1e-12)
```

### 27.2 系统接口

```python
class MechanicalSystem:
    nq: int

    def mass_matrix(self, q: np.ndarray) -> np.ndarray:
        ...

    def external_force(self, q: np.ndarray, v: np.ndarray, t: float) -> np.ndarray:
        ...

    def bias_force(self, q: np.ndarray, v: np.ndarray, t: float) -> np.ndarray:
        return np.zeros(self.nq)

    def world_point(self, q: np.ndarray, local_point: np.ndarray) -> np.ndarray:
        ...

    def world_point_jacobian(self, q: np.ndarray, local_point: np.ndarray) -> np.ndarray:
        ...
```

### 27.3 接触装配接口

```python
def assemble_contact(system, sdf, q, candidates):
    gaps = []
    normals = []
    rows = []

    for r_i in candidates:
        x_i = system.world_point(q, r_i)
        g_i = sdf.phi(x_i)
        grad_i = sdf.grad(x_i)
        normal_i = grad_i / (np.linalg.norm(grad_i) + 1e-12)
        dx_dq_i = system.world_point_jacobian(q, r_i)
        J_i = grad_i.reshape(1, -1) @ dx_dq_i

        gaps.append(g_i)
        normals.append(normal_i)
        rows.append(J_i.reshape(-1))

    g = np.asarray(gaps)
    J = np.vstack(rows) if rows else np.zeros((0, system.nq))
    normals = np.vstack(normals) if normals else np.zeros((0, system.spatial_dim))
    return g, J, normals
```

### 27.4 残差接口

```python
def residual_uncompressed(z, q_k, v_k, dt, system, sdf, candidates, eps, t_next):
    nq = system.nq
    v_next = z[:nq]
    lam = z[nq:]

    q_next = q_k + dt * v_next
    g, J, normals = assemble_contact(system, sdf, q_next, candidates)

    M = system.mass_matrix(q_next)
    Q = system.external_force(q_next, v_next, t_next)
    h = system.bias_force(q_next, v_next, t_next)

    R_dyn = M @ (v_next - v_k) - dt * (Q - h + J.T @ lam)
    R_ncp = smooth_fischer_burmeister(g, lam, eps)

    return np.concatenate([R_dyn, R_ncp])
```

---

## 28. 无候选点情况

如果

\[
N_c=0,
\]

则没有 \(\lambda\) 和 NCP 方程。

此时求解自由运动：

\[
M(v_{k+1}-v_k)=\Delta t(Q-h).
\]

点质量常质量情况下：

\[
v_{k+1}=v_k+\Delta t M^{-1}Q.
\]

\[
q_{k+1}=q_k+\Delta t v_{k+1}.
\]

代码中要避免 `np.vstack([])` 报错。

---

## 29. 候选集生成推荐实现

```python
def select_candidates(system, sdf, q_k, surface_points, d_band=None):
    candidates = []
    for r_i in surface_points:
        x_i = system.world_point(q_k, r_i)
        g_i = sdf.phi(x_i)
        if d_band is None:
            candidates.append(r_i)
        else:
            if g_i < d_band:
                candidates.append(r_i)
    return candidates
```

注意：

- 这里没有压缩；
- 没有最深点筛选；
- 没有聚类；
- 所有满足窄带条件的点都保留。

---

## 30. 推荐测试：单点接触

### 30.1 设置

点质量：

\[
q_0=[0,1]^T,
\]

\[
v_0=[0,0]^T,
\]

平面 SDF：

\[
\phi(x,y)=y.
\]

重力：

\[
g=9.81.
\]

### 30.2 预期行为

- 落体接近 \(y=0\)；
- NCP 解中最大穿透应很小；
- 接触后 \(\lambda\ge0\)；
- \(g\lambda\) 应接近零；
- contact force 应向上。

---

## 31. 推荐测试：多点无压缩

构造一个 2D 刚体矩形，下边缘采样多个点，与平面接触。

假设有 5 个采样点进入窄带，则必须满足：

```python
assert len(candidates) == 5
assert lam.shape == (5,)
assert gaps.shape == (5,)
assert J.shape == (5, 3)
assert R_ncp.shape == (5,)
assert residual.shape == (3 + 5,)
```

不能出现：

```python
lam.shape == (1,)
```

除非候选点本来只有一个。

---

## 32. 推荐测试：NCP 点态性质

对于小 \(\epsilon\)：

1. 分离状态：

\[
g=1,\quad \lambda=0.
\]

则

\[
\Phi_\epsilon(g,\lambda)\approx0.
\]

2. 接触承压状态：

\[
g=0,\quad \lambda=1.
\]

则

\[
\Phi_\epsilon(g,\lambda)\approx0.
\]

3. 非法穿透无力状态：

\[
g=-1,\quad \lambda=0.
\]

则残差应显著非零。

4. 非法拉力状态：

\[
g=1,\quad \lambda=-1.
\]

则残差应显著非零。

---

## 33. 推荐实验输出

每个时间步输出 CSV：

```text
time
q_0, q_1, ...
v_0, v_1, ...
num_candidates
max_penetration
mean_penetration
max_complementarity
ncp_residual_inf
lambda_min
lambda_max
solver_success
solver_iterations
```

对每个接触点输出长表：

```text
time
contact_id
gap
lambda
normal_x
normal_y
penetration
complementarity
ncp_residual
```

这对无压缩方法非常重要，因为论文中可以展示接触力分布。

---

## 34. 第一篇论文算法描述

可以在论文中写成如下算法。

### Algorithm: Uncompressed SDF-NCP Time Step

**Input:** \(q_k,v_k,\Delta t,\phi,\{r_i\}_{i=1}^{N_s},d_{band}\)  
**Output:** \(q_{k+1},v_{k+1},\lambda_{k+1}\)

1. Compute current world points \(x_i(q_k)\).
2. Select candidate set
   \[
   \mathcal{C}_k=\{i\mid \phi(x_i(q_k))<d_{band}\}.
   \]
3. Keep all points in \(\mathcal{C}_k\) without clustering or reduction.
4. Initialize \(v_{k+1}^{(0)}\) and \(\lambda^{(0)}\).
5. For Newton iteration \(m\):
   1. Set \(q_{k+1}^{(m)}=q_k+\Delta t v_{k+1}^{(m)}\).
   2. For every \(i\in\mathcal{C}_k\), compute \(g_i,n_i,J_i\).
   3. Assemble \(J_c\).
   4. Compute \(R_v\).
   5. Compute \(R_{\lambda,i}=\Phi_\epsilon(g_i,\lambda_i)\).
   6. Solve Newton or LM correction.
   7. Apply line search.
6. Accept converged \(q_{k+1},v_{k+1},\lambda_{k+1}\).
7. Store all per-contact-point diagnostics.

---

## 35. Codex 修正任务提示词

可以将下面提示词直接交给 Codex。

```text
Read docs/sdf_ncp_uncompressed_theory_and_codex_guide.md carefully.
Use it as the source of truth for the current first-paper implementation.

Your task is to inspect and correct the implementation so that it follows the uncompressed SDF-NCP contact formulation.

Mandatory requirements:
1. Use SDF gap g_i = phi(x_i). Do not flip the sign unless the entire codebase explicitly changes the convention.
2. Use the SDF gradient to construct each contact Jacobian row:
       J_i = grad_phi(x_i)^T @ dx_i_dq
3. Stack all contact Jacobian rows so that:
       J.shape == (num_candidates, nq)
4. Use one independent multiplier lambda_i per candidate point:
       lambda.shape == (num_candidates,)
5. Do not cluster, merge, average, or compress contact points.
6. Do not keep only the deepest penetration point.
7. A narrow-band filter is allowed, but every point inside the band must be retained.
8. Keep the candidate set fixed during Newton iterations within a time step.
9. Use the smooth Fischer-Burmeister NCP residual componentwise:
       Phi_i = sqrt(g_i^2 + lambda_i^2 + eps^2) - g_i - lambda_i
10. The time-step unknown vector must be:
       z = [v_next, lambda_1, ..., lambda_Nc]
11. The implicit Euler position update must be:
       q_next = q_k + dt * v_next
12. The dynamics residual must be:
       R_dyn = M @ (v_next - v_k) - dt * (Q - h + J.T @ lambda)
    If h is not implemented, use:
       R_dyn = M @ (v_next - v_k) - dt * (Q + J.T @ lambda)
13. The full residual must have shape:
       (nq + num_candidates,)
14. Add or fix tests for:
       - Plane SDF sign and gradient
       - Circle SDF boundary and gradient
       - SDF gap sign convention
       - contact Jacobian finite-difference check
       - smooth FB function and derivative
       - point mass plane contact
       - multi-point uncompressed contact dimensions
15. Save diagnostics for all contact points:
       gap_i, lambda_i, normal_i, penetration_i, complementarity_i, ncp_residual_i

After changes, run:
    pytest
    python scripts/run_all_experiments.py

Then report:
1. files changed
2. equations implemented
3. tests added or fixed
4. commands run
5. remaining limitations
```

---

## 36. 当前第一版应避免的复杂内容

为了保持第一篇清晰，第一版不要实现：

- AI；
- 神经网络；
- PyTorch/JAX；
- 摩擦锥；
- stick-slip；
- 自适应接触压缩；
- 接触点聚类；
- 动态 SDF 学习；
- 磨损演化；
- 大规模稀疏并行优化。

这些可以留给第二篇或后续工作。

---

## 37. 参考文献线索

以下文献方向可用于论文背景，但代码实现不依赖它们：

1. Signorini 接触条件与互补接触问题；
2. Fischer-Burmeister 函数与 NCP 重构；
3. 平滑 Fischer-Burmeister 方法与 smoothing Newton；
4. SDF 在接触检测和复杂几何接触中的应用；
5. 动态摩擦less contact 和互补条件时间离散。

建议论文中把本文方法定位为：

> SDF provides the continuous geometric gap; NCP provides the non-penetration mechanics; the uncompressed discretization provides a direct pointwise contact formulation without contact reduction.

---

## 38. 最终核心公式总结

第一篇最核心的四个公式是：

\[
g_i(q)=\phi(x_i(q)),
\]

\[
J_i(q)=\nabla\phi(x_i(q))^T\frac{\partial x_i(q)}{\partial q},
\]

\[
0\le g_i(q)\perp\lambda_i\ge0,
\]

\[
\Phi_\epsilon(g_i,\lambda_i)=\sqrt{g_i^2+\lambda_i^2+\epsilon^2}-g_i-\lambda_i=0.
\]

时间步残差：

\[
R(z)=
\begin{bmatrix}
M(v_{k+1}-v_k)-\Delta t(Q-h+J^T\lambda)\\
\Phi_\epsilon(g(q_{k+1}),\lambda)
\end{bmatrix}
=0,
\]

其中

\[
q_{k+1}=q_k+\Delta t v_{k+1}.
\]

这就是无压缩 SDF-NCP 接触动力学第一篇论文和代码实现的理论基础。
