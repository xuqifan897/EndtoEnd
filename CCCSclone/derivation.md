# Derivation of the CCCS algorithm
Before diving deeper into the existing dose calculation code, I decide to firstly re-familiarize myself with the formulations.
## DK, CK, and CCK algorithms
In this [paper](http://dx.doi.org/10.1088/0031-9155/50/4/007), Weiguo _et al._ proposed these three dose calculation algorithms, among which the CCK algorithm is the base of the GPU dose calculation [paper](http://dx.doi.org/10.1118/1.3551996). The most general convolution/superposition algorithm is defined as:
$$
\tag{1}
D(\vec{x})=\frac{\iiint T(\vec{s})\rho(\vec{s})\left[\eta(\vec{x})\bar{\eta}^2h(\vec{r}(\vec{s};\vec{x}))\right]dV_s}{\rho(\vec{x})}
$$
For a volume with $N^3$ voxels, each source (Terma) voxel will contribute energy to all $N^3$ destination (dose) voxels. Therefore, the complexity is of $O(N^6)$.

Instead of the full volume C/S, the collapsed cone C/S is more often used in radiotherapy TPS. In this regime, the energy deposition kernel is collapsed to a discrete number of directions. For a given collapsed cone direction $(\theta_m, \phi_n)$, the collapsed cone kernel is defined as:
$$
\tag{2}
k^{m,n}(r)=\int_{\Omega^{m,n}}||\vec{r}||^2h(\vec{r})\sin \theta d\theta d\phi,
$$
where $\Omega^{m,n}$ is the solid angle for the collapsed cone direction $(\theta_m, \phi_n)$. In this sense, the dose is the summation of contribution from each independent collapsed cone:
$$
\tag{3}
D(\vec{x})=\sum_{m,n}D^{m,n}(\vec{x}).
$$
We discard the superscripts $m, n$ to simplify the notation hereafter. In a single collapsed cone direction, the dose deposition can be formulated as:
$$
\tag{4}
D(x)=\int_0^xT(s)k(r(s;x))\eta(s)ds,
$$
where the 0 in the integration denotes the beam entrance. If we treat $x$ as a constant, we define $r_x(s)=r(s;x)$. Then we have
$$
\tag{5}
\eta(x)=-\frac{dr_x}{ds}.
$$
So equation 4 can be rewritten as:
$$
\tag{6}
D(x)=\int_{s=x}^0 T(s)k(r_x(s))dr_x(s).
$$
If we then further define a cumulative kernel as:
$$
\tag{7}
K(r_x)=\int_0^{r_x}k(t)dt,
$$
Then equation 6 can be rewritten as:
$$
\tag{8}
D(x)=\int_{s=x}^0 T(s)dK(r_x(s)).
$$

### Discretization
As the dose calculation is implemented on voxelated phantom, we have to convert the continuous integration to discrete sum. For this purpose, we define $r_{i,j}$ to be the radiological distance from the center of the segment $\Delta x_i$ to the center of the segment $\Delta x_j$ ($i\geq j$). Also, we define $p_i$ to be the radiological path length across voxel $i$. We thus have:
$$
\tag{9}
r_{i+1,j}=r_{i,j}+\frac{p_i+p_{i+1}}{2},
$$
$$
\tag{10}
r_{i, j+1} = r_{i,j} -\frac{p_j+p_{j+1}}{2}.
$$
We firstly use the dose to the center of segment $\Delta x_i$ to approximate the dose to the segment itself. By directly discretizing equations 6 and 8, we can get:
$$
\tag{11}
\frac{p_i}{2}k(\frac{p_i}{4})T_i+\sum_{j=i-1}^0T_jp_jk(r_{i,j})
$$
$$
\tag{12}
T_iK(\frac{p_i}{2})+\sum_{j=i-1}^0T_j\left[K(r_{i,j})-K(r_{i,j+1})\right]
$$
which are referred to as DK algorithm (differential collapsed kernel) and CK algorithm (cumulative collapsed kernel), respectively. (<span style="color: blue;">I somehow think that the term in the bracket is inaccurate. It should be </span> $K(r_{i,j}+p_j/2)-K(r_{i,j+1}+p_{j+1}/2)$)

Now, instead of using the dose to the center of the segment $\Delta x_i$, we use the average dose to it. Suppose along a certain energy transportation line, we consider two volume segments, $\Delta x_i$ and $\Delta x_j$, with Terma $T_i$ and $T_j$, density $\rho_i$ and $\rho_j$, respectively. In this form, the energy to $\Delta x_i$ from $\Delta x_j$ is defined as $E_{i,j}$. The total dose to $\Delta x_i$ is the sum of contributions from all segments:
$$
\tag{13}
D_i=\frac{1}{\rho_i\Delta S\Delta x_i}\sum_{j=i}^{0}E_{i,j}=\frac{1}{\rho_w p_i \Delta S}\sum_{j=i}^{0}E_{i,j}.
$$
For local deposition, i.e., for segment $\Delta x_i$ to deposit energy to itself, $j=i$,
$$
\begin{align}
E_{i,i}&=\int_{-p_i/2}^{p_i/2}D(s)\rho_w \Delta S ds \\
&=\rho_w\Delta S\int_{-p_i/2}^{p_i/2}D(s)ds \\
&=\rho_w\Delta S\int_{-p_i/2}^{p_i/2}ds\int_{t=s}^{p_i/2}T_idK(t-s) \\
&=\rho_w\Delta S\int_{-p_i/2}^{p_i/2}T_iK(p_i/2-s)ds \\
&=\rho_w\Delta ST_i\int_{0}^{p_i}K(q)dq.
\end{align}
\tag{14}
$$
For remote deposition ($j<i$), we have
$$
\begin{align}
E_{i,j} &= \rho_w\Delta S\int_{-p_i/2}^{p_i/2}D(s)ds \\
&= \rho_w\Delta S\int_{-p_i/2}^{p_i/2}ds \int_{r_{i,j}-p_j/2-s}^{r_{i,j}+p_j/2-s}T_jdK(q) \\
&= \rho_w \Delta S\int_{-p_i/2}^{p_i/2}dsT_j\left(K(r_{i,j}+p_j/2-s)-K(r_{i,j}-p_j/2-s)\right) \\
&=\rho_w \Delta ST_j\left(\int_{r_{i,j}+p_j/2-p_i/2}^{r_{i,j}+p_j/2+p_i/2} K(q)dq - \int_{r_{i,j}-p_j/2-p_i/2}^{r_{i,j}-p_j/2+p_i/2}K(q)dq\right).
\end{align}
\tag{15}
$$
If we define $q_{i,j}=r_{i,j}-p_i/2+p_j/2$, then
$$
\begin{split}
q_{i+1,j} &= r_{i+1,j}-p_{i+1}/2+p_j/2=r_{i,j}+p_i/2+p_j/2 \\
q_{i+1,j+1} &= r_{i+1,j+1}-p_{i+1}/2+p_{j+1}/2=r_{i,j+1}+p_i/2+p_{j+1}/2=r_{i,j}-p_j/2+p_i/2 \\
q_{i,j+1} &= r_{i,j+1}-p_i/2 + p_{j+1}/2=r_{i,j}-p_j/2-p_i/2
\end{split}
\tag{16}
$$
We then define the cumulative-cumulative collapsed cone kernel, $C(r)=\int_o^rK(t)dt$, equation 15 can be rewritten as:
$$
\begin{split}
E_{i,j} &= \rho_w \Delta ST_j\left[(C(q_{i+1,j})-C(q_{i,j}))-(C(q_{i+1,j+1})-C(q_{i,j+1}))\right] \\
&= \rho_w\Delta ST_j(c_{i,j}-c_{i,j+1}),
\end{split}
\tag{17}
$$
where $c_{i,j}=C(q_{i+1,j})-C(q_{i,j})$, $i\geq j$.
$$
\begin{align}
D_i &= \frac{1}{\rho_w p_i\Delta S}\left[ E_{i,i} + \sum_{j=i-1}^0 E_{i,j} \right] \\
&= \frac{1}{\rho_w p_i\Delta S}\left[ \rho_w \Delta S T_i C(p_i) + \rho_w\Delta S\sum_{j=i-1}^{0}T_j(c_{i,j}-c_{i,j+1})\right] \\
&= \frac{1}{p_i}T_iC(p_i) + \frac{1}{p_i}\sum_{j=i-1}^{0} T_j[c_{i,j}-c_{i,j+1}]
\end{align}
\tag{18}
$$
### Relationship between the DK, CK, and CCK algorithms
If we assume $K(r)$ is constant across one voxel, we have:
$$
\tag{19}
c_{i,j} = \int_{q_{i,j}}^{q_{i+1,j}}K(q)dq \approx K(r_{i,j}+p_j/2)p_i,
$$
$$
\tag{20}
D_i = T_iK(\frac{p_i}2)+\sum_{j=i-1}^0 T_j\left[ K(r_{i,j}+p_j/2) - K(r_{i,j}-p_j/2) \right],
$$
which means that CK (equation 12) can be viewed as an approximate of CCK. If we again assume that the differential kernel $k(r)$ is constant across one voxel, equation 19 can be further approximated as:
$$
D_i = \frac{p_i}{2}k(\frac{p_i}{4})T_i + \sum_{j=i-1}^0T_jp_jk(r_{i,j}),
$$
which is the DK algorithm (equation 11).

## Acceleration with exponential kernel
If we approximate the collapsed-cone kernel to be:
$$
\tag{21}
k^{m,n}(r)=A_{\theta}e^{-a_\theta r} + B_\theta e^{-b_\theta r}.
$$
For simplicity, we start with one component and simplify the notation: $k(r)=Ae^{-ar}$, then the cumulative kernel and the cumulative-cumulative kernel become
$$
\tag{22}
K(r)=\int_0^rk(t)dt = \frac{A}{a}(1-e^{-ar}),
$$
$$
\tag{23}
C(r)=\int_0^rK(t)dt=\frac{A}{a}r-\frac{A}{a^2}
(1-e^{-ar})$$
Substitute equation (22) to the definition of $c_{i,j}$, we get:
$$
c_{i,j} = \frac{A}{a}(q_{i+1,j}-q_{i,j})+\frac{A}{a^2}(e^{-aq_{i+1,j}}-e^{-aq_{i,j}}),
\tag{24}
$$
$$
c_{i,j+1} = \frac{A}{a}(p_{i+1}+p_i)/2+\frac{A}{a^2}(e^{-aq_{i+1,j+1}}-e^{-aq_{i,j+1}}),
\tag{25}
$$
$$
c_{i,j}-c_{i,j+1} = \frac{A}{a^2}(e^{-aq_{i+1,j}}-e^{-aq_{i,j}}-e^{-aq_{i+1,j+1}}+e^{-aq_{i,j+1}}),
\tag{26}
$$
$$
C(p_i) = \frac{A}{a}p_i - \frac{A}{a^2}(1-e^{-ap_i})=\frac{A}{a^2}(ap_i-1+e^{-ap_i}),
\tag{27}
$$
So the first term in equation 18 can be rewritten as:
$$
\frac{A}{a^2p_i}(ap_i-1+e^{-ap_i})T_i=\frac{A}{a}(1-\frac{1-e^{-ap_i}}{ap_i})T_i.
\tag{28}
$$
For the sum in equation 18,
$$
S_i=\sum_{j=i-1}^{0} T_j[c_{i,j}-c_{i,j+1}],
\tag{29}
$$
Here, we rewrite $c_{i,j}-c_{i,j+1}$ as:
$$
\begin{align}
c_{i,j}-c_{i,j+1} &= \frac{A}{a^2}(e^{-aq_{i+1,j}}-e^{-aq_{i,j}}-e^{-aq_{i+1,j+1}}+e^{-aq_{i,j+1}}) \\
&= \frac{A}{a^2}[(e^{-aq_{i,j+1}}-e^{-aq_{i,j}}) - (e^{-aq_{i+1,j+1}}-e^{-aq_{i+1,j}})] \\
&= \frac{A}{a^2}(e^{-aq_{i,j+1}}-e^{-aq_{i,j}})(1-e^{-a(p_i+p_{i+1})/2})
\end{align}
\tag{30}
$$
Substituting equation 30 into equation 29, we obtain:
$$
\begin{align}
S_i &= \sum_{j=i-1}^0 T_j\frac{A}{a^2}(1-e^{-a(p_i+p_{i+1})/2})(e^{-aq_{i,j+1}}-e^{-aq_{i,j}}) \\
&= \frac{A}{a^2} (1-e^{-a(p_i+p_{i+1})/2}) \sum_{j=i-1}^0 (e^{-aq_{i,j+1}}-e^{-aq_{i,j}})T_j
\end{align}
\tag{31}
$$
$$
\begin{align}
S_{i+1} &= \frac{A}{a^2} (1-e^{-a(p_i+p_{i+1})/2})\sum_{j=i}^0 (e^{-aq_{i+1,j+1}} - e^{-aq_{i+1,j}})T_j \\
&= \frac{A}{a^2} (1-e^{-a(p_i+p_{i+1})/2})\left((e^{-aq_{i+1,i+1}} - e^{-aq_{i+1,i}})T_i + e^{-a(p_{i+1} + p_i)/2}\sum_{j=i-1}^0 (e^{-aq_{i,j+1}}-e^{-aq_{i,j}})T_j\right) \\
&= \frac{A}{a^2} (1-e^{-a(p_i+p_{i+1})/2})\left((1 - e^{-a(p_{i+1}+p_i)/2})T_i + e^{-a(p_{i+1} + p_i)/2}\sum_{j=i-1}^0 (e^{-aq_{i,j+1}}-e^{-aq_{i,j}})T_j\right)
\end{align}
\tag{32}
$$
Here we define $X_i=\sum_{j=i-1}^0 (e^{-aq_{i,j+1}}-e^{-aq_{i,j}})T_j$. Then the $D_i$ can be calculated using the following iterations:
$$
\begin{align}
&X_i = (1 - e^{-a(p_i+p_{i-1})/2})T_i + e^{-a(p_i + p_{i-1})/2}X_{i-1}, \\
&g_i = \frac{1-e^{-ap_i}}{ap_i}, \\
&D_i = \frac{A}{a}\left((1-g_i)T_i + \frac{1-e^{-a(p_i+p_{i+1})/2}}{ap_i}X_i\right)\approx\frac{A}{a}[(1-g_i)T_i+g_iX_i]
\end{align}
\tag{33}
$$
