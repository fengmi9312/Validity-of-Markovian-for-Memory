# Validity of Markovian modeling for transient memory dependent epidemic dynamics

The extension modules, experimental codes, figure codes and application package for the paper titled "Validity of Markovian modeling for transient memory-dependent epidemic dynamics"<sup>[1]</sup> 

[Mi Feng](https://scholar.google.com/citations?user=09WNOQwAAAAJ&hl=en&oi=ao)<sup>1,2</sup>, [Liang Tian](https://physics.hkbu.edu.hk/people/tian-liang)<sup>1,3,\*</sup>, [Ying-Cheng Lai](http://chaos1.la.asu.edu/~ylai1/)<sup>4,5</sup>, [Changsong Zhou](https://physics.hkbu.edu.hk/people/zhou-chang-song)<sup>1,2,3,\*</sup> 
    
<sup>*</sup>Correspondence: [liangtian@hkbu.edu.hk](mailto:liangtian@hkbu.edu.hk), [cszhou@hkbu.edu.hk](mailto:cszhou@hkbu.edu.hk)

<sup>1</sup>[Department of Physics](https://physics.hkbu.edu.hk/), [Hong Kong Baptist University](https://www.hkbu.edu.hk/), Kowloon Tong, Hong Kong SAR 999077, China  
<sup>2</sup>[Centre for Nonlinear Studies and Beijing-Hong Kong-Singapore Joint Centre for Nonlinear and Complex Systems (Hong Kong)](https://scholars.hkbu.edu.hk/en/organisations/the-beijing-hong-kong-singapore-joint-centre-for-nonlinear-and-co), [Hong Kong Baptist University](https://www.hkbu.edu.hk/), Kowloon Tong, Hong Kong SAR 999077, China  
<sup>3</sup>[Institute of Computational and Theoretical Studies](https://www.icts.hkbu.edu.hk/), [Hong Kong Baptist University](https://www.hkbu.edu.hk/), Kowloon, Hong Kong SAR 999077, China  
<sup>4</sup>[School of Electrical, Computer and Energy Engineering](https://ecee.engineering.asu.edu/), [Arizona State University](https://www.asu.edu/), Tempe, AZ 85287, USA  
<sup>5</sup>[Department of Physics](https://physics.asu.edu/), [Arizona State University](https://www.asu.edu/), Tempe, Arizona 85287, USA

## Table of Contents
- [Paper](#Paper)
- [Dependencies](#Dependencies)
- [Experimental Code](#Experimental-Code)
- [Data Code](#Data-Code)
- [Figure Code](#Figure-Code)
- [Application](#Application)

## Paper
The main text as well as the Supplementary Information can be downloaded through <a href = "https://github.com/fengmi9312/Validity-of-Markovian-for-Memory/raw/main/Maintext-SI.pdf">Link</a>. And the main text in arXiv is available at <a href="https://arxiv.org/pdf/2306.16864">arXiv</a>. The latex version can be downloaded via <a href="https://arxiv.org/e-print/2306.16864">source</a>.

## Dependencies
Within the "<a href="https://github.com/fengmi9312/Validity-of-Markovian-for-Memory/tree/main/Dependencies">Dependencies</a>" folder, you'll find the data dependencies, code dependencies and C++ code that generates the module for conducting Monte Carlo simulations of the stratified SIR spreading dynamics model. This module plays a crucial role in the experimental and data generation processes that follow.

## Experimental Code
The "<a href="https://github.com/fengmi9312/Validity-of-Markovian-for-Memory/tree/main/ExperimentalCode">ExperimentalCode</a>" folder comprises the Python code responsible for generating all the experimental data.

Users can run the code on both Linux and Windows systems (by default) by making adjustments in the code. Please refer to the comments in the code for detailed instructions on how to make the necessary adjustments. Furthermore, The script file "validity-tasks.sh" is used to allocate computing resources with Slurm on Linux.

The folder "<a href="https://github.com/fengmi9312/Validity-of-Markovian-for-Memory/tree/main/ExperimentalData">ExperimentalData</a>" is where the experimental data can be generated.

## Data Code
The folder titled "<a href="https://github.com/fengmi9312/Validity-of-Markovian-for-Memory/tree/main/DataCode">DataCode</a>" contains the Python code that is utilized to generate the figure data based on the experimental data (in folder "<a href="https://github.com/fengmi9312/Validity-of-Markovian-for-Memory/tree/main/ExperimentalData">ExperimentalData</a>") presented in this paper.

The folder "<a href="https://github.com/fengmi9312/Validity-of-Markovian-for-Memory/tree/main/FigureData">FigureData</a>" is where the figure data can be generated.

## Figure Code
Within the "<a href="https://github.com/fengmi9312/Validity-of-Markovian-for-Memory/tree/main/FigureCode">FigureCode</a>" folder, you will find the Python code responsible for extracting the valid figure data from the experimental data and plotting all the figures showcased in both the main text and supplementary information.

The folder "<a href="https://github.com/fengmi9312/Validity-of-Markovian-for-Memory/tree/main/Figure">Figure</a>" is where the figures can be generated.

## Applications
We have developed two types of applications, namely a web-based and a Python version, to aid researchers in determining the generation time distribution, estimating the parameters of the Markovian dynamics in transient-state equivalence, and rectifying errors in $R_0$ and steady-state forecasting.

### Web-based Application
Readers can access the <a href = "https://cns.hkbu.edu.hk/toolbox/Validity-of-Markovian-for-Memory/main.html">web-based application</a> online. The source code for the web application can be found in the folder "./applications/Web-based_application".

### Python Application 
Readers can run it through source code as follows:
- Download: The corresponding package could be downloaded via <a href="https://github.com/fengmi9312/Validity-of-Markovian-for-Memory/raw/main/Applications/Python_Application.zip">Python Application</a> (5 KB, version 1.0).
- Installation: Extract the content of the enclosed Python_Application.zip file to a local directory.
- Run: Run the corresponding python script files named "Distribution.py" or "Rectification.py".
- Dependencies: numpy, scipy, and matplotlib.

### Details of Python Application 
This application demonstrates the diverse generation time distribution patterns resulting from various infection and removal time distributions, each with different parameters. Furthermore, it provides a user-friendly methodology to assess the potential underestimation or overestimation of memory-dependent epidemic dynamics through the implementation of Markovian modeling <sup>1</sup>. Here, we offer a detailed description of each input parameter used in the app and further explanation of the results.

#### Distribution
In epidemic dynamics, the infectivity of an infected individual can be characterized by the infection time distribution, denoted as $\psi_{\mathrm{inf}}(\tau)$. Here, $\tau$ represents the elapsed time between the individual's infection and the current time. The probability of an infection event occurring within the time interval $[\tau, \tau+d\tau)$ is given by $\psi_{\mathrm{inf}}(\tau)d\tau$. Similarly, the removal process is described by the removal time distribution, $\psi_{\mathrm{rem}}(\tau)$. It represents the probability of a removal event happening within the time interval $[\tau, \tau + d\tau)$. These time distributions for the infection and removal processes incorporate memory effects and are general in nature. The exponential distributions associated with memoryless processes are a special case of the broader class of memory-dependent processes.

In "Distribution.py", we visualize the infection and removal time distributions, $\psi_{\mathrm{inf}}(\tau)$ and $\psi_{\mathrm{rem}}(\tau)$,  along with their respective hazard functions, $\omega_{\mathrm{inf}}(\tau)$ and $\omega_{\mathrm{rem}}(\tau)$, and survival functions, $\Psi_{\mathrm{inf}}(\tau)$ and $\Psi_{\mathrm{rem}}(\tau)$. This provides a comprehensive overview of the statistical characteristics of the time durations involved in the epidemic dynamics. In details, the survival function can be calculated by:
$$\Psi(\tau) = \int_{\tau}^{+\infty}{\psi(\tau')d\tau'},$$
where $\Psi(\tau)$ can be either $\Psi_{\mathrm{inf}}(\tau)$ or 
$\Psi_{\mathrm{rem}}(\tau)$, while $\psi(\tau)$ can be either $\psi_{\mathrm{inf}}(\tau)$ or $\psi_{\mathrm{rem}}(\tau)$.
The hazard function can be calculated by:
$$\omega(\tau) = \frac{\psi(\tau)}{\Psi(\tau)},$$
where $\omega(\tau)$ can be either $\omega_{\mathrm{inf}}(\tau)$ or $\omega_{\mathrm{rem}}(\tau)$.
Meanwhile, the generation time distribution $\psi_{\mathrm{gen}}(\tau)$ could be calculated as follows:
$$\psi_{\mathrm{gen}}(\tau) = \frac{\omega_{\mathrm{inf}}(\tau)\Psi_{\mathrm{rem}}(\tau)}{\int_{0}^{+\infty}{\omega_{\mathrm{inf}}(\tau')\Psi_{\mathrm{rem}}(\tau')d\tau'}}.$$
$T_{\mathrm{inf}}$, $T_{\mathrm{rem}}$ and $T_{\mathrm{gen}}$ represent the the mean values of $\psi_{\mathrm{inf}}(\tau)$, $\psi_{\mathrm{rem}}(\tau)$ and $\psi_{\mathrm{gen}}(\tau)$, respectively.

The basic reproduction number $R_0$ can be calculated by:
$$R_0 = \Lambda_{\mathrm{max}}\int_{0}^{+\infty}{\omega_{\mathrm{inf}}(\tau)\Psi_{\mathrm{rem}}(\tau)d\tau}.$$
And the growth rate $g$ satisfies:
$$1 = R_0\int_{0}^{+\infty}{e^{-g\tau}\psi_{\mathrm{gen}}(\tau)d\tau}.$$
When $T_{\mathrm{gen}} = T_{\mathrm{rem}}$, the transient-state equivalence can be achieved. In this case, the infection rate ($\gamma$) and the removal rate ($\mu$), 
which are represented by the bar chart, can be calculated using the following expressions:
$$\gamma = \frac{gR_0}{\Lambda_{\mathrm{max}}(R_0 - 1)},$$
$$\mu = \frac{g}{R_0 - 1},$$
where $\Lambda_{\mathrm{max}}$ is the maximum eigenvalue of $kA\circ p$, $k$ is a parameter to adjust the overall contacts, $A$ denotes the contact matrix, and $p$ represent the population distribution.

The user has the flexibility to adjust the following input parameters in the application:
- Time Point Number: this parameter determines the number of time points in each figure. It allows the user to control the level of detail in the visual representation;
- Time Step: by adjusting the time step, the user can determine the calculation accuracy and define the total time duration covered by the data;
- Distribution Form: the application offers three distribution forms for modeling infection and removal times: Weibull, gamma, and log-normal;
- Parameters ($\alpha_{\mathrm{inf}}$, $\beta_{\mathrm{inf}}$) and  ($\alpha_{\mathrm{rem}}$, $\beta_{\mathrm{rem}}$): the subscripts "inf" and "rem" correspond to the parameters of the infection and removal time distributions, respectively; $\alpha$ and $\beta$ represent different parameters in the time distributions;
- Maximum Eigenvalue $\Lambda_{\mathrm{max}}$.
For Weibull distribution, it follows:
$$\psi(\tau) = \frac{\alpha}{\beta}(\frac{\tau}{\beta})^{\alpha-1}e^{(\frac{\tau}{\beta})^{\alpha}}.$$
For gamma distribution, it follows:
$$\psi(\tau) = \frac{1}{\Gamma(\alpha)\beta^\alpha}\tau^{\alpha-1}e^{-\frac{\tau}{\beta}}.$$
where $\Gamma(\cdot)$ denotes gamma function.
For log-normal distribution, it follow:
$$\psi(\tau) = \frac{1}{\tau\beta\sqrt{2\pi}}\exp(-\frac{(\ln\tau - \alpha)^2}{2\beta^2}).$$

#### Rectification
In "Rectification.py", we utilize the average generation and removal times from memory-dependent epidemic dynamics to adjust the estimations of various metrics derived from Markovian modeling. These metrics include the basic reproduction number, steady-state infected fraction of each age group, and total steady-state infected fraction. 

To compute the steady state using $R_0$, certain inputs are necessary, such as the contact matrix, population distribution, and initial infected fraction.
User can customize the contact matrix, age distribution, and initial infected fraction. It's important to ensure that the contact matrix must be symmetrical. 
Users can also use the contact matrix and age distribution provided in the file <sup>2, 3</sup>. Users can also make adjustments to the values of the average generation and removal times, $T_{\mathrm{gen}}$ and $T_{\mathrm{rem}}$, as well as the value of $R_0$.

In the figure, the first row displays the contact matrix, age distribution, and initial infected fraction respectively.
The second row illustrates the comparison between the Markovian estimation and the rectification results for $R_0$, steady-state infected fraction of each age group, and total steady-state infected fraction.

Based on Eq. (15), the rectified $R_0$ can be calculated by:
$$R_0 = (\hat{R}\_0)^{\eta^a},$$
where $\hat{R}\_0$ denotes the Markovian-estimated basic reproduction number, $\eta = T_{\mathrm{gen}}/T_{\mathrm{rem}}$, $a$ is a constant equal to 1.49.
Meanwhile, the steady-state infected fraction of each age group $\tilde{c}\_l$ satisfies:
$$1 - \tilde{c}\_l = \acute{s}\_l e^{-\frac{R_0}{\Lambda_{\mathrm{max}}}\sum_{m=0}^{n}{kA_{lm}p_m\[\tilde{r}\_m-\acute{r}\_m\]}}.$$
Because $\Lambda_{\mathrm{max}} / k$ represents the largest eigenvalue of $A\circ p$, we do not consider the impact of $k$.
The total cumulative infected fraction at steady state, $\tilde{c}$, can be calculated by:
$$\tilde{c} = \sum_{l=0}^{n}{\tilde{c}_lp_l}.$$

#### Reference
1. [Feng, M., Tian, L., Lai, Y.-C. & Zhou, C. Validity of Markovian modeling for transient memory-dependent epidemic dynamics. arXiv:2306.16864 (2023)](https://arxiv.org/abs/2306.16864).
2. [Prem, K., Cook, A. R. & Jit, M. Projecting social contact matrices in 152 countries using contact surveys and demographic data. PLoS Comput. Biol. 13, e1005697 (2017)](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697).
3. [United Nations, Department of Economic and Social Affairs, Population Division. World population prospects 2019 (2021). Accessed on 29 June 2021](https://population.un.org/wpp/).
