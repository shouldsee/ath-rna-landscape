### Background 

Despite the many publicly available RNA-Seq datasets, there is not yet a simple interface that displays the correlation 
between two genes from Arabidopsis thaliana, largely due to the fact that correlation-based 
hypothesis has not been adopted by plant biologists (for various reasons), and that generating 
uniformly mapped datasets on a large scale is computational-demanding.

A simple, robust interface for augmenting biology with correlation.

Here I present an online website that facilitates researching transcriptional network by
simple visualisation of the pearson-correlation between two gene expression vectors across samples.


### Pure correlation contradicts biological intuition

Though CLT applies to the esimation of most univariate statistics, the global correlation
coefficient in the limit of infinite samples is often not what's required by 
biologists. For example (blah)

Hence it is essential to select an appropriate pool of samples to enable in-context 
calculation of correlation values.


### No control means no quality

In order to control the quality of any calculated correlation, MetaAth allows users to view
the scatterplot underneath any correlation.

Hence, MetaAth allows user to browse correlation in an on-demand fashion and to 
interpret the gene-gene scatterplots manually. I found the display of meta data is useful 
in deconvolving the origin of correlation and implemented a responsive table.


### Refs

- Sun, Y., Dinneny, J.R. Q&A: How do gene regulatory networks control environmental responses in plants?. BMC Biol 16, 38 (2018). https://doi.org/10.1186/s12915-018-0506-7



TY  - JOUR
T1  - Abiotic Stress Signaling and Responses in Plants
AU  - Zhu, Jian-Kang
Y1  - 2016/10/06
PY  - 2016
N1  - doi: 10.1016/j.cell.2016.08.029
DO  - 10.1016/j.cell.2016.08.029
T2  - Cell
JF  - Cell
SP  - 313
EP  - 324
VL  - 167
IS  - 2
PB  - Elsevier
N2  - As sessile organisms, plants must cope with abiotic stress such as soil salinity, drought, and extreme temperatures. Core stress-signaling pathways involve protein kinases related to the yeast SNF1 and mammalian AMPK, suggesting that stress signaling in plants evolved from energy sensing. Stress signaling regulates proteins critical for ion and water transport and for metabolic and gene-expression reprogramming to bring about ionic and water homeostasis and cellular stability under stress conditions. Understanding stress signaling and responses will increase our ability to improve stress resistance in crops to achieve agricultural sustainability and food security for a growing world population.
AB  - As sessile organisms, plants must cope with abiotic stress such as soil salinity, drought, and extreme temperatures. Core stress-signaling pathways involve protein kinases related to the yeast SNF1 and mammalian AMPK, suggesting that stress signaling in plants evolved from energy sensing. Stress signaling regulates proteins critical for ion and water transport and for metabolic and gene-expression reprogramming to bring about ionic and water homeostasis and cellular stability under stress conditions. Understanding stress signaling and responses will increase our ability to improve stress resistance in crops to achieve agricultural sustainability and food security for a growing world population.
SN  - 0092-8674
M3  - doi: 10.1016/j.cell.2016.08.029
UR  - https://doi.org/10.1016/j.cell.2016.08.029
Y2  - 2020/04/13
ER  - 
