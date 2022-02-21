# ExAtlas
Analysis of gene expression and gene set data
Developed at the National Institute on Aging, NIA/NIH, Baltimore USA.
Programmer Alexei Sharov sharov@comcast.net

ExAtlas is a tool for meta-analysis of gene expression data. In contrast to other software, it compares multi-component data sets and generates results for all combinations (e.g., all gene expression profiles vs. all GO annotations). Main functions are: (1) standard meta-analysis (fixed and random effects, z-score, Fisher's method; analysis of (2) global correlations between gene expression data sets; (3) gene set enrichment; (4) gene set overlap; (5) gene association (Expected Proportion of False Positives, EPFP); (6) statistical analysis of gene expression (ANOVA, PCA, FDR). Results are presented graphically as heatmaps, bar-charts, or 3-D images (PCA). Gene expression data is extracted automatically from GEO/NCBI database. Several most popular public data sets (e.g., GNF, Gene Ontology, KEGG, GAD phenotypes) are pre-loaded and can be used for functional annotations.

This software is provided "AS IS".  Programmer makes no warranties, express or implied, including no representation or warranty with respect to
the performance of the software and derivatives or their safety, effectiveness, or commercial viability. Programmer does not warrant the
merchantability or fitness of the software and derivatives for any particular purpose, or that they may be exploited without infringing
the copyrights, patent rights or property rights of others. Programmer shall not be liable for any claim, demand or action for any loss, harm,
illness or other damage or injury arising from access to or use of the software or associated information, including without limitation any
direct, indirect, incidental, exemplary, special or consequential damages. This software program may not be sold, leased, transferred, exported
or otherwise disclaimed to anyone, in whole or in part, without the prior written consent of programmer.

DIRECTORY STRUSTURE

ExAtlas = exatlas.ini file

ExAtlas/exatlas = html files

ExAtlas/exatlas/bin = exatlas.cgi (CGI program)

ExAtlas/exatlas/images = images

ExAtlas/exatlas/output = writable directory for output files

ExAtlas/exatlas/download

ExAtlas/exatlasInfo/bin = other programs (C, perl), C-programs should be compiled here

ExAtlas/exatlas/info = writable directory for personal config files, login.txt

ExAtlas/exatlas/data = writable directory for data files

Compile C-programs in ExAtlas/exatlasInfo/bin:

gcc anova_oneway.c -lm -o anova_oneway

gcc togif.c -lm -o togif

gcc correlation_exatlas.c -lm -o correlation_exatlas

gcc norm_new.c -lm -o norm_new

gcc page_exatlas.c -lm -o page_exatlas

gcc pairwise.c -lm -o pairwise

gcc pca.c -lm -o pca

Configure two accounts to your name: administrator, public. Edit file login.txt and put your name and email address
Log into ExAtlas and change your passwords for both account. Create your own working account in ExAtlas
