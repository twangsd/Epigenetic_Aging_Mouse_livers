# Epigenetic_Aging_Mouse_livers
Ipython notebooks and associated scripts to generate figures for the manuscript entitled:<br> Epigenetic aging signatures in mice are slowed by dwarfism, calorie restriction and rapamycin treatment<br>
<br>
In order to run these notebooks, you must have:<br>
<br>
<b>Dependencies</b><br>
Anaconda ipython (2.7)<br>
Seaborn <br>
Anaconda R<br>
SVA (R package) <br>
Impute (R package) <br>
Minfi (R package) <br>
WGCNA (R package) <br>
Bedtools on path <br>
<br>
<b> Note: </b><br>
Many of these scripts were run in a parallel computing framework using SGE.<br> 
Thus, these notebooks represent an outline of the analyses that were conducted. They will not run out of the box<br>
After the manuscript has been accepted, data needed to regenerate the figures will be added on the Ideker lab website. These data will represent the end output of many of the scripts executed on a cluster. The zip file will contain a 'data' directory, which will need to be downloaded within the same directory as the notebooks so that the scripts displayed in the notebooks can find the data used to create each plot.<br>
<br>
<b>Figures:</b><br>
Figure 1 & Figure S1: EvolutionaryTrendsMouseHuman.ipynb<br>
Figure 2 & Figure S2: MouseEpigeneticAgingModel.ipynb & MouseEpigeneticAgingModel_RandomizationControls.ipynb<br>
Figure 3 & Figure S3: MouseEpigeneticAgingModel.ipynb & UnsupervisedAnalysis_ElasticNet_CpGs.ipynb<br>
Additional File 3 additions, for nearest genes and part of Figure S2: ModelCpGSites_gene_proximity.ipynb
