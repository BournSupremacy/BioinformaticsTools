## BioinformaticsTools

Repo with pipelines and scripts I used frequently for bioinformatics analysis of my MSc project.
Project entitled "Creating and analysing an African pan-genome".

### Analyses included:
- PythonPlots: directory containing JupyterNotebooks for all the plots presented in my dissertation.
  Analyses included are (by order of appearance in the dissertation):
    - CD-HIT cut-off and contig length distribution plots in `Plots.ipynb`
    - alignment of 1kGP data to GRCh38 in `1kGPvb38.ipynb`
    - PCA and upset plot for 1kGP data `1kGpca&upset.ipynb`
    - MultiQC analysis plots in `MultiQCplots.ipynb`
    - contamination scatter plots in `Plots.ipynb`
    - population non-ref seq bar graphs in `Plots.ipynb`
    - repeat element stacked bar graph in `Plots.ipynb`
    - pan-African PCA and upset plots in `PlotsAfricanPCA.ipynb` and `PlotsAfricanUpset.ipynb`
    - MAKER genes upset plots in `MAKERgeneUpset.ipynb`
    - all alignment plots in `PlotsAlignment.ipynb`
    - analysis of core and distributed genes at different cut offs in `PAVanalysis.ipynb`
    - core genes plots in `Plots.ipynb`
    - PAV profiles in `Plots.ipynb`

- QuastMultiQC: directory containing the MultiQC reports for the starting dataset and final dataset for the African pan-genome, and the Nextflow pipeline config files used to create those reports.

- alignment: directory containing the Python script used to analyse all alignment outputs from NUCmer and pblat.
