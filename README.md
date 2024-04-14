# Customization of bioinformatics workflows for (meta)genomics
<div aling="justify">
This project was developed by [Ricardo Oliveira(pg53501)](https://github.com/ricardofoliveira61) as part of the bioinformatics master's degree program, under the guidance of Professor Pedro Santos from the Biology Department of the University of Minho.
</div>

# Project description
<div aling="justify">
Identification and characterization of one or more organisms present in a sample through DNA-based analysis is a revolutionary approach made possible through advances in sequencing technologies and computing power. Like all sample analyses, standards are needed to benchmark performance and enable translation into applications such as clinical diagnostics, environment, and bio-surveillance. Such approach involves a series of steps, including sample collection, extraction, DNA preparation, sequencing, bioinformatics, and reporting. Bias at each step contributes to errors in the final analysis. With a wide variety of options for each step, two analyses can produce different results. Therefore, standards can help characterize the biases, enabling researchers to recognize any limitations and have confidence in the results. In our lab, we have multiple genomic, metabarcoding, and metagenomics data that has been processed through different bioinformatic workflows. However, a pertinent question that remains open is whether the elected options for such analyses were the most adequate. In this context, the current project aims at the customization of bioinformatic workflows that will include objective key checkpoints evaluating data processing and biological inference accuracy.
</div>

# Objective
<div aling="justify">
Our project aims to evaluate the performance of metabarcoding analysis tools using datasets generated with different sequencing platforms (Illumina, PacBio, and Oxford Nanopore) systematically. We intend to assess each step of the analysis workflow using various tools such as DADA2  or USEARCH/VSEARCH, and different versions of each tool to determine the most accurate approach for every step of the workflow.
For this purpose, we will resort to commercial mock community sequencing datasets (e.g. Zymobiomics and ATCC) and customized mock communities created in our lab, which facilitates bioinformatic performance assessment 

We intend to develop a set of customizable scripts that use the best tool for each phase of the analysis process by default. Nevertheless, the user will have the option to choose which tool to use and set the desired parameters for every analysis phase. Furthermore, our scripts will incorporate a set of metrics that allow measure the performance of the chosen workflow and evaluate if the results meet the established standards or if changes to the workflow are necessary to provide a more accurate result. Finally, we plan to apply the developed scripts to a real case study in our laboratory.
</div>
