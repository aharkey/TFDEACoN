# TF DEACoN 
TF DEACoN stands for Transcription Factor Discovery by Enrichment Analysis of Co-expression Networks. It is a tool for generating hypotheses about which transcription factors (TFs) may be important for observed transcriptional responses. Given a set of genes which are linked transcriptionally (e.g. because they were observed to respond to the same stimulus or across a developmental time-frame, etc.), TF DEACoN will report back which TFs show target enrichment in this set.

It can be accessed [here](https://aharkey.shinyapps.io/tfdeacon/)

This repo contains:

* app.R - the R Shiny code behind the app
* .Rmd files for the help, citations, and contact form pages
* input - the files accessed by the code containing information on the transcription factors and their targets
