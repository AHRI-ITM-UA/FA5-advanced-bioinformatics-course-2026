[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-green.svg)](https://creativecommons.org/licenses/by/4.0/) [![DOI](https://zenodo.org/badge/727651304.svg)](https://doi.org/10.5281/zenodo.18701255)

# Fighting Malaria Across Borders (FiMAB) - Bioinformatics Module

This repository houses the course notes for the Bioinformatics Module of the Fighting Malaria Across Borders (FiMAB) [International Training Programme](https://www.vliruos.be/get-funded/calls/international-training-programme-2024) that ran in April 2024. FiMAB ITP was a VLIR-UOS sponsored project created by the [Institute of Tropical Medicine Antwerp (ITM)](https://www.itg.be) in collaboration with the [University of Antwerp](https://www.uantwerpen.be/).

> The course notes for the bioinformatics module are available at the following URL: [https://pmoris.github.io/FiMAB-bioinformatics/](https://pmoris.github.io/FiMAB-bioinformatics/).

The primary goal of the overall program is to support the implementation of targeted sequencing assays (in particular, AmpliSeq) to strengthen malaria molecular surveillance and help guide national control programmes. In conjunction with laboratory training, this bioinformatics course module is intended to allow young academics around the globe to become familiar with molecular surveillance as a key activity to monitor transmission, sources of epidemics and the emergence and spread of drug resistance mutations in the _Plasmodium_ parasite.

Currently, these course notes for the bioinformatics module focus on an introduction to the Unix command line and data manipulation using R, followed by an introduction to targeted/whole genome sequencing analysis pipelines. Course materials aimed at more in-depth downstream population genetic analyses or marker discovery and interpretation were only made available on the [ITM Moodle platform](https://campus.itg.be/course/view.php?id=217), but might be added to these online notes too at a later time.

## Training tools

Accompanying GitHub Codespace (linux CLI/VSCode) environments containing all the training files and tools, can be accessed in two ways (see the course notes section 2 for more information):

1. Bioinformatics training environment (default): a VSCode environment containing a Unix terminal and several common (bioinformatics) bash tools. [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/pmoris/FiMAB-bioinformatics?template=false&devcontainer_path=.devcontainer%2Fdevcontainer.json)
2. RStudio Server (must be selected manually): an RStudio server containing tidyverse, based on the [rocker project](https://rocker-project.org/images/devcontainer/images.html). Once the environment has started (a VSCode session should appear), you can browse to the RStudio IDE URL listed under the ports tab (in the same pane as the terminal) or via the forwarded ports "radio" icon in the status bar at the bottom of the window (it will look similar to this url: `https://ubiquitous-space-eureka-5pxwr7qg96p2vr79-8787.app.github.dev/`/). [![Open in GitHub Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/pmoris/FiMAB-bioinformatics?template=false&devcontainer_path=.devcontainer%2Frstudio%2Fdevcontainer.json)

## Re-using these materials and how to cite us

Everything you find in this git repository and on the course website is licensed under CC-BY 4.0 ([Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/)), which means that we encourage you to share, re-use or modify/improve these course materials as you see fit, as long as you give us proper attribution. Ideally, you do this by citing the Zenodo DOI (_url to be included here after initial release_). Please make sure to also provide proper attribution to the primary sources that we cite ourselves, whenever this is applicable.

## Contributors and acknowledgements

The entire Fighting Malaria Across Borders (FiMAB) International Training Programme was made possible by a team of ITM instructors ([Dr. Johanna Helena Kattenberg](https://orcid.org/0000-0002-2971-5136) - malaria molecular surveillance and population genetics; [Mathijs Mutsaers](https://orcid.org/0000-0002-7633-6438) - main organiser and lab training; Ramses Deyaert - lab training) and guest lecturers ([Dr. Anita Ghansah](https://orcid.org/0000-0003-4639-1249) - applications of molecular surveillance; [Dr. Conor Meehan](https://orcid.org/0000-0003-0724-8343) - [introduction to phylogenetics](https://conmeehan.github.io/PathogenDataCourse/IntroToPhylogenetics); [Dr. Jody Phelan](https://orcid.org/0000-0001-8323-7019) - [malaria genomic analyses](https://jodyphelan.github.io/malaria-training/)), as well as the VLIR-UOS [Interational Training Progamme 2024](https://www.vliruos.be/get-funded/calls/international-training-programme-2024) funding and the [University of Antwerp](https://www.uantwerpen.be/) as an organising partner.
