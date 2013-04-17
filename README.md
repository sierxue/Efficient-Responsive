FNV-MMFE
========

    Repository:     FNV-MMFE
    Author:         Tong WANG
    Email:          tong.wang@nus.edu.sg
    Version:        v3.0 (originally written in 2011, reorganised and published on 2013-04-17)
    Copyright:      absolutely free (In case you find the code helpful in any sense and do not know how to express your appreciation, cite our paper below :)
    Description:
                    This repository hosts codes for the numerical experiments conducted in the paper "A Multiordering Newsvendor Model with Dynamic Forecast Evolution."

                    <FNV-MMFE/src> folder hosts the C++ code for the four models developed in the paper and its e-companion: additive MMFE, multiplicative MMFE, additive MMFE with fixed ordering costs, and additive MMFE with order cancelations. The codes essentially solve the Dynamic Programs developed in the paper to obtain the optimal inventory control policies and then run Monte Carlo simulations to compare the optimal multi-ordering and single-ordering strategies.

                    <FNV-MMFE/output> folder hosts the raw output data files generated by the C++ codes, the R scripts for cleaning, organizing, summarising, and plotting these data, and the final output (in the format of Tables and Figures) published in the paper.


Reference
=========

    Paper title:    "A Multiordering Newsvendor Model with Dynamic Forecast Evolution"
    Authors:        Tong Wang, Atalay Atasu and Mümin Kurtuluş
    Year:           2012
    Journal:        Manufacturing and Service Operations Management
    Volume:         14
    Number:         3
    Page:           472 -- 484
    Journal link:   http://msom.journal.informs.org/content/14/3/472.short
    BibTeX:
                    @article{WangAtasuKurtulus2012,
                    author = {Wang, Tong and Atasu, Atalay and Kurtuluş, Mümin}, 
                    title = {A Multiordering Newsvendor Model with Dynamic Forecast Evolution},
                    journal = {Manufacturing & Service Operations Management} 
                    volume = {14}, 
                    number = {3}, 
                    pages = {472-484}, 
                    year = {2012}, 
                    doi = {10.1287/msom.1120.0387}, 
                    URL = {http://msom.journal.informs.org/content/14/3/472.short}, 
                    }
    Further resources (such as the original paper in PDF format, its e-companion, and presentation slides) are available at http://bschool.nus.edu/staff/bizwt/research.html.


Instructions
============

    The C++ codes are re-written using the new C++11 standard library and the Boost C++ Library (http://www.boost.org/). In order to compile them, you will need:
        (1). An update-to-date C++ compiler (GNU GCC 4.7 or Intel C++ Compiler 13.0) supporting C++11 and OpenMP
        (2). Install the Boost Library (v1.53.0 or newer, check it out at http://www.boost.org/).
        (3). Compile by command: "g++ -std=c++11 -fopenmp -O3 -I /usr/local/boost_1_53_0 -o FNV-aMMFE.exe FNV-aMMFE.cpp", assuming the root directory of the Boost Library is "/usr/local/boost_1_53_0".

    The R scripts require a working R environment (http://www.r-project.org/). The scripts are only tested on R version 2.15, but should be running smoothly on newer versions. The scripts do not require extra R packages.

