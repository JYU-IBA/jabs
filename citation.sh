#!/bin/bash
#This script will generate CITATION.cff
today=$(date "+%Y-%m-%d")
read version < version.txt
echo "cff-version: 1.2.0
message: "If you use this software, please cite it as below."
authors:
 - family-names: "Julin"
   given-names: "Jaakko"
   orcid: "https://orcid.org/0000-0003-4376-891X"
title: "JaBS"
repository-code: "https://github.com/JYU-IBA/jabs"
license: GPL-2.0-or-later
version: $version
doi: 10.5281/zenodo.6362701
date-released: $today" > CITATION.CFF
