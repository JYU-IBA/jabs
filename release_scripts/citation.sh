#!/bin/bash
#This script will generate CITATION.cff
citationfile="../CITATION.cff"
versionfile="../version.txt"
today=$(date "+%Y-%m-%d")
read version < "$versionfile"
echo "cff-version: 1.2.0
authors:
 - family-names: "Julin"
   given-names: "Jaakko"
   orcid: "https://orcid.org/0000-0003-4376-891X"
title: "JaBS"
repository-code: "https://github.com/JYU-IBA/jabs"
license: GPL-2.0-or-later
version: $version
doi: 10.5281/zenodo.6362701
date-released: $today" > "$citationfile"
