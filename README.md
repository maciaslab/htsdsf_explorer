# HTSDSF Explorer

DSF data explorer focused in HTS.
Used to find binders to proteins which alter Tm.
Objectives: A program that is agile and allows to analyze a very high number of curves in a very short time.
The program finds wells for wwhich dTm with respect the reference is higher than a given threshold.
Then allows the user to validate the Tm with the minimum number of clicks.
The program has been tested in Mac OSX and Windows, but any modern OS with python3 and a web browser should run it.



# Installation

Just download the package and run:
python server.py
Program will start a webserver and open browserpointing at http://localhost:5555/index.html

## Windows Installation Instructions
1. Download python 3 64-bits from https://www.python.org/downloads/ and install it, checking the "Add to path"· option.

2. RUN:
python -m pip install scipy pandas xlsxwriter openpyxl matplotlib python-docx 

3. RUN
python server.py


## Requisites 
Python 3

### Python modules:
scipy pandas xlsxwriter openpyxl numpy matplotlib python-docx

### To install modules:
python3-m pip install scipy pandas xlsxwriter openpyxl numpy matplotlib python-docx

# Use
Program in python, interface through web server.

Two different modules:

* Effect of a molecule in a protein: dTm for a single concentration. HTS in which each well contains a different molecule. 
Reference values  can be in rows, columns or in individual wells.
Reports in Excel format.

* Kd calculation: Any number of rows contain a molecule at different concetrations.
Restriction: Each column has a fixed concentration 
The program is designed assuming that the dose-response dillution are in rows.
Row: molecule
Column: Concentration.
The program includes an easy to use Kd template editor.
Reports in Word format, including plots.

Tm calculations are done using first derivative, smoothed and supersampled (x5).
In case of multple Tms, the one closest to refTm is used.
Arbitrary number of reference wells, averaged and with stddev.

data:
Put data files in folders at {data_path}
Supported file types:
### txt: 
* Quantstudio 96 wells
* Quantstudio 384 wells
* Lightcycler

### xlsx:
* Biorad
### gdsf:
* Generic dsf file
A generic dsf (.gdsf) file format has been defined to allow the use of data acquired in different instruments. This format is a simple text file, with three columns separated by spaces containing in order the well, temperature and fluorescence.
``` 
A1  	20.22   2.31
A1  	20.44   2.29
A1  	20.68   2.27
A1  	20.99   2.26
A1  	21.28   2.23
...
```

# Plateinfo 
The plateinfo folder contains the information that correlates the wells in the plates with the ligands in them.
In HTS experiemnts, the sample plate templates are reused for different proteins.
You can (optionally) define these templates, and then assign a template to your plates. 
Format: 5 columns, tab separated. .txt extension. One header line.
Column 1: COMMENT
Column 2: Molecule id
Column 3: Smiles for Molecule
Column 4: Plate template name
Column 5: Well

You can have as many .txt files as needed, program will read all of them.
After matching a plate with a Plate tempalte name, the well will be used to find the molecule in the well (Molecule id and Smiles for Molecule) when generating a report.

# Settings
settings.ini contains the program settings.

Default file:

```
[Default]
data_path=data/
cache_path=cache/
persistent_path=persistent/
plateinfo_path=plateinfo/
```
# Citing this software

A paper describing this software is in process of being written. 

# License
This piece of sotware is under the GPL v3 license (https://www.gnu.org/licenses/gpl-3.0.en.html). Some parts of it are created by third parties and are under a MIT license. These parts have been unmodified and the files contain the headers corresponding to the original license and the original authors are credited below:

### Bootstrap
Bootstrap is released under the MIT license and is copyright 2018 Twitter (https://getbootstrap.com/)

### Chart.js
Chart.js is available under the MIT license. ( https://github.com/chartjs/Chart.js )

### DataTables
Bootstrap is released under the MIT license. (https://github.com/DataTables/DataTables)

### jQuery
Bootstrap is released under the MIT license. (https://jquery.org/)

