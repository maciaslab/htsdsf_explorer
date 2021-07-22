# HTSDSF Explorer

DSF data explorer focused in HTS.
Used to find binders to proteins which alter Tm.
Objectives: A program that is agile and allows to analyze a very high number of curves in a very short time.
The program finds wells for wwhich dTm with respect the reference is higher than a given threshold.
Then allows the user to validate the Tm with the minimum number of clicks.

settings.ini contains the program settings.

Default file:

[Default]

data_path=data/

cache_path=cache/

persistent_path=persistent/

plateinfo_path=plateinfo/


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
### gdsf
* Generic dsf file

plateinfo:
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




