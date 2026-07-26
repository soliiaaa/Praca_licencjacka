About this project

This tool was developed as part of a bachelor's thesis:

"Development and implementation of a tool for analyzing protein–protein interaction interfaces based on the three-dimensional structure of a complex in PyMOL" ("Opracowanie i implementacja narzędzia do analizy interfejsu oddziaływań białko-białko w oparciu o trójwymiarową strukturę kompleksu w programie PyMol")

Author: Solomiia Navachkevych

Institution: University of Warsaw, Faculty of Biology

Program: Bachelor's degree in Biotechnology

Supervisor: dr Norbert Odolczyk, Laboratory of Systems Biology (Zakład Biologii Systemów)


PDB Interface Analyzer

A PyMOL plugin for analyzing protein–protein interfaces from PDB/CIF structures.

The program automatically detects and visualizes interactions between two specified chains in a protein structure, calculating hydrogen bonds, salt bridges, hydrophobic contacts, and solvent-accessible surface area (SASA/BSA).

Installation:

Open PyMOL.

Go to Plugin → Plugin Manager → Install New Plugin.

Select the licencjat.py file.

Restart PyMOL after installation.

The plugin will appear in the menu Plugin → Analizator Interfejsu PDB.


Usage:

Launch the plugin from the Plugin menu.

Enter the identifiers of the two chains you want to compare (e.g. A and B).

Load a structure:

by typing a 4-character PDB code and clicking Fetch, or
by clicking Load file and selecting a local .pdb/.cif file.

Set the distance cutoff defining the interface (default 5.0 Å).

Click RUN ANALYSIS.

Results will appear:

directly in the PyMOL scene (coloring, sticks, interface surface, distance lines),
in the main panel (SASA, BSA),
in the interaction results window (with CSV export option),
in the interactive 2D contact map window.
