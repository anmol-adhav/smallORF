# smallORF
Small ORF Finder - Genome Analyzer
A lightweight, GUI-based Python desktop application designed to identify and extract Small Open Reading Frames (sORFs) from genomic sequences. Built with Tkinter and Biopython, this tool scans all six translation frames of a given sequence to find ORFs that match a user-defined amino acid (AA) length and tolerance.

🚀 Features
User-Friendly GUI: Simple graphical interface built with Tkinter—no command-line experience required to run analyses.

Flexible Input: Supports standard FASTA and FNA formats, including natively reading .gz compressed files.

Comprehensive Scanning: Automatically translates and scans all 6 reading frames (3 forward, 3 reverse complement).

Customizable Filtering: Define a target amino acid length and a tolerance percentage (±%) to isolate specific sORF sizes.

Interactive Results: View parsed ORFs in an integrated data table detailing sequence ID, frame, strand, coordinates, AA length, and sequence previews.

Easy Export: Save the filtered ORFs directly to a new FASTA file for downstream analysis.

📋 Prerequisites
To run this application, you need Python 3 installed on your system along with the Biopython library. The other modules (tkinter, gzip, re) are included in the Python Standard Library.

Installation
Ensure you have Python installed (Python 3.6+ recommended).

Install the required Biopython library using pip:

Bash
pip install biopython
(Note: Linux users may also need to install the python3-tk package via their system's package manager if Tkinter is not bundled with their Python installation).

💻 Usage
Launch the Application:
Run the script from your terminal or command prompt:

Bash
python small_orf_finder.py
Load a Sequence File: Click "Load FASTA/FNA File" and select your target genome file. The app will confirm how many sequences were successfully loaded in the status bar at the bottom.

Set Parameters:

Target AA length: Enter the desired size of the ORF in amino acids (default is 50).

Tolerance % (±): Enter the acceptable variance. For example, a target of 50 with a 10% tolerance will search for ORFs between 45 and 55 amino acids long.

Run the Scanner:
Click "Find ORFs". The app will populate the table with any sequences that start with an ATG codon, end with a standard stop codon (TAA, TAG, TGA), and fall within your specified length parameters.

Export Results:
Click "Save ORFs (FASTA)" to export the nucleotide sequences of your discovered ORFs.

📄 Output Format
When you save your results, the output FASTA file will contain detailed headers generated from the table data, making it easy to trace back to the source sequence.

Example output structure:

Code snippet
>SequenceID_ORF1_frame+1_len52aa_pos100-256
ATGCGT...[nucleotide sequence]...TAA
🔬 Scientific Context
Start Codon: Searches specifically for the canonical start codon ATG.

Stop Codons: Stops translation at TAA, TAG, or TGA.

Translation Table: Uses the standard genetic code.
