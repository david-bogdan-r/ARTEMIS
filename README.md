# ARTEMIS

Using [ARTEM](https://github.com/david-bogdan-r/ARTEM) to Infer Sequence alignment (ARTEMIS) is a tool for the RNA 3D structure superposition and determination of the structure-based sequence alignment

## Reference

TODO

## How ARTEMIS works

ARTEMIS reads a reference and a query structure from the specified coordinate files in PDB or in mmCIF format and finds two alignments of structure sequences.
The first alignment does not allow rearrangement and is fed to the standard output.
If verbose mode is enabled, the characteristics of the second alignment, which allows permutations, are additionally fed to the standard output.
Verbose mode is enabled automatically if the TM-score of the permutation alignment (normalised by the length of the query structure) is 10% higher than the TM-score of the alignment without permutations.
The user can choose to save the superpositions of the query structures in PDB or mmCIF format and the residue matching tables in the specified folder.
If verbose mode is disabled, only the superposition and the residue matching table for alignment without permutations are saved.
By default, ARTEMIS reads the entire first models of the input files.

### Algorithm

- perform superpositions of the query structure on the reference structure by each possible pair of residues between the sets *rseed* and *qseed* (the superposition operator is calculated using Kabsch's algorithm between 3-atom coordinate representations of residues);
- in each superposition find a set of residue pairs between structures with a distance not greater than *matchrange* angstroms between C3' atoms;
- among the *nlagest* largest sets of residue pairs, build alignments with and without permutations:
  - by the subset of mutually closest pairs superimpose the query structure on the reference structure (using 3-atom coordinate representations of residues);
  - closest pairs of residues between resuperimposed structures with the distance between C3' atoms not greater than 8 angstroms are considered to be **alignments with permutations**;
  - calculate the Score Matrix
    - multiply the matrix of distances between C3' atoms of superimposed structures by -1;
    - shift the matrix to the left by the minimum value in it (values in the matrix are not negative);
    - shift the matrix to the left by the minimum value in the cells of residue pairs from the permutation alignment (all matrix cells corresponding to pairs of residues of permutation alignment have non-negative values);
    - shift the matrix to the right by *shift* to increase the coverage in the alignment.
  - finding an **alignment** by Needleman-Wunsch Algorithm (score maximisation)
- select the better alignments by TM-score sum for the two structures.

## Installation

Clone the GitHub repository by typing

    git clone https://github.com/david-bogdan-r/ARTEMIS.git

## Dependencies

ARTEMIS requires three Python3 libraries to be installed:

- numpy
- pandas
- scipy

To install, type:

    pip install -r requirements.txt

ARTEMIS was tested with two different Python3 environments:

### Ubuntu 20.04

- python==3.8
- numpy==1.22.3
- pandas==1.4.1
- scipy==1.8.0

### MacOS Sonoma 14.0

- python==3.11.0
- numpy==1.23.4
- pandas==1.5.1
- scipy==1.9.3

## Usage

    python3 artemis.py r=FILENAME q=FILENAME [OPTIONS]

## Usage examples

    1) python3 artemis.py r=examples/6ugg/6ugg.cif  q=examples/1ivs.pdb rres=/B qres=/C saveto=result saveformat=pdb -v

    This command will superimpose chain C from 1ivs to chain B from 6ugg 
    and will save the superimposed structure into "result" sub-folder
    in PDB format.

## Options

    r=FILENAME [REQUIRED OPTION]
        Path to a reference structure in PDB/mmCIF format. For faster 
        performance, it's advised to specify the largest of the two 
        structures as the reference structure.

    q=FILENAME [REQUIRED OPTION]
        Path to a query structure, the one that ARTEMIS superimposes to 
        the reference, in PDB/mmCIF format.

    matchrange=FLOAT [DEFAULT: matchrange=3.0]
        The threshold used for searching the mutually closest residues. 
        Only those pairs of residues that have their C3' atoms at a distance 
        under the specified matchrange value can be added to the subset 
        of the mutually closest residues. The higher matchrange value 
        will produce more "noisy" matchings but won't miss anything. The lower 
        matchrange value will produce more "clean" matchings but 
        can miss something.

    rformat=KEYWORD, qformat=KEYWORD [DEFAULT: rformat=PDB,qformat=PDB] 
        The specification of the input coordinate file formats 
        (case-insensitive). By default, ARTEMIS tries to infer the format 
        from the extensions of the input filenames. ".pdb", ".cif", 
        and ".mmcif" formats can be recognized (case-insensitive). In the case 
        of any other extension ARTEMIS will treat the file as the PDB-format 
        file by default. If the "rformat" ("qformat") parameter is specified 
        and it's either "PDB", "CIF", or "MMCIF" (case-insensitive), 
        ARTEMIS will treat the reference (query) coordinate file
        as the specified format.

    rres=STRING, qres=STRING [DEFAULT: rres="#1",qres="#1"]
        The specification of the input reference (rres) and query (qres) 
        structures. Only the specified residues will considered as part 
        of the structure and all the other residues will be ignored. 
        See the format description at the end of the OPTIONS section.

    rresneg=STRING, qresneg=STRING [DEFAULT: None]
        The specification of the input reference (rresneg) and query (qresneg) 
        structures. The specified residues will be ignored and all the other 
        residues considered as part of the structure. If both "rres" 
        and "rresneg" (or "qres" and "qresneg") are specified simultaneusly, 
        ARTEMIS will ignore "rres" ("qres") and consider only "rresneg" 
        ("qresneg"). 
        See the format description at the end ofthe OPTIONS section.

    rseed=STRING, qseed=STRING [DEFAULT: rseed=rres,qseed=qres]
        The specification of the reference and query residues that ARTEMIS
        can use for single-residue matching seeds.
        See the format description at the end of the OPTIONS section.

    saveformat=KEYWORD [DEFAULT: saveformat=qformat]
        The specification of the format of the output coordinate files. 
        By default, ARTEMIS will save the coordinate files in the same format 
        as the query input file. If the "saveformat" parameter is specified 
        and it's either "PDB", "CIF", or "MMCIF" (case-insensitive), ARTEM 
        will save the output coordinate files in the specified format.

    saveres=STRING [DEFAULT: saveres=qres]
        The specification of the query structure residues that will be saved 
        in the output coordinate files.
        See the format description at the end of the OPTIONS section.

    saveto=FOLDER [DEFAULT: None]
        Path to the output folder to save the coordinate files 
        of the superimposed query structures along with the mutually 
        closest residue subsets. If the specified folder does not exist, 
        ARTEMIS will create it. If the folder is not specified, 
        nothing will be saved.

    threads=INT [DEFAULT: threads=CPU_COUNT]
        Number of CPUs to use.

    step_divider=INT [DEFAULT: step_divider=0 if len(qseed) < 500 else 100]
        To speed up the procedure of pairwise superpositions of structures 
        to find sets of mutually closest residues, ARTEMIS can skip rseed 
        residuals in steps of
        1 + number of qseed residues // step_divider
        If step_divider is 0, ARTEMIS does not skip rseed residuals.
        If the number of qseed residues exceeds 500 and step_divider
        is not specified, then step_divider automatically becomes 100.

    nlargest=INT [DEFAULT: nlargest=len(qres) if len(qseed) < 500 else 2*threads]
        Number of largest mutually nearest sets of residues for which 
        alignments are constructed.

    shift=FLOAT [DEFAULT: shift=3 if len(qseed) < 500 else 20]
        The value by which the Score Matrix is shifted for Needleman-Wunsch.
        Larger shift, greater coverage.

    -v, --verbose [DEFAULT: OFF]
        If specified, it will add the permutation alignment characteristics
        to the standard output, and save its imposition of the query structure
        with the residue matching table to the saveto folder. 
        The mode is automatically activated if the TM-score of the query 
        structure for alignment with permutation is 10% larger than
        the TM-score of alignment without permutation.

    ***************************************************************************
    ARTEMIS uses a ChimeraX-like format to specify the residues of interest 
    using the "res" parameters:

    [#[INT]][/[STRING]][:[STRING][_INT[CHAR|_INT]] - The structure specification
                                                     format. The "res" 
                                                     parameters can be defined 
                                                     with a number 
                                                     of specifications 
                                                     separated by spaces and 
                                                     enclosed in double quotes.

        #[INT]                    == Model number
        /[STRING]                 == Chain identifier
        :[STRING][_INT[CHAR|_INT] == Residue(s) specification:
            
            :STRING     == Residue type    
            :_INT[CHAR] == Residue number [with insertion code]
            :_INT_INT   == Range of residue numbers
            
    Structure specification examples:

    rres="#1/B:_-10_20 #1/A"    - Consider the entire chain A from model 1 
                                  and the range of chain B residues with 
                                  numbers from -10 to 20 from model 1 as 
                                  the reference structure.
    qres="#"                    - Consider all the residues from all 
                                  the models in the "q" file as 
                                  the query structure.
    saveres="/C:_10_20 /C:_20A" - Save the chain C residues with numbers 
                                  from 10 to 20 and the chain C residue 
                                  with number 20A (A is the insertion code).
    rseed=:A                    - Use only the model 1 adenosines as the 
                                  single-residue seeds from the reference 
                                  structure. 

## Contacts

David Bogdan, *e-mail: <dbohdan@iimcb.gov.pl>*
