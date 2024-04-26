# ARTEMIS

Using [ARTEM](https://github.com/david-bogdan-r/ARTEM) to Infer Sequence alignment (ARTEMIS) is a tool for the RNA 3D structure superposition and determination of the structure-based sequence alignment.

## Reference

[D.R. Bohdan, J.M. Bujnicki, E.F. Baulin (2024) ARTEMIS - a method for topology-independent superposition of RNA 3D structures and structure-based sequence alignment. BioRxiv. DOI: 10.1101/2024.04.06.588371](https://doi.org/10.1101/2024.04.06.588371)

## How ARTEMIS works

ARTEMIS reads a reference and a query structure from the specified coordinate files in PDB or in mmCIF format and identifies the two best superpositions.
The first superposition is strictly backbone-dependent (i.e. does not allow backbone permutations) and is fed to the standard output.
In verbose mode, the details of the second superposition, which allows permutations, are additionally fed to the standard output.
The second superposition details are printed automatically if its TM-score is at least 10% higher than the TM-score of the first alignment.
The user can choose to save the superpositions of the query structure in PDB or mmCIF format along with the residue matching tables in the specified folder.
By default, ARTEMIS saves the output files only for the first alignment, unless the second one is at least 10% better in TM-score.
By default, ARTEMIS reads the entire first models of the input files.

## Installation

Clone the GitHub repository by typing

    git clone https://github.com/david-bogdan-r/ARTEMIS.git

## Dependencies

ARTEMIS requires four Python3 libraries to be installed:

- numpy
- pandas
- scipy
- matplotlib

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
****

## Usage examples

    1) python3 artemis.py r=examples/6ugg/6ugg.cif  q=examples/1ivs.pdb rres=/B qres=/C saveto=result saveformat=pdb

    This command will superimpose chain C from 1ivs to chain B from 6ugg 
    and will save the superimposed structure into "result" sub-folder
    in PDB format. 
    
    2) python3 artemis.py r=examples/6ugg/6ugg.cif  q=examples/1ivs.pdb rres=/B qres=/C saveto=result saveformat=pdb -v

    The same command as in example 1, but with verbose mode enabled. 

## Options

    r=FILENAME [REQUIRED OPTION]
        Path to a reference structure in PDB/mmCIF format. For faster 
        performance, it's advised to specify the largest of the two 
        structures as the reference structure.

    q=FILENAME [REQUIRED OPTION]
        Path to a query structure, the one that ARTEMIS superimposes to 
        the reference, in PDB/mmCIF format.

    matchrange=FLOAT [DEFAULT: matchrange=3.5]
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
        structures. Only the specified residues will be considered as part 
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
        and it's either "PDB", "CIF", or "MMCIF" (case-insensitive), ARTEMIS 
        will save the output coordinate files in the specified format.

    saveres=STRING [DEFAULT: saveres=qres]
        The specification of the query structure residues that will be saved 
        in the output coordinate files.
        See the format description at the end of the OPTIONS section.

    saveto=FOLDER [DEFAULT: None]
        Path to the output folder to save the coordinate files 
        of the superimposed query structure along with the mutually 
        closest residue subsets. If the specified folder does not exist, 
        ARTEMIS will create it. If the folder is not specified, 
        nothing will be saved.

    threads=INT [DEFAULT: threads=CPU_COUNT]
        Number of CPUs to use.

    stepdiv=INT [DEFAULT: stepdiv=0 if len(qres) < 500 else 100]
        The step divider parameter. The parameter is used to speed up 
        the procedure for large structures. If stepdiv > 0, ARTEMIS
        will consider only each Sth reference residue as matching seed,
        where S = 1 + len(qres)//stepdiv. If the size of the query 
        structure exceeds 500 residues, stepdiv will be set to 100 
        by default.

    toplargest=INT [DEFAULT: toplargest=len(qres) if len(qres) < 500 else 2*threads]
        Number of largest mutually closest residue sets for which 
        alignments are constructed.

    shift=FLOAT [DEFAULT: shift=3 if len(qres) < 500 else 20]
        The shift2 value for the ScoreMatrix used in the Needleman-Wunsch
        algorithm. Larger shift2 provides higher coverage.

    -p, --permutation [DEFAULT: OFF]
        Permutation mode. If specified, ARTEMIS will add the topology-independent 
        alignment details to the standard output, and save its output files
        to the saveto folder (if the folder is specified). 
        The mode is automatically activated if the TM-score of the query 
        structure for topology-independent alignment is at least 10% higher 
        than the TM-score of the backbone-dependent alignment.

    -v, --verbose [DEFAULT: OFF]
        Verbose mode.

    -superonly, --superonly [DEFAULT: superonly=False]
        If specified, ARTEMIS will assume the input sequences 
        are of the same length and will try to superimpose them 
        with the perfect sequence alignment.

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

David Bogdan, *e-mail: <dav.bog.rom@gmail.com>*
