# ARTEMIS

Using [ARTEM](https://github.com/david-bogdan-r/ARTEM) to Infer Sequence alignment (ARTEMIS) is a tool for the RNA 3D structure superposition and determination of the structure-based sequence alignment

## Reference

## How ARTEMIS works

ARTEM reads a reference and a query structure from the specified coordinate files in PDB or in mmCIF format, and, by default, prints a sorted list of their local structural superpositions. The user can choose to save the superimposed versions of the query structure into a specified folder in PDB or in mmCIF format. Each of the saved files will include three models: (1) the entire (according to the "qres" parameter) superimposed query structure, (2) the subset of the reference structure residues used for the superposition, and (3) the subset of the query structure residues used for the superposition. By default, ARTEM reads the entire first models of the input files.

The ARTEM algorithm works as follows. For each possible single-residue matching seed between the reference and the query structures (as defined by "rseed" and "qseed" parameters) ARTEM superimposes the structures based on the 5-atom representations of the two residues. Then, ARTEM searches for a subset of the reference and query residues that are mutually closest to each other in the current superposition. Finally, ARTEM performs a second superposition using the subset of the mutually closest residues as the assigned residue-residue matchings. Finally, ARTEM prints a sorted list of the produced superpositions to stdout. For each superposition the output includes its ID, RMSD, SIZE, RMSD/SIZE ratio, the list of generative single-residue seeds (PRIM), and the list of the residue-residue matchings (SCND).

## Installation

Clone the GitHub repository by typing

    git clone https://github.com/david-bogdan-r/ARTEMIS.git

## Dependencies

ARTEMIS requires three Python3 libraries to be installed:

- numpy
- pandas
- scipy

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

## Time & memory usage

## Usage

    python3 artem.py r=FILENAME q=FILENAME [OPTIONS]

## Usage examples

    1) python3 artem.py r=examples/1ivs.cif  q=examples/1wz2.cif rres=/C qres=/C > output.txt

    This command will write into "output.txt" file a sorted list of all local 
    structural superpositions between the C chains of 1IVS and 1WZ2 PDB entries 
    that are two tRNAs. The user can spot the three largest mathings of size 52. 
    Then the user can save the largest mathings only into "result" folder in 
    PDB format:
    
    python3 artem.py r=examples/1ivs.cif  q=examples/1wz2.cif rres=/C qres=/C sizemin=52 saveto=result saveformat=pdb
    
    As the result three files of three different matchings of 52 residues will 
    be saved in PDB format - 1wz2_1.pdb, 1wz2_2.pdb, 1wz2_3.pdb. The latter is
    the superposition with the best RMSD. Each of the saved files lists three 
    structural models. The first model contains all "qres" residues 
    superimposed with the reference structure. The second model contains 
    the subset of the reference structure residues that were used 
    for the superposition, and the third model stores their counterpart 
    residues from the query structure. Then, the user can open the reference 
    file 1ivs.cif together with the first model of the file 1wz2_3.pdb in a 3D
    visualization tool for visual examination of the best local superposition 
    of the two structures. The user can observe a good superposition 
    of the four standard tRNA helical regions.
    
    2) python3 artem.py r=examples/motif10.pdb  q=examples/motif9.pdb rmsdsizemax=0.25
    
    This command will output a sorted list of those local structural 
    superpositions between the two topologically different motifs of 
    10 and 9 residues respectively that have a ration RMSD/SIZE under 0.25.
    The user can spot the only mathing of size 8 with RMSD of 1.713 angstroms.
    Then the user can save the superposition into "result" folder 
    in CIF format:
    
    python3 artem.py r=examples/motif10.pdb  q=examples/motif9.pdb rmsdsizemax=0.25 sizemin=8 saveto=result saveformat=cif 
    
    The only file will be saved named "motif9_1.cif". Then, the user 
    can open the reference file motif10.pdb together with the file 
    motif9_1.cif in a 3D visualization tool for visual examination. 
    The user can observe a good superposition of all the three stacked 
    base pairs. Simultaneously, two of the three A-minor-forming 
    adenosines have a counterpart residue.

# Options 

    r=FILENAME [REQUIRED OPTION]
        Path to a reference structure in PDB/mmCIF format. For faster 
        performance, it's advised to specify the largest of the two structures 
        as the reference structure.

    q=FILENAME [REQUIRED OPTION]
        Path to a query structure, the one that ARTEM superimposes to 
        the reference, in PDB/mmCIF format.

    matchrange=FLOAT [DEFAULT: matchrange=3.0]
    	The threshold used for searching the mutually closest residues. Only 
    	those pairs of residues that have their centers of mass at a distance 
    	under the specified matchrange value can be added to the subset 
    	of the mutually closest residues. The higher matchrange value 
    	will produce more "noisy" matchings but won't miss anything. The lower 
    	matchrange value will produce more "clean" matchings but 
    	can miss something.

    rformat=KEYWORD, qformat=KEYWORD [DEFAULT: rformat=PDB,qformat=PDB] 
        The specification of the input coordinate file formats 
        (case-insensitive). By default, ARTEM tries to infer the format 
        from the extensions of the input filenames. ".pdb", ".cif", 
        and ".mmcif" formats can be recognized (case-insensitive). In the case 
        of any other extension ARTEM will treat the file as the PDB-format 
        file by default. If the "rformat" ("qformat") parameter is specified 
        and it's either "PDB", "CIF", or "MMCIF" (case-insensitive), 
        ARTEM will treat the reference (query) coordinate file
        as the specified format.

    rmsdmin=FLOAT, rmsdmax=FLOAT [DEFAULT: rmsdmin=0,rmsdmax=1e10]
    	The specification of minimum and maximum RMSD thresholds. 
    	ARTEM will print and save only the superpositions that satisfy 
    	the specified thresholds. 
    	
    rmsdsizemin=FLOAT, rmsdsizemax=FLOAT [DEFAULT: rmsdsizemin=0,rmsdsizemax=1e10]
        The specification of minimum and maximum RMSD/SIZE ratio thresholds. 
        ARTEM will print and save only the superpositions that satisfy 
        the specified thresholds. 

    resrmsdmin=FLOAT, resrmsdmax=FLOAT [DEFAULT: resrmsdmin=0, resrmsdmax=1e10]
    	The specification of minimum and maximum per-residue RMSD threshold.
    	ARTEM will print and save only the superpositions that satisfy 
    	the specified thresholds.

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
    	ARTEM will ignore "rres" ("qres") and consider only "rresneg" 
    	("qresneg").
        See the format description at the end of the OPTIONS section.

    rseed=STRING, qseed=STRING [DEFAULT: rseed=rres,qseed=qres]
        The specification of the reference and query residues that ARTEM can use
        for single-residue matching seeds.
        See the format description at the end of the OPTIONS section.

    saveformat=KEYWORD [DEFAULT: saveformat=qformat]
        The specification of the format of the output coordinate files. 
        By default, ARTEM will save the coordinate files in the same format 
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
        ARTEM will create it. If the folder is not specified, 
        nothing will be saved.

    sizemin=FLOAT, sizemax=FLOAT [DEFAULT: sizemin=1,sizemax=1e10]
        The specification of minimum and maximum SIZE thresholds. 
        ARTEM will print and save only the superpositions that satisfy 
        the specified thresholds. If sizemin is set to zero, ARTEM will 
        output the dummy 0-size superpositions matching the reference 
        and query residues lacking the 5-atom representation specified. 
        The user can specify custom atom representations of any residues 
        via editing the lib/nar.py text file.

    trim=BOOL [DEFAULT: trim=None]
    	When specified, for each subset of mutually closest residues ARTEM will 
    	iteratively remove residues with the worst per-residue RMSD from the 
    	subset one by one with the following re-superpositioning based on the 
    	remaining residues of the subset until the specified thresholds for rmsdmax,
    	rmsdsizemax, resrmsdmax are reached or the subset size is less than sizemin.

    threads=INT [DEFAULT: threads=1]
        Number of CPUs to use. ARTEM multiprocessing is available only for 
        UNIX-like systems.

    ***********************************************************************
    ARTEM uses a ChimeraX-like format to specify the residues of interest 
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

# Contacts

David Bogdan, *e-mail: bogdan.d@phystech.edu* 