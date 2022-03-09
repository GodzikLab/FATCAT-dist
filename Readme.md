# FATCAT package

# (Flexible structural AlignmenT by Chaining AFPs with Twist)

## Installation 

Run 
>./Install

If Install fails, please go to the FATCATMain directory and modify the basic.h file that includes all the head files.
Then run 
>make


---

## FATCAT usages


### Program 1: FATCAT (main program)


- purpose: run Flexible alignments for a pair of proteins with given options
- command:  FATCAT [parameters...]
- for help: simply run FATCAT
- inputs:  two protein structures, and an output name (eg. pdb1.pdb2)
- avaiable outputs: 

      the short report of the results                --  to stdout (for system analysis)

      the detailed report of the results             --  to stdout (for system analysis)

      the alignment between two structures in text   --  to file or stdout (for system analysis)

      the chaining result in text file               --  pdb1.pdb2.chain.txt (for case study)

      the .ps file of all AFPs and final AFP chain   --  pdb1.pdb2.afp.ps (for case study)

      pdbfiles based on block superimpose            --  pdb1.pdb2.(block-index).pdb
                                                         pdb1.pdb2.(block-index).script (for rasmol) (for case study)

      pdbfiles based on optimized block superimpose  --  pdb1.pdb2.(block-index).pdb
                                                         pdb1.pdb2.(block-index).script (for rasmol) (for case study)

      twisted structure based on blocks              --  pdb1.pdb2.ini.twist.pdb
                                                         pdb1.pdb2.ini.twist.script (for case study)               

      twisted structure after align optimization     --  pdb1.pdb2.opt.twist.pdb
                                                         pdb1.pdb2.opt.twist.script (for case study)                



note: all the output pdb files include two chains with chain A for pdb 1 and chain B for pdb 2

----

### Program 2: FATCATQue.pl

>purpose: a script for runing FATCAT on a list of structure pairs 

>command: FATCATQue.pl log-file pair-list-file <parameters..>

----

### Program 3: FATCATSearch.pl

purpose: a script for runing FATCAT for a query against a database
command: FATCATQue.pl query target-list parameters..
note: for this purpose, you need to prepare a database of protein structures

---

##  FATCAT Tests

### 1. FATCAT for a single pair of structures

Go to ./Examples_FATCAT,

run FATCAT for structure 1a21A.pdb and 1hwgC.pdb 

    command: FATCAT -p1 1a21A.pdb -p2 1hwgC.pdb -o 1a21A_1hwgC -m -ac -t

    inputs: 
      structure 1a21A.pdb and 1hwgC.pdb

    outputs: (a twist is detected in comparing these two structures)

      1a21A_1hwgC.aln               # alignment result between the two structures
      1a21A_1hwgC.afp.color.ps      # the AFP chaining result of the FATCAT in postscript graph
      1a21A_1hwgC.opt.twist.pdb     # the pdb file in which 1a21A.pdb in chain A, twisted 1hwgC.pdb in chain B  
				                    # refer REMARK lines in the begining of .pdb file
      1a21A_1hwgC.opt.twist.script  # the rasmol script for viewing 1a21A_1hwgC.opt.twist.pdb 
				                    # blocks of 1hwgC.pdb are shown in different colors
				                    # command: rasmol -script 1a21A_1hwgC.opt.twist.script

Please compare 1a21A_1hwgC.aln with compare.aln file. 
If no difference is found between these two files, FATCAT runs ok.


### 2. FATCAT for a list of pairs

Go to ./Examples2_FATCAT

    command:
     ../FATCATMain/FATCATQue.pl timeused allpair.list -q >allpair.aln
    (note: this calculation takes about 30 min on Intel(R) XEON(TM) CPU 1.80GHz)

    input: allpair.list (the list of pairs) and the corresponding pdb files
    output: allpair.aln (the text file of all the alignments)

Please compare allpair.aln with allpair.aln.comp. If no difference is found, FATCAT runs ok.

for testing a smaller set of pairs:

    command:
    ../FATCATMain/FATCATQue.pl timeused prtpair.list -q >prtpair.aln
(note: this calculation takes about 10 min on Intel(R) XEON(TM) CPU 1.80GHz)
