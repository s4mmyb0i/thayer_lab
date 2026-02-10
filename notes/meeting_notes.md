# Notes

## February 10, 2026

Talk to Erica Taylor about exchange program in India (asian countries).
Ed-equity program about teaching 112 in India - speak with Jeffery Goetz

TODO: Speak with Pavel and Rose maybe Krizanc
TODO: Write cells of information for the paper
TODO: Re do the outline for CS department

### Outline

Important: remember that the paper is for CS department.

Write key points for discussion, but write the actual thing last

Project of networks, happens to be of a protein.

Make a thesis sentence with a title.
PDZ is a good model protein, because of extensive literature, and it happens to have some interesting things as well (allostery) so might as well look into that.

Title: Networks as a proxy for long range allosteric signalling in PDZ.

Capturing long range signalling is a problem, but go into computer science has networks.

Background:
- Networks
- Random Walks
- Markov Mathematics
  - Absorption probs (derivations)
- Bio stuff
  - Biology of PDZ (look at older papers with Lakhani and Thayer MD-MSM)
- Summary of previous work on these trajectories (bharat's paper, MD-MSM, MD-sector papers) (to figure out vocab for comparison of different methods)
- "sectors is debated phenomenon, whether it exists or not actually, we are doing something else using the entire network"
- 
  
Methods (should be citable):
- softwrae packages (python libraries)
- numpy
- plotly
- matplotlib
- bharat's simulations (MD-MSM, MD-sector papers)
- 

Results:
- the method i used to generate the random walks
- absorption probs calculations and the related figures
- random walks analyses
- bottleneck analysis
- 

## February 7, 2026

### Scripts breakdown

## In `analysis/`

1. **`node_heatmap.py`**  
   - **What it does:** Builds an interactive 2D heatmap of **node visit frequencies** from random-walk CSVs (only walks with `termination_code == 1`).  
   - **Outputs:** Writes **`heatmap.html`** (or `--output <file>`), and prints stats (total walks, unique nodes, most/least visited nodes, top 10).  
   - **Run:**  
     `python3 analysis/node_heatmap.py <csv_file> [--output heatmap.html]`

2. **`node_heatmap_3d.py`**  
   - **What it does:** Builds an interactive **3D heatmap** of atoms/nodes visited, using PDB coordinates and the same walk CSVs.  
   - **Outputs:** Writes **`heatmap_3d.html`** (or `--output`), and prints similar stats.  
   - **Run:**  
     `python3 analysis/node_heatmap_3d.py <pdb_file> <csv_file> [--output heatmap_3d.html]`

3. **`visualize_paths.py`**  
   - **What it does:** Visualizes random-walk paths in **3D** (PDB + walk CSV), with sliders for “top K paths” and max path length.  
   - **Outputs:** **Interactive only** (`fig.show()`). No file is saved unless you add a `fig.write_html(...)` call.  
   - **Run:**  
     `python3 analysis/visualize_paths.py <pdb_file> <csv_file> [max_paths]`

4. **`visualize_highways.py`**  
   - **What it does:** Aggregates paths into **“highways”** (edges weighted by how many paths use them) and shows them in 3D.  
   - **Outputs:** **Interactive only**; no file output by default.  
   - **Run:**  
     `python3 analysis/visualize_highways.py <pdb_file> <csv_file> [max_paths] [min_edge]`

5. **`bottleneck.ipynb`**  
   - **What it does:** Loads the allos→lig walk CSVs (P, PL, AP, APL), keeps successful paths, counts node frequencies, and plots **bar charts** (e.g. top nodes).  
   - **Outputs:** In-notebook figures and tables only; no written files unless you add export (e.g. `fig.write_html` or `to_csv`).

## In `linalg/absorption_probs/`

6. **`absorption_probs.py`**  
   - **What it does:** Takes a correlation/DCCM-derived matrix, builds a Markov transition matrix with absorbing nodes, and computes **absorption probabilities** (probability of being absorbed at each absorbing state from each transient state).  
   - **Outputs:** Writes **`./linalg/absorption_probs/absorption_probabilities.csv`** and prints the max value and its location.  
   - **Run:** As a module (see docstring), e.g.  
     `python3 -m linalg.absorption_probs.absorption_probs <matrix.npz> <output_dir>`  
   - The commented block in the same file (and the presence of `.fds/` and `old_data/` with `*_probs.csv`, `*_probs.npz`, `*_heatmap_absorption.html`, `*_bar_absorbing_influence.html`, `*_scatter_entropy.html`, `analysis_summary.txt`) suggests there is or was a pipeline that produces those extra artifacts (CSV, NPZ, HTML, summary); that part may live in another script or an older version.

## Upstream of “analysis” (produces data other scripts analyze)

7. **`simulate_walk/simulate_walk.py`**  
   - **What it does:** Runs **random walks** on the transition matrix (from each node, many walks, max steps).  
   - **Outputs:** Writes **CSV files** into `random_walk/` (e.g. `P_allos->lig_walks.csv`, `PL_allos->lig_walks.csv`, etc.) and `.log` files in `simulate_walk/`.  
   - **Run:**  
     `python3 -m simulate_walk.simulate_walk <matrix.npz> <max_steps> <n_walks> <random_walk_dir>`  
   - These CSVs are what `node_heatmap*.py`, `visualize_paths.py`, `visualize_highways.py`, and `bottleneck.ipynb` use.

**Summary:** The scripts that **write tangible analysis outputs** are: **`node_heatmap.py`**, **`node_heatmap_3d.py`** (HTML + stats), **`absorption_probs.py`** (CSV + console), and **`simulate_walk.py`** (walk CSVs). The path/highway visualizers and the bottleneck notebook give interactive or in-notebook results only unless you add explicit saves.


### Walk data

- There are 10,000 walks starting from each atom for both allos->lig and lig->allos walks for each construct.
- 


## February 3, 2026

both sides random walk
analyses from walks
markovian absorption prob analyses
pavel final project

make an outline, identify holes in paper and see how to fill them in
put it in a slide (powerpoint) with all results (what order do the results go in)

thesis committee members
- rose
- krizanc
- pavel
- 

see when i need to decide if i'm going in the fall

look into ba/ma program (it'll be free)


## January 30, 2026

Thesis tutorial abstract/request box thingy

Allostery is a special case of cooperativity in protein-ligand binding, where the binding of a molecule at one site on a protein affects the function at a different, distal site. Understanding the mechanism of thermocyclic signal transmission in the PDZ3 domain is essential, as PDZ3 is a crucial component in cellular signaling and the regulation of protein-protein interactions, and thus a strong candidate to be a model protein. Moreover, we know that the thermocycle directly influences allosteric regulation in PDZ3. This knowledge has significant potential for advancing our understanding of the signaling pathways involving the p53 protein, a major factor in cancer research.

This study aims to characterize the signaling pathways in the PDZ3 domain by creating a protein structure network with atoms as nodes and edges representing assumed interaction strengths. We will then perform Markov Transient Analyses and simulate random walks to map the communication channels in the network.
Through these analyses, we can understand how the allosteric and ligand binding sites and connected.


## Sometime over winter break

compare subset of nodes with each other copmare with walks from all nodes to binding site. ask questions about how interconnected these nodes are, generate datasets for null model
ask if signalling is more pronounced than background

select non paths ways and show that hte ginal is weaker (from log_prob)
see if statistically differnet
use this to normalize between dconstructrs
delta delta g
refernce to own background noise, and then ask for same pathways when A is there, is it better than if not at all

numbers wouldnt matter with A or L 

NULL HYPOTHESIS?


see reflexivity in with A and L


look at hydopohbic res with ability to flip (rings) see if the sidchains should be reoved. how does this affect the walks?


chosose non end target atoms and compare paths (is there somethign special about the successful paths?  are the sum logs less with hte non L sites?)

## December 14, 2025

Ligand construct is the 1RZX pdb.
Allosteric construct is the 1NF3 pdb.
Free Par6 PDZ is 1RY4 pdb.


Ligand-binding site (direct binding)

➡️ P38, Y42, R66, T102, M105

Allosteric effector site (Cdc42)

➡️ Thr25, Asp38, Asn39, Tyr40, Met45, Tyr64, Leu67, Leu70

Allosteric communication/interface residues

➡️ Ile3, Ser4, Pro6, Phe9, Arg10, Pro11, Ile15, Ser13, Ser14, Asp17, Val18, Asp19, Arg26, Leu78


## December 12, 2025

Get analysis on the walks (IMPORTANT!!!) - check in on email next wed, meet on this after
To what extent do these analysis work?
See if pavel's projet matches with my project, and identify bottlenecks, and cross refernce with finn to find bottlenecks
Run on sean's trajs, using Finns analyses. Could wind up having the same conclusion

## December 5, 2025


Try this on p53 and do simulation on very short 1 ns or 10 nanosecond

i have walks,
do a glass visualizainot that shows the pathways bfrom ligand pocket to allsoteric pocket
look at percetange of walks that are successful for both sets
extract the actual paths (through the atoms/residues)
find bottleneck residues, and break oit  if broken, the pathways should stop working. to break, replace bottlemext residue awith alanine (its something thats done with genetics)
arginene in p53 is usualyl the bottle nexk because its positve and so big, mayube similar in pdz
bottleneck is if almost all the walks go through one node. i could maybe find literatuer and use that info for the mutatnt, identify the bottlenecks and do the simualtions
use chatgpt deep studies to get lierature review


maybe write the paper for JOSS, so it is also good cor chem and CS so its high quality chem and high quality CS.

write email to pavel and thayer about ysug lab networks for final project


give a write up of a summary of my leetcode practice, give screenshots of questions completed, and some sample code that ive written. 1 to 2 page summary. list which algorithms, brief disceritopn, reflection on what it was like to do this. which one is my favorite, easiest and more diffcult one. how has this supplementedcourse offering at wesleyan. make a story, about how I started off with doing data structures to then using that to do algortihms (since i already did algorithms class at school). 

## November 21, 2025

Ligand construct is the 1RZX pdb.
Allosteric construct is the 1NF3 pdb.
Free Par6 PDZ is 1RY4 pdb.


From 1NF3 paper, binding pocket details are:

Cdc42 binding pocket: Thr25, Asp38, Asn39, Tyr40, Met45, Tyr64, Leu67, Leu70
Par6 binding interface: Ile133, Ser134, Pro136, Phe139, Arg140, Pro141, Ser143, Ser144, Asp147, Asp149, Ile145, Val148, Arg156, Leu208

Need to identify the Par6 binding pocket residues on the P pdb (since these matched with the clustal sequence alignment algorithm)

Binding pocket on P is the same residues as the Par6 of course since they are the same. The numbering is different though:
Ile3, Ser4, Pro6, Phe9, Arg10, Pro11, Ser13, Ser14, Asp17, Asp19, Ile15, Val18, Arg26, Leu78

From 1rzx paper:
ligand binding pocket is:
P171, Y175, R199, T235, M238

On P, this maps to:
P38, Y42, R66, T102, M105



## November 14, 2025

Look at relative differences between signalling between the successfull paths .. like gibbs free energy calc. 

what we have are baselines. 

we want to see signalls betwee nteh binding pocjet and the effector. 
(delta (delta)) 

number of steps and nuymber of times the walks go to the binding pocket.

do a scoring scheme to keep track of if it does the walks and in how many steps. 

do half of size of vertices for max steps in random walk. (this is for supplementary). Play around with walk lengths. going from allostertic effector to binding pocket. i can do walks backwards. 

look at crystal structure to find allosteric effector surface 

start from effector site, and then find the differential between paths that go to the binding pocket. is the entire surface important for allostery or just a couple of residues. find which residues on the surface have higher proiibabilitiy to hit the binding pocket. 

frequency of hits i a proxy for free energy

stat mac view (pcls class). use rules to plug into to  calculate actual energy avlues

To organize;

- markov analyses to provide sturucture for random walks
- random walk from binding pocket
- random walk from allosteritc surface


chat with nando about class offering at bigger universities in edingurgh

speak with department chairs about getting clarity on thesis. 

## November 7, 2025

Prevalence in pathways in proteins Kelly, Avram Sty, Jesse Galganoff paper has info on how the absorption capabilities are so low because there are too many pathways. As the scale of the protein gets bigger, then there are more paths

Use ROC curves. pick a delta, and put multiple thresholding values are. See how stable the cutoffs are. Jonathan Fabry did thesis on this. 
Want to check how volatile the system is. 
We have 20% from existing literature. Find a cutoff that maintains 20% of strongest signals. 
Use gephi and use the slider to see these different cutoffs.

Have a paragraph summary to detial what I'm working on with LeetCode half credit.

Look into GCC summer grant to cover summer research. 
Lookinto REUs (10-20 places). Look up old essay prompt sand write something about that.
2-3 letters of rec. Find people for that.
Try to get things in early.
After applyimg, email PIs that I am excied to be in your lab and tell them that I could be free with GCC grant. Tell them about it.
Speak with PI after applying.
Speak with Fernando about comp biophysical chem. UTMB has some cool labs.
MIT has had self funded programs.
Leticia Costa did MIT summer program
Hari's lab (MMR focused)


## October 31, 2025

We don't need 1BFE and 1BE9, it is 1NF3.
From 1NF3 paper, binding pocket details are:

Cdc42 binding pocket: Thr25, Asp38, Asn39, Tyr40, Met45, Tyr64, Leu67, Leu70
Par6 binding interface: Ile133, Ser134, Pro136, Phe139, Arg140, Pro141, Ser143, Ser144, Asp147, Asp149, Ile145, Val148, Arg156, Leu208

Need to identify the Par6 binding pocket residues on the P pdb (since these matched with the clustal sequence alignment algorithm)

Binding pocket on P is the same residues as the Par6 of course since they are the same. The numbering is different though:
Ile3, Ser4, Pro6, Phe9, Arg10, Pro11, Ser13, Ser14, Asp17, Asp19, Ile15, Val18, Arg26, Leu78






IVISMPQDFRPVSSIIDVDILPETHRRVRLCKYGTEKPLGFYIRDGSSVRVTPHGLEKVPGIFISRLPVGGLAQSTGLLAVNDEVLEVNGIEVSGKSLDQVTDMMIANSRNLIITVRPANQRN

```
CLUSTAL O(1.2.4) multiple sequence alignment


1nf3_chain_c      RKKPHIVISMPQDFRPVSSIIDVDILPETHRRVRLCKYGTEKPLGFYIRDGSSVRVTPHG	60
P                 -----IVISMPQDFRPVSSIIDVDILPETHRRVRLCKYGTEKPLGFYIRDGSSVRVTPHG	55
                       *******************************************************

1nf3_chain_c      LEKVPGIFISRLVPGGLAQSTGLLAVNDEVLEVNGIEVSGKSLDQVTDMMIANSRNLIIT	120
P                 LEKVPGIFISRLPVGGLAQSTGLLAVNDEVLEVNGIEVSGKSLDQVTDMMIANSRNLIIT	115
                  ************  **********************************************

1nf3_chain_c      VRPANQRN	128
P                 VRPANQRN	123
                  ********

```

We're not fucked!



## October 14, 2025

Uhhh we're fucked
```
$ python3 -m notes.test
✅ Absorption probability matrix saved to:
/Users/sam/Documents/Wesleyan/Thayer_Lab/2025/summer/bharat/notes/absorption_probabilities.csv
         Absorb 517  Absorb 518  Absorb 519  Absorb 520  Absorb 521  Absorb 522  Absorb 523  Absorb 524  Absorb 525  ...  Absorb 1521  Absorb 1522  Absorb 1523  Absorb 1524  Absorb 1525  Absorb 1526  Absorb 1527  Absorb 1528  Absorb 1529
Start 0    0.007840    0.007486    0.007160    0.007089    0.005490    0.005245    0.003911    0.004738    0.003745  ...     0.003138     0.002892     0.003596     0.003963     0.003618     0.003311     0.002990     0.005283     0.004641
Start 1    0.007834    0.007477    0.007157    0.007085    0.005488    0.005247    0.003908    0.004742    0.003747  ...     0.003140     0.002900     0.003592     0.003951     0.003610     0.003301     0.002997     0.005269     0.004627
Start 2    0.007845    0.007495    0.007159    0.007087    0.005487    0.005244    0.003901    0.004743    0.003749  ...     0.003134     0.002898     0.003596     0.003980     0.003637     0.003320     0.002974     0.005291     0.004645
Start 3    0.007824    0.007470    0.007146    0.007073    0.005484    0.005240    0.003914    0.004723    0.003727  ...     0.003145     0.002880     0.003601     0.003953     0.003615     0.003304     0.003000     0.005271     0.004630
Start 4    0.007863    0.007503    0.007182    0.007122    0.005503    0.005251    0.003923    0.004750    0.003758  ...     0.003127     0.002888     0.003592     0.003968     0.003615     0.003316     0.002988     0.005311     0.004670

[5 rows x 224 columns]
Largest value: 0.008670
Located at: Row = Start 515, Column = Absorb 517
```
these percentage values are way to low

## October 10, 2025

Try a toy problem (try by hand and then do with code) (2 transition and 2 absorbing)
Do it in  a way that is like code. Think about what each build in im using and its input and output
Start with mini correlation matrix (raw counts) (10x10 correlation matrix 6 nodes for things to go through, 4 absorbing)
Keep it around, use as sample data in supplementary material
Maybe co sponsor a hand calculation a workshop
email mia / dalton to askif if swap m2m


## October 3, 2025

Write some writeups fro kinds of problems I'm doing for DSA, and this can be used for background in thesis (around 3-4 of those)
Out of order execution (for making code run faster, an dusing better data strucutres for graph traversal (not transition amtrices))
Pavel bayeisian neural networks
Something for Leo with Pavel, for end of project class



## September 26, 2025

Structurla alignment over sequence alignment

clustal omega does sequence alignment
- pam amtricies (more structural) like comparing hemoglobin and myoglobin
- standard would work fine, but there are options
- will need to work out exxesive gapping, so there is more mapping

Need all reference frames for each pdb file
Choose a reference frame (this is what is used for every reference frame)

Structurally aligned mapping
Visualize the mappings are lining them up

In the project:
theres a "to" and "from", and it connects to each other

For Yehrim

Start with easier transformation:
PDB to MD numbering
i.e. do the one that is subtract 29 from all the numbers
helps to make an easy verification of the code
2 pairs (2 pdbs, 2 mds associaited with them, there is a delta betwee nthem)

Read bharat's paper, and check the alignment methodology
(its supposed to be 1be9, but could be 1bfe)
there should be an exact match with 1 of the 2

it is possible that the p/pl/ap/apl structure has some parts of the allosteric effector in the extracted pdb.
go into the pdb file, extract the  sequence chain using the ter card

use cultal omega for this

dump all the pdbs into clustal
there should be 2 sequence alignment, one is pdz, the other allosteric effector
check bharat's paper to tell what was put where

use the two crystal structures (1be9 is peptide and protein)
1bfe should be two proteins
check that hte two proteins acctually align
use unix tools to get the res sequence (3 letter code)
convert to one letter code for clustal omega
use the faster format:

> 1be9
sequence (ine one letter mapping)
>1bfe
sequence (in one letter mapping)


field width of 60 (fortran) (yehrim can do)


then when i go to MD ones, then the slices should match





## September 19, 2025

Use mapping in Pymol session
Speak with Yehrim Whang (ywhang@wes...) about remapping res num from standard pdb to md pdb standard 
Make 1-2 slides about my work on random walks, include a figure showing promise of walks for hypothesis
Provide headshot
Provide a screenshot of pymol session (to show subset of residues)

Markov:
Pseudo-counts - put a VERY tiny probability value for all transient nodes with probability < threshold
2 walkthroughs: one to set threshold and then second pass to set connectivity to all other nodes (must renormalize)
z


Used the 1BE9 paper to identify the binding pocket residues, for the peptide binding site.
Residues given in the paper:
- Carboxylate-binding loop: Arg-318, Gly-322, Leu-323, Gly-324, Phe-325
- Hydrophobic pocket: Leu-323, Phe-325, Ile-327, Leu-379
- Peptide recognition contacts: Asn-326, Gly-329, Ser-339, His-372
- β-sheet backbone stabilizers: Phe-325, Ile-327

Residues added in post:
328 (to fill out beta sheet) (it's an isoleucine) (two isoleucines next to each other)
375, 376 (to fill in the middle of the helix) (these alanines and are very small, so thats why they might have excluded it from the paper)


In the P structure (and therefore the PL, AP, APL structure (since they have the same res chain)):
32, 38-45, 98, 101-102, 105

























## September 18, 2025


Based on the 1BFE structure, by visual checking, binding pocket of:
P: K _ _ _ _ _ P L G F Y I R D (32, 38-45)

PL K _ _ _ _ _ P L G F Y I R D (32, 38-45)

AP: K _ _ _ _ _ P L G F Y I R D (32, 38-45)

APL:



1BFE paper:
- Carboxylate-binding loop: Arg-318, Gly-322, Leu-323, Gly-324, Phe-325
- Hydrophobic pocket: Leu-323, Phe-325, Ile-327, Leu-379
- Peptide recognition contacts: Asn-326, Gly-329, Ser-339, His-372

- R _ _ _ G L G F N I I G (318, 322-329)
  

1BE9 paper:
- Carboxylate-binding loop: Arg-318, Gly-322, Leu-323, Gly-324, Phe-325
- Hydrophobic pocket: Leu-323, Phe-325, Ile-327, Leu-379
- Peptide recognition contacts: Asn-326, Gly-329, Ser-339, His-372
- β-sheet backbone stabilizers: Phe-325, Ile-327




## Meeting Notes September 12, 2025

### Research & Technical Tasks
- Binding pocket identification
    - PDB & Structural Analysis
    - Read paper associated with PDB file from the web.
    - Superimpose structures of all 4 proteins and analyze binding pocket using PyMOL.
    - Use crystal structure papers to identify binding pockets.
    - Use peptide structure (with labeled residues) as a reference.
    - Superimpose proteins, align residue numbers, and ensure correct mapping.
    - Devise residue number mapping.
    - There is a standard one; need to pick a reference.
    - Use PDZ from the MD simulation (PL or AP variant used).
    - Then find the PDZ’ one.
    - Represent mapping as key–value pairs (value = list of possible places).
    - Yurm (South Korea) can help with this → acknowledge in paper.
- Alternative Approaches
    - Investigate embedding version that doesn’t rely on XYZ superposition.
    - Explore how to feed this into AI pipeline.
    - Consider pivoting to a new project senior year if needed.
- Presentations
    - Prepare a short presentation about work at Axtria.
    - Submit a paper to a journal by the first day of spring semester.

### Personal Development
- Tech
    - Coding Practice
    - Study tutorial for LeetCode with Yehor (confirm with him).
- Grad School Planning
    - Option: Go straight to grad school → higher base pay.
    - Still apply for schools when the time comes.
    - Look into research programs, labs, and PIs (not just schools).
    - Speak with Fernando → some labs he’s considering may be of interest.
    - Email PIs (though some may not respond).
    - Ask about summer internship opportunities.
    - Build relationships with labs early.
    - REU Applications
        - Apply to REUs (10 REU, 10 other programs).
        - Recycle REU applications where possible.

---

### Action Items
1.	Read the paper associated with the PDB file.
2.	Superimpose structures of all 4 proteins in PyMOL and analyze binding pockets.
3.	Label residues in peptide structure → use as reference for mapping.
4.	Devise residue number mapping (choose PDZ reference from MD simulation).
5.	Collaborate with Yurm for mapping assistance.
6.	Prepare Axtria presentation.
7.	Submit paper by spring semester start.
8.	Explore embedding/AI approach for structural analysis.
9.	Check with Yehor for LeetCode tutorial.
10.	Research grad school labs, programs, and PIs.
11.	Speak with Fernando about labs of interest.
12.	Begin contacting PIs regarding summer internships.
13.	Apply to ~20 programs (10 REU, 10 others); recycle applications where possible.

---



---

read paper asssiciated with pdb file off the web


superimose structures of all 4 and read out binding pocket. use pymol. llook at crystal papers they tell where pocket is. strucutre wih peptide use that. use the pic of the peptide with the residues labelled near it
then super impose the protein. look for the numbers in the same place. 

devise res number mapping . there is standard one. need to pick a reference. use pdz  that is in the md. one of PL or AP is the one used in the MD. Then find the pdz' one. key value pair where value is a list of possible places. yurm (a person from south korea) can do this. (ackoledgement in paper)

put together a short presetnaion about my work at axtria


submit a paper by first day of spring semester
maybe find if tehre is an embedding version of this taht is not super postion in xyz
or figure out how to feed thi sinto the ai
maybe pivot to a new project senior year


study tutorial for doing leetcode with yehor (check with him)



grad school:
go stragiht away, higher base pay
still apply for schools when the time comes
speak with fernando since some of his labs hes looking into might be interesting to me
research programs, labs, PIs. Not just schools
email PIs, but they mught nt respond
ask them about taking interns for the summer (develop relatiomship with that lab)


REU apply to them all


10 REU, 10 other
Recycle REU application
