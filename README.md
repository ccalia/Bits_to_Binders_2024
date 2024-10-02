# Team "Tumor Inhibitors"
## BioML Challenge 2024: Bits to Binders

<p align="center">
  <img src="./figures/binder_with_CD20_in_membrane.png" alt="Designed binders with CD20 in a membrane" width="700px" align="middle"/>
</p>

## Team Members

Rana Barghout, Cianna Calia, Ashish Makani, Javier Marchena-Hurtado, Matthew Williams

## Abstract

text text text...

## Workflow

<p align="center">
  <img src="./figures/pipeline_figure.png" alt="Steps of our design process" width="1100px" align="middle"/>
</p>

text text text...

## Files Included

 - **Run_ProteinMPNN_and_AlphaFold2_for_RFdiffusion_Binders_with_ColabDesign_gdrive.ipynb:** Colab notebook for running ProteinMPNN/AlphaFold2 via the ColabDesign framework for RFdiffusion binder backbones.
 - **Filter_RFD_backbones_with_membrane_check.py:** Script for filtering RFdiffusion binder backbones to remove low-quality binder structures, incorrectly located binders, and binder backbones that clash with the cell membrane. Requires a copy of 6y97_default_dppc.mpmd.finalframe.atomistic.pdb from [MemProtMD](https://memprotmd.bioch.ox.ac.uk/_ref/PDB/6y97/_sim/6y97_default_dppc/) with solvent/ions deleted.
 - **Filter_CD_AF_output_fasta_RFD_pdbs.py:** Script for identifying best designs based on AlphaFold2 metrics from outputs of Run_ProteinMPNN_and_AlphaFold2_for_RFdiffusion_Binders_with_ColabDesign_gdrive.ipynb. Basically just organizes the notebook's outputs.
 - **???** Script to convert our spreadsheet of designed sequences to the final fasta for submission, with sequences ranked based on pae_interaction and RMSD.
 - **full_protocol.txt:** Full details of our workflow, including all RFdiffusion commands.
 - **CD20_Binders.fasta:** The fasta file we submitted, with our top 500 sequences.
 - **CD20_Binders.csv:** Spreadsheet containing all sequences we obtained with pae_interaction < 22 and RMSD < 8 Å, and their AlphaFold2 metrics.

## References

text text text...
