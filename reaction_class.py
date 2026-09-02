# Reads inputed excel sheets and makes each entry an object as Reaction Class

from dataclasses import dataclass

@dataclass
class Reaction:
    # The constructor method initializes the reaction class
    doi: str                        # paper DOI
    in_paper_data_loc: str          # location of reaction in paper
    entry: str                      # entry number for specific paper, defined by inital data parsing
    sub_smiles: str                 # substrate smiles  
    sub_cpd_num_in_paper: str       # substrate compound number as defined by paper
    rh_pre_cat: str                 # rhodium precatalyst used for the reaction
    rh_pre_cat_conc: float          # mol % of rhodium catalyst used
    ligand_name: str                # name of ligand including stereochemistry
    ligand_smiles: str              # liangd smiles (from Claude)
    ligand_cpd_num: str             # ligand compound number as defined by paper
    solvent: str                    # solvent used in reaction
    reaction_time_hr: float         # reaction time
    prod_smiles: str                # product smiles (from Claude)
    prod_cpd_num: str               # product compound number as defined by paper
    byprod_smiles: str              # by product smiles (from CLaude)
    perc_yield: float               # percent product yield
    prod_yield_method: str          # experimental method for determining yield
    perc_ee: float                  # enantioselectivity of product (percent)
    dr: float                       # diastereomeric ratio
    perc_byprod_yield: float        # byproduct percent yield
    perc_byprod_ee: float           # enantioselectivity of byproduct (percent)
    byprod_yield_id_type: str       # experimental method for determining byproduct yield
    ee_method: str                  # experimental method for determining enantioselectivity
    prec_rem_sm: float              # percent remaining starting material
    notes: str                      # any additional information identified during data parsing


