from coffea import processor
from tqdm import tqdm

import csv

from uncertainties import ufloat
from uncertainties.umath import * 

from hist.intervals import ratio_uncertainty

import awkward as ak
from coffea.util import load, save

from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

import numpy as np
import config_trigger_processor as config


def build_p4(acc):
    p4 = ak.zip({'x': acc['x'].value, 
                 'y': acc['y'].value,
                 'z': acc['z'].value,
                 't': acc['t'].value}, with_name="LorentzVector")

    return p4

class TriggerProcessor(processor.ProcessorABC):
   
    '''
        
        A coffea processor class
        
     '''

    def __init__(self, analyzer_name):

        '''
        
        Initialize the processor with a name (can be any name) and an accumulator.
        
        '''
        self.analyzer_name = analyzer_name

        self._accumulator = processor.dict_accumulator({})
        
    @property
    def accumulator(self):
        return self._accumulator
     
    def process(self, file_list, pt_val_list, csv_file):

        '''

        This function is used load the summed accumulator for the desired particles. Then for each DimuStar it creates
        an awkward vector, unflatten it and apply the desired trigger.
        
        file_list (list): List with the path to accumulators


        '''

        acc = load(file_list[0])
        for i in range(1, len(file_list)):
            acc += load(file_list[i])

        # Creates the output
        output = self.accumulator.identity()
    
        # Opens the accumulator for each object (particles, vertices, trigger...

        # Takes the number of events with JpsiDstar before trigger
        DimuDstar_acc = acc['DimuDstar']
        HLT_2018_acc = acc['HLT_2018']
        DimuDstar_p4 = build_p4(DimuDstar_acc)     

        ## DimuDstar collection

        # Creates the pt, eta, phi, m lorentz vector.
        DimuDstar = ak.zip({
            'jpsi_mass' : DimuDstar_acc['Dimu']['mass'].value,
            'jpsi_pt' : DimuDstar_acc['Dimu']['pt'].value,
            'jpsi_eta' : DimuDstar_acc['Dimu']['eta'].value,
            'jpsi_phi' : DimuDstar_acc['Dimu']['phi'].value,
            'jpsi_rap' : DimuDstar_acc['Dimu']['rap'].value,
            'dstar_deltam' : DimuDstar_acc['Dstar']['deltam'].value,
            'dstar_deltamr' : DimuDstar_acc['Dstar']['deltamr'].value,
            'dstar_pt' : DimuDstar_acc['Dstar']['pt'].value,
            'dstar_eta' : DimuDstar_acc['Dstar']['eta'].value,
            'dstar_phi' : DimuDstar_acc['Dstar']['phi'].value,
            'dstar_rap' : DimuDstar_acc['Dstar']['rap'].value,
            'dimu_dstar_deltarap' : DimuDstar_acc['deltarap'].value,
            'dimu_dstar_mass' : DimuDstar_p4.mass, #is_jpsi & ~wrg_chg & dlSig & dlSig_D0Dstar
            'is_jpsi' : DimuDstar_acc['Dimu']['is_jpsi'].value,
            'wrg_chg': DimuDstar_acc['Dstar']['wrg_chg'].value,}, with_name='PtEtaPhiMCandidate')  
        
        # Unflatten it to apply trigger
        DimuDstar = ak.unflatten(DimuDstar, DimuDstar_acc['nDimuDstar'].value)
        # Takes available triggers
        hlt_filter_2018 = config.hlt_filter
        hlt_filter = hlt_filter_2018

        HLT_acc = HLT_2018_acc

        print(f"You are running with the trigger(s): {hlt_filter}")
        
        # Loop over trigger list. Applies OR condition to the trigger selection.
        trigger_cut = HLT_acc[hlt_filter[0]].value
        for i in range(0, len(hlt_filter)):
            trigger_cut |= HLT_acc[hlt_filter[i]].value

    
        with open(csv_file, 'w') as f:
            
            writer = csv.writer(f)
            header = ['First_pt[GeV/c]', 'Second_pt[GeV/c]', 'N_candidates_no_trigger', 'N_candidates_trigger', 'Efficiency_Value',
                      'Efficiency_Error_poisson', 'Efficiency_Error_up', 'Efficiency_Error_down',]
            writer.writerow(header)
            
            for c in range(len(pt_val_list)):
                
                fpt = pt_val_list[c]
                if pt_val_list[c] == pt_val_list[-1]:    
                    spt = pt_val_list[-1]
                    continue
                else:
                    spt = pt_val_list[c+1] 
                    
                print(f'First pT: {fpt} [GeV/c], Second pT: {spt} [GeV/c]')

                DimuDstar_tgcut = DimuDstar[(DimuDstar.jpsi_pt > fpt) & (DimuDstar.jpsi_pt < spt)]

                # Number of events with Jpsi Dstar before trigger
                evts_with_jpsidstar_before_trigger = len(DimuDstar_tgcut[ak.num(DimuDstar_tgcut) > 0]) 

                # Trigger cut
                DimuDstar_tgcut = DimuDstar_tgcut[trigger_cut]

                # Number of events with Jpsi Dstar after trigger
                evts_with_dimudstar_after_trigger  = len(DimuDstar_tgcut[ak.num(DimuDstar_tgcut) > 0])

                # Error calculation
                eff_error_new_down, eff_error_new_up = ratio_uncertainty(evts_with_dimudstar_after_trigger, evts_with_jpsidstar_before_trigger, uncertainty_type='efficiency')

                # Poisonian error (Not used in the final S curve)
                evts_with_jpsidstar_before_trigger_error = evts_with_jpsidstar_before_trigger**0.5
                evts_with_dimudstar_after_trigger_error = evts_with_dimudstar_after_trigger**0.5
    
                # Uses ufloat to propagate error
                x = ufloat(evts_with_dimudstar_after_trigger, evts_with_dimudstar_after_trigger_error)
                y = ufloat(evts_with_jpsidstar_before_trigger, evts_with_jpsidstar_before_trigger_error)
                eff = x/y
                eff_value = eff.nominal_value
                eff_error = eff.std_dev

                print(f"The trigger efficiency is: {eff_value:.4f} +- {eff_error:.4f}")

                line = [fpt, spt, evts_with_jpsidstar_before_trigger, evts_with_dimudstar_after_trigger, f'{eff_value:.4f}', 
                        f'{eff_error:.4f}', f'{eff_error_new_up:.4f}', f'{eff_error_new_down:.4f}']
                writer.writerow(line) 

        ## Old but gold!!!
        """ DimuDstar = DimuDstar[(DimuDstar.jpsi_pt > 25) & (DimuDstar.jpsi_pt < 25.2)]

        # Number of events with Jpsi Dstar before trigger
        #evts_with_jpsidstar_before_trigger = len(DimuDstar[ak.num(DimuDstar) > 0]) # old
        evts_with_jpsidstar_before_trigger = ak.sum(ak.num(DimuDstar))

        # Trigger cut
        DimuDstar = DimuDstar[trigger_cut]

        # Number of events with Jpsi Dstar after trigger
        #evts_with_dimudstar_after_trigger  = len(DimuDstar[ak.num(DimuDstar) > 0])
        evts_with_dimudstar_after_trigger = ak.sum(ak.num(DimuDstar))

        print(evts_with_jpsidstar_before_trigger)
        print(evts_with_dimudstar_after_trigger)

        print(f'Trigger efficiency: {(evts_with_dimudstar_after_trigger / evts_with_jpsidstar_before_trigger)*100.0} %') """


        return output

    def postprocess(self, accumulator):
        return accumulator     


if __name__ == '__main__':

    # Path to the .coffea file
    path_file = config.path_file
    files = config.files

    file_list = []
    for f in files:
        file_list.append(path_file + '/' + files[f])

    pt_val_list = config.pt_val_list

    csv_file = config.csv_file
   
    # Instatiate trigger processor, the argument is just a name, it can be anything
    p = TriggerProcessor('Trigger_Efficiency')
    #calls the function
    p.process(file_list, pt_val_list, csv_file)
