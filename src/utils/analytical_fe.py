# src/utils/analytical_fe.py
# Created on 2025-04-06

import sys
import os

os.environ["OMP_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"

sys.path.insert(0, os.path.expanduser("~/pol/Projects/Codebase/Spermatogensis/backend/NucFreeEnergy"))

from methods import nucleosome_free_energy, nucleosome_groundstate, read_nucleosome_triads, GenStiffness, soft_free_energy


class NucleosomeBreath:
    def __init__(self,  
                nuc_method='crystal', 
                hang_dna_method='md'):
        
        self.genstiff_nuc = GenStiffness(method=nuc_method)
        # self.genstiff_hang = GenStiffness(method=hang_dna_method)
        self.triadfn = os.path.expanduser('~/pol/Projects/Codebase/NucleosomeMMC/State/Nucleosome.state')
        self.nuctriads = read_nucleosome_triads(self.triadfn)

        self.fn = os.path.expanduser('~/pol/Projects/Codebase/Spermatogensis/backend/NucFreeEnergy/methods/Parametrization/Kmat_nucleosome.npy')
        self.Kmat = np.load(self.fn)

    def calculate_free_energy(self, seq601:str, left_open:int, right_open:int):
        stiff, gs = self.genstiff_nuc.gen_params(seq601)

        l_open = 2*left_open
        r_open  = 2*right_open
        
        K_resc = np.copy(self.Kmat)
        F_dict = soft_free_energy(gs,stiff,l_open,r_open, self.nuctriads, K_resc)


        F601 = F_dict['F']
        F_entrop = F_dict['F_entropy']
        F_entalap = F_dict['F_enthalpy']
        F_free = F_dict['F_free']
        F_Diff = F_dict['Fdiff']

        return F601, F_entrop, F_entalap, F_free, F_Diff
    

    def calculate_free_energy_old(self, seq147:str):
        stiff, gs = self.genstiff_nuc.gen_params(seq147)
        mid_array = self.select_phosphate_bind_sites(left=0, right=13)
        F_dict = nucleosome_free_energy(gs, stiff, mid_array, self.nuctriads, use_correction=True)


        F601 = F_dict['F']
        F_entrop = F_dict['F_entropy']
        F_entalap = F_dict['F_const']
        F_free = F_dict['F_free']
        # F_Diff = F_dict['Fdiff']

        # return F601, F_entrop, F_entalap, F_free, F_Diff
        return F601, F_entrop, F_entalap, F_free


        
    def select_phosphate_bind_sites(self, left=0, right=13):

        phosphate_bind_sites = [2, 6, 14, 17, 24, 29, 34, 38, 
                                    45, 49, 55, 59, 65, 69, 76, 
                                    80, 86, 90, 96, 100, 107, 111, 
                                    116, 121, 128, 131, 139, 143]
        
        return phosphate_bind_sites[left*2:(right*2)+2]
    

    def calculate_free_energy_slide(self, long_seq: str, step: int = 1):
        """
        Slides a window of 147 nt over the given long_seq, calculates 
        free energy for each window using calculate_free_energy_old and
        returns a dictionary mapping window start index to its energies.
        
        Parameters:
            long_seq (str): The full DNA sequence.
            step (int): The step size for the sliding window. Default is 1.
            
        Returns:
            dict: A dictionary where keys are window start indices and values are
                  tuples of (F601, F_entrop, F_entalap, F_free)
        """
        energy_dict = {}
        window_size = 147
        for i in range(0, len(long_seq) - window_size + 1, step):
            window = long_seq[i:i+window_size]
            energies = self.calculate_free_energy_old(window)
            energy_dict[i] = energies
        return energy_dict


    
#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################

if __name__=='__main__':

        
    # import multiprocessing
    # multiprocessing.set_start_method('spawn', force=True)

    import time
    from Bio import SeqIO
    from src.config import path
    import concurrent.futures
    import numpy as np
    from tqdm import tqdm


    start = time.perf_counter()

    fasta_file = path.DATA_DIR / "processed/bound_top10000.fa"
    seq_records = []
    for i, record in enumerate(SeqIO.parse(fasta_file, "fasta")):
        if i >=9:
            break
        seq_records.append(record)

    for record in seq_records:
        print(f"Sequence ID: {record.id}, Length: {len(record.seq)}")
   
    # nb = NucleosomeBreath(nuc_method='hybrid')

    # long_seq = str(seq_records[0].seq).upper()
    # energy_results = nb.calculate_free_energy_slide(long_seq, step=1)
    # # Print out energies for each window start position.
    # for start_idx, energies in energy_results.items():
    #     print(f"Window starting at {start_idx}: {energies}")

    # end = time.perf_counter()
    # print(f"Finished in {round(end - start, 2)} seconds.")



    def process_sequence(record):
        nb = NucleosomeBreath(nuc_method='hybrid')
        long_seq = str(record.seq).upper()
        return record.id, nb.calculate_free_energy_slide(long_seq, step=1)





    results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=9) as executor:

    # with ThreadPoolExecutor() as executor:
        # Submit one task per sequence.
        futures = [executor.submit(process_sequence, rec) for rec in seq_records]
        total = len(futures)

        for future in tqdm(concurrent.futures.as_completed(futures), total=total, desc="Processing sequences"):
            seq_id, energy_dict = future.result()
            print(f"Appending results for sequence {seq_id}:")
            # for start_idx, energies in energy_dict.items():
            #     print(f"  Window starting at {start_idx}: {energies}")
            results.append((seq_id, energy_dict))

    print (f"Total sequences processed: {len(results)}")
    # print (f"Results: {results}")

    end = time.perf_counter()
    print(f"Finished in {round(end - start, 2)} seconds.")




