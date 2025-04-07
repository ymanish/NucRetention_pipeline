import numpy as np
from Bio import SeqIO
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from src.config import path
from src.config.env_settings import *
from src.config.var import FreeEnergyResult
from src.utils.analytical_fe import NucleosomeBreath
from src.utils.markov_fe import calc_markov_fe
from src.utils.process_fasta import contains_non_canonical
from typing import List
from concurrent.futures import ProcessPoolExecutor


def F_weighted_markov_model(id:str, seq:str)->FreeEnergyResult:
    """
    This function calculates the weighted free energy for a given sequence.
    It uses the Markov model to compute the free energy based on the sequence.
    """
    # Placeholder for actual implementation
    P,ID = calc_markov_fe(id, 
                   seq,
                    order=3, 
                    file_l=path.PARAM_DIR / "markov_prob/MD_RoomTemp.trinucdist", 
                    file_s=path.PARAM_DIR / "markov_prob/MD_RoomTemp.dinucdist")

    # Calculate the weighted free energy
    P_min = np.min(P)
    P_rel = P/P_min
    Z_s = np.sum(P_rel)
    P_log = np.log(P)
    P_log_col = P_log.reshape(-1, 1)  # Convert to row vector
    P_rel_col = P_rel.reshape(-1, 1)  
    bF_weighted = -(P_log_col.T @ P_rel_col)[0, 0]/Z_s ### F_weighted*beta and beta is 1


    # print(P_log_col.shape)
    # print(P_rel_col.shape)
    return FreeEnergyResult(bF_weighted, ID)


def batch_calculate_free_energy(seq_records:List[SeqIO.SeqRecord]) -> List[FreeEnergyResult]:
    
    results: List[FreeEnergyResult] = []

    for record in seq_records:

        if contains_non_canonical(str(record.seq)):
            print(f"Skipping sequence {record.id} due to non-canonical characters.")
            results.append(FreeEnergyResult(-99999.0, record.id))

        elif len(record.seq) > 10000:
            print(f"Skipping sequence {record.id} because the length > 10k.")
            results.append(FreeEnergyResult(-99999.0, record.id))

        else:
            # Calculate free energy using the Markov model
            F_weight_seq = F_weighted_markov_model(record.id, str(record.seq).upper())
            results.append(F_weight_seq)
            # print(f"Free energy for sequence {F_weight_seq.id}: {F_weight_seq.energy}")


    return results



if __name__ == '__main__': 
    import time 
    start = time.perf_counter()

    # Load the sequences from the FASTA file
    # fasta_file = path.DATA_DIR / "processed/bound_top10000.fa"
    fasta_file = path.DATA_DIR / "processed/bound_seq.fa"

    seq_records = list(SeqIO.parse(fasta_file, "fasta"))

    # Define batch size
    batch_size = 100  # Adjust batch size as needed
    batches = [seq_records[i:i + batch_size] for i in range(0, len(seq_records), batch_size)]
    print(f"Total batches: {len(batches)}")

    # Use ProcessPoolExecutor for parallel processing
    all_results = []
    with ProcessPoolExecutor(max_workers=11) as executor:
        for batch_result in executor.map(batch_calculate_free_energy, batches):
            all_results.extend(batch_result)

    # Write results to a text file
    # Ensure the output directory exists
        
    output_file = path.RESULTS_DIR / "bound_free_energy_results_all.txt"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        for result in all_results:
            f.write(f"{result.id}\t{result.energy}\n")
    print(f"Results written to {output_file}")

    del all_results, seq_records, batches


     # Load the sequences from the FASTA file
    # fasta_file = path.DATA_DIR / "processed/unbound_regions.fa"
    fasta_file = path.DATA_DIR / "processed/unbound_seq.fa"

    seq_records = list(SeqIO.parse(fasta_file, "fasta"))

      # Define batch size
    batch_size = 100  # Adjust batch size as needed
    batches = [seq_records[i:i + batch_size] for i in range(0, len(seq_records), batch_size)]
    print(f"Total batches: {len(batches)}")

    # Use ProcessPoolExecutor for parallel processing
    all_results = []
    with ProcessPoolExecutor(max_workers=11) as executor:
        for batch_result in executor.map(batch_calculate_free_energy, batches):
            all_results.extend(batch_result)

    # Write results to a text file
    # Ensure the output directory exists
        
    output_file = path.RESULTS_DIR / "unbound_free_energy_results_all.txt"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    with open(output_file, "w") as f:
        for result in all_results:
            f.write(f"{result.id}\t{result.energy}\n")
    print(f"Results written to {output_file}")







    

    end = time.perf_counter()
    print(f"Finished in {round(end - start, 2)} seconds.")