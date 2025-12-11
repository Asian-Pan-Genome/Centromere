import os
import sys
import pandas as pd
from joblib import Parallel, delayed


sys.path.append(os.getcwd())
from workflow.lib.plot_func.visualizer import CentromerePairwiseVisualizer


def process_visualization_pair(row, output_dir, mask_monomer_dir, sample_data):
    """
    Worker function that processes a single pair for visualization.
    """
    try:
        source_name = row['source']
        target_name = row['target']
        chrom = row['chrom']
        
        sample_pair_dir = os.path.join(output_dir, chrom)
        os.makedirs(sample_pair_dir, exist_ok=True) 
        
        mask_file = os.path.join(mask_monomer_dir, f"{source_name}/{source_name}.{target_name}.{chrom}.monomer.mask.gz")
        source_hor_bed = sample_data.loc[source_name, 'merged_bed']
        target_hor_bed = sample_data.loc[target_name, 'merged_bed']
        
        visualizer = CentromerePairwiseVisualizer(
            source=source_name,
            target=target_name,
            chrom=chrom,
            mutation_file=mask_file,
            source_hor_bed=source_hor_bed,
            target_hor_bed=target_hor_bed,
            output_dir=sample_pair_dir,
        )
        visualizer.run()
        return f"SUCCESS: Processed {source_name} vs {target_name} on {chrom}"
    except Exception as e:
        return f"ERROR: Failed on {row.get('source', 'N/A')} vs {row.get('target', 'N/A')}: {e}"

def main():
    filter_sample_bed = snakemake.input.filter_sample
    mask_monomer_dir = snakemake.input.mask_dir
    sample_data_bed = snakemake.params.sample_bed
    output_dir = snakemake.output.plot_dir
    threads = snakemake.threads
    
    os.makedirs(output_dir, exist_ok=True)
    
    filter_sample = pd.read_csv(filter_sample_bed, sep='\t')
    sample_data = pd.read_csv(sample_data_bed, sep='\t')
    sample_data.set_index('sample', inplace=True)
    
    # filter_sample = filter_sample[filter_sample['total_snps']<=100]
    
    tasks = (
        delayed(process_visualization_pair)(row, output_dir, mask_monomer_dir, sample_data)
        for _, row in filter_sample.iterrows()
    )

    # Set n_jobs to the number of CPU cores you want to use. -1 uses all available cores.
    # The `verbose` parameter provides progress updates.
    print(f"Starting parallel processing for {len(filter_sample)} pairs...")
    results = Parallel(n_jobs=threads, verbose=10)(tasks)

    print("\nParallel processing finished.")
    # You can inspect the 'results' list to see status messages or check for errors.
    error_count = sum(1 for r in results if r.startswith("ERROR"))
    print(f"Total tasks completed. Errors: {error_count}")
    

if __name__ == '__main__':
    
    try:
        main()
    except Exception as e:
        print(f"An error occurred: {e}")