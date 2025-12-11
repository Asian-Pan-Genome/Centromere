import pandas as pd
import sys


def main():
    input_file = snakemake.input[0]
    output_file = snakemake.output[0]
    params = snakemake.params.mask_filter
    
    
    stat = pd.read_csv(input_file, sep='\t')
    sort_stat = stat.sort_values(by=['snps_outside_mask'])
    sort_stat = sort_stat[~sort_stat['source'].str.contains('NA24385|HG002') & ~sort_stat['target'].str.contains('NA24385|HG002')]
    sort_stat = sort_stat[sort_stat['masked_coverage']<=params['coverage']]
    sort_stat.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    
    # --- Snakemake 执行入口 ---
    sys.stderr = sys.stdout = open(snakemake.log[0], 'w')
    try:
        main()
    except Exception as e:
        import traceback
        print(f"A critical error occurred in the main script execution: {e}")
        traceback.print_exc()
        sys.exit(1)
    