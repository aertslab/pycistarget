# TO BE UPDATED WITH CTXCORE
import io
import logging
import numpy as np
import os
import pandas as pd
import pyranges as pr
import subprocess
import sys
import ray

from typing import Optional, Union
from typing import List, Iterable, Tuple, Dict

def cluster_buster(cbust_path: str,
                 path_to_motifs: str,
                 region_sets: Union[Dict[str, pr.PyRanges], Dict[str, List]] = None,
                 path_to_genome_fasta: str = None,
                 path_to_regions_fasta: str = None,
                 n_cpu: Optional[int] = 1,
                 motifs: Optional[List[str]] = None,
                 verbose: Optional[bool] = False,
                 **kwargs):
    """
    Add motif annotation

    Parameters
    ---------
    cluster_buster_path: str
        Path to cluster buster bin.
    path_to_motifs: str, optional.
        Path to motif collection folder (in .cb format). Only required if using a shuffled background. 
    region_sets: Dict
        A dictionary of PyRanges containing region coordinates for the regions to be analyzed. Only required
        if `path_to_regions_fasta` is not provided.
    path_to_genome_fasta: str, optional.
        Path to genome fasta file. Only required if `path_to_regions_fasta` is not provided. Default: None
    path_to_regions_fasta: str, optional.
        Path to regions fasta file. Only required if `path_to_genome_fasta` is not provided. Default: None
    n_cpu: int, optional
        Number of cores to use
    motifs: List, optional
        Names of the motif files to use (from `path_to_motifs`). Default: None (All)
    verbose: bool, optional
        Whether to print progress to screen
    **kwargs:
        Additional parameters to pass to `ray.init()`
        
    References
    ---------
    Frith, Martin C., Michael C. Li, and Zhiping Weng. "Cluster-Buster: Finding dense clusters of motifs 
    in DNA sequences." Nucleic acids research 31, no. 13 (2003): 3666-3668.
    """
    # Create logger
    level    = logging.INFO
    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    log = logging.getLogger('Cluster-Buster')
    # Generate fasta file
    if path_to_regions_fasta is None:
        path_to_regions_fasta = os.path.join(outdir,'regions.fa')
    if not os.path.exists(path_to_regions_fasta):
        log.info('Getting sequences')
        pr_regions_names_dict = {key: pyranges2names(region_sets[key]) for key in region_sets.keys()}
        pr_sequence_list = [pd.DataFrame([region_sets[key], pr.get_fasta(region_sets[key], path_to_genome_fasta).tolist()], index=['Name', 'Sequence'], columns=region_sets[key]) for key in region_sets.keys()]
        seq_df = pd.concat(pr_sequence_list, axis=1)
        seq_df = seq_df.loc[:,~seq_df.columns.duplicated()]
        seq_df.T.to_csv(path_to_regions_fasta, header=False, index=False, sep='\n')
        sequence_names =  [seq[1:] for seq in seq_df.columns]
    else:
        sequence_names = get_sequence_names_from_fasta(path_to_regions_fasta)
        
    # Get motifs and sequence name
    if motifs is None:
        motifs = os.listdir(path_to_motifs)
        motifs = grep(motifs, '.cb')
    
    log.info('Scoring sequences')
    ray.init(num_cpus=n_cpu, **kwargs)
    crm_scores = ray.get([run_cluster_buster_for_motif.remote(cbust_path, path_to_regions_fasta, path_to_motifs+motifs[i], motifs[i], i, len(motifs), verbose) for i in range(len(motifs))])
    ray.shutdown()
    crm_df = pd.concat(crm_scores, axis=1, sort=False).fillna(0).T
    # Remove .cb from motifs names
    crm_df.index = [x.replace('.cb','') for x in crm_df.index.tolist()]
    log.info('Done!')
    return crm_df

# Utils
@ray.remote
def run_cluster_buster_for_motif(cluster_buster_path: str,
                                fasta_filename: str,
                                motif_filename: str,
                                motif_name: str,
                                i: int,
                                nr_motifs: int,
                                verbose: Optional[bool] = False):
    """
    Ray method to run cluster buster for one motif
    """
    # Create logger
    level    = logging.INFO
    format   = '%(asctime)s %(name)-12s %(levelname)-8s %(message)s'
    handlers = [logging.StreamHandler(stream=sys.stdout)]
    logging.basicConfig(level = level, format = format, handlers = handlers)
    log = logging.getLogger('Cluster-Buster')
    
    if verbose == True:
        log.info('Scoring motif ' + str(i) + ' out of ' + str(nr_motifs) + ' motifs')
    # Score each region in FASTA file with Cluster-Buster
    # for motif and get top CRM score for each region.
    clusterbuster_command = [cluster_buster_path,
                             '-f', '4',
                             '-c', '0.0',
                             '-r', '10000',
                             '-t', '1',
                             '-l', #Mask repeats
                             motif_filename,
                             fasta_filename]

    try:
        pid = subprocess.Popen(args=clusterbuster_command,
                               bufsize=1,
                               executable=None,
                               stdin=None,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               preexec_fn=None,
                               close_fds=False,
                               shell=False,
                               cwd=None,
                               env=None,
                               universal_newlines=False,
                               startupinfo=None,
                               creationflags=0)
        stdout_data, stderr_data = pid.communicate()
    except OSError as msg:
        print("\nExecution error for: '" + ' '.join(clusterbuster_command) + "': " + str(msg),
              file=sys.stderr)
        sys.exit(1)

    if pid.returncode != 0:
        print("\nError: Non-zero exit status for: " + ' '.join(clusterbuster_command) + "'",
              file=sys.stderr)
        sys.exit(1)

    crm_scores_df = pd.read_csv(
        filepath_or_buffer=io.BytesIO(stdout_data),
        sep='\t',
        header=0,
        names=['motifs', 'crm_score'],
        index_col='motifs',
        usecols=['motifs','crm_score'],
        dtype={'crm_score': np.float32},
        engine='c'
    )
    
    crm_scores_df.columns=[motif_name]
    return crm_scores_df


# Utils functions for Cluster-buster
def get_sequence_names_from_fasta(fasta_filename: str):
    """
    Retrieve sequence names from fasta
    """
    sequence_names_list = list()
    sequence_names_set = set()
    duplicated_sequences = False

    with open(fasta_filename, 'r') as fh:
        for line in fh:
            if line.startswith('>'):
                # Get sequence name by getting everything after '>' up till the first whitespace.
                sequence_name = line[1:].split(maxsplit=1)[0]

                # Check if all sequence names only appear once.
                if sequence_name in sequence_names_set:
                    print(
                        'Error: Sequence name "{0:s}" is not unique in FASTA file "{1:s}".'.format(
                            sequence_name,
                            fasta_filename
                        ),
                        file=sys.stderr
                    )
                    duplicated_sequences = True

                sequence_names_list.append(sequence_name)
                sequence_names_set.add(sequence_name)

    if duplicated_sequences:
        sys.exit(1)

    return sequence_names_list

def pyranges2names(regions: pr.PyRanges):
    """
    Convert pyranges to sequence name (fasta format)
    """
    return ['>'+str(chrom) + ":" + str(start) + '-' + str(end) for chrom, start, end in zip(list(regions.Chromosome), list(regions.Start), list(regions.End))]

def grep(l: List,
         s: str):
    """
    Helper for grep
    """
    return [i for i in l if s in i]
