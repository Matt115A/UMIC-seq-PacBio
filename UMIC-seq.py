#Main script for UMI-linked consensus sequencing
#Author: Paul Jannis Zurek, pjz26@cam.ac.uk
#21/08/2019 v1.0
#07/02/2020 v1.1
## Added quick stop to the clustering, to make it more time efficient. It'll stop clustering when essentially only outliers are left (low average clustersize).
#09/04/2020 v1.1.1
## Minor changes to aesthetics
#11/02/2021 v1.1.2
## Fixed an issue where incomplete UMIs could be extracted from truncated reads


from Bio import SeqIO
from Bio.Seq import Seq
from skbio.alignment import local_pairwise_align_nucleotide
from skbio import DNA
import multiprocessing
import itertools
import os
import numpy as np
from matplotlib import pyplot as plt
import argparse
import gzip
import time
import psutil
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor

parser = argparse.ArgumentParser(description="""Main script for UMI-linked consensus sequencing.
                                 Author: Paul Zurek (pjz26@cam.ac.uk).
                                 Version 1.1.1""")
#parser.add_argument('mode', help='Select mode', choices=('UMIextract', 'clustertest', 'clusterfull'))
parser.add_argument('-T', '--threads', type=int, default=0, help='Number of threads to execute in parallel. Defaults to CPU count.')
parser.add_argument('-v', '--version', action='version', version='1.1.2')

subparsers = parser.add_subparsers(help='Select mode', dest='mode')

#Arguments for UMIextract
extract_parser = subparsers.add_parser('UMIextract', help='Extraction of UMIs from reads.')
extract_parser.add_argument('-i', '--input', help='Provide basecalled reads in fastq format.', required=True)
extract_parser.add_argument('-o', '--output', help='Specify the name of the UMI output fasta file.', required=True)
extract_parser.add_argument('--probe', help='A short sequence (eg 50 bp) adjacent to the UMI in fasta format.', required=True)
extract_parser.add_argument('--umi_loc', help='Location of UMI in reference to the probe. Upstream (up) or downstream (down).', choices=('up', 'down'), required=True)
extract_parser.add_argument('--umi_len', help='Length of the UMI to be extracted.', type=int, required=True)
extract_parser.add_argument('--min_probe_score', help='Minimal alignment score of probe for processing. Defaults to length of probe sequence.', type=int, default=0)

#Arguments for clustertest
clustest_parser = subparsers.add_parser('clustertest', help='Test suitable alignment thresholds for clustering.')
clustest_parser.add_argument('-i', '--input', help='Fasta file of extracted UMIs.', required=True)
clustest_parser.add_argument('-o', '--output', help='Prefix for output files.', required=True)
clustest_parser.add_argument('--steps', help='Accepts left border, right border and step width for sampled thresholds. Defaults to 20 70 10 (samples thresholds 20 30 40 .. 70).', nargs=3, type=int, default=[20,70,10])
clustest_parser.add_argument('--samplesize', help='Number of clusters to be sampled for threshold approximation. Defaults to 25.', type=int, default=25)

#Arguments for clusterfull
fullclus_parser = subparsers.add_parser('clusterfull', help='Full clustering of UMIs.')
fullclus_parser.add_argument('-i', '--input', help='Fasta file of extracted UMIs.', required=True)
fullclus_parser.add_argument('-o', '--output', help='Folder name for output files.', required=True)
fullclus_parser.add_argument('--reads', help='Fastq file of basecalled reads.', required=True)
fullclus_parser.add_argument('--aln_thresh', type=int, help='Alignment threshold for clustering. UMIs with alignment scores higher than aln_thresh will be clustered.', required=True)
fullclus_parser.add_argument('--size_thresh', type=int, help='Minimal size a cluster can have to be written to file.', required=True)
fullclus_parser.add_argument('--stop_thresh', type=int, default=5, required=False, help='Defaults to 5. Stops clustering if the average cluster size is smaller than this threshold. Essentially speeds up the clustering by dropping outliers. Set the threshold to 0 if you do not want the program to quit early!')
fullclus_parser.add_argument('--stop_window', type=int, default=20, required=False, help='Defaults to 20. Sets the number of clusters to be used to calculate average cluster size.')

#Arguments for clusterfull_fast (CD-HIT based)
fastclus_parser = subparsers.add_parser('clusterfull_fast', help='Fast clustering of UMIs using CD-HIT (1000x faster than clusterfull).')
fastclus_parser.add_argument('-i', '--input', help='Fasta file of extracted UMIs.', required=True)
fastclus_parser.add_argument('-o', '--output', help='Folder name for output files.', required=True)
fastclus_parser.add_argument('--reads', help='Fastq file of basecalled reads.', required=True)
fastclus_parser.add_argument('--identity', type=float, default=0.90, help='Sequence identity threshold for clustering (0-1). Default: 0.90 (90%% identity)')
fastclus_parser.add_argument('--size_thresh', type=int, help='Minimal size a cluster can have to be written to file.', required=True)

# Module-level variables that will be set when run as main
args = None
mode = None
threads = 0

#Parse arguments only when run as main script
if __name__ == '__main__':
    args = parser.parse_args()
    mode = args.mode
    threads = args.threads
    if threads == 0:
        threads = multiprocessing.cpu_count()



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#HELPER FUNCTIONS FOR COMPRESSED FILE HANDLING

def open_file_stream(filename):
    """Open a file stream, handling both compressed (.gz) and uncompressed files."""
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    else:
        return open(filename, 'r')

def parse_fastq_stream(file_handle):
    """Parse FASTQ records from a file handle (compressed or uncompressed)."""
    return SeqIO.parse(file_handle, "fastq")

def get_file_size(filename):
    """Get file size for progress estimation."""
    if filename.endswith('.gz'):
        # For compressed files, we can't easily get the uncompressed size
        # So we'll estimate based on typical compression ratios
        compressed_size = os.path.getsize(filename)
        # Typical FASTQ compression ratio is 3-4x
        return compressed_size * 3.5
    else:
        return os.path.getsize(filename)

def print_progress(count, total_estimated, start_time, success, bad_aln, short_aln):
    """Print detailed progress information."""
    elapsed_time = time.time() - start_time
    
    if count > 0:
        # Calculate rates
        seq_per_sec = count / elapsed_time
        success_rate = (success / count) * 100
        
        # Estimate remaining time
        if total_estimated > 0:
            progress_percent = (count / total_estimated) * 100
            remaining_seqs = total_estimated - count
            if seq_per_sec > 0:
                eta_seconds = remaining_seqs / seq_per_sec
                eta_minutes = eta_seconds / 60
                eta_hours = eta_minutes / 60
                
                if eta_hours >= 1:
                    eta_str = f"{eta_hours:.1f}h"
                elif eta_minutes >= 1:
                    eta_str = f"{eta_minutes:.1f}m"
                else:
                    eta_str = f"{eta_seconds:.0f}s"
            else:
                eta_str = "unknown"
        else:
            progress_percent = 0
            eta_str = "unknown"
        
        # Get memory usage
        process = psutil.Process()
        memory_mb = process.memory_info().rss / 1024 / 1024
        
        # Format elapsed time
        if elapsed_time >= 3600:
            elapsed_str = f"{elapsed_time/3600:.1f}h"
        elif elapsed_time >= 60:
            elapsed_str = f"{elapsed_time/60:.1f}m"
        else:
            elapsed_str = f"{elapsed_time:.0f}s"
        
        progress_line = (f"Processed: {count:,} sequences | "
                        f"Rate: {seq_per_sec:.1f} seq/s | "
                        f"Success: {success:,} ({success_rate:.1f}%) | "
                        f"Progress: {progress_percent:.1f}% | "
                        f"ETA: {eta_str} | "
                        f"Memory: {memory_mb:.0f}MB | "
                        f"Elapsed: {elapsed_str}")
        print(f"\r{progress_line}", end='', flush=True)

class AlignmentResult:
    """Simple class to mimic StripedSmithWaterman result interface."""
    def __init__(self, score, target_sequence, target_begin, target_end_optimal):
        self.optimal_alignment_score = score
        self.target_sequence = target_sequence
        self.target_begin = target_begin
        self.target_end_optimal = target_end_optimal

def process_sequence_batch(batch_data):
    """Process a batch of sequences for multiprocessing."""
    batch_records, probe_fwd, probe_rev, probe_minaln_score, umi_len, umi_loc = batch_data
    results = []
    
    for record in batch_records:
        rec_seq = str(record.seq)
        alnF = align_sequence_fast(probe_fwd, rec_seq)
        alnR = align_sequence_fast(probe_rev, rec_seq)
        scoreF = alnF.optimal_alignment_score
        scoreR = alnR.optimal_alignment_score
        
        if scoreF > probe_minaln_score or scoreR > probe_minaln_score:
            if scoreF > scoreR:
                # Forward orientation
                aln = alnF
                score = scoreF
            else:
                # Reverse orientation
                aln = alnR
                score = scoreR
            
            # Extract UMI based on location
            if umi_loc == 'up':
                umi_start = aln.target_begin - umi_len
                umi_end = aln.target_begin
            else:  # umi_loc == 'down'
                umi_start = aln.target_end_optimal
                umi_end = aln.target_end_optimal + umi_len
            
            if umi_start >= 0 and umi_end <= len(rec_seq):
                umi_seq = rec_seq[umi_start:umi_end]
                if len(umi_seq) == umi_len:
                    results.append((record.id, umi_seq, score))
    
    return results

def align_sequence_fast(query_seq, target_seq):
    """Ultra-fast alignment using regex patterns."""
    try:
        query_upper = query_seq.upper()
        target_upper = target_seq.upper()
        
        # Try exact match first (fastest)
        pos = target_upper.find(query_upper)
        if pos != -1:
            return AlignmentResult(
                score=len(query_seq),
                target_sequence=target_seq,
                target_begin=pos,
                target_end_optimal=pos + len(query_seq)
            )
        
        # Try fuzzy matching with regex (allow 1-2 mismatches)
        query_len = len(query_seq)
        if query_len > 20:  # Only for longer sequences
            # Create regex pattern allowing 1-2 mismatches
            pattern_parts = []
            for i, char in enumerate(query_upper):
                if i < query_len - 2:  # Allow mismatches except in last 2 positions
                    pattern_parts.append(f"[{char}]")
                else:
                    pattern_parts.append(char)
            
            pattern = ''.join(pattern_parts)
            regex = re.compile(pattern)
            
            for match in regex.finditer(target_upper):
                # Count actual matches
                matches = sum(1 for i, char in enumerate(query_upper) 
                            if target_upper[match.start() + i] == char)
                
                if matches >= query_len * 0.8:  # At least 80% match
                    return AlignmentResult(
                        score=matches * 2,
                        target_sequence=target_seq,
                        target_begin=match.start(),
                        target_end_optimal=match.end()
                    )
        
        return AlignmentResult(0, target_seq, 0, 0)
            
    except Exception:
        return AlignmentResult(0, target_seq, 0, 0)

def align_sequence(query_seq, target_seq):
    """
    Alignment function for clustering (uses skbio for accurate scoring).
    This is used for UMI-to-UMI comparisons where accuracy is more important than speed.
    """
    alignment, score, _ = local_pairwise_align_nucleotide(
        DNA(query_seq), 
        DNA(target_seq)
    )
    # Return an AlignmentResult-like object for compatibility
    return AlignmentResult(
        score=score,
        target_sequence=target_seq,
        target_begin=0,
        target_end_optimal=len(target_seq)
    )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#UMI EXTRACTION
#Extracts UMI in correct orientation
#To optimize: Multiprocessing? Should be fast enough like this...

def prnt_aln(aln):
    print('Score: %d' % aln.optimal_alignment_score)
    print(aln.aligned_target_sequence)  #READ
    print(aln.aligned_query_sequence)   #Search
    s = aln.target_begin
    d = aln.target_end_optimal
    l = len(aln.target_sequence)
    print('Left Border: %d\t Right Border: %d \n' % (s, l-d))

def extract_left(aln):
    umi_end = aln.target_begin
    umi_begin = umi_end - umi_len
    return umi_begin, umi_end

def extract_right(aln):
    umi_begin = aln.target_end_optimal + 1
    umi_end = umi_begin + umi_len
    return umi_begin, umi_end


if __name__ == '__main__' and mode == 'UMIextract':
    #Set up variables
    input_file = args.input
    probe_file = args.probe
    output_file = args.output
    umi_loc = args.umi_loc
    umi_len = args.umi_len

    proberec = SeqIO.read(probe_file, "fasta")
    probe_minaln_score = args.min_probe_score
    if probe_minaln_score == 0:
        probe_minaln_score = len(proberec.seq)

    probe_fwd = str(proberec.seq)
    probe_rev = str(proberec.reverse_complement().seq)

    count, bad_aln, short_aln, success = 0, 0, 0, 0
    
    # Open output file for streaming writes
    output_handle = open(output_file, 'w')
    
    # Estimate total sequences for progress tracking
    print("Estimating file size for progress tracking...")
    file_size = get_file_size(input_file)
    # Rough estimate: ~200 bytes per sequence (including quality scores)
    total_estimated = int(file_size / 200)
    print(f"Estimated {total_estimated:,} sequences in file")
    
    start_time = time.time()
    print("Starting UMI extraction...")
    print("Progress will be shown for every sequence...")
    
    # Use streaming approach for compressed files
    try:
        with open_file_stream(input_file) as file_handle:
            for record in parse_fastq_stream(file_handle):
                rec_seq = str(record.seq)
                alnF = align_sequence_fast(probe_fwd, rec_seq)
                alnR = align_sequence_fast(probe_rev, rec_seq)
                scoreF = alnF.optimal_alignment_score
                scoreR = alnR.optimal_alignment_score
                #Check basic alignment score
                if scoreF > probe_minaln_score or scoreR > probe_minaln_score:  
                    if scoreF > scoreR: #Target in fwd orientation
                        #Get umi location
                        if umi_loc == 'down':
                            umi_begin, umi_end = extract_right(alnF)   #FWD of downstream is right
                        elif umi_loc == 'up':
                            umi_begin, umi_end = extract_left(alnF)   #FWD of upstream is left
                        #append to UMI to record list
                        if umi_end < len(alnF.target_sequence) and umi_begin > 0: #UMI could be out of bounds
                            umi = alnF.target_sequence[umi_begin:umi_end]
                            # Write directly to output file instead of storing in memory
                            output_handle.write(f">{record.id}\n{umi}\n")
                            success += 1
                        else:
                            short_aln += 1
                    else: #Target in rev orientation
                        #Get umi location
                        if umi_loc == 'down':
                            umi_begin, umi_end = extract_left(alnR)   #REV of downstream is left
                        elif umi_loc == 'up':
                            umi_begin, umi_end = extract_right(alnR)   #REV of upstream is right
                        #append to UMI to record list
                        if umi_begin > 0 and umi_end < len(alnR.target_sequence): #UMI could be out of bounds
                            umiR = alnR.target_sequence[umi_begin:umi_end]
                            # Write directly to output file instead of storing in memory
                            umi_revcomp = str(Seq(umiR).reverse_complement())
                            output_handle.write(f">{record.id}\n{umi_revcomp}\n")
                            success += 1
                        else:
                            short_aln += 1
                else:
                    bad_aln += 1
                count += 1

                # Print progress for every sequence
                elapsed_time = time.time() - start_time
                seq_per_sec = count / elapsed_time if elapsed_time > 0 else 0
                success_rate = (success / count) * 100 if count > 0 else 0
                
                print(f"\r[{time.strftime('%H:%M:%S')}] Seq: {count:,} | Rate: {seq_per_sec:.1f}/s | Success: {success:,} ({success_rate:.1f}%) | Mem: {psutil.Process().memory_info().rss / 1024 / 1024:.0f}MB", end='', flush=True)
    finally:
        # Ensure output file is closed even if there's an error
        output_handle.close()
    
    # Final progress update
    print_progress(count, total_estimated, start_time, success, bad_aln, short_aln)
    print()  # New line after progress
    print('\nUMIs extracted: %d' % success)
    print("Discarded: %.2f%%:" % ((bad_aln+short_aln)/count * 100))
    print("Bad alignment: %d" % bad_aln)
    print("Incomplete UMI: %d" % short_aln)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#GENERAL CLLUSTERING FUNCTIONS

#Calculates SW alignment score between sequences i and j
def aln_score(query_seq, umi):
    aln = align_sequence(query_seq, umi)
    score = aln.optimal_alignment_score
    return score

def simplesim_cluster(umis, thresh, max_clusters=0, save_scores=False, clussize_thresh=0, clussize_window=20):
    N_umis = len(umis)
    #Generate query list for alignment
    query_lst = [str(umi.seq) for umi in umis]

    clusnr = 0
    lab = [None for i in range(N_umis)]
    seq_index = list(range(N_umis))
    remaining_umis = [str(umi.seq) for umi in umis]
    clussize = []

    score_lst_lst = []
    pool = multiprocessing.Pool(processes=threads)
    #Goes through seq_index and only compares those indices and removes from there!
    while len(seq_index) > 0:
        #Calculate a list of similar umis to the first one left in seq_index 
        # Pass the actual query sequence instead of index
        query_seq = query_lst[seq_index[0]]
        score_lst = pool.starmap(aln_score, zip(itertools.repeat(query_seq), remaining_umis))
        if save_scores:
            score_lst_lst.append(score_lst)

        #Only propagate what is similar!
        #Go through the index list with similar sequences
        sim_lst = []    
        for i in reversed(range(len(seq_index))):
            if score_lst[i] > thresh:
                sim_lst.append(seq_index[i])
                lab[seq_index[i]] = clusnr
                del seq_index[i]
                del remaining_umis[i]
        
        #Threshold must be smaller than similarity to itself for the simplesim clusterer to work!! (could be a problem for non-relative metric, just the aln score)
        assert len(sim_lst) >= 1,  'Threshold too high'

        clussize.append(len(sim_lst))
        clusnr += 1
        
        if clussize_thresh > 0:
            #Check last few clusters to see if average is small
            if len(clussize) >= clussize_window:
                av_size = sum(clussize[-clussize_window:]) / clussize_window
                if av_size < clussize_thresh:
                    print("Ended early due to small average cluster size.")
                    break
                print("Cluster %d: %d entries.    \tSeq remaining: %d    \t Average size: %.2f" % (clusnr, len(sim_lst), len(seq_index), av_size), end='\r')
            else:
                print("Cluster %d: %d entries.    \tSeq remaining: %d    " % (clusnr, len(sim_lst), len(seq_index)), end='\r')
        else:
            print("Cluster %d: %d entries.    \tSeq remaining: %d    " % (clusnr, len(sim_lst), len(seq_index)), end='\r')
        
        if (clusnr >= max_clusters) and (max_clusters > 0):  #Sampling mode!
            break

    pool.close()
    if save_scores:
        return clusnr, clussize, lab, score_lst_lst
    return clusnr, clussize, lab



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#CLUSTER THRESHOLD APPROXIMATION

#Analyse similarity within one cluster (used in averages_withincluster)
def within_cluster_analysis(clus_id, labels, umis, maxiter=200):
    #Get a list of the actual read ID for reads belonging to the given cluster
    read_ids = list(range(len(labels)))
    cluster_members = [i for i in read_ids if labels[i] == clus_id]
    if len(cluster_members) <= 1:  #if the threshold was so high that no clusters are formed I'll return a 0 as score for now.
        return 0.0, len(cluster_members)
    else:
        #Calculate alignment score for all combinations within this cluster
        #but not all all combinations, maximally for the first 100 sequences in the cluster, otherwise the calculations are too long
        if len(cluster_members) < maxiter:
            maxiter = len(cluster_members)
        count = 0
        total_score = 0
        for i in range(maxiter):
            #Align: Doesn't need to be multicore, as there are only few alignments to be done (nr = maxiter)
            #Main calculation burden is in the clustering, which already is multiprocessed
            query_seq = str(umis[cluster_members[i]].seq)
            for j in range(i+1, maxiter):
                aln = align_sequence(query_seq, str(umis[cluster_members[j]].seq))
                score = aln.optimal_alignment_score
                total_score += score
                count += 1
        #print("Score: %.2f" % (total_score/count))
        return total_score/count, len(cluster_members)

#Calculates the within cluster similarity average for the first n clusters
def averages_withincluster(n_samples, labels, umis):
    similarity_lst = []
    length_lst = []
    for i in range(n_samples):
        sim, l = within_cluster_analysis(i, labels, umis)
        similarity_lst.append(sim)
        length_lst.append(l)
    return np.average(similarity_lst), np.average(length_lst), np.median(length_lst)

#Plots similarity histogram
def similarity_histogram(score_lst_lst, outname):
    #Remove self-alignment
    for i in range(len(score_lst_lst)):
        maxind = score_lst_lst[i].index(max(score_lst_lst[i]))
        del score_lst_lst[i][maxind]
    #Plotting
    binwidth = 2
    f, axarr = plt.subplots(5, 5, sharex=True, sharey=True, figsize=(12,8))
    f.subplots_adjust(hspace = 0.2, wspace = 0.1)
    f.suptitle("Individual similarity histograms", fontsize=12)
    f.text(0.5, 0.04, 'Alignment score', ha='center')
    f.text(0.04, 0.5, 'Count', va='center', rotation='vertical')
    i = 0
    for a in range(5):
        for b in range(5):
            axarr[a,b].hist(score_lst_lst[i], bins=np.arange(min(score_lst_lst[i]), max(score_lst_lst[i])+binwidth, binwidth))
            axarr[a,b].set_yscale('log')
            i += 1
    plt.savefig(outname + "_similarityscores_hist.pdf", bbox_inches='tight')


#Function performing a threshold approximation for clustering
def threshold_approx(umis, ssize, left, right, step, outname):
    threshholds = list(range(left, right+step, step))
    similarities_lst, sizes_lst = [], []
    print(f"Starting threshold approximation ({threads} threads)")
    for thresh in threshholds:
        _, _, read_labels, scores_lst_lst = simplesim_cluster(umis, thresh, max_clusters=ssize, save_scores=True)
        average_sim, _, median_size = averages_withincluster(ssize, read_labels, umis)
        similarities_lst.append(average_sim)
        sizes_lst.append(median_size)
        print(f"Threshold {thresh}: Average similarity {average_sim:.2f} and median size {median_size:.0f}                ")

    if ssize >= 25:
        #Plot a similarity histogram
        similarity_histogram(scores_lst_lst, outname)

    #Plot threshold approximation
    fig, ax1 = plt.subplots(figsize=(6,6))
    ax2 = ax1.twinx()
    ax1.plot(threshholds, similarities_lst, color="k", marker="o", linestyle="-", linewidth=2)
    ax2.plot(threshholds, sizes_lst, color="b", marker="o", linestyle="-", linewidth=2)
    #ax1.errorbar(threshholds, similarities_lst, yerr=sim_std_lst, color="k", marker="o", linestyle="-", linewidth=2, capsize=3, elinewidth=0.5)
    #ax2.errorbar(threshholds, sizes_lst, yerr=size_std_lst, color="b", marker="o", linestyle="-", linewidth=2, capsize=3, elinewidth=0.5)    ax1.set_xlabel("Alignment score threshold", fontsize=14)
    ax1.set_ylabel("Average cluster similarity", color="k", fontsize=14)
    ax2.set_ylabel("Median cluster size", color="b", fontsize=14)
    ax2.set_ylim(0, max(sizes_lst)+10)
    plt.savefig(outname + "_thresholdapproximation.pdf", bbox_inches='tight')


if __name__ == '__main__' and mode == 'clustertest':
    #Setup
    input_file = args.input
    sample_size = args.samplesize
    output_name = args.output
    left, right, step = args.steps  #Probably a good idea to have sensible default values based on the UMI length. Maybe a range from 0.5 x len to 1.5 x len.

    #Load umis and start approximation
    with open_file_stream(input_file) as file_handle:
        umis = list(SeqIO.parse(file_handle, "fasta"))
    threshold_approx(umis, sample_size, left, right, step, output_name)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FULL CLUSTERING

##Change max_clusters to eg 10 for testing, leave out/set to 0 for production
def cluster_sequences(umis, reads_file_path, aln_thresh, size_thresh, max_clusters=0, clussize_thresh=0, clussize_window=20):
    print("Beginning clustering...")
    clus_N, cluster_sizes, labels = simplesim_cluster(umis, aln_thresh, max_clusters=max_clusters, clussize_thresh=clussize_thresh, clussize_window=clussize_window) 

    #Printing some metrics
    print("\nClustering done!")
    print(f"Number of clusters: {clus_N}")
    #print(f"Number of clusters (TEST): {max(labels) + 1}")
    large_clus = [1 for size in cluster_sizes if size >= size_thresh]
    print(f"Clusters with >= {size_thresh} members: {sum(large_clus)}")
    seq_large_clus = sum([size for size in cluster_sizes if size >= size_thresh])
    rel_seq_large_clus = seq_large_clus / len(umis) * 100
    print(f"Total number of sequences in clusters with >= {size_thresh} members: {seq_large_clus} ({rel_seq_large_clus:.2f}%)")

    #Plot clustersizes histogram
    binwidth = 1
    """
    plt.figure()
    plt.hist(cluster_sizes, bins=np.arange(min(cluster_sizes), max(cluster_sizes)+binwidth, binwidth))
    plt.xlabel("Cluster size")
    plt.ylabel("Number of clusters")
    plt.savefig(output_folder + "clustersizes_clusters.pdf", bbox_inches='tight')
    """
    fig = plt.figure()
    hist_clussize = np.repeat(cluster_sizes, cluster_sizes)
    med = np.median(hist_clussize)
    plt.hist(hist_clussize, bins=np.arange(min(hist_clussize), max(hist_clussize)+binwidth, binwidth))
    textstr = f"Median: {med}"
    print(f"Median number of sequences per cluster: {med}")
    plt.text(0.85, 0.85, textstr, transform=fig.transFigure, ha='right')
    plt.xlabel("Cluster size")
    plt.ylabel("Number of sequences")
    plt.savefig(output_folder + "_clustersizes_sequences.pdf", bbox_inches='tight')
    print("Clustersize distributions plotted")

    count = 0
    for i in range(clus_N):
        if cluster_sizes[i] >= size_thresh:
            #Get reads according to UMI label using streaming
            count += 1
            clus_member_ids = [rec.id for rec, lab in zip(umis, labels) if lab == i]
            
            # Stream through reads file to find matching sequences
            clus_members = []
            with open_file_stream(reads_file_path) as file_handle:
                for record in parse_fastq_stream(file_handle):
                    if record.id in clus_member_ids:
                        clus_members.append(record)
                        if len(clus_members) == len(clus_member_ids):
                            break  # Found all members
            
            fname = output_folder + "/cluster_" + str(count) + ".fasta"
            SeqIO.write(clus_members, fname, "fasta")
    print("Cluster files written")


if __name__ == '__main__' and mode == 'clusterfull':
    #Setup
    input_UMIfile = args.input   #'./2_UMIextract/UMIexRS_BC03.fasta'
    input_READSfile = args.reads  #'./1_demultiplex/UMIrandomsample_BC03.fastq'
    aln_thresh = args.aln_thresh  #50
    size_thresh = args.size_thresh  #5  
    output_folder = args.output    #'./4_clustering/BC03_RS-clusters/'
    clussize_thresh = args.stop_thresh  #Defaults to 5. Set to 0 if you don't want it to stop clustering!
    clussize_window = args.stop_window  #defaults to 20. 
    
    #Make output directory
    if not os.path.exists(output_folder):
        os.makedirs(output_folder) 

    # Use streaming approach for compressed files
    # SeqIO.index doesn't work with gzipped files, so we'll use streaming
    print(f"Using streaming approach for reads file: {input_READSfile}")
    
    # Handle compressed files for UMIs
    with open_file_stream(input_UMIfile) as file_handle:
        umis = list(SeqIO.parse(file_handle, "fasta"))

    #NOW ENDS EARLY IF CLUSSIZE_WINDOW AND THRESH ARE SET!
    #Calculates average clusterisze over the last X clusters and stops clustering if this is below thresh!
    cluster_sequences(umis, input_READSfile, aln_thresh, size_thresh, clussize_thresh=clussize_thresh, clussize_window=clussize_window)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#FAST CLUSTERING USING CD-HIT

def parse_cdhit_clusters(clstr_file):
    """Parse CD-HIT .clstr file and return cluster assignments."""
    clusters = {}  # cluster_id -> list of sequence IDs
    current_cluster = -1
    
    with open(clstr_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>Cluster'):
                current_cluster = int(line.split()[1])
                clusters[current_cluster] = []
            elif line:
                # Parse sequence ID from line like: "0	52nt, >sequence_id... at 100%"
                seq_id = line.split('>')[1].split('...')[0]
                clusters[current_cluster].append(seq_id)
    
    return clusters

def run_cdhit_clustering(input_umi_file, output_folder, identity_threshold, size_thresh, reads_file_path):
    """Run CD-HIT clustering on UMI sequences."""
    print(f"Running CD-HIT clustering with {identity_threshold:.2f} identity threshold...", flush=True)
    
    # Run CD-HIT
    cdhit_output = os.path.join(output_folder, "cdhit_output")
    cdhit_cmd = [
        'cd-hit-est',
        '-i', input_umi_file,
        '-o', cdhit_output,
        '-c', str(identity_threshold),
        '-n', '8',  # word size for 90% identity
        '-d', '0',  # full sequence ID in output
        '-M', '0',  # unlimited memory
        '-T', str(threads)
    ]
    
    print(f"Running: {' '.join(cdhit_cmd)}")
    result = subprocess.run(cdhit_cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"CD-HIT failed: {result.stderr}")
        return False
    
    print("CD-HIT clustering completed!")
    
    # Parse cluster file
    clstr_file = cdhit_output + '.clstr'
    clusters = parse_cdhit_clusters(clstr_file)
    
    print(f"Total clusters: {len(clusters)}")
    
    # Filter by size and write cluster files
    cluster_sizes = [len(members) for members in clusters.values()]
    large_clusters = {cid: members for cid, members in clusters.items() if len(members) >= size_thresh}
    
    print(f"Clusters with >= {size_thresh} members: {len(large_clusters)}")
    
    if len(large_clusters) == 0:
        print("No clusters meet size threshold!")
        return True
    
    # SINGLE-PASS extraction: Stream once, write to temp batch files, then split
    print(f"Extracting cluster members via single-pass streaming...", flush=True)
    print(f"Processing {len(large_clusters):,} clusters with {sum(len(m) for m in large_clusters.values()):,} total sequences", flush=True)
    
    # Build read-to-cluster mapping
    print(f"Building read-to-cluster mapping...", flush=True)
    read_to_cluster_file = {}
    cluster_to_file = {}
    
    # Assign each cluster to a batch file (max 1000 clusters per batch file to manage memory)
    BATCH_SIZE = 1000
    cluster_list = list(large_clusters.items())
    num_batch_files = (len(cluster_list) + BATCH_SIZE - 1) // BATCH_SIZE
    
    for idx, (cluster_id, member_ids) in enumerate(cluster_list):
        batch_file_idx = idx // BATCH_SIZE
        cluster_to_file[cluster_id] = batch_file_idx
        for read_id in member_ids:
            read_to_cluster_file[read_id] = (cluster_id, batch_file_idx)
    
    print(f"Will use {num_batch_files} temporary batch files", flush=True)
    
    # Create batch file handles
    import tempfile
    temp_dir = tempfile.mkdtemp(dir=output_folder)
    batch_files = []
    for i in range(num_batch_files):
        fname = os.path.join(temp_dir, f"batch_{i}.txt")
        batch_files.append(open(fname, 'w'))
    
    # SINGLE PASS: Stream through reads and write to batch files
    print(f"Streaming through reads file (single pass)...", flush=True)
    total_found = 0
    total_needed = len(read_to_cluster_file)
    
    try:
        with open_file_stream(reads_file_path) as file_handle:
            for i, record in enumerate(parse_fastq_stream(file_handle)):
                if record.id in read_to_cluster_file:
                    cluster_id, batch_idx = read_to_cluster_file[record.id]
                    # Write to batch file: cluster_id|seq_id|sequence
                    batch_files[batch_idx].write(f"{cluster_id}|{record.id}|{str(record.seq)}\n")
                    total_found += 1
                    
                    if total_found % 50000 == 0:
                        print(f"Progress: {total_found:,}/{total_needed:,} sequences ({100*total_found/total_needed:.1f}%)", flush=True)
                    
                    if total_found == total_needed:
                        break
                
                if (i + 1) % 2000000 == 0:
                    print(f"Scanned {(i+1):,} reads, found {total_found:,}/{total_needed:,}", flush=True)
    finally:
        # Close batch files
        for f in batch_files:
            f.close()
    
    print(f"✓ Single pass complete: found {total_found:,}/{total_needed:,} sequences", flush=True)
    
    # Now split batch files into individual cluster files
    print(f"Writing individual cluster files...", flush=True)
    written = 0
    
    for batch_idx in range(num_batch_files):
        batch_fname = os.path.join(temp_dir, f"batch_{batch_idx}.txt")
        
        # Read batch file and group by cluster
        cluster_seqs = {}
        with open(batch_fname, 'r') as f:
            for line in f:
                cluster_id, seq_id, sequence = line.strip().split('|', 2)
                if cluster_id not in cluster_seqs:
                    cluster_seqs[cluster_id] = []
                cluster_seqs[cluster_id].append((seq_id, sequence))
        
        # Write cluster files
        for cluster_id, seqs in cluster_seqs.items():
            written += 1
            fname = os.path.join(output_folder, f"cluster_{written}.fasta")
            with open(fname, 'w') as f:
                for seq_id, sequence in seqs:
                    f.write(f">{seq_id}\n{sequence}\n")
        
        # Clean up batch file
        os.remove(batch_fname)
        
        if (batch_idx + 1) % 10 == 0:
            print(f"Processed {batch_idx+1}/{num_batch_files} batch files, written {written} clusters", flush=True)
    
    # Clean up temp directory
    os.rmdir(temp_dir)
    
    print(f"\n✓ All cluster files written: {written:,}", flush=True)
    
    # Plot cluster size distribution
    fig = plt.figure()
    hist_clussize = []
    for size in cluster_sizes:
        hist_clussize.extend([size] * size)
    
    if hist_clussize:
        med = np.median(hist_clussize)
        plt.hist(hist_clussize, bins=np.arange(1, max(hist_clussize)+2, 1))
        textstr = f"Median: {med}"
        print(f"Median number of sequences per cluster: {med}")
        plt.text(0.85, 0.85, textstr, transform=fig.transFigure, ha='right')
        plt.xlabel("Cluster size")
        plt.ylabel("Number of sequences")
        plt.savefig(output_folder + "_clustersizes_sequences.pdf", bbox_inches='tight')
        print("Clustersize distribution plotted")
    
    return True

if __name__ == '__main__' and mode == 'clusterfull_fast':
    #Setup
    input_UMIfile = args.input
    input_READSfile = args.reads
    identity_thresh = args.identity
    size_thresh = args.size_thresh
    output_folder = args.output
    
    #Make output directory
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    print(f"Using CD-HIT fast clustering...")
    print(f"Using streaming approach for reads file: {input_READSfile}")
    
    # Run CD-HIT clustering
    success = run_cdhit_clustering(input_UMIfile, output_folder, identity_thresh, size_thresh, input_READSfile)
    
    if success:
        print("✓ Fast clustering completed successfully!")
    else:
        print("✗ Fast clustering failed!")


