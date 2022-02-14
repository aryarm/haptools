import numpy as np
import time
from .admix_storage import GeneticMarker, HaplotypeSegment

# TODO at a certain point we are going to need to ensure populations in model file are also in invcf files
#      This is only required for outputting the haplotypes in a vcf 

def simulate_gt(model_file, coords_file, seed=None):
    """
    Simulate admixed genotypes based on the parameters of model_file. 
    Arguments
        model
            File with the following structure. (Must be tab delimited)
            Header = # samples, Admixed, {all pop labels}
            Below  = generation#, frac, frac
            ex: 40    Admixed    CEU   YRI
                1       0        0.05  0.95
                2       0.20     0.05  0.75
        coords_file
            File containing genetic map coords in cM used for recombination points
        seed
            Seed used for randomization.
    Return

    """
    # initialize seed used for breakpoints
    if seed:
        np.random.seed(seed)
        print(f"Using seed {seed}")

    # load population samples and labels to be simulated 
    mfile = open(model_file, 'r')
    num_samples, not_used, *pops = mfile.readline().strip().split()
    num_samples = int(num_samples)

    # coord file structure chr variant cMcoord bpcoord
    coords = []
    with open(coords_file, 'r') as cfile:
        for line in cfile:
            # create marker from each line and append to coords
            data = line.strip().split()
            gen_mark = GeneticMarker(int(data[0]), float(data[2]), int(data[3])) 
            coords.append(gen_mark)

    # number of haplotypes simulated per generation
    haps_per_gen = max(5000, 20 * num_samples)

    # starting generation is 0
    prev_gen = 0
    next_gen_samples = []

    # Time code
    start = time.time()

    # iterate over generations in model file
    for gen in mfile:
        # setup population proportions and generations to simulate
        cur_gen, *pop_fracs = gen.strip().split()
        cur_gen = int(cur_gen)
        pop_fracs = np.array(pop_fracs).astype(np.float) 
        sim_gens = cur_gen - prev_gen
        
        assert sim_gens > 0
        assert np.sum(pop_fracs) == 1

        # sim generation
        print(f"Simulating generation {prev_gen+1}")
        next_gen_samples = _simulate(haps_per_gen, pop_fracs, prev_gen, coords, next_gen_samples)

        # simulate remaining generations
        for i in range(1, sim_gens):
            print(f"Simulating generation {prev_gen+i+1}")
            # simulate next generations using previous generations to sample from for admixture
            next_gen_samples = _simulate(haps_per_gen, pop_fracs, prev_gen+i, coords, next_gen_samples)

        prev_gen = cur_gen 

    end = time.time()
    print(f"Time elapsed for simulation: {end - start}")

    mfile.close()
    return num_samples, next_gen_samples

def write_breakpoints(samples, breakpoints, out):
    breakpt_file = out + '.bp'
    print(f"Outputting breakpoint file {breakpt_file}")

    # randomly sample breakpoints to get the correct amount of samples to output
    breakpoints = np.random.choice(breakpoints, size=2*samples, replace=False)

    with open(breakpt_file, 'a') as output:
        for ind, sample in enumerate(breakpoints):
            # Get sample number and haplotype number
            haplotype = (ind)%2 + 1
            sample_num = ind//2 + 1

            # write header for sample
            output.write(f"Sample_{sample_num}_{haplotype}\n")

            # write all segments for sample
            for segment in sample:
                # write all segments for current sample
                pop = segment.get_pop()
                chrom = segment.get_chrom()
                end_coord = segment.get_end_coord()
                output.write(f"{pop}\t{chrom}\t{end_coord}\n")
    return

def _simulate(samples, pop_fracs, pop_gen, coords, prev_gen_samples=None):
    # generate all samples
    hap_samples = []
    
    # pre compute haplotypes and parent population 
    # if there is no previous generation randomly choose population based on frac
    parent_pop = np.zeros(samples)
    if not prev_gen_samples:
        parent_pop = np.random.choice(np.arange(len(pop_fracs)), 
                                      size=samples, p=pop_fracs)

    # choose a haplotype to copy
    # if previous generation of samples randomize two haplotypes
    # otherwise only choose a single haplotype
    haplotypes = -1*np.ones(2*samples)
    if not parent_pop[0]:
        # choose two of previous samples as haplotypes to choose from
        haplotypes = np.random.randint(samples, size=2*samples)

        # ensure no two haplotypes are the same
        for hap in range(samples):
            while haplotypes[2*hap] == haplotypes[2*hap+1]:
                haplotypes[2*hap+1] = np.random.randint(samples)

    # store all probabilities to compare for recombination
    prob_vals = np.random.rand(samples*len(coords))

    # generate all samples
    for sample in range(samples):
        segments = []
        prev_chrom = -1
        prev_map_pos = 0.0
        start_bp = 0
        end_bp = 0
        p_pop = parent_pop[sample]
        haps = haplotypes[2*sample:2*sample+2]

        # haplotypes only matters when pop == 0 otherwise we are sampling from a population which has a label 
        # choose our starting chromosome
        homolog = np.random.randint(2)

        # iterate over all coords in the coords file
        for cind, coord in enumerate(coords):
            # get current chromosome and position
            cur_chrom = coord.get_chrom()
            cur_map_pos = coord.get_map_pos()

            # start of new segment used for storing previous generation samples  
            if segments:
                start_bp = segments[-1].get_end_coord() + 1

            # first segment so update chrom and starting map pos
            if prev_chrom < 0:
                prev_chrom = cur_chrom
                prev_map_pos = cur_map_pos

            # Ended chromosome write segment
            if cur_chrom != prev_chrom:
                # Store haplotype
                segments.extend(get_segment(p_pop, haps[homolog], prev_chrom, 
                                            start_bp, end_bp, prev_gen_samples))
                
                # update homolog, previous chrom, and previous position
                prev_chrom = cur_chrom
                homolog = np.random.randint(2)
                prev_map_pos = cur_map_pos
                continue
        
            # check if a recombination event occurs and if so write segment and swap homolog
            dist = cur_map_pos - prev_map_pos
            recomb_prob = 1-np.exp(-dist)
            if prob_vals[sample*samples+cind] < recomb_prob:
                # Store haplotype segments switching homologs
                segments.extend(get_segment(p_pop, haps[homolog], prev_chrom,
                                            start_bp, end_bp, prev_gen_samples))
                homolog = 1-homolog

            # update positions and chromosome
            prev_map_pos = cur_map_pos
            prev_chrom = cur_chrom
            end_bp = coord.get_bp_pos()

        # Record remaining segments
        if segments: start_bp = segments[-1].get_end_coord() + 1
        segments.extend(get_segment(p_pop, haps[homolog], prev_chrom,
                                    start_bp, end_bp, prev_gen_samples))

        # append segments of 1 sample
        hap_samples.append(segments)
    return hap_samples

def get_segment(pop, haplotype, chrom, start_coord, end_coord, prev_gen_samples):
    """
    Create a segment or segments for an individual of the current generation
    using either a population label (>0) or the previous generation's samples if 
    the admix pop type (0) is used. 
    Arguments
        pop - population to generate haplotype segments
        haplotype - index of range [0, len(prev_gen_samples)] to identify
                    the parent haplotype to copy segments from
        chrom - chromosome the haplotype segment lies on
        start_coord - starting coordinate from where to begin taking segments
                      from previous generation samples
        end_coord - ending coordinate of haplotype segment
        prev_gen_samples - the previous generation simulated used as the parents
                           for the current generation
    Returns
        A list of HaplotypeSegments storing the population type and end coordinate
    """
    # Take from population data not admixed data
    if pop:
        return [HaplotypeSegment(pop, chrom, end_coord)]
    
    # Take from previous admixed data
    else:
        segments = []
        prev_gen_segments = prev_gen_samples[haplotype]

        # iterate over haplotype segments to collect relevant ones
        # use binary search to find starting segment to collect information
        start_seg = start_segment(start_coord, chrom, prev_gen_segments)
        for prev_segment in prev_gen_segments[start_seg:]:
            if prev_segment.get_end_coord() > end_coord:
                break

            segments.append(prev_segment)
            
        # If there is not a segment within the boundary use the segment
        #    that spans the boundary respective population
        if not segments:
            out_pop = prev_segment.get_pop()
        else:
            out_pop = segments[-1].get_pop()

        # Append last segment using previous segments population
        segments.append(HaplotypeSegment(out_pop, chrom, end_coord))
        return segments

def start_segment(start, chrom, segments):
    """
    Find first segment that is on chrom and its end coordinate is > start via binary search.
    """
    low = 0
    high = len(segments)-1
    mid = 0

    # first segment > start implies segment prior end coord < start
    while low <= high:
        mid = (high+low) // 2
       
        # collect coordinate and chrom information
        cur_coord = segments[mid].get_end_coord()
        cur_chrom = segments[mid].get_chrom()
        if mid == 0:
            prev_coord = -1
            prev_chrom = -1
        else:
            prev_coord = segments[mid-1].get_end_coord()
            prev_chrom = segments[mid-1].get_chrom()
    
        # check if chromosomes match otherwise update
        if chrom == cur_chrom:
            # check for current coords loc
            if cur_coord < start:
                low = mid + 1

            elif cur_coord >= start:
                if prev_chrom < cur_chrom:
                    return mid

                if prev_chrom == cur_chrom and prev_coord < start:
                    return mid
                else:
                    high = mid - 1
                    
            else:
                return len(segments)

        elif chrom < segments[mid].get_chrom():
            high = mid - 1

        else:
            low = mid + 1

    return len(segments)


if __name__ == '__main__':
    test()