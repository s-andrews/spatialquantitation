#!/usr/bin/env python3
import gzip
import argparse
import sys
import pysam

def main():
    global options
    options = get_options()
    genes = parse_gtf(options.gtf_file)

    counts = count_overlaps(options.bam_file,genes)
    write_output(counts,options.outfile)

def write_output(counts, outfile):
    if not options.quiet:
        print("Writing counts to",outfile, file=sys.stderr)

    # Collect the full list of genes
    genes = set()
    for cell in counts.keys():
        for gene in counts[cell].keys():
            genes.add(gene)

    genes = sorted(list(genes))

    with open(outfile,"wt") as out:

        # Print the header first
        cell_ids = sorted(list(counts.keys()))

        line=["Gene"]
        line.extend(cell_ids)
        print("\t".join(line), file=out)


        # Print the data
        for gene in genes:
            line = [gene]

            for cell in counts.keys():
                if gene in counts[cell]:
                    line.append(str(len(counts[cell][gene])))
                else:
                    line.append("0")

            print("\t".join(line), file=out)





def count_overlaps(file,genes):
    if not options.quiet:
        print("Counting "+file)

    # We save and restore the verbosity so we can suppress an otherwise
    # annoying warning that we don't have an index (which we don't expect)
    save = pysam.set_verbosity(0)
    bamfile = pysam.AlignmentFile(file,"rb")
    pysam.set_verbosity(save)

    # We're going to keep track of what we've seen in this structure
    # we index by cell identity then record the UMIs we've seen for
    # each gene so we can deduplicate
    #
    # TODO: We should check whether we need to store position as well
    # as UMI
    counts = {}

    current_chromosome = ""
    chr_genes = []
    start_index = 0

    cell_count = 0

    for x,read in enumerate(bamfile.fetch(until_eof=True)):
        chromosome = read.reference_name
        if chromosome.lower().startswith("chr"):
            chromosome = chromosome[3:]

        if chromosome != current_chromosome:
            if not options.quiet:
                print("Searching chromosome",chromosome)
            current_chromosome = chromosome
            if chromosome in genes:
                chr_genes = genes[chromosome]
            else:
                chr_genes = []

            start_index = 0

        for i in range(start_index,len(chr_genes)):
            # See if this is a feature we never need to look at again
            if i==start_index and read.reference_start > chr_genes[i]["end"]:
                start_index += 1
                continue
        
            # Can we skip based on strand
            if read.is_reverse != chr_genes[i]["is_reverse"]:
                continue

            # Can we stop looking
            if chr_genes[i]["start"] > read.reference_end:
                break

            # Do we overlap
            if read.reference_start < chr_genes[i]["end"] and read.reference_end > chr_genes[i]["start"]:
                # We need to check whether an alignment segment actually overlaps, rather
                # than just an intron
                found_block = False
                for rs,re in read.blocks:
                    if rs < chr_genes[i]["end"] and re > chr_genes[i]["start"]:
                        found_block = True
                        break

                if not found_block:
                    # None of our blocks actually overlap this gene
                    continue

                # We need to extract the UMI and cell identity from the read
                sections = read.query_name.split(":")
                umi = sections[-1]

                # We add the read start position to the umi for completeness
                if read.is_reverse:
                    umi += str(read.reference_end)
                else:
                    umi += str(read.reference_start)

                cell_id = ":".join(sections[-3:-1])
                if not cell_id in counts:
                    cell_count += 1
                    counts[cell_id] = {}

                if not chr_genes[i]["name"] in counts[cell_id]:
                    counts[cell_id][chr_genes[i]["name"]] = set()


                if umi in counts[cell_id][chr_genes[i]["name"]]:
                    continue

                counts[cell_id][chr_genes[i]["name"]].add(umi)               
 
    if not options.quiet:
        print("Found",cell_count,"cells", file=sys.stderr)
    return counts

def parse_gtf(file):
    if not options.quiet:
        print(f"Parsing genes from {file}", file=sys.stderr)

    genes = {}

    gene_count = 0
    failed_count = 0

    infh = None

    if file.lower().endswith(".gz"):
        infh = gzip.open(file,"rt", encoding="utf8")

    else:
        infh = open(file,"rt",encoding="utf-8")

    for line in infh:
        if line.startswith("#"):
            continue

        sections = line.strip().split("\t")

        if sections[2] != "gene":
            continue

        chromosome = sections[0]
        start = int(sections[3])
        end = int(sections[4])
        is_reverse = sections[6] == "-"
        name = None

        ##### TESTING ONLY #####
#            if chromosome != "1":
#                break

        annotations = sections[8].split(";")
        for annot in annotations:
            if annot.strip().startswith("gene_name"):
                name = annot.strip().split(" ",2)[1].replace("\"","").strip()
                break

        if name is None:
            failed_count += 1
            continue
        else:
            gene_count += 1

        if not chromosome in genes:
            genes[chromosome] = []

        genes[chromosome].append({"start":start, "end":end, "is_reverse": is_reverse, "name":name})

    infh.close()

    # It looks like the genes aren't listed in order in the gtf file so
    # we need to fix that.

    for n in genes.keys():
        g = genes[n]
        g.sort(key=lambda x: x["start"])


    if gene_count == 0:
        print(f"No genes found.  Rejected {failed_count} genes for having no gene_name attribute", file=sys.stderr)

    else:
        if not options.quiet:
            print(f"Found {gene_count} genes and rejected {failed_count} genes with no name", file=sys.stderr)

    return genes


            


def get_options():

    parser = argparse.ArgumentParser(description="Quantitate single cell spatial transcriptomic data")

    parser.add_argument("gtf_file", type=str, help="The GTF file of features")
    parser.add_argument("bam_file", type=str, help="The BAM file to quantitate")
    parser.add_argument("outfile", type=str, help="The output file to write to")

    parser.add_argument("--quiet",action="store_true", default=False, help="Suppress progress messages")

    return(parser.parse_args())


if __name__ == "__main__":
    main()