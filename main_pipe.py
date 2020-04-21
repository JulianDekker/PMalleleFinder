import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from BCBio import GFF
import sys, traceback
import pprint
import subprocess
import pandas as pd
import scripts.gene as gene
from datetime import datetime

now = datetime.now()
outdir = 'output/'+now.strftime("%Y%m%d%H%M%S")


def get_rss(genelist):
    """
    Creates position objects for rss sequences
    :param genelist:
    :return:
    """
    rsslist = []
    for gene in genelist:
        boundaries = gene.get_gene_boundaries()
        if gene.get_name()[3] == 'V':
            if boundaries[2] == 1:
                rsslist = rsslist + [(gene.get_name() + '_f_rss', (int(boundaries[1]), int(boundaries[1]) + 60,
                                                      gene.get_exons()[list(gene.get_exons().keys())[0]][2]),
                                 'RSS', gene.get_chro(), boundaries[2])]
            else:
                rsslist = rsslist + [(gene.get_name() + '_f_rss', (int(boundaries[0]) - 60, int(boundaries[0]),
                                                                   gene.get_exons()[list(gene.get_exons().keys())[0]][
                                                                       2]),
                                      'RSS', gene.get_chro(), boundaries[2])]
        elif gene.get_name()[3] == 'J':
            if boundaries[2] == 1:
                rsslist = rsslist + [(gene.get_name()+'_r_rss', (int(boundaries[0])-60, int(boundaries[0]),
                                                    gene.get_exons()[list(gene.get_exons().keys())[0]][2]),
                                 'RSS', gene.get_chro(), "1")]
            else:
                rsslist = rsslist + [(gene.get_name() + '_r_rss', (int(boundaries[1]), int(boundaries[1]) + 60,
                                                                   gene.get_exons()[list(gene.get_exons().keys())[0]][
                                                                       2]),
                                      'RSS', gene.get_chro(), boundaries[2])]
        elif gene.get_name()[3] == 'D':
            if boundaries[2] == 1:
                rsslist = rsslist + [(gene.get_name()+'_f_rss', (int(boundaries[0])-60, int(boundaries[0]),
                                                    gene.get_exons()[list(gene.get_exons().keys())[0]][2]), 'RSS',
                                  gene.get_chro(), "1"),
                    (gene.get_name()+'_r_rss', (int(boundaries[1]), int(boundaries[1])+60,
                                                gene.get_exons()[list(gene.get_exons().keys())[0]][2]), 'RSS',
                     gene.get_chro(), boundaries[2])]
            else:
                rsslist = rsslist + [(gene.get_name() + '_f_rss', (int(boundaries[1]), int(boundaries[1])+60,
                                                                   gene.get_exons()[
                                                                       list(gene.get_exons().keys())[0]][2]), 'RSS',
                                      gene.get_chro(), "1"),
                                     (gene.get_name() + '_r_rss', (int(boundaries[0])-60, int(boundaries[0]),
                                                                   gene.get_exons()[
                                                                       list(gene.get_exons().keys())[0]][2]), 'RSS',
                                      gene.get_chro(), boundaries[2])]
    return rsslist


def make_all_fasta():
    """
    reads all fasta files in output directory and compiles them into a multi fasta file
    :return:
    """
    import glob
    read_files = glob.glob(outdir + '/sequences/*.fasta')
    with open(outdir + '/sequences/'+ outdir.split('/')[len(outdir.split('/'))-1] + '_all.fasta', "wb") as outfile:
        for f in read_files:
            with open(f, "rb") as infile:
                outfile.write(infile.read())


def report(gff, vcf, ref, populationfile, genelist, nopedlist, filter, threshold):
    """
    writes a report from the pipeline analysis
    :return:
    """
    print("Reporting results..")
    with open(outdir + '/report_' + outdir.split('/')[len(outdir.split('/'))-1] + '.txt', 'a') as report:
        sequences_total = 0
        reportlist = []
        for gene in genelist:
            seqs = [seq for seq in SeqIO.parse(open(outdir + '/sequences/' + gene.get_name() + '.fasta'),
                                               'fasta') if 'rss' not in gene.get_name()]
            sequences_total += len(seqs)
            try:
                with open(outdir + '/vcf/' + gene.get_name() + '-py.csv', 'r') as pos:
                    line = pos.readline().strip()
                    mutated_pos = len(line.split(','))
            except:
                mutated_pos = 0
            if 'rss' in gene.get_name():
                seqs = [seq for seq in
                        SeqIO.parse(open(outdir + '/sequences/' + gene.get_name() + '.fasta'), 'fasta')]
                reportlist.append([gene.get_name(), len(seqs), mutated_pos])
            else:
                reportlist.append([gene.get_name(), len(seqs), mutated_pos])
        genecount = len([i for i in genelist if 'rss' not in i.get_name()])
        report.write("--Report save location: " + outdir + "---\n\nGff file: " + gff + "\nVcf files: " + vcf + "\nReference fasta: " + ref +
                     "\nPopulation file: " + populationfile +"\n\nTotal gene count: " + str(genecount) +
                     "\nTotal allele count: " + str(sequences_total)
                     + "\nSkipped exons (no variation): " + str(len(nopedlist)) +"\nThreshold: " + str(threshold) +"\n")
        if filter:
            report.write("Applied filters: ["+",".join(filter)+"]\n")
        report.write('\nGene\talleles\tmutated_positions\n')
        for i in reportlist:
            report.write(str(i[0])+'\t'+str(i[1])+'\t'+ str(i[2])+'\n')


def hap_seq(gene):
    """
    Creates alleles for all variants in the haplotype files.
    input: genelist
    """
    keys = []
    values = []
    ref = SeqIO.parse(open(outdir + '/reffasta/' + gene.get_name() + '.fasta'), 'fasta')
    ref = [seq for seq in ref]
    for seq in ref:
        start, stop = seq.id.split(':')[1].split('-')
        values = values + list(seq.seq)
        for i in range(int(start), int(stop) + 1):
            keys.append(str(i))
    try:
        with open(outdir + '/vcf/' + gene.get_name() + '-py.csv', 'r') as n:  # Open the file
            lines = n.readlines()   # Read all the lines in the CSV file, with each line as a
        h = []
        for j in lines:
            h.append(j.split()[0].split(","))  # Split each line in into strings based on delimiters
        indices = [keys.index(i.split('.')[0]) for i in h[0]]
        h1 = []
        for y in range(1, len(h)):
            l1new = values
            for i, j in zip(indices, h[y]):
                l1new[i] = j
            h1.append(''.join(l1new))
        h1.append(''.join(values))
        h1 = list(sorted(set(h1), key=h1.index))
        open(outdir + '/sequences/' + gene.get_name() + '.fasta', "w").write('')
        with open(outdir + '/sequences/' + gene.get_name() + '.fasta', "w") as output:
            for i in range(len(h1)): #
                if gene.get_exons()[list(gene.get_exons().keys())[0]][2] == -1:
                    output.write(">" + gene.get_name() + "_" + str(i + 1) + "\n" +
                                 str(Seq(str(h1[i])).reverse_complement() + "\n"))
                else:
                    output.write(">" + gene.get_name() + "_" + str(i + 1) + "\n" + str(h1[i]) + "\n")
    except Exception as e:
        open(outdir + '/sequences/' + gene.get_name() + '.fasta', "w").write('')
        with open(outdir + '/sequences/' + gene.get_name() + '.fasta', "w") as output:
            seq = ''
            for sequence in ref:
                seq += str(sequence.seq)
            if gene.get_exons()[list(gene.get_exons().keys())[0]][2] == -1:
                seq = str(Seq(seq).reverse_complement())
            output.write(">" + gene.get_name() + "_" + str(1) + "\n" + (seq) + "\n")


def ped2hap(gene, popfile, nopedlist, threshold):
    """
    Creates haplotype files from ped files
    :param genelist:
    :return:
    """
    prefix = outdir + '/vcf/'
    exons = sorted([int(k) for k in gene.get_exons().keys()])
    exons = [prefix+gene.get_name()+'-exon'+str(exon) for exon in exons if gene.get_name()+'-exon'+str(exon) not in nopedlist]
    if len(exons) >= 1:
        command('Rscript scripts/Hapmerge-exall.R '+popfile+' '+threshold+' '+' '.join(exons))
        command('Rscript scripts/CSV-Py-v3.R ' + outdir + '/vcf/ ' + gene.get_name())


def vcf2ped(genpos, popfile):
    """Creates PED files from VCF files
        input is genpos and population annotations file in sample\tlocation format.
    """
    print("Generating PED files..")
    popfile = pd.read_csv(popfile, sep='\t', header=0)
    nopedlist=[]
    for file in genpos:
        filename = outdir + '/vcf/' + file[0] + "-exon" + str(file[4]) + '.vcf'
        try:
            df1 = pd.read_table(filename, comment='#', header=None)
            df = df1.drop(columns=[0, 1, 2, 5, 6, 7, 8])
            colsf = popfile['sample'].to_list()
            df.columns = ['REF', 'ALT'] + colsf
            # df = df[df.ALT != '<CN0>']
            # df.index = list(range(len(df)))
            cols = list(df.columns)
            seperator = '    '
            # %%
            master_list = []
            for row in df.index:
                ref = df.loc[row, 'REF']  # ref for each row
                alt = df.loc[row, 'ALT'].split(',')  # Alt for each row
                aa = []
                for column in cols[2:]:  # Spilt each column and search
                    p = df.loc[row, column].split('|')
                    a_0 = [n1 for n1, n2 in enumerate(p) if n2 == '0']  # INdex of 0
                    a_1 = [n1 for n1, n2 in enumerate(p) if n2 == '1']  # INdex of 1
                    a_2 = [n1 for n1, n2 in enumerate(p) if n2 == '2']  # INdex of 2
                    a_3 = [n1 for n1, n2 in enumerate(p) if n2 == '3']  # INdex of 3
                    # Now there are 2 possibilities
                    # 1) There is atleast one '0' present --> reference + alt
                    # 2) There is no '0' present ---> alt
                    # Case 1
                    if '0' in p:
                        # Now find the index of '0'. Also, there can either be one 0 or two 0s
                        # Now either both are 0 or only one is 0
                        # case 1_1 both are zero
                        if len(a_0) == 2:
                            p[0] = ref
                            p[1] = ref
                        # Case 1_2 only one is zero, the other is any one of 1, 2, 3
                        elif len(a_0) == 1:
                            p[a_0[0]] = ref
                            if len(a_1) == 1:  # If the other number is 1
                                p[a_1[0]] = alt[0]
                            elif len(a_2) == 1:  # If the other number is 2
                                p[a_2[0]] = alt[1]
                            elif len(a_3) == 1:  # If the other number is 3
                                p[a_3[0]] = alt[2]
                        aa.append(p[0] + seperator + p[1])
                    # Case 2 ---> 0 is not present at all, alt
                    elif '0' not in p:
                        # Case 2_1, both are 1
                        if len(a_1) == 2:
                            p[0] = alt[0]
                            p[1] = alt[0]
                        # Case 2_2, only 1 is one
                        elif len(a_1) == 1:
                            p[a_1[0]] = alt[0]
                            if len(a_2) == 1:  # Case 2_2_1, the other place is 2
                                p[a_2[0]] = alt[1]
                            elif len(a_3) == 1:  # Case 2_2_3, the other place is 3
                                p[a_3[0]] = alt[2]
                        elif len(a_2) == 2:  # Both places are 2
                            p[0] = alt[1]
                            p[1] = alt[1]
                        elif len(a_2) == 1:  # Only one place is 2
                            p[a_2[0]] = alt[1]
                            if len(a_1) == 1:  # Case 2_2_1, the other place is 1
                                p[a_1[0]] = alt[0]
                            elif len(a_3) == 1:  # Case 2_2_3, the other place is 3
                                p[a_3[0]] = alt[2]
                        elif len(a_3) == 2:  # Both places are 3
                            p[0] = alt[2]
                            p[1] = alt[2]
                        elif len(a_3) == 1:  # Only one place is 3
                            p[a_3[0]] = alt[2]
                            if len(a_1) == 1:  # Case 2_2_1, the other place is 1
                                p[a_1[0]] = alt[0]
                            elif len(a_2) == 1:  # Case 2_2_3, the other place is 2
                                p[a_2[0]] = alt[1]
                        aa.append(p[0] + seperator + p[1])
                master_list.append(aa)
            ReF = df['REF']
            AlT = df['ALT']
            final_data = pd.DataFrame(master_list)
            final_data.columns = colsf
            dd = final_data.T
            dd.insert(0, 'index2', dd.index)
            dd.to_csv(filename + '.ped', header=False, sep='\t')
        except Exception as e:
            nopedlist.append(file[0] + "-exon" + str(file[4]))
    return nopedlist


def cutVCF(item, vcf):
    """
    cut's the vcf in pieces defined by genpos items
    :param genpos:
    :param vcf: vcf files for chr7 and 14 seperated by ;
    :return:
    """
    name, loc, type, chro, ex = item
    output, error = '', ''
    if 'chr' in chro:
        chrom = chro[3:]
    else:
        chrom = chro
    if len(vcf.split(';')) == 2:
        vcf7, vcf14 = vcf.split(';')
        if int(loc[0]) > int(loc[1]):
            loc = (loc[1], loc[0])
        if chro == 'chr7':
            output, error = command("tabix -fh "+vcf7+" "+chrom+":"+str(loc[0])+"-"+str(loc[1]))
        elif chro == 'chr14':
            output, error = command("tabix -fh " + vcf14 + " " + chrom + ":" + str(loc[0]) + "-" + str(loc[1]))
    else:
        output, error = command("tabix -fh " + vcf + " " + chrom + ":" + str(loc[0]) + "-" + str(loc[1]))
    if output is not '':
        open(outdir + '/vcf/' + name + "-exon" + str(ex) + ".vcf", 'wb').write(output)


def cutFasta(item, fasta):
    """
    cut's the fasta in pieces defined by genpos items
    :param genpos:
    :param fasta: fasta reference file
    :return:
    """
    name, loc, type, chro, ex = item
    if int(loc[0]) > int(loc[1]):
        loc = (loc[1], loc[0])
    output, error = command("samtools faidx "+fasta+" "+chro+":"+str(loc[0])+"-"+str(loc[1]))
    open(outdir + '/reffasta/' + name + ".fasta", 'ab').write(output)


def extract_pos(gff, filter=None):
    """
    retrieves information from gff files, also is able to filter on certain genes.
    :param gff:
    :param filter:
    :return:
    """
    print("Processing gencode annotation..")
    in_handle = open(gff)
    limit_info = dict(
        gff_type=['exon', 'CDS']
    )
    if filter:
        print('Filtering on: ', filter)
        if filter[0].startswith('#'):
            chrfilter = filter[0].split()[1].strip().split(',')
            limit_info['gff_id'] = ['chr'+(str(x)) for x in chrfilter]
            filter = filter[1::]
    genlist = set()
    for rec in GFF.parse(in_handle, limit_info=limit_info):
        for feature in rec.features:
            if feature.type == 'inferred_parent':
                for f in feature.sub_features:
                    if containsTR(f.qualifiers['gene_name'][0], filter):
                        if f.type == 'CDS' and 'chr' in rec.id:
                            if containCDSGENE(genlist, f.qualifiers['gene_name'][0], 'gene'):
                                genlist = discarditem(genlist, f.qualifiers['gene_name'][0])
                            genlist.add((f.qualifiers['gene_name'][0], (f.location.start+1, f.location.end,
                                                                        f.location.strand), f.type, rec.id,
                                         f.qualifiers['exon_number'][0]))
                            entries = [genname for genname in genlist if genname[0] == f.qualifiers['gene_name'][0]]
                            for entry in entries:
                                if entry[2] == 'exon':
                                    genlist.remove(entry)
                        elif f.type == 'exon' and 'chr' in rec.id:
                            if not containCDSGENE(genlist, f.qualifiers['gene_name'][0], 'CDS'):
                                genlist.add((f.qualifiers['gene_name'][0], (f.location.start+1, f.location.end,
                                                                            f.location.strand), f.type, rec.id,
                                             f.qualifiers['exon_number'][0]))
            else:
                if containsTR(feature.qualifiers['gene_name'][0], filter):
                    try:
                        ex_number = feature.qualifiers['exon_number'][0]
                    except:
                        ex_number = 1
                    if feature.type == 'CDS' and 'chr' in rec.id:
                        if containCDSGENE(genlist, feature.qualifiers['gene_name'][0], 'gene'):
                            genlist = discarditem(genlist, feature.qualifiers['gene_name'][0])
                        genlist.add((feature.qualifiers['gene_name'][0], (feature.location.start+1, feature.location.end,
                                                                          feature.location.strand), feature.type,
                                     rec.id, ex_number))
                        entries = [genname for genname in genlist if genname[0] == feature.qualifiers['gene_name'][0]]
                        for entry in entries:
                            if entry[2] == 'exon':
                                genlist.remove(entry)
                    elif feature.type == 'exon' and 'chr' in rec.id:
                        if not containCDSGENE(genlist, feature.qualifiers['gene_name'][0], 'CDS'):
                            genlist.add((feature.qualifiers['gene_name'][0], (feature.location.start+1,
                                                                              feature.location.end,
                                                                              feature.location.strand), feature.type,
                                         rec.id, ex_number))
        #for item in genlist:
        #    if "chr" not in item[3]:
        #        genlist.remove(item)
    return sorted(genlist)


def containsTR(str, list=None):
    """
    checks if str contains item in list
    :param str:
    :return:
    """
    if not list:
        list = ['TRAV', 'TRAC', 'TRAJ',
                'TRBV', 'TRBC', 'TRBJ', 'TRBD',
                'TRGV', 'TRGC', 'TRGJ',
                'TRDV', 'TRDC', 'TRDJ', 'TRDD']
    for item in list:
        if item in str:
            return True
    return False


def command(command):
    """
    Creates a bash command
    :param command:
    :return:
    """
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output, error


def containCDSGENE(genlist, name, type):
    """
    checks if [CDS/GENE] exists in set for matching name
    :param genlist:
    :param name:
    :param type:
    :return:
    """
    for i in genlist:
        if i[0] == name and i[2] == type:
            return True
    return False


def discarditem(genlist, item):
    """
    Removes item from set if first element matches item
    :param genlist:
    :param item:
    :return:
    """
    for g in genlist:
        if g[0] == item:
            genlist.remove(g)
            return genlist
    return genlist


def in_genelist(genelist, name):
    """
    checks if a name is in a list of gene objects
    :param genelist:
    :param name:
    :return:
    """
    for i in genelist:
        if i.get_name() == name:
            return True
    return False

def printhelp():
    print("""
    TR allele pipeline:

    --gff=input gff <file paths>
    --vcf="inputvcf1;inputvcf2" <file paths>
    --popfile=inputpopfile <file path>
    --ref=ref.fasta <file path>
    --filter=filterlist <file path>
    --threshold=value <integer>
    --rss=TRUE|FALSE <boolean>
    --outdir=directory <string>

    gff is a valid annotation file in gencode format
    vcf requires a vcf.gz and requires a vcf.tbi to be present. needs to be two filenames in quotes.
    popfile is a tab seperated file containing sample and populaltion information: SampleID Pop Superpop
    ref is the reference genome in fasta format.
    filter is an optional argument to filter on a set of genes.
    threshold is an integer defining the minimum support a allele must have to be saved in the HAP file.
    rss is the option to enable generating rss sequences as well.
    outdir is the desired output directory name.
                """)


def run_normal(gff, genepos, vcf, ref, populationfile, filterfile, threshold):
    """
    Run TR diversity normally
    :param genepos:
    :param vcf:
    :param ref:
    :param populationfile:
    :param filterfile:
    :param threshold:
    :return:
    """
    genelist = []
    for pos in genepos:
        if in_genelist(genelist, pos[0]):
            for i in genelist:
                if i.get_name() == pos[0]:
                    i.add_exon(pos[4], pos[1])
                    break
        else:
            genelist.append(gene.Gene(pos[0], {str(pos[4]): pos[1]}, chro=pos[3]))
    print("Slicing VCF/reference fasta..")
    command('mkdir ' + outdir + '/')
    command('mkdir ' + outdir + '/vcf/')
    command('mkdir ' + outdir + '/reffasta/')
    command('mkdir ' + outdir + '/sequences/')
    for item in genepos:
        cutVCF(item, vcf)
        cutFasta(item, ref)
    nopedlist = vcf2ped(genepos, populationfile)
    print('generating HAP files and allele sequences..')
    for genes in genelist:
        ped2hap(genes, populationfile, nopedlist, threshold)
        hap_seq(genes)
    make_all_fasta()
    report(gff, vcf, ref, populationfile, genelist, nopedlist, filterfile, threshold)  # add filtergenelist = []


def run_rss(gff, genepos, vcf, ref, populationfile, filterfile, threshold):
    """
    Run TR diversity with RSS sequences included
    :param genepos:
    :param vcf:
    :param ref:
    :param populationfile:
    :param filterfile:
    :param threshold:
    :return:
    """
    genelist = []
    for pos in genepos:
        if in_genelist(genelist, pos[0]):
            for i in genelist:
                if i.get_name() == pos[0]:
                    i.add_exon(pos[4], pos[1])
                    break
        else:
            genelist.append(gene.Gene(pos[0], {str(pos[4]): pos[1]}, chro=pos[3]))
    rsspos = get_rss(genelist)
    rsslist = []
    for pos in rsspos:
        if in_genelist(rsslist, pos[0]):
            for i in rsslist:
                if i.get_name() == pos[0]:
                    i.add_exon(pos[4], pos[1])
                    break
        else:
            rsslist.append(gene.Gene(pos[0], {str(pos[4]): pos[1]}))
    genelist = genelist + rsslist
    print("Slicing VCF/reference fasta..")
    command('mkdir ' + outdir + '/')
    command('mkdir ' + outdir + '/vcf/')
    command('mkdir ' + outdir + '/reffasta/')
    command('mkdir ' + outdir + '/sequences/')
    #command('mkdir output/vcf/' + date_time)
    #command('mkdir output/reffasta/' + date_time)
    for item in genepos+rsspos:
        cutVCF(item, vcf)
        cutFasta(item, ref)
    get_rss(genelist)
    nopedlist = vcf2ped(genepos+rsspos, populationfile)
    #command('mkdir output/sequences/' + date_time)
    print('generating HAP files and allele sequences..')
    for genes in genelist:
        ped2hap(genes, populationfile, nopedlist, threshold)
        hap_seq(genes)
    make_all_fasta()
    report(gff, vcf, ref, populationfile, genelist, nopedlist, filterfile, threshold)  # add filtergenelist = []


def main():
    """
    --gff=input gff <file paths>
    --vcf="inputvcf1;inputvcf2" <file paths>
    --popfile=inputpopfile <file path>
    --ref=ref.fasta <file path>
    --filter=filterlist <file path>
    --threshold=value <integer>
    --rss=TRUE|FALSE <boolean>
    --outdir=directory <string>
    gff is a valid annotation file in gencode format
    vcf requires a vcf.gz and requires a vcf.tbi to be present. needs to be two filenames in quotes.
    popfile is a tab seperated file containing sample and populaltion information: SampleID Pop Superpop
    ref is the reference genome in fasta format.
    filter is an optional argument to filter on a set of genes.
    threshold is an integer defining the minimum support a allele must have to be saved in the HAP file.
    :return:
    """
    global outdir
    gff, vcf, populationfile, ref = '', '', '', ''
    filterfile = None
    rss = False
    threshold = '4'
    for arguments in sys.argv[1::]:
        if arguments.startswith('--gff'):
            gff = arguments.split('=')[1]
        elif arguments.startswith('--vcf'):
            vcf = arguments.split('=')[1]
        elif arguments.startswith('--popfile'):
            populationfile = arguments.split('=')[1]
        elif arguments.startswith('--ref'):
            ref = arguments.split('=')[1]
        elif arguments.startswith('--filter'):
            filterfile = [x.strip() for x in open(arguments.split('=')[1], 'r').readlines()]
        elif arguments.startswith('--threshold'):
            threshold = str(arguments.split('=')[1])
        elif arguments.startswith('--rss'):
            rss = arguments.split('=')[1]
            if rss.upper() == "TRUE":
                rss = True
        elif arguments.startswith('--outdir'):
            outdir = str(arguments.split('=')[1])
        elif arguments.startswith('--help'):
            printhelp()
        elif '=' not in arguments:
            printhelp()
            raise ValueError(
                'Argument "' + arguments + '" recognised.. Please use --help to list the arguments available.\n')
        elif arguments.split('=')[0] not in ['--filter', '--vcf', '--popfile', '--ref', '--gff', '--threshold', '--rss',
                                             '--outdir', '--help']:
            printhelp()
            raise ValueError(
                '\nArgument "'+arguments+'" recognised.. Please use --help to list the arguments available.\n')

    if gff == '' or vcf == '' or populationfile == '' or ref == '':
        print("Not enough arguments to run.\n\nPlease make sure to specify:\n\t--vcf=\"input1.vcf;input2.vcf\"\n\t--gff=input.gff\n\t--popfile=inputpopfile\n\t--ref=ref.fasta\nIn the command line options.")
    else:
        genepos = extract_pos(gff, filterfile)
        if rss:
            run_rss(gff, genepos, vcf, ref, populationfile, filterfile, threshold)
        else:
            run_normal(gff, genepos, vcf, ref, populationfile, filterfile, threshold)
        print("""
========================================
  ______ _       _     _              _ 
 |  ____(_)     (_)   | |            | |
 | |__   _ _ __  _ ___| |__   ___  __| |
 |  __| | | '_ \| / __| '_ \ / _ \/ _` |
 | |    | | | | | \__ \ | | |  __/ (_| |
 |_|    |_|_| |_|_|___/_| |_|\___|\__,_|
 =======================================                               
        """)
        print("Output ID: " + outdir)


main()

