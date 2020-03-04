import Bio
from Bio import SeqIO
from BCBio import GFF
import sys
import pprint
import subprocess
import pandas as pd
import scripts.gene as gene
from datetime import datetime

now = datetime.now()
date_time = now.strftime("%Y%m%d%H%M%S")

def hap_seq(genelist):
    """
    Creates alleles for all variants in the haplotype files.
    input: genelist
    """
    print("Generating allele sequence files..")
    bashCommand = 'mkdir output/sequences/'+date_time
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    for gene in genelist:
        keys = []
        values = []
        ref = SeqIO.parse(open('output/reffasta/'+date_time+'/'+gene.get_name()+'.fasta'), 'fasta')
        ref = [seq for seq in ref]
        for seq in ref:
            start, stop = seq.id.split(':')[1].split('-')
            values = values + list(seq.seq)
            for i in range(int(start), int(stop) + 1):
                keys.append(str(i))
        try:
            with open('output/vcf/'+ date_time + '/' + gene.get_name() + '-py.csv', 'r') as n:  # Open the file
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
            h1 = list(set(h1))
            open('output/sequences/'+date_time+'/' + gene.get_name() + '.fasta', "w").write('')
            with open('output/sequences/'+date_time+'/' + gene.get_name() + '.fasta', "w") as output:
                for i in range(len(h1)): #
                    if gene.get_exons()[list(gene.get_exons().keys())[0]] == -1:
                        output.write(">" + gene.get_name() + "_" + str(i + 1) + "\n" + str(h1[i])[::-1] + "\n")
                    else:
                        output.write(">" + gene.get_name() + "_" + str(i + 1) + "\n" + str(h1[i]) + "\n")
                    #print('writing result to: output/reffasta/seq/' + gene.get_name())
        except:
            open('output/sequences/'+date_time+'/' + gene.get_name() + '.fasta', "w").write('')
            with open('output/sequences/'+date_time+'/' + gene.get_name() + '.fasta', "w") as output:
                seq = ''
                for sequence in ref:
                    seq += str(sequence.seq)
                if gene.get_exons()[list(gene.get_exons().keys())[0]] == -1:
                    seq = seq[::-1]
                output.write(">" + gene.get_name() + "_" + str(1) + "\n" + (seq) + "\n")
                #print('writing result to: output/reffasta/seq/' + gene.get_name())


def ped2hap(genelist, popfile, nopedlist, threshold):
    """
    Creates haplotype files from ped files
    :param genelist:
    :return:
    """
    print("Generating HAP files..")
    prefix = 'output/vcf/'+date_time+'/'
    for genes in genelist:
        exons = sorted(genes.get_exons().keys())
        exons = [exon for exon in exons if genes.get_name()+'-exon'+str(exon) not in nopedlist]
        nexon = len(exons)

        if nexon == 1:
            bashCommand = 'Rscript scripts/Ped2Hap_shark.R '+popfile+' '+prefix+genes.get_name()+'-exon'+exons[0]+' '+threshold
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
        elif nexon == 2:
            ex1, ex2, = prefix+genes.get_name() + '-exon' + exons[0], prefix+genes.get_name()+'-exon'+exons[1]
            bashCommand = 'Rscript scripts/Hapmerge.R '+popfile+' '+ex1+' '+ex2+' '+threshold
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
        elif nexon == 3:
            ex1, ex2, ex3 = prefix+genes.get_name() + '-exon' + exons[0], \
                            prefix+genes.get_name() + '-exon' + exons[1], \
                            prefix+genes.get_name() + '-exon' + exons[2]
            bashCommand = 'Rscript scripts/Hapmerge-ex3.R '+popfile+' '+ex1+' '+ex2+' '+ex3+' '+threshold
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
        elif nexon == 4:
            ex1, ex2, ex3, ex4 = prefix+genes.get_name() + '-exon' + exons[0], \
                                 prefix+genes.get_name() + '-exon' + exons[1], \
                                 prefix+genes.get_name() + '-exon' + exons[2], \
                                 prefix+genes.get_name() + '-exon' + exons[3]
            bashCommand = 'Rscript scripts/Hapmerge-ex4.R '+popfile+' '+ex1+' '+ex2+' '+ex3+' '+ex4+' '+threshold
            process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
            output, error = process.communicate()
        elif nexon > 4:
            print(len(exons), exons)
            raise ValueError("Too many exons to process")
    bashCommand = 'Rscript scripts/CSV-Py-v2.R output/vcf/'+date_time+'/'
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()


def vcf2ped(genpos, popfile):
    """Creates PED files from VCF files
        input is genpos and population annotations file in sample\tlocation format.
    """
    print("Generating PED files..")
    popfile = pd.read_csv(popfile, sep='\t', header=0)
    nopedlist=[]
    for file in genpos:
        filename = 'output/vcf/'+ date_time + '/' + file[0] + "-exon" + str(file[4]) + '.vcf'
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
            #print("No variation found in: "+filename, e)
    return nopedlist


def cutVCF(genpos, vcf):
    """
    cut's the vcf in pieces defined by genpos items
    :param genpos:
    :param vcf: vcf files for chr7 and 14 seperated by ;
    :return:
    """
    print("Slicing VCF files..")
    bashCommand = 'mkdir output/vcf/' + date_time
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    vcf7, vcf14 = vcf.split(';')
    for item in genpos:
        name, loc, type, chro, ex = item
        if int(loc[0]) > int(loc[1]):
            loc = (loc[1], loc[0])
        bashCommand = ''
        if chro == 'chr7':
            bashCommand = "tabix -fh "+vcf7+" "+chro[3:]+":"+str(loc[0])+"-"+str(loc[1])
            #print('Getting '+name+' from: '+vcf7+" "+chro[3:]+":"+str(loc[0])+"-"+str(loc[1]))
        elif chro == 'chr14':
            bashCommand = "tabix -fh " + vcf14 + " " + chro[3:] + ":" + str(loc[0]) + "-" + str(loc[1])
            #print('Getting ' + name + ' from: ' + vcf14 + " " + chro[3:] + ":" + str(loc[0]) + "-" + str(loc[1]))
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        open('output/vcf/'+date_time+'/'+name+"-exon"+str(ex)+".vcf", 'wb').write(output)


def cutFasta(genpos, fasta):
    """
    cut's the fasta in pieces defined by genpos items
    :param genpos:
    :param fasta: fasta reference file
    :return:
    """
    print("Slicing reference fasta..")
    bashCommand = 'mkdir output/reffasta/' + date_time
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    for item in genpos:
        name, loc, type, chro, ex = item
        open('output/reffasta/'+date_time+'/' + name + ".fasta", 'w').write('')
    for item in genpos:
        name, loc, type, chro, ex = item
        if int(loc[0]) > int(loc[1]):
            loc = (loc[1], loc[0])
        bashCommand = "samtools faidx "+fasta+" "+chro+":"+str(loc[0])+"-"+str(loc[1])
        #print('Getting fasta for: '+name+' from: '+fasta+" "+chro[3:]+":"+str(loc[0])+"-"+str(loc[1]))
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        open('output/reffasta/'+date_time+'/'+name+".fasta", 'ab').write(output)


def extract_pos(gff, filter=None):
    print("Processing gencode annotation..")
    if filter:
        print('Filtering on: ', filter)
    in_handle = open(gff)
    limit_info = dict(
        gff_id=["chr7", "chr14"],
        gff_type=['gene', 'CDS']
    )
    genlist = set()
    for rec in GFF.parse(in_handle, limit_info=limit_info):
        for feature in rec.features:
            if feature.type == 'inferred_parent':
                for f in feature.sub_features:
                    if containsTR(f.qualifiers['gene_name'][0], filter):
                        if f.type == 'CDS':
                            if containCDSGENE(genlist, f.qualifiers['gene_name'][0], 'gene'):
                                genlist = discarditem(genlist, f.qualifiers['gene_name'][0])
                            genlist.add((f.qualifiers['gene_name'][0], (f.location.start, f.location.end,
                                                                        f.location.strand), f.type, rec.id,
                                         f.qualifiers['exon_number'][0]))
                        elif f.type == 'gene':
                            if not containCDSGENE(genlist, f.qualifiers['gene_name'][0], 'CDS'):
                                genlist.add((f.qualifiers['gene_name'][0], (f.location.start, f.location.end,
                                                                            f.location.strand), f.type, rec.id,
                                             f.qualifiers['exon_number'][0]))
            else:
                if containsTR(feature.qualifiers['gene_name'][0], filter):
                    try:
                        ex_number = feature.qualifiers['exon_number'][0]
                    except:
                        ex_number = 1
                    if feature.type == 'CDS':
                        if containCDSGENE(genlist, feature.qualifiers['gene_name'][0], 'gene'):
                            genlist = discarditem(genlist, feature.qualifiers['gene_name'][0])
                        genlist.add((feature.qualifiers['gene_name'][0], (feature.location.start, feature.location.end,
                                                                          feature.location.strand), feature.type,
                                     rec.id, ex_number))
                    elif feature.type == 'gene':
                        if not containCDSGENE(genlist, feature.qualifiers['gene_name'][0], 'CDS'):
                            genlist.add((feature.qualifiers['gene_name'][0], (feature.location.start,
                                                                              feature.location.end,
                                                                              feature.location.strand), feature.type,
                                         rec.id, ex_number))
    return sorted(genlist)


def containsTR(str, list=None):
    """
    checks if str contains item in list
    :param str:
    :return:
    """
    if not list:
        list = ['TRAV', 'TRAC', 'TRAJ',
                'TRBV', 'TRBC', 'TRBJ', 'TRBC', 'TRBD',
                'TRGV', 'TRGC', 'TRGJ',
                'TRDV', 'TRDC', 'TRDJ', 'TRDD']
    for item in list:
        if item in str:
            return True
    return False


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
    for i in genelist:
        if i.get_name() == name:
            return True
    return False

def main():
    """
    --gff=input gff
    --vcf="inputvcf1;inputvcf2"
    --popfile=inputpopfile
    --ref=ref.fasta
    --filter=filterlist
    --threshold=integer

    gff is a valid annotation file in gencode format
    vcf requires a vcf.gz and requires a vcf.tbi to be present. needs to be two filenames in quotes.
    popfile is a tab seperated file containing sample and populaltion information: SampleID Pop Superpop
    ref is the reference genome in fasta format.
    filter is an optional argument to filter on a set of genes.
    threshold is an integer defining the minimum support a allele must have to be saved in the HAP file.
    :return:
    """
    gff, vcf, populationfile, ref = '', '', '', ''
    filterfile = None
    threshold = '4'
    for arguments in sys.argv:
        if arguments.startswith('--gff'):
            gff = arguments.split('=')[1]
        if arguments.startswith('--vcf'):
            vcf = arguments.split('=')[1]
        if arguments.startswith('--popfile'):
            populationfile = arguments.split('=')[1]
        if arguments.startswith('--ref'):
            ref = arguments.split('=')[1]
        if arguments.startswith('--filter'):
            filterfile = [x.strip() for x in open(arguments.split('=')[1], 'r').readlines()]
        if arguments.startswith('--treshold'):
            threshold = str(arguments.split('=')[1])
        if arguments.startswith('--help'):
            print("""
TR allele pipeline:

--gff=input gff
--vcf="inputvcf1;inputvcf2"
--popfile=inputpopfile
--ref=ref.fasta
--filter=filterlist
--threshold=integer

gff is a valid annotation file in gencode format
vcf requires a vcf.gz and requires a vcf.tbi to be present. needs to be two filenames in quotes.
popfile is a tab seperated file containing sample and populaltion information: SampleID Pop Superpop
ref is the reference genome in fasta format.
filter is an optional argument to filter on a set of genes.
threshold is an integer defining the minimum support a allele must have to be saved in the HAP file.
            """)

    if gff == '' or vcf == '' or len(vcf.split(';')) != 2 or populationfile == '' or ref == '':
        print("Not enough arguments to run.\n\nPlease make sure to specify:\n\t--vcf=\"input1.vcf;input2.vcf\"\n\t--gff=input.gff\n\t--popfile=inputpopfile\n\t--ref=ref.fasta\nIn the command line options.")
    else:
        genepos = extract_pos(gff, filterfile)
        genelist = []
        for pos in genepos:
            if in_genelist(genelist, pos[0]):
                for i in genelist:
                    if i.get_name() == pos[0]:
                        i.add_exon(pos[4], pos[1])
                        break
            else:
                genelist.append(gene.Gene(pos[0], {str(pos[4]): pos[1]}))
        cutVCF(genepos, vcf)
        cutFasta(genepos, ref)
        nopedlist = vcf2ped(genepos, populationfile)
        ped2hap(genelist, populationfile, nopedlist, threshold)
        hap_seq(genelist)
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
        print("Output ID: "+ date_time)
    #print(genepos)


main()

