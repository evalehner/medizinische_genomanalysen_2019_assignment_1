import mysql.connector
import os
import pysam
import subprocess
import shlex


__author__ = 'Eva V. Lehner'

##
## Concept:
## TODO
##

currentGenome = "hg38"
currentUCSC = "coordinates"
currentBam = "../data/chr21.bam"
currentGenename = "ETS2"


class Assignment1:
    
    def __init__(self, reference_genome, ucscFile2Write, bam, genename):
        self.gene = genename
        self.reference = reference_genome
        self.ucsc = ucscFile2Write
        self.bam = bam

    def download_gene_coordinates(self):
        # Only tests if file exists in current working directory
        DataFiles = os.listdir('.')
        if not self.ucsc in DataFiles:
            print("Connecting to UCSC to fetch data")
            cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=self.reference) # connect to DB
            cursor = cnx.cursor() # set cursor

            ## Build query fields
            query_fields = ["refGene.name2", "refGene.name", "refGene.chrom", "refGene.txStart",
                            "refGene.txEnd", "refGene.strand", "refGene.exonCount", "refGene.exonStarts",
                            "refGene.exonEnds"]
            ## Build query
            query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)
            ## Execute query
            cursor.execute(query)

            ## Write result to file, this takes time
            with open(self.ucsc, "w") as fh:
                for row in cursor:
                    fh.write(str(row) + "\n")

            ## Close cursor & connection
            cursor.close()
            cnx.close()

            print("Done fetching data")

        else:
            print("Coordinates were already downloaded")
        
    def get_coordinates_of_gene(self):
        gene_coordinates = []
        with open(self.ucsc, "r") as fh:
            for line in fh.readlines():
                if line.find(self.gene) != -1:
                    # could be simplyfied
                    line = line.replace("'", "")
                    line = line.replace("(", "")
                    line = line.replace(")", "")
                    line = line.replace(",", "")
                    line = line.split()
                    gene_coordinates = line[2:5]
                    break # take only first entry in file

        return gene_coordinates

    def get_gene_symbol(self):
        return self.gene

    def get_sam_header(self):
        currentFile = pysam.AlignmentFile(self.bam, "rb")
        header = currentFile.header
        print(type(header))
        currentFile.close()
        return header

    def get_properly_paired_reads_of_gene(self):
        currentFile = pysam.AlignmentFile(self.bam, "rb")
        coordinates = self.get_coordinates_of_gene()

        paired_read_count = 0
        for read in currentFile.fetch(contig=coordinates[0], start=int(coordinates[1]), stop=int(coordinates[2])):
            if read.is_paired:
                paired_read_count = paired_read_count + 1
        currentFile.close()

        return paired_read_count

    def get_gene_reads_with_indels(self):
        # extract reads where there is I od D in cigar string (col6)
        # memory inefficient :(

        view = shlex.split("samtools view "+ self.bam)
        callView = subprocess.Popen(view, stdout=subprocess.PIPE)
        cut = shlex.split("cut -f 6")
        callCut = subprocess.Popen(cut, stdin=callView.stdout, stdout=subprocess.PIPE)
        countCigar = shlex.split("grep -c -E 'I|D'")
        nInDels = subprocess.Popen(countCigar, stdin=callCut.stdout, stdout=subprocess.PIPE)
        out = nInDels.stdout.read()

        return out

    def calculate_total_average_coverage(self):
        # Coordinates in pysam are 0 based!
        # Coordinates of UCSC File are 0 based!
        # --> No conversion conflicts or sth!
        header = self.get_sam_header()
        totalRefLength=0
        for item in header.references:
            totalRefLength += header.get_reference_length(item)

        samfile = pysam.AlignmentFile(self.bam, "rb")
        coverage = 0

        for pileupcolumn in samfile.pileup():
            coverage += pileupcolumn.n
        samfile.close()

        averageCoverageOverGenome = coverage / totalRefLength
        return round(averageCoverageOverGenome, 2)
        
    def calculate_gene_average_coverage(self):
        # Coordinates in pysam are 0 based!
        # Coordinates of UCSC File are 0 based!
        # --> No conversion conflicts or sth!
        samfile = pysam.AlignmentFile(self.bam, "rb")
        coordinates = self.get_coordinates_of_gene()
        chrom = coordinates[0]

        coverage = 0
        positionsInGeneCovered = 0
        for pileupcolumn in samfile.pileup(chrom, start = int(coordinates[1]) , stop = int(coordinates[2])):
            coverage += pileupcolumn.n
            positionsInGeneCovered +=1
        samfile.close()

        genelength = int(coordinates[2]) - int(coordinates[1])
        averageCoverageOverGene = coverage / genelength
        return round(averageCoverageOverGene, 2)

    def get_number_mapped_reads(self):
        currentFile = pysam.AlignmentFile(self.bam, "rb")
        mappedReads = currentFile.mapped
        currentFile.close()
        return mappedReads

    def get_region_of_gene(self):
        exon_coordinates = []
        with open(self.ucsc, "r") as fh:
            for line in fh.readlines():
                if line.find(self.gene) != -1:
                    line = line.split(',')
                    exon_coordinates.append(line[7:-1])
                    break # If more than one entries of gene coordinates exist in file, process only first one.

        # output could be nicer formated
        return exon_coordinates

        
    def get_number_of_exons(self):
        # Gene is contained twice in ucsc file. Here, only exon count of first record is included.
        exon_coordinates = []
        with open(self.ucsc, "r") as fh:
            for line in fh.readlines():
                if line.find(self.gene) != -1:
                    line = line.split(',')
                    exon_coordinates.append(line[7:-1])

        exonCount = 0
        # use only first record
        exon_coordinates = exon_coordinates[0]
        #turn list around and itterate over it till 'b' is found
        for item in reversed(exon_coordinates):
            if item.find('b') == -1:
                exonCount = exonCount + 1
            else:
                exonCount = exonCount + 1
                break
        return exonCount

    def print_summary(self):
        self.download_gene_coordinates()

        print("Coordinates of Gene %s: \n \t %s " % (self.gene, self.get_coordinates_of_gene()))
        print("Number of mapped reads: \t%s " % (self.get_number_mapped_reads()))
        print("Number of Exons: \t %s  "% (self.get_number_of_exons()))
        print("Number of properly paired reads: \t %s  " % (self.get_properly_paired_reads_of_gene()))
        print("Gene Symbol: \t %s " % (self.get_gene_symbol()))
        print("Exon Coordinates:  \n \t %s" % (self.get_region_of_gene())) # exon coordinates

        print("Average Coverage over gene: \t%s  " % (self.calculate_gene_average_coverage()))
        print("Average Coverage over genome: \t%s  " % (self.calculate_total_average_coverage()))

        print("Number of Alignments with Indesl:\t%s" %(self.get_gene_reads_with_indels()))


        # Fix this
        header = self.get_sam_header()
        print("Header of file %s" % (self.bam))
        # for line in header.iteritems():
        #    print(line)
def main():
    print("Assignment 1")
    assignment1 = Assignment1(currentGenome, currentUCSC, currentBam, currentGenename)
    assignment1.print_summary()
    print("Done with assignment 1")
    
        
if __name__ == '__main__':
    main()
    
    
