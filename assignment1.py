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


class Assignment1:
    
    def __init__(self, reference_geneome, ucscFile2Write, bam):
        ## Your gene of interest
        # insert gen name
        self.gene = "ETS2"
        print("Query Gene: " + self.gene)

        self.reference = reference_geneome
        self.ucsc = ucscFile2Write
        self.bam = bam

    
    def download_gene_coordinates(self):
        # Only tests if file exists in current working directory
        DataFiles = os.listdir('.')
        if not self.ucsc in DataFiles:
            print("Connecting to UCSC to fetch data")

            cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
            ## Get cursor
            cursor = cnx.cursor()

            ## Build query fields
            query_fields = ["refGene.name2",
                            "refGene.name",
                            "refGene.chrom",
                            "refGene.txStart",
                            "refGene.txEnd",
                            "refGene.strand",
                            "refGene.exonCount",
                            "refGene.exonStarts",
                            "refGene.exonEnds"]
            ## Build query
            query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields)

            ## Execute query
            cursor.execute(query)

            ## Write to file
            ## TODO this may need some work
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
                    # could be done with regex
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
        currentFile.close()

        return header

    def get_properly_paired_reads_of_gene(self):
        currentFile = pysam.AlignmentFile(self.bam, "rb")
        coordinates = self.get_coordinates_of_gene(self.ucsc, self.gene)

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

        averageCoverageOverGene = coverage / totalRefLength
        print("average genome coverage: %s \ntotal genome length: %s \n" % (averageCoverageOverGene,  totalRefLength))

        
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
        print("average gene coverage: %s \ngenelength: %s \nbases in gene covered: %s " % (averageCoverageOverGene, genelength, positionsInGeneCovered))
        print("UNFINISHED")


    def get_number_mapped_reads(self):
        currentFile = pysam.AlignmentFile(self.bam, "rb")
        mappedReads = currentFile.mapped
        currentFile.close()
        return mappedReads

    def get_region_of_gene(self, ucsc_file, gene_name):
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
        #self.download_gene_coordinates()
        #coordinates= self.get_coordinates_of_gene()
        #print("Coordinates: ", coordinates)
        #header = self.get_sam_header()

        #mapped = self.get_number_mapped_reads()

        #print("Number of mapped reads:", mapped)
        #print("Number of Exons: ", self.get_number_of_exons())

        # change to count
        #print("Number of properly paired reads :", self.get_properly_paired_reads_of_gene())
        #print("Gene Symbol: ", self.get_gene_symbol())
        #print("Exon Coordinates: ", self.get_region_of_gene()) # exon coordinates


        self.calculate_gene_average_coverage()
        self.calculate_total_average_coverage()

        #nReadsWithIndels = self.get_gene_reads_with_indels(inputBam='../data/chr21.bam') # get counts
        #print("Number of Reads with Indels: " , nReadsWithIndels)
    
    
def main():
    print("Assignment 1")
    assignment1 = Assignment1(currentGenome, currentUCSC, currentBam)

    assignment1.print_summary()

    print("Done with assignment 1")
    
        
if __name__ == '__main__':
    main()
    
    
