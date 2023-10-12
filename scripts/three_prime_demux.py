import sys
import os

# All possible three prime barcodes
barcode_data = {
    "ATCACG": "1.fastq",
    "CGATGT": "2.fastq",
    "TTAGGC": "3.fastq",
    "TGACCA": "4.fastq",
    "ACAGTG": "5.fastq",
    "GCCAAT": "6.fastq",
    "CAGATC": "7.fastq",
    "ACTTGA": "8.fastq",
    "GATCAG": "9.fastq",
    "TAGCTT": "10.fastq",
    "GGCTAC": "11.fastq",
    "CTTGTA": "12.fastq",
    "AGTTCC": "14.fastq",
    "ATGTCA": "15.fastq",
    "GTCCGC": "18.fastq",
    "GTTTCG": "21.fastq",
    "CGTACG": "22.fastq",
    "GAGTGG": "23.fastq",
    "GGTAGC": "24.fastq",
    "ACTGAT": "25.fastq"
}

if len(sys.argv) != 3:
    print("Usage: python three_prime_demux.py five_prime_trim.fastq output_directory")
    sys.exit(1)

fastq_path = sys.argv[1]
output_directory = sys.argv[2]

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)

fastq = open(fastq_path)

# Associate barcodes and the output file names
fileDict = {}
barcodeDict = {}
for barcode in barcode_data:
    fileDict[barcode] = os.path.join(output_directory, barcode_data[barcode])
    barcodeDict[barcode] = 0

files = fileDict.values()
files = set(files)
fileOpen = {}
for file in files:
    fileOpen[file] = open(file, 'w') 

barcodes = barcodeDict.keys()

counter = 1

for line in fastq:
    if counter == 1:
        ID = line.strip()
        counter = 2

    elif counter == 2:
        seq = line.strip()
        counter = 3

    elif counter == 3:
        strand = line.strip()
        counter = 4
                
    elif counter == 4:
        score = line.strip()
        counter = 5
                
    if counter == 5:
        # Obtain LAST 6 base-pairs of read
        seq6bp = seq[-6:]
                
        # Obtain number of matches for each barcode
        for barcode in barcodes:
            matchCount = 0

            for i in range(6):
                if seq6bp[i] == barcode[i]:
                    matchCount += 1

            if matchCount >= 5:
                fileOpen[fileDict[barcode]].write(ID + '\n' + seq + '\n' + strand + '\n' + score + '\n')
        counter = 1

sys.exit(0)

