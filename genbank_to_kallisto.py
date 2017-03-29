#python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

import csv

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("genbankfile")
parser.add_argument("output_transcripts")
parser.add_argument("output_table")


args = parser.parse_args()
print("You provided: " + args.genbankfile)
print("Output will be: " + args.output_transcripts + " and " + args.output_table)

def genbank_to_kallisto(file):

    #list of keys
    valid_types = ["CDS","rRNA","tmRNA","tRNA","misc_RNA"]
    list_transcripts = []

    fp = open(args.output_table, "w")
    tsv = csv.writer(fp, delimiter = '\t')

    for seq_record in SeqIO.parse(file, "genbank"):
    #there could be multiple seq records if lots of contigs

        print("now processing seq record " + seq_record.name)
        for feature in seq_record.features: #these are now SeqFeature objects
            if feature.type in valid_types:
                print("now at ", feature.type, feature.location)
                transcript = feature.extract(seq_record)  #extract from seq_record, rtns SeqObject
                #change transcript name
                transcript.id = get_qualifier(feature, "locus_tag")
                #print("id is " + transcript.id)
                #change transcript description
                if "product" not in feature.qualifiers:
                    transcript.description = feature.qualifiers["note"][0]

                else:
                    transcript.description = feature.qualifiers["product"][0]
                #print("description is : " + transcript.description)

                list_transcripts.append(transcript)

                tsv.writerow([get_qualifier(feature, "locus_tag"), feature.type, get_qualifier(feature, "gene"), get_qualifier(feature, "EC_number"), get_qualifier(feature, "product")])


    SeqIO.write(list_transcripts, args.output_transcripts, "fasta")

    print("finished")


def get_qualifier(feature, qualifier):
    if qualifier not in feature.qualifiers:
        return ""
    else:
        return feature.qualifiers[qualifier][0]

genbank_to_kallisto(args.genbankfile)
