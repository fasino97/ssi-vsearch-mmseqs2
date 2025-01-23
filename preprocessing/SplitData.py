from Bio import SeqIO
import random

#Specify what percentage of data should be training data
percent=0.25

#Open full FASTA file and randomly assign a percentage of the sequences to
#a training FASTA file and the rest to a testing FASTA file
with open("silvaFullSet.fasta") as input_file:
  sequences = list(SeqIO.parse(input_file, "fasta"))
  training_count = int(len(sequences) * percent)
  
#Pick sequences randomly
  random.shuffle(sequences)
  training_sequences = sequences[:training_count]
  testing_sequences = sequences[training_count:]

#Write the files to the folder so they can be used elsewhere
  SeqIO.write(training_sequences, "silvatraining2024.fasta", "fasta")
  SeqIO.write(testing_sequences, "silvatesting2024.fasta", "fasta")
