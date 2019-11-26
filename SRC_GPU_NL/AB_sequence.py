# -*- coding: utf-8 -*-
import sys
#sample:
#python AB_sequence.py 2GB1
#python AB_sequence.py <pdb id of the fasta file>
#OBS> make sure that the file fasta is in the same directory as the program AB_sequence.py


#This function classifie each amino acid in hydrophobic (A) or polar (B) with Alberts classification
def convert(amino):
    if amino == 'G':
        return 'A'
    elif amino == 'A':
        return 'A'
    elif amino == 'L':
        return 'A'
    elif amino == 'M':
        return 'A'
    elif amino == 'F':
        return 'A'
    elif amino == 'W':
        return 'A'
    elif amino == 'K':
        return 'B'
    elif amino == 'Q':
        return 'B'
    elif amino == 'E':
        return 'B'
    elif amino == 'S':
        return 'B'
    elif amino == 'P':
        return 'A'
    elif amino == 'V':
        return 'A'
    elif amino == 'I':
        return 'A'
    elif amino == 'C':
        return 'A'
    elif amino == 'Y':
        return 'B'
    elif amino == 'H':
        return 'B'
    elif amino == 'R':
        return 'B'
    elif amino == 'N':
        return 'B'
    elif amino == 'D':
        return 'B'
    elif amino == 'T':
        return 'B'
    else:
        print 'Error in the amino acid ', amino
        return 0
    
def main(protein_name):
    fasta_name = '../INPUT/'+ protein_name + '.fasta.txt'

    with open(fasta_name, "r") as f:
        sequence = f.readlines()[1]
    f.close()

    print sequence
        
    name = '../INPUT/'+  protein_name.upper() + '.txt'
    sequence_AB = []
    
    sequence = list(sequence)
    tamanho = len(sequence)-1
    
    for i in xrange(0,tamanho):
        if sequence_AB.append(convert(sequence[i]) ) == 0 :
            
            break

    sequence_AB = ''.join(sequence_AB)
    print sequence_AB

    arquivo = open(name, 'w')
    arquivo.write(sequence_AB)
    arquivo.close()

    return 'O.K :)'
    
if __name__ == '__main__':

    main(sys.argv[1].lower())