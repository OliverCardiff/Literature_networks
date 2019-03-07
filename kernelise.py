# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 08:38:07 2019

@author: Oliver
"""

import sys

#These are the default file names used in the example
file = 'pubmed_result.txt'

#You need to have the sleepy children text file in the same folder
#as the script when you execute it
file2 = 'sleepy_children.txt'

#This is the name of the file that will be created when this script is run
outfile = 'frequency_kernel.txt'

#You can include a positional argument which specifies
#name of your main abstracts file if it is different to the default
if len(sys.argv) > 2:
    file = sys.argv[1]
    
minfreq = 200
wordlen = 4

#if you have specified the name of your abstract file, you may also
#add two more positional arguments to control the minimum word length
#and the the minimum frequency of interest
if len(sys.argv) > 3:
    minfreq = sys.argv[2]
    wordlen = sys.argv[3]

#This is a function which decomposes a text file into a dictionary of
#word frequencies. It takes the filename as the only argument
def Decomp(filename):
    #This line opens the file for reading
    abstracts = open(filename, 'r', encoding="utf8")
    #this line instantiates an empty dictionary we can use to store words
    term_freq = {}
    
    #This for loop scans each line of the text file, 'x' being the string
    #variable in which the next line of text is stored
    for x in abstracts:
        
        #Here we replace common punctuation marks with spaces, just so our
        #words don't get them mixed in
        x = x.replace(",", " ")
        x = x.replace(".", " ")
        x = x.replace(":", " ")
        x = x.replace("-", " ")
        
        #This function 'rstrip' removes all trailing newline/whitespace chars
        x = x.rstrip()
        #this function converts all characters to lower case, so we can
        #avoid case based mismatches
        x = x.lower()
        #the split function creates a list of strings that were separated by
        #the specified delimiter, in this case a space
        words = x.split(' ')
        
        #This loop iterates through the list of words
        for w in words:
            #If the word 'w' already exists in the dictionary, increment its
            #value by one
            if w in term_freq:
                term_freq[w] += 1
            #Else, if there is currently no entry in the dictionary, create
            #an entry with a paired value of 1
            else:
                term_freq[w] = 1
    #The function then 'returns' to the caller a dictionary filled with
    #word frequency values
    return term_freq

#Here we run the above function on the two files of abstracts
term_freq = Decomp(file)
ref_freq = Decomp(file2)

#Here we create a new dictionary to store the final filtered values
t_freq2 = {}

#This for loop iterates through all the keys in the dictionary of our
#abstracts of interest, and filters out words which do not match our
#criteria
for x in term_freq:
    #This line checks if the frequency is above the minimum
    #AND that the word is above the minimum length
    #Both of these values can be specified above or as runtime arguments
    if(term_freq[x] > minfreq) & (len(x) > wordlen):
        #If the word passes the above check, then we check to see if it is
        #also present in the irrelevant dictionary
        if x not in ref_freq: 
            #In this case the word was not found
            t_freq2[x] = term_freq[x]
        #This last check evaluates whether, if the word IS present in the
        #irrelevant dictionary BUT at a very low frequency, then we should
        #still include it
        elif ref_freq[x] < minfreq:
            t_freq2[x] = term_freq[x]

#this rather compact line of code sorts the dictionary entries by the
#frequencies of the words used and outputs a list of tuples
sortword = sorted(t_freq2.items(), key=lambda kv: kv[1])

#Here a new file is created with the 'w+' argument to 'open'
fout = open(outfile, 'w+')

#iterate through the list of tuples and write them to the newly created file
for x in sortword:  
    fout.write(str(x[0]) + "\t" + str(x[1]) + "\n")
    
#close the file!
fout.close()