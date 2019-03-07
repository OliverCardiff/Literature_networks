# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 10:40:41 2019

@author: Oliver
"""
import sys
import pandas as pd

#These are the default names of files expected
file1 = 'pubmed_result.txt'
file2 = 'gene_list.txt'

#These are the defaul names of the files created by this script. Change these
#if you are running the script multiple times and don't want to overwrite
#your previous results!
outAttributes = 'Attribute_file.txt'
outNetwork = 'Network_file.txt'

#This variable controls the output filtering based on the occurence of a given
#term. This will limit the output to all terms with a frequency above this
#number
node_freq_filter = 10
#similarly to above, all connections in the network which are evidenced less
#then this number will be excluded from the output
edge_freq_filter = 5

#if the user wishes, they can pass positional arguments to the function
#to specify the names of the abstracts and the gene list
if len(sys.argv) > 1:
    file1 = sys.argv[1]
if len(sys.argv) > 2:
    file2 = sys.argv[2]
    
def DecompMini(lines):
    term_freq = {}
    #This for loop scans each line of the text file, 'x' being the string
    #variable in which the next line of text is stored
    for x in lines:
        #Here we replace common punctuation marks with spaces, just so our
        #words don't get them mixed in
        x = x.replace(",", " ")
        x = x.replace(".", " ")
        x = x.replace(":", " ")
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
                
    return term_freq

#This function converts a file into a list of strings, one for each line
def ToLineList(abstracts):
    #creates an empty list to store the strings
    lines = []
    #iterates through each line of the file
    for x in abstracts:
        #adds the line, as a string, to the list
        lines.append(x)
    #returns the list of strings to the function caller
    return lines


#This function finds the frequencies of phrases in a set of lines, and adds
#the frequencies to a pre-existing frequency dictionary
def DecompPhrases(fdict, phrases, lines):
    #merge all lines to search into a single string for efficiency
    all_lines = ' '.join(lines)
    #then for each phrase to search..
    for p in phrases:
        #check if the phrase is a substring within the line
        if p in all_lines:
            frq = all_lines.count(p)
            if p in fdict:
                fdict[p] += frq
            else:
                fdict[p] = frq
        

#This is a function which decomposes a text file into a dictionary of
#word frequencies. It takes the filename as the only argument
def Decomp(filename, phrases = []):
    #This line opens the file for reading
    abstracts = open(filename, 'r', encoding="utf8")
    #this line instantiates an empty dictionary we can use to store words
    term_freq = DecompMini(abstracts)
    #close the abstracts file
    abstracts.close()
    #if there are any phrases to obtain frequencies for..
    if len(phrases) > 0:
        #re-open the abstracts file
        abstracts = open(filename, 'r', encoding="utf8")
        #convert the file into a list of strings to search
        lines = ToLineList(abstracts)
        #then run the phrase decomposition function
        DecompPhrases(term_freq, phrases, lines)
    #The function then 'returns' to the caller a dictionary filled with
    #word frequency values
    return term_freq

#This function filters a dictionary of word frequencies based on a prior
#set of terms
def FilterDict(terms, miniD):
    #A list of 'keys' (words) from the dictionary is extracted
    kys = pd.Series(list(miniD.keys()))
    #the list of keys is filtered for entries which match terms in the
    # 'terms' list, using the pandas 'isin' function 
    fkys = kys.loc[kys.isin(terms)]
    #A new dictionary is created to store the filtered info
    dict2 = {}
    
    #assign an entry in the new dictionary for every non-filtered key
    for k in fkys:
        dict2[k] = 1
    
    #return the new filtered dictionary to the function caller
    return dict2

#This function merges a dictionary of word frequencies 'miniD' from a single 
#abstract with the network dictionary 'mainD', creating connections between
#all entries found in the same abstract frequency dictionary
def MergeToNetwork(miniD, mainD):
    #first a list of terms is extracted from the abstract dict, as keys
    kys = list(miniD.keys())
    
    #This nested for loop iterates over the set of terms in the the mini
    #dictionary. for each term it creates an entry in the network between
    #that term and all subsequent terms in the single-abstract list
    for i in range(0, (len(kys) - 1)):
        for j in range(i+1, len(kys)):
            
            #this bit of code ensures that network entries aren't duplicated
            #by always arranging node-pairs by their lexicographical order
            llw = kys[i]
            mlw = kys[j]
            if llw > mlw:
                mlw = llw
                llw = kys[j]
            
            #if the least lexicographical word is not in the network as a
            #primary key, add it
            if llw not in mainD:
                mainD[llw] = {}
            #if the other term is not a sub-key of the first, add it
            if mlw not in mainD[llw]:
                mainD[llw][mlw] = 1
            #if the connection already exists in the network, simply increment
            #its frequency
            else:
                mainD[llw][mlw] += 1
    #This function does not need to return anything because the network
    #dictionary is passed by reference as a default in python, meaning that 
    #the original version of the dictionary supplied to the function will
    #be modified in memory instead of a copy being modified

#This function simply combines the above three functions into a process which
#turns a list of strings representing an abstract into a set of network entries
def DecompFilterMerge(terms, lines, net_dict):
    #If the list of lines is at least 5 long (this is a check against) a
    #poorly formatted abstract or other erronous file reading!
    if len(lines) > 5:
        #Frequency decompose the abstract
        mini_dict = DecompMini(lines)
        #Filter for only the relevant terms
        mini_dict = FilterDict(terms[0], mini_dict)
        #Check if there are phrases to be searched for as well
        if len(terms[1]) > 0:
            #if so, perform frequency searches for the set phrase list
            DecompPhrases(mini_dict, terms[1], lines)
        #If relevant terms are found, and there is more than one of them..
        if len(mini_dict) > 1:
            #...then merge those entries with the current network
            MergeToNetwork(mini_dict, net_dict)
    #Similar to the above function, this function does not need to return
    #anything as all dictionaries are passed by reference

#This function simply converts the network dictionary to a tabluar format
def PrintNetwork(netD, freqs):
    #Create a file for the network output
    fout = open(outNetwork, 'w+')
    #write a header for the three columns
    fout.write("Node1\tConnections\tNode2\n")
    #nested for loop which iterates through all the primary keys in the
    #network dictionary, and they iterates through all the keys of each
    #sub-dictionary, writing their associated value (the connection freq.)
    #to the network file
    for k1 in netD:
        for k2 in netD[k1]:
            frq = netD[k1][k2]
            f1 = freqs[k1]
            f2 = freqs[k2]
            #make sure both node frequencies pass the filter threshold
            if (f1 > node_freq_filter) & (f2 > node_freq_filter):
                #only output connections which occur more than the filter level
                if frq > edge_freq_filter:
                    fout.write(str(k1) + "\t" + str(frq) + "\t" + str(k2) + "\n")
    #close the file for writing
    fout.close()

#This script appends the frequencies of the words to the annotated gene symbol
#list that you created
def PrintAttributes(terms, freqs):
    #create a file for node arrtibutes
    fout = open(outAttributes, 'w+')
    #write a header to the node attribute file
    fout.write("Node\tAnnotation_Type\tFrequency\n")
    
    #iterate through all the input single word terms
    for i in range(0,len(terms[0])):
        #if the term was found in the abstracts, output it, with the frequency
        if terms[0][i] in freqs:
            frq = freqs[terms[0][i]]
            #check that the node is high-enough frequency to include
            if frq > node_freq_filter:
                fout.write(str(terms[0][i]) + "\t" + str(terms[2][i]) + "\t" + str(frq) + "\n")
            
    #iterate through all the input phrases
    for i in range(0,len(terms[1])):
        #if the phrase was found in the abstracts, output it, with the frequency
        if terms[1][i] in freqs:
            frq = freqs[terms[1][i]]
            #check that the node is high-enough frequency to include
            if frq > node_freq_filter:
                fout.write(str(terms[1][i]) + "\t" + str(terms[3][i]) + "\t" + str(frq) + "\n")
    #closes the file
    fout.close()

#This function separates a list of search terms into single words and phrases
#This is needed because we want to search for the phrases differenently to the
#single word terms.
def FindPhrases(searchtms):
    #Extract a list of search terms from the discrete annotation gene list
    trms = searchtms.iloc[:,0].tolist()
    #extract a list of discrete annotations from the gene list
    annots = searchtms.iloc[:,1].tolist()
    
    #Here we prepare four empty lists to store the segregated terms and phrases
    singletons = []
    single_annots = []
    phrases = []
    phrase_annots = []
    
    #Iterate through all input terms
    for i in range(0, len(trms)):
        #if a whitespace exists in the term, classify it as a phrase
        if ' ' in trms[i]:
            #append the phrase to the phrase list
            phrases.append(trms[i].lower())
            #append the phrase annotation to its own list
            phrase_annots.append(annots[i])
        else:
            #otherwise append single words to the singletons list
            singletons.append(trms[i].lower())
            #append the singleton annotations to their own list
            single_annots.append(annots[i])
    #this return creates a single object containing the four lists created
    return [singletons, phrases, single_annots, phrase_annots]
    
'''
vvvvv Execution of the defined functions occurs below vvvvv
'''

#Read in the gene symbol and terminology file you created
search_terms = pd.read_csv(file2, delimiter="\t")
#divide term list into a list of lists of phrases and words
processed_terms = FindPhrases(search_terms)
#create an empty dictionary used to store the network dictionary
network_dict = {}
#a reference to the list used to store lines of text from an abstract 
abs_lines = []
#Decompose the abstracts file to get the term frequencies
freqs = Decomp(file1, processed_terms[1])
#Re-open the abstracts file to iterate through
fin = open(file1, 'r', encoding='utf8')
#for each abstract in the file, gather its lines into a list, to frequency
#decompose, filter, and add to the nework
for x in fin:    
    test = x.split('.')
    #The only lines in the abstract file with more than four periods are the
    #title lines. We use this information to reset the collection of abstract
    #lines and process them in batches
    if (len(test) > 4) & (len(test[0]) < 6):
        #The main function is run on the set of 
        DecompFilterMerge(processed_terms, abs_lines, network_dict)
        #reset the list of lines to empty
        abs_lines = []
    
    #each line is added to the list of abstract lines
    abs_lines.append(x)

#Run the main function a final time to process the last abstract
DecompFilterMerge(processed_terms, abs_lines, network_dict)

#Print the network file for Cytoscape
PrintNetwork(network_dict, freqs)
#Print a node attribute file for Cytoscape
PrintAttributes(processed_terms, freqs)
