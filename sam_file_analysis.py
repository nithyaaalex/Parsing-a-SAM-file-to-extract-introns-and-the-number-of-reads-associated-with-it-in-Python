#IMPORTING MODULES/LIBRARIES
import sys
import re
#MAIN (equivalent) starts from line 110


#TAKING FILES AS INPUT FROM COMMAND LINE

#using sys.argv to take user input from command line, first file is sam format, second is text
sam_file = sys.argv[1] 

try: #checking if the sam file name is valid
    file = open(sam_file)
except FileNotFoundError:
    sys.exit(f'File {sam_file} was not found. Please check the sam file name and try again!')

gene_file = sys.argv[2] 

try: #checking if the text file name is valid
    file = open(gene_file)
except FileNotFoundError:
    sys.exit(f'File {gene_file} was not found. Please check the text file name and try again!')



#FUNCTION DEFINITIONS 

#function for parsing cigar string and splitting it, returns a match as a list of all matches as tuples 
def cigar_parse(cigar_string):
    #decided to use findall instead of finditer because finditer returns a iterable object which I wasn't able check using assert
    matches = re.findall(r"(\d+)([MIDSN])", cigar_string)

    #checking if the cigar string is valid and doesn't contain anything other than the letters MIDSN and numbers 0-9
    try:
        assert re.findall(r"[^MIDSN\d]", cigar_string) == []
    except AssertionError:
        print("Invalid cigar string, has letter other than MIDSN please check your sam file.")
    
    return matches

#function for parsing the match object from cigar_parse and determining location of junction, returns list of tuples of chromosome name and the location
def find_junction(chromosome_name, match_object, start_position):
    start = int(start_position)
    end = 0
    list_of_junctions = []

    #iterating through list from cigar parse one tuple at a time to evaluate each part of cigar string
    for match in match_object:
        #unpacking the tuple into the letter and number and making the number an integer
        number_associated, last_letter = match
        number_associated = int(number_associated)
            
        #using flow control to check what kind of alignment the read has to the reference
        if last_letter == "M" or last_letter == "D": #moving start location for junction whenever matches or deletions are found
            start += number_associated
        elif last_letter == "I" or last_letter == "S": #ignoring insertions or splits in the read
            continue 
        elif last_letter == "N": #calculating the end of the junction, storing it in the list, and continue to look for extra introns
            end = start + number_associated 
            list_of_junctions.append((chromosome_name, start, end))
            start = start + number_associated

    return list_of_junctions

#function for extracting gene location from the gene file, it returns a tuple of the chromosome it's on and the start and end location of the gene
def find_gene_location(third_column_in_file):      
        match_location = re.search(r"(\:.+\.)(\..+\()", third_column_in_file) #this match finds the start and end (between the special characters)
        match_chromosome = re.search(r"(^.+\:)", third_column_in_file) #this match finds the chromosome name

        #saving the chromosome name and removing the colon at the end
        chromosome = match_chromosome.group()[0:len(match_chromosome.group())-1] 
        
        #extracting start and end from the location and removing the special characters 
        start,end = match_location.groups()

        #first it removes the characters at index 0 and the end of the string
        start = start[1:len(start)-1]
        end = end[1:len(end)-1]

        #then it removes the commas between the numbers and turns the string into int
        start = int(start.replace(",",""))
        end = int(end.replace(",",""))
        return chromosome, start, end



#ASSERTIONS FOR FUNCTIONS

#checking if cigar parse splits the cigar string appropriately
assert cigar_parse("63M146N37M") == [('63', 'M'), ('146', 'N'), ('37', 'M')]
assert cigar_parse("32434243") == [] #empty match object 
assert cigar_parse("MINSD") == []

#checking if find junction gives appropriate output for matches, insertions, deletions, soft clips and skipped regions (aka junctions)
assert find_junction("TGME49_chrVIII", [('63', 'M'), ('146', 'N'), ('37', 'M')], 3000000) == [("TGME49_chrVIII", 3000063,3000209)]
assert find_junction("TGME49_chrVIII", [('63', 'M'), ('5', 'I'), ('146', 'N'), ('37', 'M')], 3000000) == [("TGME49_chrVIII", 3000063,3000209)]
assert find_junction("TGME49_chrVIII", [('63', 'M'), ('5', 'S'), ('146', 'N'), ('37', 'M')], 3000000) == [("TGME49_chrVIII", 3000063,3000209)]
assert find_junction("TGME49_chrVIII", [('63', 'M'), ('5', 'D'), ('146', 'N'), ('37', 'M')], 3000000) == [("TGME49_chrVIII", 3000068,3000214)]

#checking if find junction finds multiple introns or no introns
assert find_junction("TGME49_chrVIII", [('63', 'M'), ('146', 'N'), ('37', 'M'), ('23', 'N'), ('34', 'M')], 3000000) == [("TGME49_chrVIII", 3000063,3000209),('TGME49_chrVIII', 3000246, 3000269)]
assert find_junction("TGME49_chrVIII", [('222', 'M')], 3000000) == []

#checking if find gene location splits the third column appropriately
assert find_gene_location("TGME49_chrVIII:3,000,000..3,000,900(+)") == ("TGME49_chrVIII", 3000000, 3000900)
assert find_gene_location("TGME49_chrVIII:6,631,349..6,636,865(+)") == ("TGME49_chrVIII", 6631349, 6636865) 



#CODE FOR OPENING SAM FILE AND STORING JUNCTIONS AND READS

#dictionary with key = chromosome and location of junction(as a tuple), value = number of reads associated with it (int)
junction_dict = {} 

#opening sam file & iterating line by line
with open(sam_file) as my_input_file: 
    for line in my_input_file: 
        #for skipping the headers of the sam file
        if line[0] == "@":
            next(my_input_file)
            continue

        #removing newline character and spaces & storing each column as element in the list 
        line = line.rstrip("\n").strip(" ")
        list_of_all = line.rsplit("\t") 
    
        #variables for saving the last column and columns 3,4 & 6 respectively (-1 to the actual column number because of zero indexing)
        read_align_check = list_of_all[len(list_of_all)-1]  
        rname = list_of_all[2]
        pos = list_of_all[3]
        cigar = list_of_all[5]

        #making sure that the last column was retrieved properly by minusing 1 to the length of line
        try:
            assert read_align_check[0:len(read_align_check)-1] == "NH:i:"
        except AssertionError:
            print("Last column isn't in the format NH:i:x, please check your file and try again")

    
        #deleting the line to save memory, the other columns are not needed
        del list_of_all 

        if read_align_check == "NH:i:1": #checks if end of "NH:i:x" is 1
            #checking if the line has a junction
            check_nmatch = re.findall(r"N",cigar)
            if check_nmatch:
                #parsing and finding junction from the cigar string
                all_matches = cigar_parse(cigar)
                temp = (find_junction(rname, all_matches, pos)) #temporary variable for control statements

                #making sure that the junction found is not empty
                try:
                    assert temp != []
                except AssertionError:
                    print("Invalid junction found")

                for junction in temp:
                    #storing the junction in the junction dictionary
                    if junction in junction_dict:
                        junction_dict[junction] += 1
                    
                    elif junction not in junction_dict:
                        #if it already exists as a key: add a count to the read, if it doesn't exist: add it to the dictionary 
                        junction_dict[junction] = 1



#CODE FOR OPENING GENE FILE AND OUTPUT FILE TO WRITE EACH GENE AND THE JUNCTIONS ASSOCIATED WITH IT

#writing the header for my output file
with open("2933044.txt", 'w') as my_output_file:
    my_output_file.write(f"Gene ID\tJunction Start\tJunction End\tNumber of reads supporting the junction\n")

#opening the gene file
with open(gene_file) as gene_input:
    next(gene_input) #skipping header

    #opening output file in append mode
    with open("2933044.txt", 'a') as my_output_file: 
        #iterating over each line in the gene file
        for line in gene_input: 
            #removing the spaces and newline character and assigning each column to a variable
            line = line.rstrip(" ").rstrip("\n")
            gene_ID, transcript_ID, genomic_location = line.rsplit("\t")

            #finding location from the third column of the gene file
            location_data = find_gene_location(genomic_location)

            #creating a flag for keeping track of when the gene actually has junctions, and is written to the file
            flag = 0

            #iterating over the junction dictionary for each gene
            for junction_data,number_of_reads in junction_dict.items():
                #unpacking both tuples from the junction dictionary and from the find gene location data
                chromosome_name, start_junction, end_junction = junction_data
                chromosome_ref, start_gene, end_gene = location_data
                
                #if the junction is lying between the gene coordinates and if it's on the same chromosome, it writes this data out to a output file
                if start_junction>=start_gene and end_junction<=end_gene and chromosome_name == chromosome_ref:
                    my_output_file.write(f"{gene_ID}\t{start_junction}\t{end_junction}\t{number_of_reads}\n") 
                    flag = 1 #setting the flag to one since we saved this gene

            #if the gene was written to the file, we leave a line before the for loop goes to the next gene in gene file
            if flag ==1: 
                #without using this flag if a gene didn't have any junctions, there were multiple newlines in the output
                my_output_file.write("\n") 