'''
Title:   Manipulate Files

Author:  Gubbin Eel (Samuel Lindsay, East Carolina University)

Purpose: Defines functions used to read, write, and modify files, 
         generally.

NOTE:    Intended improvements are denoted by the #%# flag.
         Searching for that flag will identify code that needs to be
         updated or will later be improved.
'''

whitespace = "\t\n\x0b\x0c\r "

# Return the contents of a file as a list of strings.
''' Purpose:   Return the contents of a file as a list of strings.
	Arguments: (1) The name of a file with its full path.
'''
def import_file_contents(file):
    with open(file, 'r') as f:
        contents = f.readlines()
    return contents

''' Purpose: Write a list of strings to a file.
    Arguments: (1) The name of a file with its full path.
               (2) A list containing strings representing the lines
                   of each file to be written.
'''
def write_file(file, file_list):
	with open(file, 'w') as f:
		for line in file_list:
			f.writelines(line)
