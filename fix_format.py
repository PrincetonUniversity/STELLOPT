import sys
import subprocess

filename = sys.argv[1]
print('processing '+filename)

old_list = ['[math](math)', r'\\ '[:-1], r'\_ '[:-1], r'\^', r'\|\| \|\|', r'\|\|', r"\`", r"\'", r'\$', r'\*', r'\[', r'\]', r'\<', r'\>']
new_list = ['$$', r'\ '[:-1], r"_ "[:-1], r"^", r'| \n |', '|', "'", "'", '$$', '*', '[', ']', '<', '>']

with open(filename, 'r') as f:
    contents = f.readlines()
    for i, line in enumerate(contents):
        for old, new in zip(old_list, new_list):           
            contents[i] = contents[i].replace(old, new)
            
with open(filename, 'w') as f:
    f.writelines(contents)

exit
