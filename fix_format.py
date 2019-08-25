
import sys

filename = sys.argv[1]
print(filename)

a =  r'\\ '[:-1]
b =  r'\ '[:-1]

old_list = ['[math](math)', a, r'\_ '[:-1], r'\^']
new_list = ['$$', b, r"_ "[:-1], r"^"]


with open(filename, 'r') as f:
    contents = f.readlines()
    for i, line in enumerate(contents):
        for old, new in zip(old_list, new_list):
            contents[i] = line.replace(old, new)

with open(filename, 'w') as f:
    f.writelines(contents)

exit
