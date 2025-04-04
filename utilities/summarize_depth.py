import os

# Author: Julius Falck
# change directory to the bams_indels_realigned folder in the results

directory= ""

dir_list = os.listdir(directory)

depth_table = "depth_table.txt"


for f in dir_list:
    if f.endswith('dpstats.txt'):
        with open(directory + '/' + f, 'r') as file:
            line = file.readlines()[0]
            with open(depth_table, 'a') as file:
                file.write(f.split('.')[0] + ' ' + line + '\n')
            file.close()
        file.close()

