import csv
import numpy as np

"""
Write out a file in place/lon/lat for place names
extracted from the
"""

# list of Australian in column 1, state in column 2
cities_file = "nsha_scenario_locations_v2.csv"
# name of outfile
outfile = "New_localities_NSHA.csv"

def read_loc_file():
    simple_locs = []
    with open("cities1000_au.txt", 'r', encoding='utf8') as locfile:
        locations_all = locfile.readlines()
        for line in locations_all:
            line = line.split("\t")
            if line[-2] == "Australia/Sydney":
                state = "NSW"
            elif line[-2] == "Australia/Brisbane":
                state = "QLD"
            elif line[-2] == "Australia/Melbourne":
                state = "VIC"
            elif line[-2] == "Australia/Adelaide":
                state = "SA"
            elif line[-2] == "Australia/Perth":
                state = "WA"
            elif line[-2] == "Australia/Hobart":
                state = "TAS"
            elif line[-2] == "Australia/Darwin":
                state = "NT"
            else:
                state = "Unknown"
            line_new = (line[1], line[5], line[4], state)
            simple_locs.append(line_new)
    return simple_locs


# output simplified file with name/lon/lat
simple_locs = read_loc_file()

complete_file = []
outF = open(outfile, "w")

with open(cities_file) as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',')
    for row1 in csvreader:
        place = row1[0]
        for row2 in simple_locs:
            #print(row2[0], place) # Check logic
            if row2[0] == place and row2[3] == row1[1]:
                new_line = row2[1] + ",  " + row2[2] + ",  " + place
                outF.write(new_line)
                outF.write("\n")
outF.close()


