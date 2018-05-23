""" script to update source ids in case multiple sources have the same id
"""
import os, sys

def number_ids(source_file, prefix = ''):
    """Main function
    """
    f_in = open(source_file, 'r')
    new_source_file = source_file[:-4] + '_new_ids.xml'
    f_out = open(new_source_file, 'w')
    index = 1
    for line in f_in.readlines():
        if line.lstrip().startswith(r'id="'):
            new_line = '        id="%s_%06d"\n' % (prefix, index)
            index += 1
            f_out.write(new_line)
        else:
            f_out.write(line)
        
    f_in.close()
    f_out.close()

if __name__ == "__main__":
    pathname = '../zones/2018_mw/AUS6/input/collapsed/'
    source_filename = 'AUS6_collapsed_pts_geom_weighted_merged_parallel.xml'
    source_file = os.path.join(pathname, source_filename)
    number_ids(source_file, prefix = 'AUS6')
