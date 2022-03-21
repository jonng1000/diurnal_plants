# -*- coding: utf-8 -*-
"""
Spyder Editor

Editing sequence IDs as spe aa file has gene desc which should not be
the case.
"""

import csv

orig = 'SequenceIDs.txt'
new = 'SequenceIDs_edit.txt'

with open(orig, newline='') as origf, open(new, 'w', newline='') as newf:
    reader = csv.reader(origf, delimiter=' ')
    writer = csv.writer(newf, delimiter=' ')
    for row in reader:
        if len(row) == 2:
            writer.writerow(row)
        elif len(row) > 2:
            writer.writerow([row[0],row[1]])
        else:
            raise ValueError('Row length doesnt make sense')