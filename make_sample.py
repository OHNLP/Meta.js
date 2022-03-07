import random
import string
import csv

# the number of outcomes
N_OCS = 10000

# the number of max studies
N_STUS = 50

# the output file name
fn = 'sample.csv'

# make the studies
study_dict = {}
for i in range(N_STUS):
    # make the study name
    study = 'S%02d' % i

    # make the study treatment and control
    t = random.choice(string.ascii_uppercase)
    c = random.choice(string.ascii_uppercase.replace(t, ''))

    study_dict[study] = {
        'phase': random.choice([2, 3]),
        'year': random.randint(2000, 2021),
        'type_of_therapy': random.choice([1, 2]),
        'type_of_cancer': random.randint(0, 9),
    }

# for saving the generated records
rs = []
for i in range(N_OCS):
    # make an outcome name
    oc_name = 'M%04d' % i

    # generate the number of studies to be generated
    n = random.randint(1, N_STUS)

    # generate those studies
    for j in range(n):
        # generate a study name
        study = 'S%02d' % random.randint(0, N_STUS-1)

        # generate the Et, Ec
        Et = random.randint(0, 50)
        Ec = random.randint(0, 50)

        # generate the Nt, Nc
        Nt = Et + random.randint(1, 50)
        Nc = Ec + random.randint(1, 50)

        t = 'T'
        c = 'C'
        # add this record to rs
        rs.append(
            [oc_name, study, Et, Nt, Ec, Nc, t, c]
        )

print('* generated all records %s' % (len(rs)))
print('* writing recrods to %s' % fn)

# the columns
cols = [
    'outcome', 'study', 'Et', 'Nt', 'Ec', 'Nc', 
    'treatment', 'control'
]

# write the records
with open(fn, 'w') as csvfile:
    # create a csv writer
    cw = csv.writer(csvfile)

    # write columns
    cw.writerow(cols)

    # write the records
    cw.writerows(rs)

print("* done!")