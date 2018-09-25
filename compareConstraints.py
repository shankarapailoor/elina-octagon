import sys

def addConstraints(file):
    constraints = []
    with open(file, 'rb') as f:
        inConstraint = False
        constraint = None
        for line in f:
            if not inConstraint:
                if "PRINTING" in str(line.strip()):
                    inConstraint = True
                    constraint = []
            elif inConstraint and ("PRINTED" in str(line.strip())):
                inConstraint = False
                constraints.append(constraint)
            elif inConstraint:
                constraint.append(line.strip())
            else:
                continue
    return constraints



constraints1 = addConstraints(sys.argv[1])
constraints2 = addConstraints(sys.argv[2])
#print(constraints2)
#print(constraints1)
#print(len(constraints1), len(constraints2))
n = min(len(constraints1), len(constraints2))
for i in range(n):
    c1 = constraints1[i]
    c2 = constraints2[i]
    for j in range(len(c1)):
        #print(c1[j], c2[j])
        if c1[j] != c2[j]:
            print("NOT EQUAL CONSTRAINTS: ", c1[j], c2[j])


