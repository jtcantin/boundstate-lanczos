from sys import argv
import string

filename = argv[1];

file = open(filename, 'r')

xDisp = float(argv[2])
yDisp = float(argv[3])
zDisp = float(argv[4])

outfilename = argv[5]
outfile = open(outfilename, 'w')

title = file.readline() #ignore Title
outfile.write(title)

header = file.readline() #ignore column header
outfile.write(header)

while 1:
        line = file.readline()
        if not line: break
        x, phi, y, theta, z, chi = map(float, string.split(line)[1:]) #ignore first value (water #)
        new_id = string.split(line)[0]

        x -= xDisp;
        y -= yDisp;
        z -= zDisp;

        outfile.write(new_id + " %13.5e %13.5e %13.5e %13.5e %13.5e %13.5e" % (x, phi, y, theta, z, chi) +"\n")

file.close()
outfile.close()
    
