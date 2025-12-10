import sys

if len(sys.argv) != 2:
    print("Need filename")
    
file = sys.argv[1]
lines = list()
with open(file, "r") as f:
    for line in f:
        lines.append(line.strip())

with open(file, "w") as f:
    i = 0
    for line in lines:
        if i < 2 or i == len(lines) - 1:
            f.write(f"{line}\n")
        else:
            l = line.split()
            l[0] = "(" + l[0] + ","
            l[1] = l[1] + ")"
            res = " ".join(l)
            f.write(f"{res}\n")
        i += 1
    
