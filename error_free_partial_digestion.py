import random
import os
import linecache
import sys
import errno
def find_all(a_str, sub):
    start = 0
    while True:
        start = a_str.find(sub, start)
        if start == -1: return
        yield start
        start += len(sub) 

if __name__ == "__main__":
    if len(sys.argv)<3:
        print "Please add the pattern and the fileName to the command line.\nUsage: python error_free_partial_digestion.py patern fName"
    else:
        pattern = sys.argv[1]
        fName = sys.argv[2]
        file = open(fName,'r')
        line = file.read();
        s=line.strip()
        file.close()
        """
        array_line=[]
        If you download genome from ncbi, htis steps will make the genome in one line
        for line in file:
            s=line.replace("\n",'')
            array_line.append(s)
        """
        positions=list(find_all(s, pattern))
        positions.append(0)
        positions.append(len(s))
        Distances = []
        for i in range(0, len(positions)):
            for j in range(i+1, len(positions)):
                Distances.append(abs(positions[i]-positions[j]))
        Distances.sort()
        f_out = open("distancses_"+pattern+".txt",'w')
        for x in Distances:
            f_out.write(str(x)+"\t")
        f_out.close()
        f_out = open("positions_"+pattern+".txt",'w')
        for x in positions:
            f_out.write(str(x)+"\t")
        f_out.close()
        print ("No. of distancses (N) is "+str(len(Distances)))
        #print Distances
        #print positions
        print ("No. of restriction sites (n) is "+str(len(positions)))
   
