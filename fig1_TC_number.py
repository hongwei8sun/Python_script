import numpy as np
import math

f = open ("TC1988.txt","r")
c = f.readlines()

iyear = np.linspace(1988, 2017, num=30)

for iyr in iyear:
	i = 0
	for line in c:
	
		if line[4:8]==str(int(iyr)) :
			i = i + 1
	print(int(iyr),i)


