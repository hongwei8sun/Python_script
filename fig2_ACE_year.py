import numpy as np
import math

f = open ("TC1988.txt","r")
c = f.readlines()

iyear = np.linspace(1988, 2017, num=30)

ACE_year = iyear*0.0

for line in c:
	for iyr in iyear:
		i = int(iyr - 1988)
		if line[0:4]==str(int(iyr)):
			#print(line[0:8],line[39:41])
			ACE_year[i] = ACE_year[i] + 10**(-4)*float(line[39:41])**2

for iyr in iyear:
	i = int(iyr - 1988)
	print(int(iyr),ACE_year[i])
