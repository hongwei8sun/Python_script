import numpy as np
import math

f = open ('TC1988.txt','r')
c = f.readlines()
print(c[1][4:6])
imonth = ['01','02','03','04','05','06','07','08','09','10','11','12']

TC_number_month = np.linspace(0.0, 0.0, num=12)

num = 0
for line in c:
	num = num + 1
	if line[0:2]=='EP' or line[0:2]=='CP':
		i = 0
		for imon in imonth:
			print(imon)
			if c[num][4:6]==imon:
				TC_number_month[i] = TC_number_month[i] + 1.0
			i = i + 1

print(TC_number_month)
