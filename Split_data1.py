

f = open ("TC1988.txt","r")

c = f.readlines()
i = 1

for line in c:
	
	if line[0:2]=='EP' or line[0:2]=='CP':
		print(line[0:2])
		f.close()
		f = open(str(i)+'_'+str(line[0:7])+'.txt',"w+")
		i = i + 1
	
	f.write(line)
