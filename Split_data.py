

f = open ("TC1988.txt","r")

c = f.readlines()
i = 0
tc = []

for line in c:
	
	if line[0:2]=='EP' or line[0:2]=='CP':
		tc.append([])
		i = i + 1
	
	tc[i-1].append(line)
print(tc[0][1])
