#!/usr/bin/env python

import random

names = ['Paul','Ivan','Eric','Nate','Bryon','Jacob','Joe','Aaron','Geoff','Chris','Bradford','Anders','Tenzing','Carlos','Arlo','Erick']

team1 = ['Ivan']
team2 = ['Jacob']

n = names
for i in xrange(len(names)):
	num = random.randint(0,len(n)-1)
	if n[num] in team1 or n[num] in team2:
		n.remove(n[num])
		continue
	if divmod(i,2)[1] == 0:
		team1.append(n[num])
	else:
		team2.append(n[num])

	n.remove(n[num])

print team1
print team2