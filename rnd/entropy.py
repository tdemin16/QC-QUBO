import statistics as st
f = open("nums.txt", "r")

l = [0 for i in range(128)]
max = 0
for n in f:
    l[int(n)] += 1

print(str(st.mean(l)))
print(str(st.stdev(l)))

f.close()