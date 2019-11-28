from math import sqrt
from bisect import bisect_left

def getSums(n):
	if n== 0:
		return [0]
	prev = getSums(n-1)
	return prev + [sqrt(n) + i for i in prev]

if __name__ == '__main__':
	N = 20
	sums = getSums(N)
	target = sum(sqrt(n) for n in range(1, N+1))/2
	print(sums[-1] - target*2, target)
	sums.sort()
	print(sums[bisect_left(sums, target)] - target)