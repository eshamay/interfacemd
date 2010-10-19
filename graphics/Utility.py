import os
import fnmatch

# an iterator that grabs all matches to a filename patter in a dir-tree. Much like `find . -name 'pattern'`
def SearchDirectoryTree(srcpath, pattern):
	for root, dirnames, filenames in os.walk(srcpath):
		for filename in fnmatch.filter(filenames, pattern):
			yield os.path.join(root, filename)

group = lambda t, n: zip(*[t[i::n] for i in range(n)])

