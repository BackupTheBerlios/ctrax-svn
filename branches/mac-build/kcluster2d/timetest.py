import timeit

t = timeit.Timer('test_kcluster.test0()','import test_kcluster')
res = t.timeit(10)

print "elapsed time for test0 = " + str(res)

t = timeit.Timer('test_kcluster.test()','import test_kcluster')
res = t.timeit(10)

print "elapsed time for test = " + str(res)
