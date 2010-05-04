from hotshot import stats
s = stats.load("hotshot_ctrax_stats")
s.sort_stats("time").print_stats()

