import learn_bg_prior
import numpy as num

movienames = ['/media/data/data3/pGAWBCyO-GAL4_TrpA_Rig1Plate01Bowl1_20100820T140408/movie.fmf',]
annnames = ['/media/data/data3/pGAWBCyO-GAL4_TrpA_Rig1Plate01Bowl1_20100820T140408/movie.fmf.ann',]

BGModelPrior = learn_bg_prior.BGModelPrior(movienames,annnames)

BGModelPrior.est_prior()
