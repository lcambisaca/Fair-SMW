## FAIR-SMW

https://arxiv.org/abs/2510.20220

### code structure
```
-- src
 -- alg1.m: SC 
 -- alg2.m: FairSC
 -- alg3.m: s-FairSC
 -- alg5.m: AFF-SMW-SC
 -- alg6Sym.m: SYM-SMW-SC
 -- alg6Rw.m: RW-SMW-SC

 -- Afun.m: spmv function to eig in alg3.m
 -- Afun2.m: spmv function to eig in alg1.m and alg2.m
 -- SMWAfun.m: spmv function to eig in alg5.m, alg6SYM.m, alg6RW.m

 -- Note alg5n.m, alg4Rw, alg4Sym were the started point for our algorithm and were used as a baseline for their development
-- test
 -- m-SBM: experiment on the Modified Stochastic Block Model 
 -- lastfm: experiment on Last.fm dataset
 -- friendship: experiment on FacebookNet dataset
 -- ranLap: experiment on random Laplacian
 -- German: experiment on German dataset
 -- Deezer: experiment on Deezer dataset
 -- PersonalTest: experimented on dense and sparse SBM matrix, as well as simulated unique matrix that pose problems for traditonal Fair-SC inorder to test robustness

alg1.m, alg2.m, and part of m-SBM are credited to https://github.com/matthklein/fair_spectral_clustering

alg3.m, our test code format and cleaned data set for lastfm and freindship are credited to https://github.com/jiiwang/scalable_fair_spectral_clustering

German, Deezer cleaned data sets were gathered from https://github.com/JiaLi2000/FNM 

Real dataset: [Last.fm](http://snap.stanford.edu/data/feather-lastfm-social.html) | [FacebookNet](http://www.sociopatterns.org/datasets/high-school-contact-and-friendship-networks/) | [Deezer](https://snap.stanford.edu/data/feather-deezer-social.html) | [German] (https://github.com/yushundong/Graph-Mining-Fairness-Data/tree/main/dataset/german )

