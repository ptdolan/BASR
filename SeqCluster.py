#Usage: python seqCluster.py <myfasta.fa>
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sklearn as skl
from multiprocessing import Pool
from sklearn.cluster import AgglomerativeClustering
from Bio import SeqIO
import sys

inFile=sys.argv[1]
if len(sys.argv)>2:
    K=int(sys.argv[2])
else:
    K=9

ClusterN=K

def distMatch(packet):
    (pos,dims,fragList,kmerL)=packet
    matchList=[[(np.sum(np.array(fragList[pos][taxA] != fragList[pos][taxB]))) for taxA in range(dims[1])] for taxB in range(dims[1])]
    return(np.array(matchList))

def SeqClust(fastaAln,kmerL=21):
    inF=fastaAln
    records = list(SeqIO.parse(inF, "fasta"))

    fragList=np.array([[i.seq[j:j+kmerL] for i in records ] for j in range(0,len(records[0].seq)-kmerL,kmerL)])

    dims=np.shape(fragList)

    print("Computing K-mer Distances, with K-mer length:"+str(kmerL))
    with Pool(10) as p:
        distList=p.map(distMatch, [(i,dims,fragList,kmerL) for i in range(dims[0])])
    print("Done.")
    npDL=np.array(distList)
    al=np.zeros([dims[1],dims[1]])
    n=0

    for i in npDL:
        al+=i
        n+=1
    dists=al/n
    print(pd.DataFrame(dists))
    #mpl.pyplot.imshow(dists)

    ac = AgglomerativeClustering(ClusterN, affinity='precomputed',linkage="single")
    clusters=ac.fit_predict(np.matrix(dists))

    scaled=skl.manifold.MDS(n_components=2)
    Emb=scaled.fit(dists).embedding_
    #print(Emb)

    plt.scatter([i[0] for i in Emb],[i[1] for i in Emb],c=clusters,cmap="Spectral")
    plt.show()

    return()

if __name__ == "__main__":
    SeqClust(inFile)
