#Usage: python3 seqCluster.py <myfasta.fa> <K> <Clustering Method: "AC" or "BIRCH">
#levenshtein
#MODULES
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from multiprocessing import Pool
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import Birch
from Bio import SeqIO
import sys
import random
import Levenshtein
#FUNCTIONS

def distMatch(packet):#Finds the match counts for for each pair of sequences.
    (pos,dims,fragList,kmerL)=packet  #Break apart packet of input data, necessary for Pool().
    matchList=[[(np.sum(np.array(fragList[pos][taxA] != fragList[pos][taxB]))) for taxA in range(0,dims[1])] for taxB in range(0,dims[1])]
    return(np.array(matchList))

def SeqClust(fastaAln,K,kmerL):# function that takes in a Fasta Alignment, breaks into Kmers and then performs parallelized comparison of the kmers with Pool()
    inF=fastaAln
    records = list(SeqIO.parse(inF, "fasta"))
    subsampRecords=[records[i] for i in random.sample(range(len(records)),int((len(records))/10))]
    print("Before sub-sampling:")
    print(len(records))
    print("After sub-sampling:")
    print(len(subsampRecords))
    #make fragment list
    records=subsampRecords
    print("computing")
    print(records[0])
    dists=[[Levenshtein.distance(str(i),str(j)) for i in records] for j in records]

    print(pd.DataFrame(dists))
    #mpl.pyplot.imshow(dists)
    if clusterMethod=="AC":
        ac = AgglomerativeClustering(ClusterN, affinity='precomputed',linkage="single")
        clusters=ac.fit_predict(np.matrix(dists))
        print(model)
        scaled=MDS(n_components=2)
        Emb=scaled.fit(dists).embedding_
        #print(Emb)

        plt.scatter([i[0] for i in Emb],[i[1] for i in Emb],c=clusters,cmap="Spectral")
        plt.savefig(inFile+clusterMethod+"_SeqCluster.png")
    elif clusterMethod=="BIRCH":
        model = Birch(threshold=0.01, n_clusters=K)
        clusters=model.fit_predict(np.matrix(dists))
        print(model)
        scaled=MDS(n_components=2)
        Emb=scaled.fit(dists).embedding_
        #print(Emb)

        plt.scatter([i[0] for i in Emb],[i[1] for i in Emb],c=clusters,cmap="Spectral")
        plt.savefig(inFile+clusterMethod+"_SeqCluster.png")
    selectRecords=[[records[i] for i in range(len(clusters)) if clusters[i]==cluster] for cluster in np.unique(clusters)]
    allDF=pd.DataFrame()
    for cluster in np.unique(clusters):
        print(cluster)
        clustDict=zip([cluster]*len(selectRecords[cluster]),[selectRecords[cluster][i].description for i in range(len(selectRecords[cluster]))])
        print(clustDict)
        DF=pd.DataFrame(clustDict,columns=["cluster","info"])
        print(DF)
        allDF=pd.concat([allDF,DF])
        firsts=[i[0] for i in selectRecords]
    return(allDF,firsts)

#MAIN
if __name__ == "__main__":
    kmerL=21
    inFile=sys.argv[1]
    if len(sys.argv)>=2: #Read desired cluster number
        K=int(sys.argv[2])
    else:
        K=9
    if len(sys.argv)>=3: #Read desired Cluster method
        clusterMethod=sys.argv[3]
        print("Cluster Method: "+clusterMethod)
    if clusterMethod not in ["BIRCH","AC"]:
            print("Clustering method must be Birch: 'BIRCH' or Agglomerative Clustering: 'AC'. Defaulting to Birch. ")
    else:
        clusterMethod="BIRCH"
    outputclusters, reps = SeqClust(inFile,K,kmerL=kmerL)
    outputclusters.to_csv(inFile+clusterMethod+"_outputClusters.csv")
    with open(inFile+clusterMethod+"_outputSeqs.fasta", "w") as output_handle:#write out representative samples.
        SeqIO.write(reps, output_handle, "fasta")
