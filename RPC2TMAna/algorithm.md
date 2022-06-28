# Algorithm.

I think the RPCHitCleaner strores only RPC Digis which satisfy all the following conditions:

1. Occur in the Barrel region.

2. Occur at bx window such that: |bx| <= 3.

3. Belong to clusters such that: cluster_size <= 3.
- So we constructed the first loop over chambers and filled the (vcluster_size) vector which assigns a cluster_id for each cluster and stores its cluster_size in the corresponding index in the vector.

4. Belong to the cluster such that: the cluster has the min_bx in its corresponding chamber.
- So we constructed the second loop over chambers and first inner loop over digis to fill the (std::map<RPCDetId, int> bx_hits) map whcih stores the min_bx for a cluster (value) in each chamber (key) to use it for checking later.

5. Store all digis which belong to clusters of cluster_size < 3. And for digis which belong to clusters of clu_size = 3, store only the second digi of the cluster in (m_outrpcDigis) RPCDigiCollection.
- So we constructed the second loop over digis in the second loop over chambers to fill (m_outrpcDigis) RPCDigiCollection.
