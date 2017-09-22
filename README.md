# sparse-permanents

The functions here allow to compute permanents of dense or sparse matrices.
The function sparsePerm uses the tree decomposition algorithm from

- [Cifuentes, Parrilo (2016)](https://arxiv.org/abs/1507.03046), "An efficient tree decomposition method for permanents and mixed discriminants", *Linear Algebra and its Applications*, 493:45--81

The complexity of the algorithm is *O(n 2^k)*, where *k* is the width of the decomposition.
For a band matrix the complexity is *O(n 2^{w1+w2})*, where *w1*, *w2* are the lower and upper bandwidths of the matrix.
