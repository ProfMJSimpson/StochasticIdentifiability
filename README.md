# StochasticIdentifiability
Matlab and Maple code for "Profile likelihood analysis for a stochastic model of diffusion in heterogeneous media" preprint, available at https://arxiv.org/abs/2011.03638
Authors: Matthew J Simpson, Alexander P Browning, Christopher Drovandi, Elliot J Carr, Oliver J Maclaren, Ruth E Baker.

Maple code is provided to solve the system of boundary value problems (2.5)-(2.7).  This applies to an arbitrary system with m layers and n moments.  Results with m=2,3,4 and n=2 are used to generate results in Figure 3.

Several folders of Matlab code are provided:

Folder Onelayer provides Matlab code that produces results in Section S3, Figure S3.

Folder Twolayer provides Matlab code that produces results in Section 2(f), Figure 4 (separate codes for each part of Figure 4 are listed).

Folder Threelayer provides Matlab code that produces results in Section 2(g), Figure 5. 

Folder Fourlayer provides Matlab code that produces results in Section S4(a), Figure S4. 

Folder Threelayer_homogenisation provides Matlab code that produces results in Section 2(h), Figure 6.

Folder Fourlayer_homogenisation provides Matlab code that produces results in Section S4(b), Figure S5. 

Folder Moments_Benchmarking provides Matlab code that evaluates the Maple code and compares exact and stochastic results in Section 2(c), Figure 3.
