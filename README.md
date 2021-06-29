# lsrm_comparison
Code for  comparing hidden structure of cbcl and ysr using LSRM


- CBCL : 

  - Contains cbcl binay matrix in 'data' folder. 
  - Use lsrm_DA.ipynb for fitting lsrm with cbcl data. 
- YSR : 

  - Contains ysr binay matrix in 'data' folder. 
  - Use lsrm_DA.ipynb for fitting lsrm with ysr data.

- notebook : 
  - Code for comparing lsrm result of YSR with lsrm result of CBCL.
  - Code contains 
    1. defining functions like calculating overlapped portion, KL-divergence, Jaccard similarity
    2. loading and  post processing of each result
    3. Point estimate 
    4. Calculate basic matrix using pre-defined function
    5. Basic visualization using basic matrix
   
- Supplementary : 

  - Contains supplementary information about data, model, and analysis