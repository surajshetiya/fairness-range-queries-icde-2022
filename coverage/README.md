# Coverage-based-Query-Rewriting

This repository is modified from https://github.com/chia-a/Coverage-based-Query-Rewriting

The changes include

1. db\_enter.py : can be used enter a CSV into a postgres database
2. The coverage paper has been modified to work with demographic disparity measure that we propose in https://www.cs.uic.edu/∼indexlab/assets/fairRQ.pdf

We consider nondiscrimination expressed in terms of coverage constraints and its impact on data transformation, proposing solutions for detecting the bias, mitigating it, and checking whether the mitigation was effective.

*Coverage constraints* guarantee that the dataset includes enough examples for each (protected) category of interest, defined in terms of sensitive attribute values.

The proposed Coverage-based Query Rewriting approach aims at guaranteeing that the result of a transformation, defined as a Select-Project-Join query, over a tabular dataset, satisfies a set of given coverage constraints.
Specifically, given an SPJ query over a tabular dataset and a set of coverage constraints, the algorithm produces a rewritten query (that is the "closest" one to the input query) satisfying those constraints.


The proposed approach is approximate and relies on a sample-based cardinality estimation, thus it introduces a trade-off between the accuracy and the efficiency of the process.
For evaluating and quantifying the error that can be generated three groups of measures have been introduced: grid-based, solution-based and sample-based accuracy measures.

This repository contains both the code for running the coverage-based query rewriting algorithm and computing the proposed accuracy measures.
The code is written in Python3. PostgreSQL is required to run this code (but you can easily replace it with another relational DBMS)


### Organization of the repository
This repository is organized as follows:

- *1_coverage_rewriting_sql.py* : code for running the three main variations of the rewriting algorithms (CRBase, CRBaseP, CRBaseI, CRBaseIP) (See [1] for further details)
- *2_calc_measures_sol.py* : code for computing grid-based and solution-based measures (See [2] for further details)
- *2_calc_measures_sample.py* : code for computing sample-based measures (See [2] for further details)




### References

[1] *Coverage-based Rewriting for Data Preparation*. C. Accinelli, S. Minisi and B. Catania. EDBT/ICDT Workshops 2020

[2] *The impact of rewriting on coverage constraint satisfaction*. C. Accinelli, B. Catania, G. Guerrini and S. Minisi. EDBT/ICDT Workshops 2021.

[3] *covRew: a Python Toolkit for Pre-Processing Pipeline Rewriting Ensuring Coverage Constraint Satisfaction*. C. Accinelli, B. Catania, G. Guerrini and S. Minisi. EDBT 2021.
