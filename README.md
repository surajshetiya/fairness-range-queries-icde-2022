# Fairness in range queries

This repostory is contains code to test the techniques proposed in paper[1].

As a prerequisite, the C++ code requires boost libraries. The python code requires numpy, pandas and other minor packages.

To compile the IBFSMP code, compile the file using gcc.

> g++ AST\_uniform.cpp -o AST\_uniform

To run the program with x >= 10 AND x <= 20 AND y >=30 AND y <= 40, use below format.

> ./AST\_uniform 10 30 20 40

To run the web interface, naviagate to the folder web.

Install numpy and pandas to run the script, server.py

> python server.py


To execute coverage based paper[2] to compare our measure with naviagate to the folder coverage and execute the script 1\_coverage\_rewriting\_sql.py

> python 1\_coverage\_rewriting\_sql.py

Change the postgres credentials in the script to work with the database.

[1] Fairness-Aware Range Queries for Selecting Unbiased Data
[2] Coverage-based Rewriting for Data Preparation

