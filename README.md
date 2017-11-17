# Pearson
Calculates all pairwise correlations (Pearson's and Spearman's) between pairs of rows in an input matrix

usage: pearson3 filename -flag1 -flag2 . . .

filename - Tab-separated file of expression values. The first row should contain column labels and the first column should contain labels for each row. The upper-left cell is treated as a dummy label and can be blank. Any non-numeric entry, including an empty string, is interpreted as a missing value. Blank lines or lines with fewer than 2 values are ignored.

--- Flags ---
-h - output a frequency histogram
-w - output a weighted graph
-e - output an unweighted graph
-m - output a correlation matrix
-s - use Spearman's correlation (default is Pearson's correlation"
-n - output integer vertex labels instead of the original labels
-pos - only output edges with positive correlations
-neg - only output edges with negative correlations
-pvalue - use correlation p-values instead of correlations (suitable for data with missing values)

Exactly one of -h, -w, -e, or -m must be specified. If -w or -e is specified, then an additional threshold parameter between 0 and 1 must also be specified. Any edges at or above the threshold will be retained in the output graph. If -pvalue is specified, then the two-tailed correlation p-value is used instead of the correlation; values at or below the threshold will be retained. In output p-values, 0 means <0.0001.;

Rows with ultra-low-variance (< .000001) are assumed to have 0 correlation with all other rows.

Example 1
Output an unweighted graph from input file test.txt at a threshold of |0.8|.
> pearson3 test.txt -e .8

Example 2
Output a weighted graph from input file test1.txt at a threshold of |0.7|.
> pearson3 test1.txt -w .7

Example 3
Output a correlation histogram for input file test.txt
> pearson3 test.txt -h

Example 4
Output a correlation matrix for input file test.txt
> pearson3 test.txt -m

Example 5
Output an unweighted graph for input file test4.txt using a p-value of 0.01 as a threshold
> pearson3 test4.txt -e .8 -pvalue

Example 6
Output a weighted graph for input file test4.txt, retaining only positive correlations above 0.6
> pearson3 test4.txt -w .6 -pos

Example 7
Output an unweighted graph for input file test3.txt, using Spearman's correlation and a threshold of |0.6|
> pearson3 test3.txt -e 0.6 -s
