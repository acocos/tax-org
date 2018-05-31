## Predicted relations

This sqlite3 database file contains predicted hypernym, synonym[equiv], and entailment[entmax=max(hypernym,equiv)] edges.

If you would like to experiment with your own predicted relations, you can either use this database schema, or provide a relation .csv file with the following headers:

w1    w2    PPDB20Score    hypernym    equiv    entmax

The PPDB20Score is used for filtering in some of the code; To ignore this, simply set all rows' PPDB20Score to some high number (e.g. 100).

Values for the hypernym, equiv, and entmax fields should be floats.

w1 and w2 give the ordered pair of entities for which a relation is predicted.
