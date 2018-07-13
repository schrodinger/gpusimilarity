# gpusimilarity postgres foreign data wrapper

This code is a [multicorn-based](https://multicorn.org/) foreign data wrapper for gpusimilarity.  Its intention is to be used with a postgres cartridge like the [RDKit catridge](http://www.rdkit.org/docs/Cartridge.html) to allow for super-fast similarity searching from a database context.

## Dependencies
* GPUSimilarity web endpoint available
* Postgres:  A new enough version to support multicorn
* [multicorn](https://multicorn.org/)

## Installation
python2 setup.py install

Add a foreign data wrapper 'server' once on your database:

```
CREATE SERVER gpusim_srv foreign data wrapper multicorn options (
    wrapper 'gpusim_fdw.GPUSimilarityFDW'
);
```

## Example usage
Perform a search by creating a foreign table:

```
CREATE FOREIGN TABLE similarity_search (
    SMILES character varying,
    id character varying,
    similarity real
) server gpusim_srv options (
    server 'http://localhost',
    port '8080',
    db_name 'all',
    SMILES 'CCCC',
    max_results '35'
);
```

Now you have a table that contains the results of that search:

`select * from similarity_search;`

```
             smiles              |      id       | similarity
---------------------------------+---------------+------------
 CCCC                            | CHEMBL134702  |          1
 CCCCC                           | ZINC01698513  |        0.5
 CCCCCC                          | ZINC01532209  |        0.5
 CCCCC                           | CHEMBL16102   |        0.5
 CCCCCC                          | CHEMBL15939   |        0.5
 CCC                             | CHEMBL135416  |        0.5
......
```
