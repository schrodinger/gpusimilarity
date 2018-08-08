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
First create a foreign table that relates to a set of search options and a server:
```
CREATE FOREIGN TABLE similarity_search (
    SMILES character varying,
    id character varying,
    similarity real,
    query character varying
) server gpusim_srv options (
    server 'http://localhost',
    port '8080',
    db_name 'all',
    max_results '35'
);
```

Now search against that table, varying the search smiles by providing 'query':

`select * from similarity_search where query='CCOCCC';`

```
gpusimilarity=# select * from similarity_search where query='CCOCCC';
             smiles             |      id       | similarity | query
--------------------------------+---------------+------------+--------
 CCCOCC                         | ZINC02031623  |          1 | CCOCCC
 CCCOCCC                        | ZINC02041060  |   0.692308 | CCOCCC
 CCOCCOCC                       | ZINC02031604  |   0.692308 | CCOCCC
 CCCOCCC                        | CHEMBL3187166 |   0.692308 | CCOCCC
 CCOCCOCC                       | CHEMBL1877517 |   0.692308 | CCOCCC
 CCOCCOCCOCC                    | ZINC02041052  |   0.642857 | CCOCCC
 CCOCCOCCOCC                    | CHEMBL1235106 |   0.642857 | CCOCCC
 CCCCOCC                        | ZINC02031618  |      0.625 | CCOCCC
 CCCCOCC                        | CHEMBL3561108 |      0.625 | CCOCCC
 CCCOCOCCC                      | ZINC01680429  |        0.6 | CCOCCC
 CCOCCCCOCC                     | ZINC01656310  |        0.6 | CCOCCC
 CCCCCOCC                       | ZINC01555811  |   0.588235 | CCOCCC
 CCOCCCCCCCCCOCC                | ZINC02173592  |     0.5625 | CCOCCC
 CCCOC[C@H](COCC)O              | ZINC86522641  |       0.55 | CCOCCC
 CCCOC[C@@H](COCC)O             | ZINC86522640  |       0.55 | CCOCCC
 CCOCC                          | ZINC01657408  |   0.538462 | CCOCCC
 CCOCC                          | CHEMBL16264   |   0.538462 | CCOCCC
 CCOCCOCCOCOCCOCCOCC            | ZINC04974293  |   0.529412 | CCOCCC
 CCOCCOC                        | ZINC02563430  |   0.529412 | CCOCCC
 CCOCC[NH3+]                    | ZINC01847864  |   0.529412 | CCOCCC
 CCOCCO                         | ZINC01648262  |   0.529412 | CCOCCC
 CCOCCO                         | CHEMBL119596  |   0.529412 | CCOCCC
 CCCOC[C@H](COCC)[NH3+]         | ZINC86524603  |    0.52381 | CCOCCC
 CCCOC[C@@H](COCC)[NH3+]        | ZINC86524602  |    0.52381 | CCOCCC
 CCCOCC(=O)COCC                 | ZINC86521471  |    0.52381 | CCOCCC
 CCCOCCO                        | ZINC04283944  |        0.5 | CCOCCC
 CCCOCC[NH3+]                   | ZINC02563311  |        0.5 | CCOCCC
 CCOCC[NH2+]CCOCC               | ZINC02389863  |        0.5 | CCOCCC
 CCOCCOCCOC                     | ZINC02040081  |        0.5 | CCOCCC
 CCCOCCO                        | CHEMBL3189002 |        0.5 | CCOCCC
 CCOCCP(CCOCC)CCP(CCOCC)CCOCC   | CHEMBL1615784 |        0.5 | CCOCCC
 CCCOC[C@H](COCC)[NH2+]C        | ZINC86524608  |   0.478261 | CCOCCC
 CCCOC[C@@H](COCC)[NH2+]C       | ZINC86524607  |   0.478261 | CCOCCC
 CCOCCOCCOCCOCCOCCOCCOCCOCCOCCO | ZINC61389413  |   0.473684 | CCOCCC
 CCCOCCOCCOCCO                  | ZINC19075767  |   0.473684 | CCOCCC
```
