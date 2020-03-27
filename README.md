# lumpy2site

These two scripts were used in the publication (Helou _et al_, 2020) [work in progress].

They allow post-processing the results of the Lumpy structural variant detection tool (RM Layer et al, 2014).

## Extraction of interchromosomal events (mapping to transposon and genome)

```
usage : extract_bnd_with_evid.py [-in IN] [-out OUT] [-name NAME]
```

**Required arguments**
```
  -in IN      VCF file output from lumpy
  -out OUT    output file with results
  -name NAME  name of the transposon (from your bwa bank) for boundary to filter ex: -name IFP2_neo
```

## Gff file grouping one line per event

```
usage: lumpy2gff.py [-in IN] [-out OUT] [-size SIZE]
```

**Required arguments**
```
  -in IN      file from extract_bnd_with_evid.py
  -out OUT    output gff file with results
  -size SIZE  size of the transposon ex: -size 2133
```

