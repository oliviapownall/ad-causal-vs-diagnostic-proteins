# GO enrichment analysis

Gene Ontology (GO) enrichment of the 24 MR-significant proteins was run using
the PANTHER Overrepresentation Test via the Gene Ontology website
(geneontology.org). This step was performed through the web interface, so no
analysis code is included here. The settings below are sufficient to reproduce it.

## Settings

| Setting | Value |
|---|---|
| Tool | PANTHER Overrepresentation Test |
| PANTHER version | [e.g. 19.0, released YYYY-MM-DD] |
| GO release | [GO Ontology database DOI and release date] |
| Organism | *Homo sapiens* |
| Annotation datasets | GO biological process complete; GO molecular function complete; GO cellular component complete |
| Test | Fisher's exact test |
| Correction | False discovery rate (FDR) |
| Reference list | Custom: all proteins tested in the MR analysis |
| Identifier type | [gene symbols / UniProt IDs] |

## Inputs

- Analysed list: the 24 MR-significant proteins
- Reference list: all 1875 proteins included in the MR analysis

## Steps to reproduce

1. Go to geneontology.org and paste the analysed list into the
   "Gene Ontology enrichment analysis" box.
2. Select the annotation dataset (biological process, molecular function or
   cellular component) and *Homo sapiens*, then launch the analysis.
3. On the PANTHER page, change the reference list to the custom list above
   and select the matching identifier type.
4. Set the test to Fisher's exact and the correction to FDR, then relaunch.
5. Repeat for each of the three GO aspects.