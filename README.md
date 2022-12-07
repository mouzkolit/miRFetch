<h1> miRFetch Package to provide easy access to the DIANA microT and microCDS Webserver </h1>
<p> This package allows to submit RNA sequences provided from tRNA fragments ect and determine the Target Spaces using Selenium WebScraping and
Automatisation using an easy accessible API</p>


<h2> miRT Fetching Segment: </h2>

<p> To start analysis first a dictionary of sequences must generated, consisting of a key (the name of the RNA sequence) and a list harboring the individual RNA
sequences like shown below. In Future we will provide also support of pd.DataFrame as well as output from MintMap a Pipeline to annotate 
tRNA fragments</p>


```
from rnaFetch.mirTFetch import mirTFetch \n
from rnaFetch.FetchBioMart import GeneBiomart\n
from rnaFetch.mirCDSFetch import microTCDS\n

RNA = {"GlyCCC": ["GCATTGGTGGTTCAGTGGTAGAATTCTC", 
                  "GCATTGGTGGTTCAGTGGTAGAATTCTCGCC", 
                  "GCATTGGTGGTTCAGTGGTAGAATTCT"],
       "LysTTT": ["GGGAGCGCCCGGATAGCTCAGTCGGTAGAGCATCAGACTTTT",
                  "TCGGGCGGGAGTGGTGGCTTTT",
                  "TCGGGCGGGAGTGGTGGCTTT"],
       "ThrAGT": ["TCGAATCCCAGCGGTGCCTCCA",
                  "ATCCCAGCGGTGCCTCCA",
                  "ATCCCAGCGGTGCCTCCG"]
      }
```

<p> Then you can initialize the using "Chrome", "Firefox" or "Edge" </p>

```
# Change to Firefox or Edge if you prefer
# Selenium Driver is initialized in headless mode but you can ask for Browser Window setting headless = None
fetcht = mirTFetch("Chrome")
```

<p> We can then set the threshold to consider a target; Can also be added manually via pandas when table is generated
And then we can run the Pipeline to let the miRWebserver determine the Target Spaces </p>

```
fetcht.threshold = 0.95
# this will return a table or save a table in self.prediction_data
# In addition UTR sequence Table will be provided in self.utr_table
# data table will be also returned
data = fetcht.run_miRNA_analysis(dictionary)
```

<h2> BioMart Segment: </h2>
<p> Since the Table generate only contains Transcript IDs we also ofer direct conversion using the BioMart ID via multithreading
</p>
<p> This will generate a Table with the ensembl_gene_id and the external_gene_name which can be used for further analysis and furthermore we merge the 
prediction table with the biomart table for direct interaction and downstream analysis </p>

```
# Provide the column of the Table harboring the Transcript ID
# Currently only Transcript ID is supported
biomart_data = GeneBiomart(fetcht.prediction_data["Transcript_ID"].unique().tolist())
return_table = biomart_data.query_data()
final_table = biomart_data.merge_annotation(fetcht.prediction_data, return_table)

```

<h2> Get RNA miRNA overlap </h2>

<p> We also provided overlapping target spaces between miRNAs and queried sequences using the mirCDSFetch module </p>
<p> Input the final_table generate after biomart annotation into the following code snippet
The ouptut is a table of miRNA:sequence prediction partner shared in the grouped table
</p>

```
# This will connect to the microTCDS webpage via a Selenium Driver
# 500 Genes per run will be supplied in chunks
# Threshold can be set to a float between 0-1 and will be automatically set 
fetchcds = microTCDS(final_table, threshold = 0.95)
new_table = fetchcds.run_miRNA_analysis()
overlap, grouped = fetchcds.get_mt_cds_overlap(final_table, new_table)
```





