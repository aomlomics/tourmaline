# tourmaline

<center><img src="https://upload.wikimedia.org/wikipedia/commons/0/00/Tourmaline-121240.jpg" width=200></center>

This workflow uses [QIIME 2](https://qiime2.org) for most tasks. Two methods of amplicon data processing are supported, both of which generate ASVs (amplicon sequence variants, approximating the "true"/"exact" sequences in the sample) rather than OTUs (operational taxonomic units, which blur sequencing errors and microdiversity through clustering):

* [Deblur](https://github.com/biocore/deblur)
* [DADA2](https://github.com/benjjneb/dada2)

