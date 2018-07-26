# tourmaline

<img src="https://en.wikipedia.org/wiki/Tourmaline#/media/File:Tourmaline-121240.jpg" width=100px>

This workflow uses [QIIME 2](https://qiime2.org) for most tasks. Two methods of amplicon data processing are supported, both of which generate ASVs (amplicon sequence variants, approximating the "true"/"exact" sequences in the sample) rather than OTUs (operational taxonomic units, which blur sequencing errors and microdiversity through clustering):

* [Deblur](https://github.com/biocore/deblur)
* [DADA2](https://github.com/benjjneb/dada2)

