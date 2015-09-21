#GOrichr - a command line GO term enrichment function
This tool was designed to perform a GO term enrichment.
It implements a modified Fisher's exact test: Fisher's non-central hypergeometric distribution with an odd's ratio bias of w (w defined by user).
This tool calls significant an enrichment of at least w between
(the number of genes present in gene set and GO term) / (number of genes in gene set but not in GO term)


#Installation
No installation necessary. Everything you need is in the folder. Just pull and use.

#How to use
In order to use the tool, you must launch from your command line the tool, followed by a .txt file containing the geneset you want to analyse (no quotation marks, one gene per line)
as simple as this:
```python GOrichr.py geneset.txt```

The tool has been optimised to be used with PyPy You need to download and install it. The tool might be a lot slower (10 fold) without PyPy.
FOr more information look here: http://pypy.org/


#Motivation
This tool has been used to perform GO term enrichment for the following paper:

Acute monocytic leukemia cells display non-oncogene addiction to immunoproteasomes.
by Alexandre Rouette, Assya Trofimov, David Haberl, Geneviève Boucher, Josée Hébert, Guy Sauvageau, Sébastien Lemieux and Claude Perreault

Manuscript under review.
