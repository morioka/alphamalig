# ALPHAMALIG(ALPHAbet Multiple ALIGnment)

※これは、自分で新たに作成したコードではない

## 何であるか

[ALPHAMALIG - Alignment of sequences of a finite alphabet](https://alggen.cs.upc.edu/recerca/align/alphamalig/intro-alphamalig.html)のページから:
```
ALPHAMALIG(ALPHAbet Multiple ALIGnment) is a tool which allows you to align many but no more than 200 sequences. These sequences will be formed by characters of a finite alphabet.

All you have to do is upload a file with the sequences in FASTA format and another file with the information of
the alphabet. A limitation of ALPHAMALIG is that the original sequences must not be longer than 2000 characters.
In fact, sequences shorter than 1800 chracters are strongly recommended.
The web service will be no longer available but a tar file with the source code and the running files for linux can be downloded here .

Please report any bug or comment to Xavier Messeguer.
Last updated: 17 June 2003
This software was made by Jordi Escribano.
```

## 出典

Universitat Politècnica de Catalunya / BarcelonaTech (UPC)のNLPグループの2002年頃の成果物から。

- [Multiple Sequence Alignment for Linguistics](https://www.cs.upc.edu/~nlp/msa.html)
- [ALPHAMALIG - Alignment of sequences of a finite alphabet](https://alggen.cs.upc.edu/recerca/align/alphamalig/intro-alphamalig.html)
  - https://alggen.cs.upc.edu/recerca/align/alphamalig/alphamalig.tar
- [Example application of Multiple Sequence Alignment (MSA) to linguistic phenomena](https://www.cs.upc.edu/~nlp/exampleMSA.html)
- [Tools, demos and resources](https://www.cs.upc.edu/~nlp/tools.html)

古いコードなせいか、プロトタイプ宣言など修正。権利関係は不明。

## インストール方法

```bash
apt install libgd-dev libjpeg-dev
make
./alfm alphabetexample.txt sequencesexample
```
