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

[ALGEN - framealign;](https://alggen.lsi.upc.es/recerca/align/frame-align.html)のページから:
```
AlphaMALIG 1.1
Alpha Multiple ALIGnment tool is a collaborative research project with Laura Alonso. The software has been designed by J. Escribano of the Universitat Politècnica de Catalunya.
Given a set of sequences of any alphabet and the parameters of the alignment (the score of the match, mismatch, insertion and deletion) this tools builds the multialignment of sequences.
```

## 出典

Universitat Politècnica de Catalunya / BarcelonaTech (UPC)のNLPグループの2002年頃の成果物から。

- [Multiple Sequence Alignment for Linguistics](https://www.cs.upc.edu/~nlp/msa.html)
- [ALPHAMALIG - Alignment of sequences of a finite alphabet](https://alggen.cs.upc.edu/recerca/align/alphamalig/intro-alphamalig.html)
  - https://alggen.cs.upc.edu/recerca/align/alphamalig/alphamalig.tar
- [Example application of Multiple Sequence Alignment (MSA) to linguistic phenomena](https://www.cs.upc.edu/~nlp/exampleMSA.html)
- [Tools, demos and resources](https://www.cs.upc.edu/~nlp/tools.html)
- [ALGGEN - RECERCA](https://alggen.lsi.upc.es/recerca/frame-recerca.html)
- [Multiple sequence alignments in linguistics](https://dl.acm.org/doi/pdf/10.5555/1642049.1642052)

古いコードなせいか、プロトタイプ宣言など要修正。権利記載がない。

## インストール方法

```bash
apt install libgd-dev libjpeg-dev
make
./alfm alphabetexample.txt sequencesexample
```

## 修正

- プロトタイプ宣言
- 関数名、変数名及びコメントの翻訳 (カタロニア語->英語)
- 英大文字のみの対応から、英小文字を含めた対応。hex表記でのアルファベット定義対応
  - non-printable-character のFASTAファイルについてはMAFFTのmaffttext2hex, hex2maffttext との組み合わせが前提
    - [Systems Immunology Lab / mafft · GitLab](https://gitlab.com/sysimm/mafft)
    - [MAFFT - a multiple sequence alignment program](https://mafft.cbrc.jp/alignment/software/)
    - [Non-biological sequences : MAFFT - a multiple sequence alignment program](https://mafft.cbrc.jp/alignment/software/textcomparison.html)
  - アルファベット＋コスト定義ファイルの作成補助ツール(python)
  
## 参考

- [Publications Section.](https://www.cs.upc.edu/~nlp/papers.html) 
- [SemEval-2007 - UPC Universitat Politècnica de Catalunya](https://www.cs.upc.edu/~nlp/semeval/msacs_home.html)
- [SemEval-2007 - UPC Universitat Politècnica de Catalunya](https://www.cs.upc.edu/~nlp/semeval/msacs_download.html)

- [Multiple sequence alignment as a sequence-to-sequence learning problem | OpenReview](https://openreview.net/forum?id=8efJYMBrNb)

- [Using deep reinforcement learning approach for solving the multiple sequence alignment problem | SpringerLink](https://link.springer.com/article/10.1007/s42452-019-0611-4)
- [End-to-end learning of multiple sequence alignments with differentiable Smith–Waterman | Bioinformatics | Oxford Academic](https://academic.oup.com/bioinformatics/article/39/1/btac724/6820925?login=false)
- [DeepMSA: constructing deep multiple sequence alignment to improve contact prediction and fold-recognition for distant-homology proteins - PubMed](https://pubmed.ncbi.nlm.nih.gov/31738385/)

- [Multiple Sequence Alignment | Papers With Code](https://paperswithcode.com/task/multiple-sequence-alignment)
- [MSA Transformer | Papers With Code](https://paperswithcode.com/paper/msa-transformer)
- [MSA Transformer | bioRxiv](https://www.biorxiv.org/content/10.1101/2021.02.12.430858v3)

- [Protein language models trained on multiple sequence alignments learn phylogenetic relationships | Nature Communications](https://www.nature.com/articles/s41467-022-34032-y)
- [[2110.07609] Application of Sequence Embedding in Protein Sequence-Based Predictions](https://arxiv.org/abs/2110.07609)
- [nlp - Building NER using Sequence Alignment algorithms - Stack Overflow](https://stackoverflow.com/questions/34365621/building-ner-using-sequence-alignment-algorithms)
- [sarahamick/NLP_sequence_alignment: An algorithm for calculating the minimum cost of aligning a sequence n with a sequence m given gap and mismatch penalties](https://github.com/sarahamick/NLP_sequence_alignment)
- [The module for multiple sequence alignments, AlignIO · Biopython](https://biopython.org/wiki/AlignIO)
- [MSA-Multiple-sequence-alignment-の作成 / スッキリわかるAlphaFold2 - どこから見てもメンダコ](https://horomary.hatenablog.com/entry/2021/10/01/194825#MSA-Multiple-sequence-alignment-%E3%81%AE%E4%BD%9C%E6%88%90)
- [HMMER: biosequence analysis using profile hidden Markov models](http://hmmer.org/)
- [Generative power of a protein language model trained on multiple sequence alignments | bioRxiv](https://www.biorxiv.org/content/10.1101/2022.04.14.488405v2)
- [encounter1997/SFA: Official Implementation of "Exploring Sequence Feature Alignment for Domain Adaptive Detection Transformers"](https://github.com/encounter1997/SFA)
- [OpenFold2](https://lupoglaz.github.io/OpenFold2/msa.html)

- [シーケンスアラインメント - Wikipedia](https://ja.wikipedia.org/wiki/%E3%82%B7%E3%83%BC%E3%82%B1%E3%83%B3%E3%82%B9%E3%82%A2%E3%83%A9%E3%82%A4%E3%83%B3%E3%83%A1%E3%83%B3%E3%83%88)
