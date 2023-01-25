# ALPHAMALIG(ALPHAbet Multiple ALIGnment)

※これは、自分で新たに作成したコードではない

## 概要

ALPHAMALIGは、テキスト列に対するマルチプルアライメントツールである。

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

## 出所と修正内容

Universitat Politècnica de Catalunya / BarcelonaTech (UPC)のNLPグループとゲノム系グループの2002年頃の成果物に基づく。

- [ALPHAMALIG - Alignment of sequences of a finite alphabet](https://alggen.cs.upc.edu/recerca/align/alphamalig/intro-alphamalig.html)
  - https://alggen.cs.upc.edu/recerca/align/alphamalig/alphamalig.tar
  - [Multiple Sequence Alignment for Linguistics](https://www.cs.upc.edu/~nlp/msa.html)
  - [Example application of Multiple Sequence Alignment (MSA) to linguistic phenomena](https://www.cs.upc.edu/~nlp/exampleMSA.html)
  - [Tools, demos and resources](https://www.cs.upc.edu/~nlp/tools.html)
- [ALGGEN - RECERCA](https://alggen.lsi.upc.es/recerca/frame-recerca.html)
- paper: [Multiple sequence alignments in linguistics](https://dl.acm.org/doi/pdf/10.5555/1642049.1642052) (LaTeCH - SHELT&R@EACL 2009)


古いコードのせいかプロトタイプ宣言の修正が必要で、また権利記載がない。

以下の修正を行った。

- プロトタイプ宣言
- 関数名、変数名及びコメントの翻訳 (カタロニア語->英語)
- アルファベット定義の拡張
  - 英大文字のみから、英小文字を含めた対応
  - あわせて hex表記での non-printable アルファベット定義対応
    - non-printable-character のFASTAファイルの取り扱いは、`MAFFT` の `maffttext2hex`, `hex2maffttext` と組み合わせて実施する前提
      - [Systems Immunology Lab / mafft · GitLab](https://gitlab.com/sysimm/mafft)
      - [MAFFT - a multiple sequence alignment program](https://mafft.cbrc.jp/alignment/software/)
      - [Non-biological sequences : MAFFT - a multiple sequence alignment program](https://mafft.cbrc.jp/alignment/software/textcomparison.html)
  - アルファベット＋コスト定義ファイルの作成補助ツール `create_alphabet_matrix.py` の作成


## インストール

```bash
apt install libgd-dev libjpeg-dev
make
./alfm alphabetexample.txt sequencesexample
```

## 利用方法

```bash
alfm alphabetexample.txt sequencesexample
```

### アルファベットとコスト定義ファイル

- アルファベット文字数( Gap文字を含む)
- アルファベット文字の羅列。空白区切り。末尾はGap ('-')文字。
- コスト行列。一致コスト。不一致コスト(文字間で対称)。最終行はGap挿入コスト

```text
6 
o p s c n -
2 
-1 15
-2 -2 1
-2 -2 0 1
-2 -2 -1 -1 1
-2 -2 0 0 0 0
```

non-printable-characterを含む場合は hex 表記とすること。(下記の内容は上記と同じ)

```text
6 
6f 70 73 63 6e 2d
2 
-1 15
-2 -2 1
-2 -2 0 1
-2 -2 -1 -1 1
-2 -2 0 0 0 0
```

次の文字は含まないこと。
- NUL (0x00)
- '>' (0x3e)
- '=' (0x3d)
- '<' (0x3c)
- Space (0x20)
- Carriage Return (0x0d) 
- Line Feed (0x0a)

以下の文字をGap文字として、アルファベット定義の末尾に含めること。
- '-' (0x2d)

上記以外の 248文字を有効なアルファベットとして利用できる。


典型例は以下のとおり。
- 文字ごとの一致コストは文字によらず同じ 100
- 文字ごとの不一致コストは文字の組み合わせによらず同じ -10
- Gap挿入コストは 0

```
6 
a b c d e -
100 
-10 100
-10 -10 100
-10 -10 -10 100
-10 -10 -10 -10 100
0 0 0 0 0 0
```

non-printable-characterとして alphabets = [0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x2d] を定義した場合(末尾は Gap '-' (0x2d))。

```
6 
18 19 1a 1b 1c 2d
100 
-10 100
-10 -10 100
-10 -10 -10 100
-10 -10 -10 -10 100
0 0 0 0 0 0
```


アルファベットとコスト定義ファイルを作成する補助ツールとして `create_alphabet_matrix.py` を用意している。
コストとアルファベット定義を修正して実行すると、上記の形式でアルファベットとコスト定義を出力する。

```python
# create_alphabet_matrix.py

    costs = {
        'match': 100.0,     # match はアルファベットによらず固定値。
        'mismatch': -10.0,  # mismatch はアルファベットによらず固定値。かつ対称
        'gap_penalty': 0.0, # gap_penalty は前後のアルファベットや継続長によらず固定値。
    }

    alphabets =  [ 'a', 'b', 'c', '-']
```

```bash
python3 create_alphabet_matrix.py > alphabet_matrix
```

```text
4
a b c -
100.0
-10.0 100.0
-10.0 -10.0 100.0
0.0 0.0 0.0 0.0
```

### テキストシーケンス列のファイル

入力はFASTA形式のファイル。

```text
>1
ppposnonoccsnopoosoononoscs (以降、略)
>2
osspoooosososoposopop (以降、略)
>3
nospococooospoosposnopossososososooponos (以降、略)
>4
spopoopopocpocnpocopsossonpo (以降、略)
>5
ooosposonospposoononopocososoponocsnop (以降、略)
 (以降、略)
```

アルファベットを拡張した際にnon-printable-characterを含む場合の扱いは、"[MAFFT の Non-biological sequences の扱い](https://mafft.cbrc.jp/alignment/software/textcomparison.html)" にならう。

```text
>sequence1
01 02 03 4e 6f 72 74 68 65 72 6e 5f 70 61 ... (以降、略)
>sequence2
01 02 03 4e 6f 72 74 68 65 72 6e 5f 70 61 ... (以降、略)
>sequence3
a3 6f 5f 47 6f 6d 65 73 ... (以降、略)
>sequence4
01 02 03 46 c3 b3 67 6f 3a 5f 6e 6f 72 74 68 65 72 6e 5f 70 61 ... (以降、略)
 (以降、略)
```

non-printable-characterを含むアルファベットとシーケンスに対する作業パイプラインは次のとおり。
```bash
$ /usr/local/libexec/mafft/hex2maffttext input.hex > input.ASCII
$ alfm alphabet_matrix input.ASCII | grep -A 2000 "Number of sequnces=" | tail -n +2 | sort -g | sed '/^[[:blank:]]*$/d' > output.ASCII
$ /usr/local/libexec/mafft/maffttext2hex output.ASCII > output.hex
```

### 出力

出力のうち、アライメント部分はCLUSTALフォーマットか。後処理のために当該部分を抽出するには、次のようにするとよい。

```bash
alfm alphabetexample.txt sequencesexample | grep -A 2000 "Number of sequnces=" | tail -n +2 | sort -g | sed '/^[[:blank:]]*$/d'
```

(出力のうちCLUSTALフォーマットのアライメント部分を切り出して、シーケンス番号順に並べ替え、空行を削除。クラスタの維持を優先するならば、シーケンス番号順の並び替えは不要だろう)

結果: 
```
1          ------------pp-po---s-no-no--ccsno-p-o--o-so---o-no------no---s-cspo-- (以降、略)
2          -----os---s-po--o-oos--os-o----s-o-p-o-s---o------p-------o-------po-- (以降、略)
3          ---n-os-----poc-o----c-o--o------osp-o-----o--s---p---o-sno-------poss (以降、略)
4          ------s-----p---o------------------p-o-----o------p-------o-------p--- (以降、略)
5          -----o-------o--o---s--------------p-o-son-o--s---p---------------po-- (以降、略)
 (以降、略)
```


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
