#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""dataset_maker  is a module for creating datasets from a file or a string

This module comprises several different classes, FastaBuilder, FastaReader,
and DataSetMaker.
FastaReader is a module for reading objects from a file or a string,
FastaBuilder is a module for creating data objects.
DataSetMaker is a high-level module that creates dataset from its source.
"""
import os.path
import re

from tappm import dataset
from tappm import fasta


class FastaReader(object):
    """FastaUtil A utility class for fasta sequences.

    hogehoge
    """

    def __init__(self, builder, protein=True):
        """Constructor.
        """
        self.builder = builder

    def parse_file(self, filename, protein=True):
        """複数のfasta配列が入ってるファイルをパースして、
        Fastaオブジェクトのリストにして返す。とりあえず
        アミノ酸配列だけ対応したよ。"""
        if os.path.exists(filename):
            f = open(filename, 'r')
        else:
            raise ValueError(filename + " not found.")

        seq = ''
        fasta_list = []
        # seq_num = 0
        for line in f:
            if line == '' or line == "\n" or line[0] == '#':
                continue
            elif len(seq) > 0 and line[0] == '>':
                # seq_num += 1
                fasta_list.append(self.builder.create(seq))
                seq = line
            elif len(seq) == 0 and line[0] == '>':
                seq = line
            else:
                seq += line.strip()
        # Add the last line
        fasta_list.append(self.builder.create(seq))
        return fasta_list

    def parse_string(self, s, protein=True):
        """
        TODO: codes below are duplicated as parse_file(), so needs refactoring
        """
        seq = ''
        fasta_list = []
        line_no = 0
        for line in s.splitlines():
            line_no += 1
            # print(line_no, line)
            if line == '' or line == "\n" or line[0] == '#':
                continue
            elif len(seq) > 0 and line[0] == '>':
                # seq_num += 1
                # print(len(fasta_list), seq)
                fasta_list.append(self.builder.create(seq))
                seq = line + "\n"
            elif line[0] == '>':
                # the first line
                seq += line + "\n"
            else:
                seq += line.strip()
        # Add the last line
        fasta_list.append(self.builder.create(seq))
        return fasta_list


class FastaBuilder(object):
    """FastaBuilder  is a module for creating data objects

    値の登録、登録された値からの適切なプロトタイプの選択などの
    基本的な機能を実装しておく。"""

    def __init__(self, default):
        """コンストラクタ。
        """
        self.regex = {}  # 正規表現を持つ。
        self.prototypes = {}  #
        self.default_prototype = default

    def register(self, prototype_name, prototype, regex):
        """prototype_nameに対応するクラスと正規表現を
        登録する。すでに登録してある場合は、エラー。"""
        if prototype_name in self.prototypes or prototype_name in self.regex:
            raise ValueError(prototype_name + " already registered.")
        if isinstance(regex, str):
            regex = re.compile(regex)
        if not isinstance(prototype, type):
            raise TypeError(prototype + " is not a Fasta object.")
        self.prototypes[prototype_name] = prototype
        self.regex[prototype_name] = regex

    def guess_database(self, unknown_str):
        """unknown_strがどのプロトタイプに適合するのかを
        判定し、適切なプロトタイプを返す。

        @unknown_str Fasta配列そのものを想定。"""
        # 最初の行のみを取得
        unknown_str = unknown_str.split("\n", 1)[0]
        if len(unknown_str) == 0:
            raise ValueError("Empty string.")
        elif unknown_str[0] != '>':
            raise ValueError("")
        for prototype_name, regex in list(self.regex.items()):
            result = regex.match(unknown_str)
            if result:  # 正規表現がマッチした場合
                return self.prototypes[prototype_name]
        # ここまで処理が流れてくると、どのプロトタイプにもマッチしなかった
        # ということ。デフォルトのプロトタイプを返す。
        return self.default_prototype

    def create(self, unknown_str):
        """インスタンスを生成するメソッド。
        """
        prototype = self.guess_database(unknown_str)
        return prototype(unknown_str)

    def new_class(self, classname, baseclass, attr, regex):
        """動的に新しいクラスを生成して、登録します。

        @classname 新たなクラスの名前。既存のものと衝突しないものである
        必要があります。
        ＠baseclass 基底クラス。通常はFastaまたはProteinFastaあたりを指定
        すると良いと思います。
        @attr そのクラスが持つ属性値。re_identifierなどの正規表現を登録する
        必要があります。
        @regex 新たなクラスであるかどうかを判別するための正規表現です。
        """
        # 名前の衝突をチェック
        if classname in self.prototypes:
            raise ValueError(classname + "は既に登録されています。")
        # 正規表現は新たに登録する必要がある。
        required_attrs = ['re_identifier', 're_accession', 're_organism']
        for required_attr in required_attrs:
            if not required_attr in attr:
                raise ValueError(required_attr + "は必須です。")
        new_class = type(classname, baseclass, attr)
        self.register(classname, new_class, regex)


class ProteinFastaBuilder(FastaBuilder):
    """ProteinFastaManager アミノ酸配列のFastaオブジェクトを作成する。

    あらかじめよく使うクラスを登録しておく。"""

    def __init__(self):
        """コンストラクタ"""
        self.regex = {}
        self.prototypes = {}
        self.default_prototype = fasta.BasicProteinFasta

        # nested classes cannot be pickled.
        # SwissProt
        # self.new_class('SwissProt', (fasta.BasicProteinFasta,), {
        #    're_identifier' : re.compile("^>sp\|[^\|]+\|(\S+) .*"),
        #    're_accession'  : re.compile("^>sp\|([^\|]+)\|.*"),
        #    're_organism'   : re.compile(".* OS=(\w+ \w+(?: [^=]+)?) .*"),
        #    }, re.compile("^>sp") )
        # # SwissProt
        # self.new_class('TrEMBLE', (fasta.BasicProteinFasta,), {
        #    're_identifier' : re.compile("^>tr\|[^\|]+\|(\S+) .*"),
        #    're_accession'  : re.compile("^>tr\|([^\|]+)\|.*"),
        #    're_organism'   : re.compile(".* OS=(\w+ \w+(?: [^=]+)?) .*"),
        #    }, re.compile("^>tr") )
        # # GenBank / RefSeq
        # self.new_class('GenBank_refseq', (fasta.BasicProteinFasta, ), {
        #    're_identifier' : re.compile("^>gi\|\d+\|ref\|([^\|]+)\|? .*"),
        #    're_accession'  : re.compile("^>gi\|(\d+)\|ref.*"),
        #    're_organism'   : re.compile(".* \[(\w+ \w+(?: \w+)?)\].*"),
        #    }, re.compile("^>gi\|\d+\|ref") )


class DataSetMaker(object):
    """DataSetMaker データセットを作るためのクラス。

    このクラスのサブクラスで実際に実装するよ！"""

    def __init__(self, builder, reader, data_type):
        """builder 個々の配列をどういう方法で構築するか、具体的に実装してあるクラス。
        現状だと、FastaManagerを使うかFastaMakerを使うか、選択できるようにする。
        reader ファイルなりデータベースなりを読み込むためのクラス。

        流れは、readerの属性値としてbuilderをセットする→readerはその
        builderを使って、データのリストを返す→このクラスがそれを
        data_typeに従って、DataSetクラスのインスタンスにして返す。

        両クラスとも、createメソッドで実際のFastaオブジェクトを生成している。
        今後新たなクラスを追加するときも、このインターフェースを
        継承すればおk。
        """
        self.builder = builder()
        self.reader = reader(builder)

    def read_from_file(self, filename, name='', label=None):
        """filenameのファイルから読み込む。

        readerがデータのリストを返してくれるので、それを使って
        各DataSetクラスは新しいインスタンスを生成する。
        """
        if not name:
            name, ext = os.path.splitext(os.path.basename(filename))
        data_list = self.reader.parse_file(filename)
        return self.data_type(data_list, name=name, origin=name, labels=label)

    def read_from_string(self, s, name='', label=None):
        """read sequences from a given string"""
        if not name:
            name = 'test'
        data_list = self.reader.parse_string(s)
        return self.data_type(data_list, name=name, origin=name, labels=label)

    def read_from_sqlite(self, dbname, query, name=""):
        """SQLite3データベースから、データセットを作成する。
        dbnameとqueryは必須。"""
        if isinstance(self.builder, type):
            self.builder = self.builder()
        data_list = []
        for r in self.reader.parse_sqlite(dbname, query):
            data_list.append(r)
        return self.data_type(data_list, name=name, origin=name)

    def create_new_maker(self, classname, baseclasses, attr):
        """自分で好きなDataSetMakerを作る。
        このとき、attrにはbuilderとreaderをセットする。"""
        return type(classname, baseclasses, attr)

    def import_from_db(self, dbname, username='', password=''):
        """データベースからデータセットをインポートする。"""
        data_list = self.reader.parse_db(dbname, username, password)
        return self.data_type(data_list)


class FastaDataSetMaker(DataSetMaker):
    """FastaDataSetMaker Fastaファイルから新しいデータセットを作る。
    """

    def __init__(self, builder=ProteinFastaBuilder,
                 reader=FastaReader,
                 data_type=dataset.FastaDataSet):
        """builderはFastaMakerあるいはFastaManager。
        """
        self.builder = builder
        self.reader = reader(builder())
        self.data_type = data_type
