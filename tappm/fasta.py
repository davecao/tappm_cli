#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""Fasta  FASTA形式のデータを扱うためのクラス群。

Fastaクラスのインスタンス一つ一つが個々の配列を表す。
Fastaクラス自体は抽象クラスで、実際の配列はサブクラスの
どれかになる(予定)。
具体的には、その配列が格納されているデータベース毎
(e.g. SwissProt, GenBank, ...etc)に個々のサブクラスを
作り、サブクラス固有のメソッドとしてヘッダー行の解析や
データベースへのアクセスなどを考えている。"""

import re
import abc


class Fasta(object):
    __metaclass__ = abc.ABCMeta
    """Fasta Fasta形式の配列データを扱うための抽象クラス。

    ヘッダー行のパース、データベースへの接続などは
    サブクラスによって実装される。抽象クラスは、インスタンスを
    作れないクラスのことなので、メソッドを実装するのは問題ない
    ・・・はず。"""
    # ヘッダーをパースする際に用いる、正規表現オブジェクト。
    # 実際はサブクラスで定義されます。
    re_identifier = re.compile("")
    re_accession = re.compile("")
    re_organism = re.compile("")
    # ヘッダーをパースした際に、値が見つからない場合の
    # デフォルトの値。辞書形式にしておく。
    parser_default_values = {
            'identifier': 'NotAvailable',
            'accession': 'NotAvailable',
            'organism': 'NotAvailable',
            }

    @abc.abstractmethod
    def parse_header(self, header=''):
        """ヘッダー行をパースして、情報を取り出すためのメソッドです。
        抽象メソッドです。"""
        return

    @abc.abstractmethod
    def connect_db(self, query=''):
        """元のデータベースにこの配列に関するクエリーを投げます。
        抽象メソッドです。"""
        return

    @abc.abstractmethod
    def is_valid_char(self, char):
        """ある文字がそのFasta形式で許容される文字かどうかを判定します。
        抽象メソッドです。"""
        return

    @abc.abstractmethod
    def delete_invalid_char(self):
        """is_valid_charによって不適当と判断された文字を消去する。
        """
        return

    def __init__(self, sequence, header=""):
        """コンストラクタです。
        Fasta形式そのままのテキストか、ヘッダーと配列部分を別々に
        要求します。"""
        # print("sequence: {0}\nheader: {1}".format(sequence, header))
        # 引数の最初が'>'で始まっている＝ヘッダー行が含まれる
        if sequence[0] == '>':
            (self.header, self.sequence) = sequence.splitlines()
        elif header == "":
            raise ValueError("No header line found.")
        else:
            self.sequence = sequence
            self.header = header

        # (DNA|アミノ酸)配列中の改行文字を削除
        self.sequence = self.sequence.replace("\n", "")
        # ヘッダー行をパース
        self.parse_header()
        # 配列の長さもセットしておく。
        self.seqlen = len(self.sequence)

    def __repr__(self):
        """公式の文字列を返します。実態はFASTA形式の文字列です"""
        return self.identifier

    def __str__(self):
        """ほとんど__repr__と同じですが、適宜改行を含めます"""
        return self.header + "\n" + self.sequence

    def __lt__(self, other):
        """比較演算子で呼び出されるメソッドです。
        比較演算子で呼び出された場合、(DNA|アミノ酸)配列の長さを
        比べます。"""
        return len(self) < len(other)

    def __le__(self, other):
        """__lt__を参照。"""
        return len(self) <= len(other)

    def __gt__(self, other):
        """__lt__を参照。"""
        return len(self) > len(other)

    def __ge__(self, other):
        """__lt__を参照。"""
        return len(self) >= len(other)

    def __eq__(self, other):
        """__lt__などとは異なり、selfとotherが同じ配列かどうかを
        判定します。判定に際しては、配列が同じであるかどうかだけが
        考慮され、ヘッダー行の中身は考慮されません。"""
        return self.sequence == other.sequence

    def __ne__(self, other):
        """__eq__と同様です。"""
        return self.sequence != other.sequence

    def __len__(self):
        """配列の長さ(整数)を返します。"""
        return len(self.sequence)

    def __getitem__(self, key):
        """key番目の文字を返します。"""
        # if key > len(self):
        #    raise IndexError(key + "は大きすぎます。")
        # elif key < -1 * len(self):
        #    raise IndexError(key + "は小さすぎます。")
        # else:
        #    return self.sequence[key]
        return self.sequence[key]

    def __setitem__(self, key, value):
        """配列中のkey番目の文字を変更します。正しくない文字に
        変更しようとした場合、エラーを返します。"""
        if not self.is_valid_char(value):
            raise ValueError(value + "は不正な値です。")
        elif key > len(self):
            raise IndexError(key + "は大きすぎます。")
        elif key < -1 * len(self):
            raise IndexError(key + "は小さすぎます。")
        else:
            self.sequence[key] = value

    def __delitem__(self, key):
        """配列中のkey番目の文字を削除します。"""
        if key == 0:
            self.sequence = self.sequence[1:]
        elif key == -0:
            self.sequence = self.sequence[:-0]
        else:
            self.sequence = self.sequence[:key - 1] + self.sequence[key + 1:-1]
        # 配列の長さは１文字短くなる。
        self.seqlen -= 1

    def __iter__(self):
        """Iteratorオブジェクトを返します。反復処理は配列に対して
        行われます。"""
        self.current_index = 0
        return self

    def __next__(self):
        """Iterator処理に使われます。"""
        if self.current_index > self.seqlen:
            self.current_index += 1
            return self.sequence[self.current_index - 1]
        else:
            raise StopIteration

    def __contains__(self, regex):
        """特定の正規表現で表されるモチーフが配列中に含まれているか
        どうかを返します。"""
        if not isinstance(regex, type(self.re_accession)):
            regex = re.compile(regex)
        return regex.match(self.sequence)

    def get_seq(self):
        """配列を返す。"""
        return self.sequence

    def register_regex(self, target, regex):
        """ヘッダー行をパースするための正規表現を登録する。
        この部分を変更することで、個々のデータベースに対応した
        パースができるようになる。

        @target どのフィールドに対する正規表現なのかを指定する
        @regex  正規表現オブジェクトまたは文字列"""
        # targetが正しい値かチェック
        if not target in ['identifier', 'accession', 'organism']:
            raise ValueError("target がありません。" + target)
        # regexが正しい型がチェック
        if isinstance(regex, str):
            regex = re.compile(regex)
        elif isinstance(regex, type(re_identifier)):
            pass
        else:
            raise TypeError(regex + "の型が不正です。")
        if target == 'identifier':
            re_identifier = regex
        elif target == 'accession':
            re_accession = regex
        elif target == 'organism':
            re_organism = regex

    def show(self):
        """fasta形式の文字列を返す。"""
        return str(self.header + "\n" + self.sequence)

# ##XXX###XXX###XXX###XXX###XXX###XXX###XXX###XXX###XXX###XXX###XXX###XXX###
# 次の二つのクラスについては、基底クラスと言うよりインターフェースに
# した方がいいかも。インスタンスを生成するときに、どっちを継承するか
# 決められた方が、後々柔軟性があってやりやすそう。
# ##XXX###XXX###XXX###XXX###XXX###XXX###XXX###XXX###XXX###XXX###XXX###XXX###


class ProteinFasta(Fasta):
    """ProteinFasta アミノ酸配列を扱うためのクラスです。
    このクラスも抽象クラスで、個々のサブクラスで具体的な実装が
    期待されています。
    また、このクラスはBやZのような曖昧な文字を許可しません。"""
    __metaclass__ = abc.ABCMeta

    valid_chars = "ACDEFGHIKLMNPQRSTVWY"
    invalid_chars = "BJOUXZ"

    def is_valid_char(self, char):
        """charがアミノ酸配列として正しい文字であるかどうかを
        判定し、booleanを返します。"""
        return char in "ACDEFGHIKLMNPQRSTVWY"

    def delete_invalid_char(self):
        """不適当な文字を削除する。この実装は重いかもしれん
        高速化するなら、最初から文字のリストを決めうちにして、
        replaceを使えばいいかな"""
        for c in self.invalid_chars:
            self.sequence = self.sequence.replace(c, '')


class ExProteinFasta(ProteinFasta):
    """ExProteinFasta アミノ酸配列を扱うためのクラスです。

    このクラスも抽象クラスです。ProteinFastaクラスと異なり、
    このクラスのインスタンスはBやZのような曖昧な文字を配列中に
    持つことが出来ます。ただしJおよびO、Uは許可されません。"""
    __metaclass__ = abc.ABCMeta
    valid_chars = "ABCDEFGHIKLMNPQRSTVWXYZ"

    def is_valid_char(self, char):
        """charがアミノ酸配列(拡張版)として正しい文字であるかどうかを
        判定します。"""
        return char in "ABCDEFGHIKLMNPQRSTVWXYZ"


class DNAFasta(Fasta):
    """DNAFasta DNA配列を扱うためのクラス。

    このクラスも抽象クラス。"""
    __metaclass__ = abc.ABCMeta

    def is_valid_char(self, char):
        """charがDNA配列として正しい文字であるかどうかを判定する。
        Nのような拡張的な文字も受け付ける。"""
        pass


class BasicProteinFasta(ProteinFasta):
    """BasicProteinFasta 基本的なアミノ酸配列を扱うためのクラスです。

    このクラスは、特定のフォーマットに対応していないFasta配列を扱います。
    そのため、connect_dbおよびparse_headerの一部はサポートされていません。"""

    re_identifier = re.compile("^>(.+)$")
    re_accession = re.compile("^>([^\s^\|]+)[\s\|]+")
    re_organism = re.compile("^>([^\s^\|]+)[\s\|]+")

    def parse_header(self, header=''):
        """ヘッダー行をパースして情報を取り出します。
        このクラスでは、識別のための文字列しか取り出しません。"""
        if header == "":
            header = self.header
        matches = {}
        matches['identifier'] = self.re_identifier.findall(header)
        matches['accession'] = self.re_accession.findall(header)
        matches['organism'] = self.re_organism.findall(header)
        # 目的のものが見つからなかった場合に、デフォルトの値をセットする。
        for key in ('identifier', 'accession', 'organism'):
            if len(matches[key]) == 0:
                matches[key] = [self.parser_default_values[key]]
                if key == 'identifier':
                    raise InvalidValueWarning(
                        "No " + key + " found in the header line: " + header)
        self.identifier = matches['identifier'][0]
        self.accession = matches['accession'][0]
        self.organism = matches['organism'][0]

    def connect_db(self, query=''):
        """データベースに接続するためのメソッドです。
        まだ実装されていません。デフォルトでどこか有名なDBへ
        検索クエリーを投げるようにしてもいいかもです。"""
        pass


class SwissProt(BasicProteinFasta):
    """SwissProt SwissProtデータベースのfastaファイルを扱う。

    SwissProtのFastaファイルを扱うクラス。とりあえずタンパク質の方に
    対応している。動的にタンパク質とDNAを切り替えられるようにした方が
    いいのかもしれない。どうせDNA使わないからあれなんだけど。"""

    re_identifier = re.compile("^>sp\|[^\|]+\|(\S+) ")
    re_accession = re.compile("^>sp\|([^\|]+)\|.*")
    re_organism = re.compile(".* OS=(\w+ \w+(?: \w+)?) ")


class TrEMBLE(BasicProteinFasta):
    """TrEMBLE TrEMBLEデータベースのfastaファイルを扱う。

    ほとんどSwissProtと一緒だけど、こちらは機械的にアノテーション
    されたもの。まーフォーマットは最初の二文字が違うってだけなんだけど。"""

    re_identifier = re.compile("^>tr\|[^\|]+\|(\S+) ")
    re_accession = re.compile("^>tr\|([^\|]+)\|.*")
    re_organism = re.compile(".* OS=(\w+ \w+(?: \w+)?) ")


class GenBank_refseq(BasicProteinFasta):
    """GenBank_refseq GenBankフォーマットのもののうち、ReferenceSequenceの
    データベースの配列のためのクラス。

    他にもGenBank_gbやらGenBank_eucみたいなのがあるはず。"""

    re_identifier = re.compile("^>gi\|\d+\|ref\|([^\|]+)\|? ")
    re_accession = re.compile("^>gi\|(\d+)\|ref.*")
    re_organism = re.compile(".* \[(\w+ \w+(?: \w+)?)\].*")


class InvalidValueWarning(Warning):
    """InvalidValueWarning 不正な値を検出したが、処理を続行したいときの警告。

    ValueErrorをraiseするほどでもないが、不正な値を検出したという時に
    使うための警告。"""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
