#!/usr/bin/env python
# -*- coding:utf-8 -*-

import random
import copy
from tappm import fasta


class DataSet(object):
    """DataSet タンパク質/DNAの配列のデータセットを表すクラス。

    このクラスはabstractなクラス。実装はサブクラスに。
    ProteinDataSet、DNADataSet、RNADataSetなどなど。。。
    また、インスタンスの生成に関してはDataSetMakerに丸投げする。
    Factory Methodパターンを採用。もっといいのがあれば変えるかも。

    やっぱりただのTemplate Methodパターンにする。
    「どこからデータが来たのか」をこのクラスで管理するのは、
    お門違いな感じがするため。そういうのはdataset_makerに任せてしまう。"""

    def __init__(self):
        """コンストラクタ。データのリストが来たら、それをデータセット
        クラスのインスタンスによしなにしてやる。
        将来的には、ストリームデータとかにも対応できるとよい。
        """
        self.current_index = 0
        self.seqnum = 0
        # self.container   = {} # 実際にデータを格納する場所.一応辞書を想定。
        # self.seqnum      = 0# 格納されているデータの数
        # self.data_type   = SomeType   # どんなデータを格納しているか
        # self.identifiers = [] # 識別用。
        # self.name        = "" # 名前。
        # self.labels       = {} # ラベル
        # self.attributes  = [self.seqnum, self.data_type, self.identifiers,
        #                     self.name, self.labels]
        # リストを想定している(removeメソッドが使える必要あり)

    def __len__(self):
        """データの数を返す。"""
        return self.seqnum

    def __contains__(self, key):
        """keyが存在しているかを判定する。具体的にはidentifierで
        判断する。"""
        return key in self.identifiers

    def __getitem__(self, key):
        """keyに対応するデータを取り出す。keyとして文字列(ID)または
        整数値(インデクス)を受け付ける。両方あり得る場合、数字の方が
        優先される。"""
        if isinstance(key, list) or isinstance(key, tuple):
            return [self[k] for k in key]
        if isinstance(key, int):
            return self.container[self.identifiers[key]]
        elif key in self.identifiers:
            return self.container[key]
        else:
            raise IndexError(key + "に対応する値は存在しません。")

    def __setitem__(self, key, val):
        """keyに対応するデータを登録する。
        既に登録済みの場合はエラーを出力し、keyはIDとしてあつかう。"""
        if isinstance(key, list) and isinstance(val, list):
            for i in range(len(key)):
                self[key[i]] = val[i]
        if key in self.identifiers:
            raise IndexError(key + "は既に存在しています。")
        elif not isinstance(val, self.data_type):
            print(type(val))
            raise TypeError(str(val) + " の型が不正です。")
        else:
            self.container[key] = val
            self.identifiers.append(key)
            self.seqnum += 1

    def __delitem__(self, key):
        """keyに対応するアイテムを削除する。"""
        if isinstance(key, list) or isinstance(key, tuple):
            for k in key:
                del(self[k])
        if key in self.identifiers:
            self.identifiers.remove(key)
            del(self.container[key])
            self.seqnum -= 1
        elif key < self.seqnum:
            del(self.container[self.identifiers[key]])
            del(self.identifiers[key])
            self.seqnum -= 1
        else:
            raise IndexError(key + "がありません。")

    def __iter__(self):
        """Iterationを作成する。個々のデータを取り出す。
        """
        self.current_index = 0
        return self

    def __next__(self):
        """イテレータの一部。次の要素を返す。
        """
        if self.current_index < self.seqnum:
            self.current_index += 1
            return self.get_by_index(self.current_index - 1)
        else:
            raise StopIteration

    def __str__(self):
        '''REturn string'''
        txt = 'Name: ' + self.name + "\n"
        txt += '\tDataType: ' + str(self.data_type) + "\n"
        txt += '\tSeqNum  : ' + str(self.seqnum) + "\n"
        txt += '\tID 00000: ' + self.identifiers[0] + "\n"
        txt += '\t...\n'
        txt += "\t...{0} sequence(s)".format(len(self))
        return txt

    def __repr__(self):
        txt = '<DataSet>\n' + self.__str__()
        return txt

    def _calc_stat(self):
        """コピーなどでseqnumなどの内部情報と実際のデータの数との
        間に齟齬が生じた場合に、適切な値を計算し直す。"""
        self.seqnum = len(self.identifiers)

    def get_by_index(self, i):
        """内部的な配列内の順番で取り出す。
        """
        return self.container[self.identifiers[i]]

    def copy(self, idlist, do_deepcopy=False):
        """idlist内にあるデータをコピーして、新たなデータセットを
        返す。新たなリストを作って、そこからインスタンスを生成
        している。"""
        seq_list = self.copy_container(idlist, do_deepcopy)
        labels = self.copy_labels(idlist)
        return self.__class__(seq_list, name=self.name, labels=labels)
        # for identifier in idlist:
        #    if do_deepcopy:
        #        seq_list.append(copy.deepcopy(self[identifier]))
        #    else:
        #        seq_list.append(self[identifier])
        # return self.__class__(seq_list, name=self.name)

    def copy_container(self, idlist, do_deepcopy=False):
        '''Return the copy of self.container'''
        if do_deepcopy:
            return [copy.deepcopy(self.container[id]) for id in idlist]
        else:
            return [self.container[id] for id in idlist]

    def copy_labels(self, idlist):
        '''Return copy version of self.labels'''
        # return {id: self.get_label(id) for id in idlist}
        label = {}
        for key in idlist:
            label[key] = self.get_label(key)
        return label

    def sample(self, num):
        """データセットからnum個のデータをサンプリングする。"""
        l = random.shuffle(list(range(len(self))))
        for i in range(num):
            print(l)
            pass

    def split_to(self, num, is_random=True, do_copy=False):
        """データセットをnum個に分けて、同じくDatasetクラスの
        インスタンスのリスト(あるいはタプル)として返す。
        randomがTrueの時は、一回順番をシャッフルしてから
        splitする。"""
        if is_random:
            idlist = self.identifiers[:]
            random.shuffle(idlist)
        else:
            idlist = self.identifiers[:]
        generated_datasets = []
        size_per_subset = self.seqnum / num
        fraction = self.seqnum % num
        splitted = 0
        for i in range(num):
            start = splitted
            end = splitted + size_per_subset
            if fraction > 0:
                end += 1
                fraction -= 1
            splitted = end
            # generated_datasets.append(self.copy(range(start,end), do_copy))
            generated_datasets.append(self.copy(idlist[start:end], do_copy))
        return generated_datasets

    def split_by(self, num, is_random=True, do_copy=False):
        """データセットをnum個ずつのサブ・データセットに
        分けて、そのリストあるいはタプルを返す。
        split_toと同じように、random=Trueならシャッフルしてから
        分割を行う。"""
        if is_random:
            idlist = self.identifiers[:]
            random.shuffle(idlist)
        else:
            idlist = self.identifiers[:]
        generated_datasets = []
        splitted = 0
        while splitted < self.seqnum:
            start = splitted
            end = splitted + num
            splitted = end
            generated_datasets.append(self.copy(idlist[start:end], do_copy))
        return generated_datasets

    def prepare_cross_validation(self, fold=5, is_random=True,
                                 do_deepcopy=False):
        """上のsplit_*系の応用となるメソッドで、
        何フォールドかを指定することで、それに適したデータセットの
        サブセットを返してくれる。辞書のリストを返す。
        [{'test':dataset_1, 'train':dataset_2}, {'test':dataset_3, ..
        のような感じ。"""
        splitted = self.split_to(fold, is_random, do_deepcopy)
        cv_data = []
        for i in range(fold):
            train = None
            for j in range(fold):
                if i == j:
                    next
                elif train is None:
                    # ここでcopyを使わないと、splitted[j]のサイズが
                    # どんどん増えてしまう。
                    train = copy.deepcopy(splitted[j])
                else:
                    train.merge(splitted[j], do_deepcopy)
            cv_data.append({'test': splitted[i], 'train': train})
        return cv_data

    def cv(self, fold=5, is_random=True):
        """prepare_cross_validationへのエイリアス。名前長い(；´Д｀)"""
        return self.prepare_cross_validation(fold, is_random)

    def merge(self, other, do_copy=False, verbose=False):
        """他のデータセット(一つ)とマージする。
        データの種類が違うとき(e.g.アミノ酸配列と塩基配列)は
        エラーを出力する。
        複数のデータセットをマージしたい時はmerge_allを参照
        IDは、self[key] = ...のところで自動的に追加されている。
        たぶん配列の数も同じだと思う。。。"""
        if self.data_type != other.data_type:
            raise TypeError(other + "のデータの種類が違います。")
        if verbose:
            print(self, other)
        for key in other.identifiers:
            if key in self.identifiers:
                continue
            else:
                if do_copy:
                    self[key] = copy.copy(other[key])
                else:
                    self[key] = other[key]
                if key in other.labels:
                    self.labels[key] = other.get_label(key)
                # self.identifiers.append(key)
        self._calc_stat()
        if verbose:
            print(self)

    def merge_all(self, others, return_new=False):
        """複数のデータセットを一つにまとめる。

        return_new Trueの場合、新たなインスタンスとして返す。
        """
        if not isinstance(others, list) or not isinstance(others, tuple):
            raise TypeError(others + "はリストである必要があります。")
        if return_new:
            original = copy.deepcopy(self)
        else:
            original = self
        for dataset in others:
            original.merge(dataset)
        if return_new:
            original

    def set_name(self, name):
        """データセットに名前をつける。"""
        self.name = name

    def set_label(self, id, label):
        """データセット中に、ラベルを設定する。"""
        if not self.labels:
            self.labels = {}
        if id in self.labels:
            self.labels[id] = label
        elif isinstance(id, int) and id < self.seqnum:
            self.labels[self.identifiers[id]] = label
        else:
            self.labels[id] = label

    def set_labels(self, label):
        """データセット中のデータすべてに、labelを設定する。"""
        for id in self.identifiers:
            self.labels[id] = label

    def get_label(self, index):
        """index番目のデータのラベルを取得する。"""
        if index == None:
            return self.labels
        elif index in self.labels:
            return self.labels[index]
        elif isinstance(index, int):
            return self.labels[self.identifiers[index]]

    def get_labels(self, type='list'):
        """すべてのラベルをリストとして返す。
        順番は、self.identifiersに入っているものと同様。"""
        return [self.get_label(id) for id in self.identifiers]

    def convert2num(self, charlist):
        """文字列を数字にして返す。"""
        char_dic = {charlist[i]: i for i in range(len(charlist))}
        return [[char_dic[c] for c in seq] for seq in self]


class FastaDataSet(DataSet):
    """FastaDataSet Fasta配列のデータセット。

    データセットを二つに分けたり、
    クロスバリデーション用にデータを分割して提供したり
    とかそんな感じ。"""

    def __init__(self, seq_list, name='', origin='', labels=None):
        """seq_listから、新しいデータセットを作る。
        seq_listはリスト(あるいはタプルでもおk？)で、各要素は
        self.data_typeと一致していなければならない。"""
        self.identifiers = [seq.identifier for seq in seq_list]
        self.container = {}
        self.data_type = fasta.Fasta
        self.seqnum = 0
        for seq in seq_list:
            self.container[seq.identifier] = seq
            self.seqnum += 1
        if name:
            self.name = name
        if origin:
            self.origin = origin
        elif name:
            self.origin = name
        if labels != None and not type(labels) in (list, tuple, dict):
            self.labels = {}
            for key in self.identifiers:
                self.labels[key] = labels
        elif isinstance(labels, dict):
            self.labels = labels
        else:
            self.labels = {}
            self.set_labels(0)


class DNADataSet(DataSet):
    """DNADataSet DNA配列を扱うデータセット。

    出来ることはタンパク質版のと大体一緒。
    翻訳してProteinDataSetのインスタンスを作れるようにしても
    おもしろいかも(たぶんそこまで出来ないけど)"""
    pass


class RNADataSet(DataSet):
    """RNADataSet RNA配列を扱うデータセット。

    ほぼDNADataSetと同じ。たぶん使わないし作らない。"""
    pass
