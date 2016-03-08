#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""Method HMMやSVMなどのラッパー。正確にはラッパーのラッパー。

目的は、Modelクラスから共通のAPIで操作を可能にすること。"""


class Method(object):
    """Method 予測手法の抽象クラス。

    train DataSetオブジェクトを引数にとり、モデルを訓練する。
    predict DataSetオブジェクトを引数にとり、予測する。結果を
    表すResultクラスのオブジェクトみたいなのを返す。
    cross_valid 同じくDataSetオブジェクトを引数にとり、クロス
    バリデーションを行う。返す値はpredictと同じ。
    save モデルをセーブする(ファイルに書き出す)。
    load モデルをロードする(ファイルから読み込む)。"""

    def __init__(self):
        """コンストラクタ。
        実装はサブクラスに。"""
        self.method = None

    def initialize(self):
        """モデルを初期化する。たとえば、クロスバリデーションを
        行う際に、一回目が終わった後にモデルを初期化しておいて、
        二個目のデータセットで再訓練する、というような使い方を
        想定している。実際にどう初期化するかはサブクラスに任せる。
        自分でモデルを作っているクラスならそうするし、そうじゃない
        ならセーブしてあるファイルからロードする。"""
        pass

    def train(self, dataset):
        """dataset DataSetオブジェクト。今のところ、このメソッド内で
        手法特有のデータセットに変換する方式を取る予定。あまりに遅い
        などの弊害があれば、あらかじめ作っておいたデータを用いるよう
        にする"""
        # データセットを、適用可能な形に変換する。
        dataset_tmp = self.convert_dataset(dataset)
        # 変換したデータセットを用いて訓練する。このメソッド名は
        # 変更する必要があるだろう。
        result_tmp = self.method.train(dataset_tmp)
        # (あれば)帰ってきた結果を、わかりやすい形に変換する
        result = self.convert_result(result_tmp)
        return result

    def convert_dataset(self, dataset):
        """dataset DataSetオブジェクト。一般的なDataSetオブジェクトを
        その手法に適用できる形に整える。"""
        pass

    def convert_result(self, result):
        """帰ってきた結果を、共通のAPIで扱えるように調整する。
        """
        pass

    def predict(self, dataset):
        """datasetを用いて予測を行う。
        """
        dataset_tmp = self.convert_dataset(dataset)
        result_tmp = self.method.predict(dataset_tmp)
        return self.convert_result(result_tmp)

    def cross_valid(self, dataset, fold=5, is_random=True, **args):
        """datasetを用いてクロスバリデーションを行う。
        """
        datasets_for_cv = dataset.cv(fold=fold, is_random=is_random)
        result_of_cv = None
        for dic in datasets_for_cv:
            self.initialize()
            self.train(dic['train'], **args)
            r = self.predict(dic['test'], **args)
            if result_of_cv == None:
                result_of_cv = r
            else:
                result_of_cv.merge(r)
        del(datasets_for_cv)  # 後片付け。fold数を増やしたときのために。
        return result_of_cv

    def save(self, filename):
        """ファイルにモデルを保存する。"""
        pass

    def load(self, filename):
        """ファイルに保存してあるモデルを読み込む。
        """
        pass

    def grid(self, dataset, prange=[]):
        '''datasetに対して、グリッドサーチを行う。
        グリッドサーチ(Grid search)は、ハイパーパラメータを網羅的に探索
        して最適化する手法'''
        pass
