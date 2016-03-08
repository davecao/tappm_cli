# -*- coding:utf-8 -*-

import tappm.hmm.hmm as hmm
import tappm.hmm.hmm_mp as hmm_mp
import tappm.hmm.util as hmmutil
import tappm.dataset
import numpy as np
import tappm.method as method


class MyHmmPredictor(method.Method):
    """MyHmmPredictor  A wrapper of my implementation of HMM.

    It offers converting models from and to those of GHMM,
    make the datasets into numerical form that suit my implementation."""

    def __init__(self, filename='', cpus=1,
                 valid_chars="ACDEFGHIKLMNPQRSTVWY"):
        '''Read an XML file of GHMM and convert it.'''
        self.method_name = 'hmm'
        self.model_file = filename
        self.method = None
        self.valid_chars = valid_chars
        self.valid_char_dic = {
            self.valid_chars[i]: i for i in range(len(self.valid_chars))}
        self.decoder = ""
        self.load(filename, cpus)

    def load(self, filename, cpus=1):
        """Read an XML file of GHMM."""
        (t, e, i) = hmmutil.load_ghmmxml(filename)
        if cpus == 1:
            self.method = hmm.HMM(t, e, i)
        elif cpus > 1:
            self.method = hmm_mp.MultiProcessHMM(t, e, i, worker_num=cpus)

    def initialize(self, cpus=1):
        """Reload hmm files"""
        if self.model_file:
            self.load(self.model_file, cpus)

    def predict(self, dataset, reverse=False, **args):
        """Predict (or Decode) a sequence by Viterbi algorithm."""
        dataset_tmp = self.convert_dataset(dataset, reverse)
        # i: identifier
        # d: (converted) data
        result_tmp = {i: self.method.viterbi(d, return_omega=True)
                      for i, d in list(dataset_tmp.items())}
        return self.convert_result(result_tmp, reverse=reverse)

    def train(self, dataset, reverse=False, if_debug=False, **args):
        """Train sequences using Baum-Welch algorithm."""
        dataset_tmp = self.convert_dataset(dataset, reverse)
        if if_debug:
            return self.method.baum_welch(
                list(dataset_tmp.values()), do_debug=True, **args)
        else:
            self.method.baum_welch(list(dataset_tmp.values()), **args)

    def convert_dataset(self, dataset, reverse=False, missing='ignore'):
        """Convert DataSet objects into numerical form.

        @param dataset  is a DataSet object.
        @param reverse  is a boolean"""
        converted = {}
        for seq in dataset:
            converted_tmp = []
            for c in seq.sequence:
                try:
                    converted_tmp.append(self.valid_char_dic[c])
                except KeyError:
                    if missing == 'ignore':
                        print("invalid character %s found." % c)
                    elif missing == 'error':
                        raise ValueError("Invalid character: " + c)
            if reverse:
                converted_tmp = converted_tmp[::-1]
            converted[seq.identifier] = converted_tmp
        return converted

    def convert_result(self, results, reverse=False):
        """Convert numerical representation into more readable form."""
        converted = {}
        for i, result in list(results.items()):
            converted_tmp = ""
            for n in result[0]:  # result[1] is a likelihood
                try:
                    converted_tmp += self.decoder[int(n)]
                except IndexError:
                    print("%d is out of range (only %d states registered)." %
                          (n, len(self.decoder)))
            if reverse:
                converted_tmp = converted_tmp[::-1]
            converted[i] = {'path': converted_tmp,
                            'pathnum': result[0],
                            'likelihood': result[1]}
            if len(result) > 2:
                converted[i]['omega'] = result[2]
        return converted

    def reset_valid_chars(self, chars=""):
        """Reset valid characters, which are used in convert_dataset method."""
        if len(chars) > 0:
            self.valid_chars = chars
        else:
            raise ValueError("chars must contain at least one character.")

    def decode(self, char):
        """Decode a numerical state into a character."""
        if self.decoder == "":
            raise ValueError("set the decoder list.")
        return self.decoder[char]

    def set_decoder(self, charlist):
        """Set decoder which is used in converting numerical states."""
        self.decoder = charlist

    def cross_valid(self, dataset, fold=5, is_random=True,
                    pseudocounts=[0, 0, 0], cpus=1, **args):
        """datasetを用いてクロスバリデーションを行う。
        """
        datasets_for_cv = dataset.cv(fold=fold, is_random=is_random)
        result_of_cv = {}
        for dic in datasets_for_cv:
            self.initialize(cpus)
            self.method.add_pseudocounts(pseudocounts)
            self.train(dic['train'], pseudocounts=pseudocounts, **args)
            result_of_cv.update(self.predict(dic['test'], **args))
        del(datasets_for_cv)  # 後片付け。fold数を増やしたときのために。
        return result_of_cv


class HMMResultSet(tappm.dataset.DataSet):
    """HMMResultSet  is a class that concatenates several results of viterbi.

    hoge"""
    def __init__(self):
        """Constructer method. This method does nothing.
        Use add_dataset method instead."""
        self.data_type = dict
        self.container = {}
        self.identifiers = []
        self.models = []
        self.seqnum = 0
        self.dataset_num = 0
        self.origins = []
        self.labels = {}

    def add_dataset(self, dataset, test='', model='', label=1):
        """Add a dataset to this object.

        @param dataset  is a dictionary"""
        self.dataset_num += 1
        self.models.append(model)
        for name, result in list(dataset.items()):
            # result is a dictionary which has two keys, likelihood and path.
            if name not in self.identifiers:
                self.identifiers.append(name)
                self.seqnum += 1
            if name not in self.container:
                self.container[name] = {model: result, 'origin': test}
                if test not in self.origins:
                    self.origins.append(test)
            else:
                self.container[name][model] = result

    def get_likelihood(self, name, model=''):
        """Return likelihood of a sequence with 'name'.

        if model is an empty string or a list, it returns likelihood values
        as a dictionary (if empty, values of all models will be returned)."""
        if model == '':
            l = {}
            for m in self.models:
                l[m] = self[name][m]['likelihood']
        elif isinstance(model, str):
            l = self[name][model]['likelihood']
        elif isinstance(model, (list, tuple, set)):
            l = {}
            for m in model:
                try:
                    l[m] = self[name][m]['likelihood']
                except KeyError:
                    print("No model %s is registered." % m)
                    continue
        return l

    def get_seqlen(self, name):
        """Return sequence length of a specified sequence."""
        return len(self[name][self.models[0]]['path'])

    def iter_by_origin(self, origin):
        """Return a generator for each sequence with 'origin'."""
        for i in self.identifiers:
            if self[i]['origin'] == origin:
                yield self[i]
            else:
                continue

    def iter_names_by_origin(self, origin):
        """REturn a generator that yields identifiers."""
        for i in self.identifiers:
            if self[i]['origin'] == origin:
                yield i
            else:
                continue

    def likelihood_dif(self):
        """Calculate differences of log likelihood values.

        Each calculated value is stored in
            self.container[id]['diff']['minuend-substrahend']
        The values that already exist are ignored when Calculating."""
        for i in range(len(self.models)):
            for j in range(i + 1, len(self.models)):
                m1 = self.models[i]
                m2 = self.models[j]
                keyname = m1 + "-" + m2
                keyname_rev = m2 + "-" + m1
                for name in self.identifiers:
                    if 'diff' not in self[name]:
                        self[name]['diff'] = {}
                    if keyname not in self[name]['diff'] and\
                       keyname_rev not in self[name]['diff']:
                        self[name]['diff'][keyname] = (
                            self.get_likelihood(name, m1) -
                            self.get_likelihood(name, m2)
                            ) / self.get_seqlen(name)

    def get_ldiff(self, name, minuend, substrahend):
        """Retrieve the difference of log likelihood"""
        keynames = [minuend + "-" + substrahend, substrahend + "-" + minuend]
        if keynames[0] in self[name]['diff']:
            return self[name]['diff'][keynames[0]]
        elif keynames[1] in self[name]['diff']:
            return self[name]['diff'][keynames[1]] * -1
        elif minuend in self.models and substrahend in self.models:
            self.likelihood_dif()
            self.get_ldiff(name, minuend, substrahend)
        elif minuend not in self.models:
            raise KeyError("there is no model %s" % minuend)
        elif substrahend not in self.models:
            raise KeyError("there is no model %s" % substrahend)

    def get_ldiffs(self, minuend, substrahend, origin=None):
        """Returns differences along all identifiers"""
        if origin is None:
            return [self.get_ldiff(name, minuend, substrahend)
                    for name in self.identifiers]
        elif isinstance(origin, str) and origin in self.origins:
            return [self.get_ldiff(name, minuend, substrahend)
                    for name in self.iter_names_by_origin(origin)]

    def roc(self, positive, negative, by='likelihood'):
        """Calculate likelihood by likelihood or differences of it.

        @param pos is a string or list that indicates its origin or IDs.
        @param neg is similar to pos."""
        # switch what to compare
        if by == 'likelihood':
            choose_seq = lambda name: self.get_likelihood(name)
        elif isinstance(by, (list, tuple)):
            choose_seq = lambda name: self.get_ldiff(name, by[0], by[1])
        # Positive data
        pos = self.get_subset(positive, choose_seq)
        neg = self.get_subset(negative, choose_seq)
        return roc(pos, neg)

    def get_subset(self, what, how=None):
        """Retrieve a subset specified by the variable.

        @param what  can be a string, list, or tuple.
        If being string, it should indicate the origin.
        List or tuple can specify multiple origins, as well as directly
        dictate identifiers to choose."""
        if how is None:
            how = lambda name: self.get_likelihood(name)
        if isinstance(what, str):  # Assume it indicates an origin.
            return [how(name) for name in self.iter_names_by_origin(what)]
        elif isinstance(what, (tuple, list)):
            if all([elem in self.identifiers for elem in what]):
                # Assume 'what' is a list of the identifiers
                return [how(name) for name in what]
            elif all([elem in self.origins for elem in what]):
                # Assume 'what' is a list of the origins
                l = []
                for origin in what:
                    l.extend([
                        how(name) for name in
                        self.iter_names_by_origin(origin)])
                return l


def roc(pos, neg):
    """Given positive and negative list, calculate and return
    an array that consists of true positive rate and false positive
    rate.
    This function returns (roc, auc, ber, thr), where roc is a tuple
    that contains two numpy arrays, true-positive rates and false-positive
    rates. auc is the AUC value itself. ber is the balanced error rate.
    thr is the threshold value that maximizes ber."""
    true = len(pos) + 0.0
    pos = np.sort(pos)[::-1]  # Descent order
    false = len(neg) + 0.0
    neg = np.sort(neg)[::-1]  # Descent order
    all_sorted = np.concatenate((pos, neg))
    all_sorted.sort()
    mds = np.concatenate(
        (np.array([all_sorted.min() - 1]),
         np.array(
            [(i + j) / 2 for i, j in
             zip(all_sorted[:-1], all_sorted[1:])]),
         np.array([all_sorted.max() + 1])))  # Ascent order
    tp = fp = 0
    ber = 1.0
    thr = 0
    ans = np.zeros((2, len(mds)))
    for x, n in zip(mds[::-1], list(range(len(mds)))):
        while tp < true and pos[tp] >= x:
            tp += 1
        while fp < false and neg[fp] >= x:
            fp += 1
        ber_new = (1 - tp / true + fp / false) / 2
        if ber_new < ber:
            ber = ber_new
            thr = x
        ans[:, n] = np.array([tp / true, fp / false])
    auc = np.array(
        [(ans[1][i] - ans[1][i - 1]) *
         (ans[0][i] + ans[0][i - 1]) / 2
         for i in range(1, len(ans[0]))]).sum()
    return {'roc': ans, 'auc': auc, 'ber': ber, 'thr': thr,
            'tp': sum(pos >= thr), 'fn': sum(pos < thr),
            'fp': sum(neg >= thr), 'tn': sum(neg < thr)}
