# -*- coding: utf-8 -*-
# @Author: davecao
# @Date:   2016-01-18 15:42:33
# @Last Modified by:   davecao
# @Last Modified time: 2016-01-31 12:34:23
import re
from itertools import chain
from tappm.io.report import ResultItems, TPL_ITEM_MAP


def convert_numpy_types(predicted, predicted_mp, threshold, fmt, fastalist):
    """ Prepare for render """
    resultItemsList = []
    tmd_15H = 'HHHHHHHHHHHHHHH'
    append = resultItemsList.append
    for seq_id, dic in predicted.items():
        identifier = seq_id[:seq_id.find(' ')]
        description = seq_id
        sequence = ''
        for x in chain(fastalist):
            if x.identifier == seq_id:
                sequence = x.sequence
                break
        vpath = dic['path']
        likelihood = dic['likelihood'].item()
        likelihood_mp = predicted_mp[seq_id]['likelihood'].item()
        score = (likelihood - likelihood_mp) / len(vpath)
        omega = [i.item() for i in dic['omega']]
        pathnum = [i.item() for i in dic['pathnum']]
        has_tmd = tmd_15H in vpath  # has >=15 continuous H
        tmd_position = None
        if has_tmd:
            tmd_position = [(a.start(), a.end())
                            for a in list(re.finditer(tmd_15H, vpath))]
        isTA = score >= threshold

        resultItems = ResultItems(
            identifier=identifier,
            description=description,
            sequence=sequence,
            vpath=vpath,
            score=score,
            omega=omega,
            likelihood=likelihood,
            likelihood_mp=likelihood_mp,
            pathnum=pathnum,
            has_tmd=has_tmd,
            tmd_position=tmd_position,
            isTA=isTA,
            threshold=threshold,
            tplName=TPL_ITEM_MAP[fmt])
        append(resultItems)
    return resultItemsList
