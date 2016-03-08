# -*- coding: utf-8 -*-
# @Author: davecao
# @Date:   2016-01-18 17:15:06
# @Last Modified by:   davecao
# @Last Modified time: 2016-01-20 20:25:15

import os
import numpy as np

from time import localtime, strftime
from tappm import jinja2_ENV, VERSION
from tappm.utils import wrapText

TPL_MAP = {
    'text': 'report.txt',
    'tabular': 'report.tab',
    'xml': 'report.xml'
}

TPL_ITEM_MAP = {
    'text': 'item.txt',
    'tabular': 'item.tab',
    'xml': 'item.xml'
}


class ResultItems(object):
    def __init__(self, identifier=None, description=None, sequence=None,
                 vpath={}, score=-np.inf, omega=[], likelihood=None,
                 pathnum=[], likelihood_mp=None, has_tmd=False, isTA=False,
                 threshold=-0.016722, tmd_position=None,
                 hasTable=False, tableColName=[], hasImg=False,
                 img_path=None, img_name=None,
                 tplName="item.html"):

        super(ResultItems, self).__init__()
        CterDist = 50
        self.identifier = identifier if identifier else ''
        self.description = description if description else ''
        self.sequence = sequence
        self.omega = omega
        self.vpath = vpath
        self.score = score
        self.omega = omega
        self.likelihood = likelihood
        self.likelihood_mp = likelihood_mp
        self.pathnum = pathnum
        self.has_tmd = has_tmd
        self.isTA = isTA
        self.tmd_position = tmd_position
        self.TAprotein = False
        # "".join(str(tmd_position)).strip('[]')
        # str(tmd_position).strip('[]')
        self.tmd_POS = None
        self.NumOfTMD = len(tmd_position) if tmd_position else 0
        self.CterTMDPos = None
        self.threshold = threshold
        self.hasTable = hasTable
        self.tableColName = tableColName
        self.hasImg = hasImg
        self.img_path = img_path if img_path else ''
        self.img_name = img_name if img_name else ''

        if tmd_position:
            self.tmd_POS = ';'.join(map(str, tmd_position))

        if has_tmd:
            # find the TMD segments closed to C terminus
            Cter = len(self.sequence) - CterDist
            tmp = [a for a in tmd_position if a[0] >= Cter]
            if tmp:
                self.CterTMDPos = ', '.join(map(str, tmp))
        if isTA and has_tmd and\
           (self.NumOfTMD == 1) and (self.CterTMDPos is not None):
            self.TAprotein = True
        self.template = jinja2_ENV.get_template(tplName)

    def __get_elements__(self):
        elems = {
            'name': self.identifier,
            'description': self.description,
            'sequence': self.sequence,
            'wrapCotent': zip(wrapText(self.sequence, width=60).split(),
                              wrapText(self.vpath, width=60).split()),
            'path': self.vpath,
            'likelihood': self.likelihood,
            'likelihood_mp': self.likelihood_mp,
            'TAprotein': self.TAprotein,
            'score': self.score,
            'isTA': self.isTA,
            'seqLen': len(self.sequence),
            'threshold': self.threshold,
            'hasTMD': self.has_tmd,
            'NumOfTMD': self.NumOfTMD,
            'tmd_POS': self.tmd_position,
            'CterTMDPos': self.CterTMDPos,
            'hasImg': self.hasImg,
            'hasTable': self.hasTable,
            'id': self.identifier,
            'tableColName': self.tableColName
        }
        return elems

    def render(self, format='html'):
        return self.template.render(self.__get_elements__())


def do_rollover(outfile, backupCount=5):
    """ Determine if rollover should occur """
    # odir = os.path.dirname(os.path.abspath(outfile))
    # fname = os.path.basename(outfile)
    if backupCount > 0:
        for i in range(backupCount - 1, -1, -1):
            if i == 0:
                sfn = outfile
            else:
                sfn = "%s.%d" % (outfile, i)
            dfn = "%s.%d" % (outfile, i + 1)
            if os.path.exists(sfn):
                os.rename(sfn, dfn)


def report(resultItemsList, outfile, fmt='html', backupCount=5,
           ncpu=4, threshold=-0.0167222981, totalSeq=0):
    tplName = TPL_MAP[fmt]
    tpl = jinja2_ENV.get_template(tplName)
    # do roll-over
    do_rollover(outfile, backupCount=backupCount)
    # Add resultItems
    block = ''
    for item in resultItemsList:
        block += item.render()

    with open(outfile, 'w') as fout:
        fout.write(
            tpl.render({
                'package': "TAPPM ver. " + VERSION,
                'version': VERSION,
                'rptTime': strftime("%d %b %Y %H:%M:%S", localtime()),
                'threshold': threshold,
                'ncpu': ncpu,
                'totalSeq': totalSeq,
                'body_content': block,
            }))
