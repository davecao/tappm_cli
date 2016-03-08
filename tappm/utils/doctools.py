# -*- coding: utf-8 -*-
# @Author: davecao
# @Date:   2016-01-19 11:32:47
# @Last Modified by:   davecao
# @Last Modified time: 2016-01-19 11:35:16

from textwrap import wrap

__all__ = ['wrapText']


def wrapText(text, width=70, join='\n', **kwargs):
    """Return wrapped lines from :func:`textwrap.wrap` after *join*\ing them.
    """

    try:
        indent = kwargs.pop('indent')
    except KeyError:
        pass
    else:
        kwargs['initial_indent'] = kwargs['subsequent_indent'] = ' ' * indent
    if join:
        return join.join(wrap(text, width, **kwargs))
    else:
        return wrap(text, width, **kwargs)
