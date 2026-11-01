# based on /lmfdb/elliptic_curves/ec_stats.py
# Authors: David Neill Asanza, Albert Ford, Ngi Nho, Jen Paulhus, Kevin Wang

from flask import url_for
from lmfdb import db
from psycodict import SQL
from sage.all import lazy_attribute

from lmfdb.utils import comma, display_knowl, StatsDisplay


def _key(k):
    # count(..., groupby=[col]) is indexed by tuples
    return k[0] if isinstance(k, tuple) else k


class HGCWAstats(StatsDisplay):
    """
    Statistics for families of higher genus curves with automorphisms.

    The distributions below are computed from hgcwa_topological by StatsDisplay;
    the table by genus (``by_genus_data``) sums the counts stored in hgcwa_families.
    """
    table = db.hgcwa_topological
    baseurl_func = '.index'

    short_display = {'genus': 'genus',
                     'group_order': 'group order',
                     'dim': 'dimension',
                     'r': 'number of branch points',
                     'full_index': r'index $[\mathrm{Aut}:G]$',
                     'num_braid': 'braid classes',
                     'num_vectors': 'generating vectors'}
    top_titles = {'group_order': 'group orders',
                  'dim': 'dimensions',
                  'r': 'numbers of branch points',
                  'full_index': r'indices $[\mathrm{Aut}:G]$ of $G$ in the full automorphism group',
                  'num_braid': 'numbers of braid classes',
                  'num_vectors': 'numbers of generating vectors'}
    knowls = {'genus': 'ag.curve.genus',
              'group_order': 'group.order',
              'dim': 'curve.highergenus.aut.dimension',
              'r': 'curve.highergenus.aut.branchpoints',
              'full_index': 'curve.highergenus.aut.full_index',
              'num_braid': 'curve.highergenus.aut.braid_equivalence',
              'num_vectors': 'curve.highergenus.aut.generatingvector'}
    buckets = {'group_order': ['2-9', '10-19', '20-49', '50-199', '200-999', '1000-'],
               'num_braid': ['1', '2', '3-5', '6-10', '11-20', '21-100', '101-1000', '1001-'],
               'num_vectors': ['1', '2', '3-5', '6-10', '11-20', '21-100', '101-1000', '1001-']}
    stat_list = [
        {'cols': ['r'],
         'intro': ["The number of topological classes of actions by the number $r$ of branch points; the family has dimension $r-3$."]},
        {'cols': ['dim']},
        {'cols': ['group_order']},
        {'cols': ['full_index'],
         'intro': [r"The index is $1$ exactly when $G$ is the full automorphism group of the generic curve in the family."]},
        {'cols': ['num_braid'],
         'intro': ["The number of topological classes with a given number of braid classes."]},
        {'cols': ['num_vectors'],
         'intro': ["The number of topological classes with a given number of generating vectors."]},
    ]

    def __init__(self):
        self.topological_knowl = display_knowl('curve.highergenus.aut.topological_equivalence', title='topological classes')
        self.braid_knowl = display_knowl('curve.highergenus.aut.braid_equivalence', title='braid classes')
        self.vectors_knowl = display_knowl('curve.highergenus.aut.generatingvector', title='generating vectors')
        self.families_knowl = display_knowl('curve.highergenus.aut.family', title='families')

    @lazy_attribute
    def by_genus_data(self):
        """
        One row per genus: the number of groups, families, topological classes
        (all, full, hyperelliptic, cyclic trigonal), braid classes and generating vectors.
        """
        cur = db._execute(SQL(
            'SELECT genus, COUNT(DISTINCT "group"), COUNT(*), SUM(num_topological), '
            'SUM(num_topological_full), SUM(num_braid), SUM(num_vectors) '
            'FROM hgcwa_families GROUP BY genus ORDER BY genus'))
        hyp = {_key(k): v for k, v in db.hgcwa_topological.count({'hyperelliptic': True}, groupby=['genus']).items()}
        trig = {_key(k): v for k, v in db.hgcwa_topological.count({'cyclic_trigonal': True}, groupby=['genus']).items()}
        rows = []
        for genus, groups, families, top, full, braid, vectors in cur:
            rows.append({'genus': genus, 'groups': groups, 'families': families,
                         'topological': int(top), 'full': int(full),
                         'hyperelliptic': hyp.get(genus, 0), 'cyclic_trigonal': trig.get(genus, 0),
                         'braid': int(braid), 'vectors': int(vectors)})
        return rows

    @lazy_attribute
    def totals(self):
        keys = ['families', 'topological', 'full', 'hyperelliptic', 'cyclic_trigonal', 'braid', 'vectors']
        tot = {key: sum(row[key] for row in self.by_genus_data) for key in keys}
        tot['groups'] = db.hgcwa_families.count_distinct('group')
        return tot

    @lazy_attribute
    def genus_list(self):
        return [row['genus'] for row in self.by_genus_data]

    @lazy_attribute
    def hyperelliptic_genera(self):
        return [row['genus'] for row in self.by_genus_data if row['hyperelliptic']]

    @lazy_attribute
    def cyclic_trigonal_genera(self):
        return [row['genus'] for row in self.by_genus_data if row['cyclic_trigonal']]

    def _counts(self):
        t = self.totals
        return (comma(t['topological']), self.topological_knowl, comma(t['full']),
                comma(t['families']), self.families_knowl,
                comma(t['braid']), self.braid_knowl, comma(t['vectors']), self.vectors_knowl,
                min(self.genus_list), max(self.genus_list))

    @property
    def short_summary(self):
        return (
            r'The database currently contains %s %s of group actions on curves (%s of them with $G$ the full '
            r'automorphism group) in %s %s, made of %s %s and %s %s, for genus $%s \le g \le %s$ '
            r'and quotient genus $0$. Here are some <a href="%s">further statistics</a>.'
            % (self._counts() + (url_for('.statistics'),)))

    @property
    def summary(self):
        return (
            r'The database currently contains %s %s of group actions on curves (%s of them with $G$ the full '
            r'automorphism group) in %s %s, made of %s %s and %s %s, for genus $%s \le g \le %s$ '
            r'and quotient genus $0$.' % self._counts())
