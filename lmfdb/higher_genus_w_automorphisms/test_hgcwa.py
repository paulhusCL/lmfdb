from lmfdb.tests import LmfdbTest
from lmfdb.higher_genus_w_automorphisms.hgcwa_helpers import (
    label_kind, split_label, parent_label, relative_label, label_sort_key, vector_labels,
    family_label_parts, unrank_perm, rank_perm, cycle_string, perm_order, compose_lr,
    is_generating_vector_product_one, signature_display, jacobian_tex, class_tuple_summary,
    parse_signature, parse_periods, parse_group_labels)

# The page tests use labels from the test sample (HGCWA_redo/runs/sample40), which
# is a subset of the full data.

BASE = '/HigherGenus/C/Aut/'


class HigherGenusWithAutomorphismsTest(LmfdbTest):

    # Pure Python helpers (no database needed)

    def test_labels(self):
        assert label_kind("4.5-1.0.5-5-5-5") == "family"
        assert label_kind("4.5-1.0.5-5-5-5.2") == "topological"
        assert label_kind("21.32-47.0.2-2-4-4-4.1.17") == "braid"
        assert label_kind("10.6-1.0.2-2-2-2-2-2-2-2-2-2.1.1.3280") == "vector"
        assert label_kind("junk") is None
        assert label_kind("4.5-1.0.5-5-5-5.x") is None
        assert split_label("10.6-1.0.2-2-2-2-2-2-2-2-2-2.1.1.3280") == ("10.6-1.0.2-2-2-2-2-2-2-2-2-2", [1, 1, 3280])
        assert parent_label("4.5-1.0.5-5-5-5.2.1.1") == "4.5-1.0.5-5-5-5.2.1"
        assert parent_label("4.5-1.0.5-5-5-5") is None
        assert relative_label("21.32-47.0.2-2-4-4-4.1.17") == "1.17"
        assert family_label_parts("46.1080-260.0.3-3-4") == (46, "1080.260", [0, 3, 3, 4])
        labels = ["x.1-1.0.2-2.1.10", "x.1-1.0.2-2.1.2", "x.1-1.0.2-2.1.1"]
        assert sorted(labels, key=label_sort_key) == ["x.1-1.0.2-2.1.1", "x.1-1.0.2-2.1.2", "x.1-1.0.2-2.1.10"]
        assert vector_labels([{"label": "a.b", "num_vectors": 2}, {"label": "a.c", "num_vectors": 1}]) == ["a.b.1", "a.b.2", "a.c.1"]

    def test_permutations(self):
        # The LMFDB encoding: rank in the lexicographic order of Sym(d), from 0
        assert unrank_perm(0, 4) == [1, 2, 3, 4]
        assert unrank_perm(23, 4) == [4, 3, 2, 1]
        assert all(rank_perm(unrank_perm(n, 5)) == n for n in range(120))
        # A permutation code larger than 2^63 (numeric column)
        assert rank_perm(unrank_perm(10**20, 25)) == 10**20
        assert cycle_string([2, 1, 4, 5, 3]) == "(1,2)(3,4,5)"
        assert cycle_string([1, 2, 3]) == "()"
        assert perm_order([2, 1, 4, 5, 3]) == 6
        # Products are read from left to right: apply (1,2) first, then (2,3)
        assert compose_lr([[2, 1, 3], [1, 3, 2]]) == [3, 1, 2]
        # The representative vector of 4.5-1.0.5-5-5-5.2: (1,5,4,3,2)^2 (1,2,3,4,5)^2 = 1
        c = rank_perm([5, 1, 2, 3, 4])
        cinv = rank_perm([2, 3, 4, 5, 1])
        assert cycle_string(unrank_perm(c, 5)) == "(1,5,4,3,2)"
        assert is_generating_vector_product_one([c, c, cinv, cinv], 5)
        assert not is_generating_vector_product_one([c, c, c, cinv], 5)

    def test_display_and_parsing(self):
        assert signature_display([0, 2, 3, 7]) == "[0;2,3,7]"
        assert signature_display([2]) == "[2;-]"
        assert jacobian_tex([[1, 2], [3, 1]]) == r"E^{2} \times A_{3}"
        assert parse_signature("[0;7,3,2]") == [0, 2, 3, 7]
        assert parse_signature("[0,2,3,7]") == [0, 2, 3, 7]
        assert parse_periods("4,2,4") == [2, 4]
        assert parse_group_labels("4.2,[8,3]") == ["4.2", "8.3"]
        summary = class_tuple_summary([{"braid": "b1", "con": ["2A", "2A"]},
                                       {"braid": "b1", "con": ["2A", "2A"]},
                                       {"braid": "b2", "con": ["2B", "2A"]}])
        assert summary == [{"con": ["2A", "2A"], "num_vectors": 2, "num_braid": 1},
                           {"con": ["2B", "2A"], "num_vectors": 1, "num_braid": 1}]

    # Pages

    def test_index_page(self):
        L = self.tc.get(BASE).get_data(as_text=True)
        assert 'Cyclic trigonal' in L
        assert 'topological classes' in L

    def test_family_page(self):
        L = self.tc.get(BASE + '4.5-1.0.5-5-5-5').get_data(as_text=True)
        assert '[0;5,5,5,5]' in L
        for T in [1, 2, 3]:
            assert '4.5-1.0.5-5-5-5.%s' % T in L
        assert 'Underlying data' in L

    def test_topological_page(self):
        L = self.tc.get(BASE + '4.5-1.0.5-5-5-5.2').get_data(as_text=True)
        assert '(1,5,4,3,2)' in L            # an entry in cycle notation
        assert '5A1' in L                    # its conjugacy class label
        assert '20.4' in L                   # the full automorphism group
        assert '4.20-4.0.2-2-2-5' in L       # the family of the full action
        assert 'not the LMFDB permutation representation' in L

    def test_unlabelled_full_group(self):
        L = self.tc.get(BASE + '46.1080-260.0.3-3-4.1').get_data(as_text=True)
        assert '2160' in L and 'not in the LMFDB' in L

    def test_many_braid_classes(self):
        L = self.tc.get(BASE + '21.32-47.0.2-2-4-4-4.1').get_data(as_text=True)
        assert 'show all' in L
        assert 'id="21.32-47.0.2-2-4-4-4.1.20"' in L
        assert 'id="21.32-47.0.2-2-4-4-4.1.21"' not in L
        L = self.tc.get(BASE + '21.32-47.0.2-2-4-4-4.1?braids=all').get_data(as_text=True)
        assert 'id="21.32-47.0.2-2-4-4-4.1.192"' in L

    def test_braid_and_vector_labels(self):
        L = self.tc.get(BASE + '21.32-47.0.2-2-4-4-4.1.17')
        assert L.status_code == 302 and L.location.endswith('21.32-47.0.2-2-4-4-4.1#21.32-47.0.2-2-4-4-4.1.17')
        L = self.tc.get(BASE + '10.6-1.0.2-2-2-2-2-2-2-2-2-2.1.1.3280')
        assert L.status_code == 302 and 'start=3250' in L.location

    def test_vector_listing(self):
        L = self.tc.get(BASE + '10.6-1.0.2-2-2-2-2-2-2-2-2-2.1.1/vectors').get_data(as_text=True)
        assert '1-50 of 3280' in L
        L = self.tc.get(BASE + '10.6-1.0.2-2-2-2-2-2-2-2-2-2.1.1/vectors?start=3250').get_data(as_text=True)
        assert 'id="10.6-1.0.2-2-2-2-2-2-2-2-2-2.1.1.3280"' in L

    def test_bad_labels(self):
        L = self.tc.get(BASE + 'junk', follow_redirects=True)
        assert 'No family with label' in L.get_data(as_text=True)
        L = self.tc.get(BASE + '4.5-1.0.5-5-5-5.99', follow_redirects=True)
        assert 'No topological class with label' in L.get_data(as_text=True)

    def test_search_topological(self):
        L = self.tc.get(BASE + '?genus=4&group=5.1').get_data(as_text=True)
        assert '4.5-1.0.5-5-5-5.2' in L
        # the group format of the previous version of this section
        L = self.tc.get(BASE + '?group=%5B5%2C1%5D').get_data(as_text=True)
        assert '4.5-1.0.5-5-5-5.2' in L

    def test_search_families(self):
        L = self.tc.get(BASE + '?genus=4&search_type=Families').get_data(as_text=True)
        assert 'family search results' in L and '4.5-1.0.5-5-5-5' in L
        L = self.tc.get(BASE + '?hyperelliptic=yes&search_type=Families', follow_redirects=True).get_data(as_text=True)
        assert 'can only be searched for topological classes' in L

    def test_search_aut(self):
        # 4.5-1.0.5-5-5-5.2 is not full; its full automorphism group 20.4 has order 20
        L = self.tc.get(BASE + '?aut_order=20').get_data(as_text=True)
        assert '4.5-1.0.5-5-5-5.2' in L
        L = self.tc.get(BASE + '?aut_group=20.4').get_data(as_text=True)
        assert '4.5-1.0.5-5-5-5.2' in L
        L = self.tc.get(BASE + '?genus=2&hyperelliptic=yes').get_data(as_text=True)
        assert '2.10-2.0.2-5-10.1' in L

    def test_search_jacobian_and_signature(self):
        L = self.tc.get(BASE + '?jac_dim=1&genus=2').get_data(as_text=True)
        assert '2.16-8.0.2-4-8.1' in L
        L = self.tc.get(BASE + '?signature=[0;5,5,5,5]').get_data(as_text=True)
        assert '4.5-1.0.5-5-5-5.3' in L
        L = self.tc.get(BASE + '?periods=5&r=4').get_data(as_text=True)
        assert '4.5-1.0.5-5-5-5.3' in L

    def test_downloads(self):
        L = self.tc.get(BASE + '4.5-1.0.5-5-5-5.2/download/magma').get_data(as_text=True)
        assert 'vector_data' in L and 'Sym(d)' in L
        L = self.tc.get(BASE + '4.5-1.0.5-5-5-5/download/gap').get_data(as_text=True)
        assert 'PermList' in L
        L = self.tc.get(BASE + '4.5-1.0.5-5-5-5.2/download/sage').get_data(as_text=True)
        assert 'SymmetricGroup(d)' in L
        L = self.tc.get(BASE + '?genus=4&download=1&query=%7B%27genus%27%3A+4%7D&Submit=magma')
        assert L.status_code == 200

    def test_data_page(self):
        L = self.tc.get(BASE + 'data/4.5-1.0.5-5-5-5.2').get_data(as_text=True)
        assert 'hgcwa_topological' in L and 'hgcwa_braid' in L

    def test_stats_page(self):
        L = self.tc.get(BASE + 'stats').get_data(as_text=True)
        assert 'Counts by' in L

    def test_random(self):
        L = self.tc.get(BASE + 'random')
        assert L.status_code == 307
