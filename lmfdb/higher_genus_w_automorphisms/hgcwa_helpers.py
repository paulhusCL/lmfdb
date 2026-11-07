"""
Pure Python helpers for families of higher genus curves with automorphisms.

Nothing in this file imports Sage, Flask or the LMFDB, so the functions can be
tested on their own (see test_hgcwa.py).  They handle:

* labels: a family ``g.n-i.g0.m1-...-mr``, a topological class ``family.T``,
  a braid class ``family.T.B`` and a generating vector ``family.T.B.V``;
* permutations: the LMFDB permutation codes stored in ``hgcwa_vectors.gen_vectors``
  (the rank of a permutation in the lexicographic order of Sym(d), counting from 0,
  as in Sage's ``Permutations(d).unrank`` and ``decode_as_perm`` in
  ``lmfdb/groups/abstract/web_groups.py``), cycle notation and products;
* display strings for signatures and Jacobian decompositions;
* parsing of the search inputs for signatures and group labels.

Permutations are lists of images ``[p(1), ..., p(d)]``.  Products are read from
left to right, as in Magma and GAP: in ``c_1 c_2 ... c_r`` the permutation
``c_1`` is applied first.  The generating vectors in the database satisfy
``c_1 c_2 ... c_r = 1`` in this convention.
"""

import re
from math import factorial, gcd

# Labels ####################################################################

# The group part of a family label is the LMFDB group label with "-" for ".";
# its second part can be letters for groups not identified by counter (e.g. 2304-ld).
_FAMILY = r"\d+\.\d+-[0-9a-z]+\.\d+\.\d+(?:-\d+)*"
FAMILY_RE = re.compile(r"^(%s)$" % _FAMILY)
TOPOLOGICAL_RE = re.compile(r"^(%s)\.(\d+)$" % _FAMILY)
BRAID_RE = re.compile(r"^(%s)\.(\d+)\.(\d+)$" % _FAMILY)
VECTOR_RE = re.compile(r"^(%s)\.(\d+)\.(\d+)\.(\d+)$" % _FAMILY)

GROUP_LABEL_RE = re.compile(r"^\d+\.(\d+|[a-z]+)$")
OLD_GROUP_LABEL_RE = re.compile(r"^\[(\d+),(\d+)\]$")


def label_kind(label):
    """
    The kind of object a label names: 'family', 'topological', 'braid', 'vector',
    or None if it is not a valid label.

    >>> label_kind("4.5-1.0.5-5-5-5"), label_kind("4.5-1.0.5-5-5-5.2")
    ('family', 'topological')
    >>> label_kind("4.5-1.0.5-5-5-5.2.1"), label_kind("4.5-1.0.5-5-5-5.2.1.1"), label_kind("junk")
    ('braid', 'vector', None)
    """
    if not isinstance(label, str):
        return None
    for kind, regex in [("family", FAMILY_RE), ("topological", TOPOLOGICAL_RE),
                        ("braid", BRAID_RE), ("vector", VECTOR_RE)]:
        if regex.match(label):
            return kind
    return None


def split_label(label):
    """
    Split a label into the family label and the list of integers following it.

    >>> split_label("10.6-1.0.2-2-2-2-2-2-2-2-2-2.1.1.3280")
    ('10.6-1.0.2-2-2-2-2-2-2-2-2-2', [1, 1, 3280])
    """
    for regex in [FAMILY_RE, TOPOLOGICAL_RE, BRAID_RE, VECTOR_RE]:
        m = regex.match(label)
        if m:
            return m.group(1), [int(x) for x in m.groups()[1:]]
    raise ValueError("Invalid label %s" % label)


def family_label_parts(family):
    """
    The genus, the LMFDB group label and the signature encoded in a family label.

    >>> family_label_parts("4.5-1.0.5-5-5-5")
    (4, '5.1', [0, 5, 5, 5, 5])
    """
    g, grp, g0, periods = family.split(".")
    sig = [int(g0)] + [int(m) for m in periods.split("-") if m != "0"]
    return int(g), grp.replace("-", "."), sig


def parent_label(label):
    """
    The label of the object one level up (vector -> braid class -> topological
    class -> family), or None for a family.
    """
    family, nums = split_label(label)
    if not nums:
        return None
    return ".".join([family] + [str(n) for n in nums[:-1]])


def relative_label(label):
    """
    The part of a label after the family label, as displayed in tables of a
    single family ("T", "T.B" or "T.B.V").

    >>> relative_label("21.32-47.0.2-2-4-4-4.1.17")
    '1.17'
    """
    family, nums = split_label(label)
    return ".".join(str(n) for n in nums) if nums else family


def label_sort_key(label):
    """
    Sort key for labels of one kind: numeric in T, B and V (the labels sort as
    text in the database, where "x.10" comes before "x.2").
    """
    family, nums = split_label(label)
    return (family, nums)


def vector_labels(braid_rows):
    """
    The labels of all generating vectors in the given braid classes, in order.

    INPUT: a list of dictionaries with keys ``label`` and ``num_vectors``, sorted
    by braid class.  Vectors are numbered 1, ..., num_vectors within their braid
    class, so the labels are determined without a database query.
    """
    return ["%s.%s" % (rec["label"], v)
            for rec in braid_rows for v in range(1, rec["num_vectors"] + 1)]


# Permutations ##############################################################

def unrank_perm(code, d):
    """
    The permutation of {1, ..., d} with the given LMFDB code: its rank, counting
    from 0, in the lexicographic order of Sym(d).  Agrees with Sage's
    ``Permutations(d).unrank(code)``.  Returns the list of images.

    >>> unrank_perm(0, 3), unrank_perm(5, 3), unrank_perm(1, 3)
    ([1, 2, 3], [3, 2, 1], [1, 3, 2])
    """
    code = int(code)  # numeric columns come back as Decimal
    if d < 0 or code < 0 or code >= factorial(d):
        raise ValueError("%s is not the code of a permutation of degree %s" % (code, d))
    elems = list(range(1, d + 1))
    out = []
    for i in range(d - 1, -1, -1):
        q, code = divmod(code, factorial(i))
        out.append(elems.pop(q))
    return out


def rank_perm(images):
    """
    The LMFDB code of a permutation given by its images (inverse of ``unrank_perm``).

    >>> rank_perm([3, 2, 1])
    5
    """
    d = len(images)
    elems = list(range(1, d + 1))
    code = 0
    for i, x in enumerate(images):
        j = elems.index(x)
        code += j * factorial(d - 1 - i)
        elems.pop(j)
    return code


def cycles(images):
    """
    The cycles of length at least 2, each starting with its smallest element,
    ordered by their smallest elements.

    >>> cycles([2, 1, 4, 5, 3])
    [[1, 2], [3, 4, 5]]
    """
    seen = set()
    out = []
    for i in range(1, len(images) + 1):
        if i in seen:
            continue
        cyc = []
        j = i
        while j not in seen:
            seen.add(j)
            cyc.append(j)
            j = images[j - 1]
        if len(cyc) > 1:
            out.append(cyc)
    return out


def cycle_string(images):
    """
    Cycle notation, as Sage prints elements of permutation groups.

    >>> cycle_string([1, 5, 2, 3, 4]), cycle_string([1, 2, 3])
    ('(2,5,4,3)', '()')
    """
    cs = cycles(images)
    if not cs:
        return "()"
    return "".join("(" + ",".join(str(x) for x in c) + ")" for c in cs)


def perm_order(images):
    """
    The order of a permutation.

    >>> perm_order([2, 1, 4, 5, 3])
    6
    """
    n = 1
    for c in cycles(images):
        n = n * len(c) // gcd(n, len(c))
    return n


def compose_lr(perms):
    """
    The product c_1 c_2 ... c_r read from left to right (apply c_1 first),
    as the list of images.

    >>> compose_lr([[2, 1, 3], [1, 3, 2]])
    [3, 1, 2]
    """
    d = len(perms[0])
    x = list(range(1, d + 1))
    for p in perms:
        x = [p[v - 1] for v in x]
    return x


def decode_vector(codes, d):
    """
    The entries of a generating vector, stored as LMFDB codes, as lists of images.
    """
    return [unrank_perm(c, d) for c in codes]


def is_generating_vector_product_one(codes, d):
    """
    Whether c_1 c_2 ... c_r = 1 (products from left to right) for the
    generating vector with the given LMFDB codes.
    """
    return compose_lr(decode_vector(codes, d)) == list(range(1, d + 1))


# Display ###################################################################

def signature_display(sig):
    """
    The signature [g0; m_1, ..., m_r] for display in math mode, from the stored
    list [g0, m_1, ..., m_r].  Returns "" for a missing value.

    >>> signature_display([0, 5, 5, 5, 5]), signature_display([2])
    ('[0;5,5,5,5]', '[2;-]')
    """
    if not sig:
        return ""
    if len(sig) == 1:
        return "[%s;-]" % sig[0]
    return "[%s;%s]" % (sig[0], ",".join(str(m) for m in sig[1:]))


def jacobian_factor_tex(dim, mult):
    """
    >>> jacobian_factor_tex(1, 2), jacobian_factor_tex(3, 1)
    ('E^{2}', 'A_{3}')
    """
    s = "E" if dim == 1 else "A_{%s}" % dim
    if mult != 1:
        s += "^{%s}" % mult
    return s


def jacobian_tex(dims):
    """
    The isogeny decomposition of the Jacobian from the pairs [dimension, multiplicity],
    for math mode.  Returns "" for a missing value.

    >>> jacobian_tex([[1, 2], [3, 1]])
    'E^{2} \\\\times A_{3}'
    """
    if not dims:
        return ""
    return r" \times ".join(jacobian_factor_tex(d, m) for d, m in dims)


def class_tuple_summary(rows):
    """
    The distinct tuples of conjugacy class labels among the given vectors, with the
    number of vectors and of braid classes having each, in order of first appearance.

    INPUT: an iterable of dictionaries with keys ``braid`` and ``con``.
    """
    order = []
    nvec = {}
    braids = {}
    for rec in rows:
        key = tuple(rec["con"])
        if key not in nvec:
            order.append(key)
            nvec[key] = 0
            braids[key] = set()
        nvec[key] += 1
        braids[key].add(rec["braid"])
    return [{"con": list(key), "num_vectors": nvec[key], "num_braid": len(braids[key])}
            for key in order]


# Parsing ###################################################################

def parse_signature(inp):
    """
    Parse a signature typed by a user: [0;2,3,7], [0,2,3,7], 0;2,3,7 or (0,2,3,7).
    The periods are sorted, as in the database.

    >>> parse_signature("[0;7,2,3]"), parse_signature("0,2,3,7")
    ([0, 2, 3, 7], [0, 2, 3, 7])
    """
    s = inp.replace(" ", "").replace(";", ",")
    if s[:1] in "[(" and s[-1:] in "])":
        s = s[1:-1]
    if not re.fullmatch(r"\d+(,\d+)*", s):
        raise ValueError("It needs to be a signature [g0;m1,...,mr] such as [0;2,3,7] or [0,2,3,7]")
    L = [int(x) for x in s.split(",")]
    if any(m < 2 for m in L[1:]):
        raise ValueError("The periods m1, ..., mr must be at least 2")
    return [L[0]] + sorted(L[1:])


def parse_periods(inp):
    """
    Parse a comma separated list of periods.

    >>> parse_periods("4,2")
    [2, 4]
    """
    s = inp.replace(" ", "").strip("[]()")
    if not re.fullmatch(r"\d+(,\d+)*", s):
        raise ValueError("It needs to be a period such as 5, or a list of periods such as 2,4")
    L = sorted(set(int(x) for x in s.split(",")))
    if any(m < 2 for m in L):
        raise ValueError("Periods must be at least 2")
    return L


def parse_group_labels(inp):
    """
    Parse a list of LMFDB group labels.  The old form [n,i] of the previous
    version of this section (still used in links from other pages) is accepted too.

    >>> parse_group_labels("4.2,[8,3]")
    ['4.2', '8.3']
    """
    s = inp.replace(" ", "")
    pieces = re.findall(r"\[\d+,\d+\]|[^,]+", s)
    if not pieces:
        raise ValueError("It needs to be an LMFDB group label such as 4.2")
    out = []
    for piece in pieces:
        m = OLD_GROUP_LABEL_RE.match(piece)
        if m:
            piece = "%s.%s" % m.groups()
        if not GROUP_LABEL_RE.match(piece):
            raise ValueError("%s is not an LMFDB group label, such as 4.2 or 32.47" % piece)
        out.append(piece)
    return out
