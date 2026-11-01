# This Blueprint is about higher genus curves with automorphisms
# Authors: Jen Paulhus, Lex Martin, David Neill Asanza, Nhi Ngo, Albert Ford
# (initial code copied from John Jones Local Fields)
#
# The data lives in four tables:
#   hgcwa_families     one row per (genus, group, signature)
#   hgcwa_topological  one row per topological class of actions (the main objects)
#   hgcwa_braid        one row per braid class
#   hgcwa_vectors      one row per generating vector, up to simultaneous conjugation
# Labels are family = g.n-i.g0.m1-...-mr, then family.T, family.T.B, family.T.B.V.

import os
import re

import yaml
from flask import render_template, request, url_for, redirect, abort
from psycodict.utils import SearchParsingError

from lmfdb import db
from lmfdb.utils import (
    flash_error, to_dict, comma, Downloader,
    SearchArray, EmbeddedSearchArray, TextBox, YesNoBox, CountBox, SneakyTextBox,
    parse_ints, parse_bool, parse_element_of, clean_input,
    search_wrap, embed_wrap, redirect_no_cache)
from lmfdb.utils.interesting import interesting_knowls
from lmfdb.utils.search_parsing import search_parser, collapse_ors, parse_range2, LIST_RE
from lmfdb.utils.search_columns import (
    SearchColumns, LinkCol, MathCol, ProcessedCol, MultiProcessedCol, CheckCol)
from lmfdb.api import datapage
from lmfdb.groups.abstract.main import abstract_group_display_knowl
from lmfdb.higher_genus_w_automorphisms import higher_genus_w_automorphisms_page
from lmfdb.higher_genus_w_automorphisms.hgcwa_stats import HGCWAstats
from lmfdb.higher_genus_w_automorphisms.hgcwa_helpers import (
    label_kind, split_label, parent_label, relative_label, label_sort_key, vector_labels,
    family_label_parts,
    unrank_perm, cycle_string, signature_display, jacobian_tex, jacobian_factor_tex,
    class_tuple_summary, parse_signature, parse_periods, parse_group_labels)

# Number of braid classes shown on a topological class page before "show all"
BRAID_ROWS_SHOWN = 20
# Never show more braid classes than this on one page (use the vector listing instead)
BRAID_ROWS_MAX = 2000
# Above this many vectors, the class tuple summary is not computed on the page
CLASS_TUPLE_MAX_VECTORS = 10000
# Default number of vectors per page in the vector listing
VECTORS_PER_PAGE = 50

# Groups whose genus 2 curves can be found in the genus 2 section
# (the automorphism groups of genus 2 curves)
G2_AUT_GROUPS = ["2.1", "4.2", "8.3", "10.2", "12.4", "24.8", "48.29"]


def get_bread(tail=[]):
    base = [("Higher genus", url_for(".index")),
            ("C", url_for(".index")),
            ("Aut", url_for(".index"))]
    if not isinstance(tail, list):
        tail = [(tail, " ")]
    return base + tail


def learnmore_list():
    return [('Source and acknowledgments', url_for(".how_computed_page")),
            ('Completeness of the data', url_for(".completeness_page")),
            ('Reliability of the data', url_for(".reliability_page")),
            ('Labeling convention', url_for(".labels_page"))]


def url_for_label(label):
    if label == "random":
        return url_for(".random_topological")
    return url_for(".by_label", label=label)


def family_bread(family, group_tex=None, tail=[]):
    """
    Bread crumbs genus -> group -> signature (the family) -> tail.
    """
    genus, group, sig = family_label_parts(family)
    name = "$%s$" % group_tex if group_tex else group
    return get_bread([(str(genus), url_for(".index", genus=genus)),
                      (name, url_for(".index", genus=genus, group=group)),
                      ("$%s$" % signature_display(sig), url_for(".by_label", label=family))] + tail)


# Groups ####################################################################
# Names, links and permutation representations come from the LMFDB groups tables.

def group_cache(labels, extra=[]):
    """
    One query for the records of the given groups in gps_groups, as a dictionary
    indexed by label.
    """
    labels = sorted({lab for lab in labels if lab})
    if not labels:
        return {}
    proj = ["label", "tex_name"] + list(extra)
    return {rec["label"]: rec for rec in db.gps_groups.search({"label": {"$in": labels}}, proj)}


def group_tex(label, cache):
    rec = (cache or {}).get(label)
    if rec is None:
        return None
    return rec.get("tex_name")


def group_knowl(label, cache):
    """
    A group as a knowl showing its name; the label is shown if the group has no name.
    """
    if not label:
        return ""
    tex = group_tex(label, cache)
    return abstract_group_display_knowl(label, name=("$%s$" % tex) if tex else label)


def group_link(label, cache):
    """
    The name of a group, linked to its page, followed by its label.
    """
    tex = group_tex(label, cache)
    url = url_for("abstract.by_label", label=label)
    if tex:
        return '<a href="%s">$%s$</a> (%s)' % (url, tex, label)
    return '<a href="%s">%s</a>' % (url, label)


def perm_rep(rec):
    """
    The degree and the generators (as lists of images) of the LMFDB permutation
    representation of a group, from its gps_groups record; (None, []) if missing.
    """
    try:
        rep = rec["representations"]["Perm"]
        d = rep["d"]
        return d, [unrank_perm(g, d) for g in rep["gens"]]
    except (KeyError, TypeError, ValueError):
        return None, []


def cc_knowl(group, label):
    """
    A conjugacy class label of a group as a dynamic knowl (as on abstract group pages).
    """
    if not label:
        return ""
    return ('<a title="%s [lmfdb.object_information]" knowl="lmfdb.object_information" '
            'kwargs="func=cc_data&args=%s%%7C%s%%7Ccomplex">%s</a>' % (label, group, label, label))


def class_tuple_display(group, con):
    if not con:
        return ""
    return ", ".join(cc_knowl(group, cl) for cl in con)


def qchar_link(group, char):
    """
    A rational character label, linked to the rational character table of the group.
    """
    url = url_for("abstract.Qchar_table", label=group, char_highlight=char)
    return '<a href="%s">%s</a>' % (url, char)


def perm_display(perm):
    """
    A permutation (list of images, or the string from ``entry_data``) in cycle notation, in math mode.
    """
    if perm is None or perm == "":
        return ""
    if isinstance(perm, list):
        perm = cycle_string(perm)
    return "$%s$" % perm


def entry_data(code, d):
    """
    The cycle notation of a stored permutation code, or the code itself if the
    degree of the permutation representation is unknown.
    """
    if d is None:
        return "code %s" % code
    return cycle_string(unrank_perm(code, d))


# Search columns ############################################################
# Derived values are computed by search columns, as in the rest of the LMFDB.

def group_col(group, tex_cache):
    return group_knowl(group, tex_cache)


def aut_col(is_full, group, full_auto, full_auto_order, tex_cache):
    """
    The full automorphism group of the generic curve: G itself when the action is full.
    """
    if is_full:
        return group_knowl(group, tex_cache)
    if full_auto:
        return group_knowl(full_auto, tex_cache)
    if full_auto_order:
        return "order $%s$ (not in the LMFDB)" % full_auto_order
    return ""


def aut_download(is_full, group, full_auto, full_auto_order, tex_cache):
    if is_full:
        return group
    return full_auto


def sig_col(sig):
    return signature_display(sig) if sig else ""


def jac_col(dims):
    return jacobian_tex(dims) if dims else ""


def tex_cache_postprocess(res, info, query):
    """
    Adds the names of all groups in the results with one gps_groups query.
    """
    cache = group_cache([rec.get("group") for rec in res] + [rec.get("full_auto") for rec in res])
    for rec in res:
        rec["tex_cache"] = cache
    return res


def family_top_postprocess(res, info, query):
    """
    As tex_cache_postprocess, and adds the class tuple of the representative vector
    of each topological class (one hgcwa_vectors query).
    """
    res = tex_cache_postprocess(res, info, query)
    reps = [rec["rep_vector"] for rec in res if rec.get("rep_vector")]
    cons = {}
    if reps:
        cons = {rec["label"]: rec["con"] for rec in db.hgcwa_vectors.search({"label": {"$in": reps}}, ["label", "con"])}
    for rec in res:
        rec["rep_con"] = cons.get(rec.get("rep_vector"))
    return res


def _db_cols(columns, extra=[]):
    return sorted(set(sum([C.orig for C in columns], [])) | set(extra))


label_col = LinkCol("label", "dq.curve.highergenus.aut.label", "Label", url_for_label)
genus_col = MathCol("genus", "ag.curve.genus", "Genus")
group_name_col = MultiProcessedCol("group", "group.name", "Group", ["group", "tex_cache"], group_col,
                                   apply_download=lambda group, tex_cache: group)
group_order_col = MathCol("group_order", "group.order", "Group order", default=False)
signature_col = ProcessedCol("signature", "curve.highergenus.aut.signature", "Signature", sig_col,
                             mathmode=True, align="center")
r_col = MathCol("r", "curve.highergenus.aut.branchpoints", "$r$", short_title="branch points", default=False)
dim_col = MathCol("dim", "curve.highergenus.aut.dimension", "Dimension")
is_full_col = CheckCol("is_full", "curve.highergenus.aut.full", "Full", short_title="full")
aut_group_col = MultiProcessedCol("aut", "curve.highergenus.aut.full", "Aut",
                                  ["is_full", "group", "full_auto", "full_auto_order", "tex_cache"],
                                  aut_col, apply_download=aut_download, short_title="full automorphism group")
signH_col = ProcessedCol("signH", "curve.highergenus.aut.signature", "Aut signature", sig_col,
                         mathmode=True, align="center", default=False)
full_index_col = MathCol("full_index", "curve.highergenus.aut.full_index", r"$[\mathrm{Aut}:G]$",
                         short_title="index [Aut:G]", default=False)
hyp_col = CheckCol("hyperelliptic", "ag.hyperelliptic_curve", "Hyperelliptic", unknown="")
trig_col = CheckCol("cyclic_trigonal", "ag.cyclic_trigonal", "Cyclic trigonal", unknown="")
jac_dims_col = ProcessedCol("jacobian_dims", "curve.highergenus.aut.groupalgebradecomp", "Jacobian",
                            jac_col, mathmode=True, align="center", short_title="Jacobian decomposition")
num_braid_col = MathCol("num_braid", "curve.highergenus.aut.braid_equivalence", "Braid classes")
num_vectors_col = MathCol("num_vectors", "curve.highergenus.aut.generatingvector", "Vectors",
                          short_title="generating vectors")

# Search results: topological classes
hgcwa_columns = SearchColumns([
    label_col, genus_col, group_name_col, group_order_col, signature_col, r_col, dim_col,
    is_full_col, aut_group_col, signH_col, full_index_col, hyp_col, trig_col,
    jac_dims_col, num_braid_col, num_vectors_col,
    LinkCol("family", "curve.highergenus.aut.family", "Family", url_for_label, default=False)])
hgcwa_columns.languages = ['text', 'gap', 'magma', 'sage', 'csv']

# Search results: families
family_columns = SearchColumns([
    label_col, genus_col, group_name_col, group_order_col, signature_col, r_col, dim_col,
    MathCol("num_topological", "curve.highergenus.aut.topological_equivalence", "Topological classes"),
    MathCol("num_topological_full", "curve.highergenus.aut.full", "Full", short_title="full topological classes"),
    num_braid_col, num_vectors_col])
family_columns.languages = ['text', 'gap', 'magma', 'sage', 'csv']

# The topological classes of one family, on the family page
_family_top_cols = [
    label_col, is_full_col, aut_group_col,
    ProcessedCol("signH", "curve.highergenus.aut.signature", "Aut signature", sig_col, mathmode=True, align="center"),
    MathCol("full_index", "curve.highergenus.aut.full_index", r"$[\mathrm{Aut}:G]$", short_title="index [Aut:G]"),
    hyp_col, trig_col, jac_dims_col, num_braid_col, num_vectors_col,
    MultiProcessedCol("rep_con", "curve.highergenus.aut.conjugacyclasses", "Class tuple", ["group", "rep_con"],
                      class_tuple_display, apply_download=lambda group, con: con,
                      short_title="class tuple of the representative vector")]
family_top_columns = SearchColumns(_family_top_cols, db_cols=_db_cols(_family_top_cols, ["rep_vector"]))
family_top_columns.languages = ['text', 'gap', 'magma', 'sage', 'csv']


def entry_col(i, period):
    """
    The column for the entry c_i of the generating vectors in a table of vectors:
    its conjugacy class label above the permutation.
    """
    def disp(group, entry):
        if not entry:
            return ""
        cl, perm = entry
        return "%s<br>%s" % (cc_knowl(group, cl), perm_display(perm))
    return MultiProcessedCol("c%s" % i, "curve.highergenus.aut.generatingvector",
                             "$c_{%s}$ (order $%s$)" % (i, period), ["group", "c%s" % i], disp,
                             apply_download=lambda group, entry: entry, short_title="c%s" % i)


def braid_label_display(label):
    return '<span id="%s">%s</span>' % (label, relative_label(label))


def vector_count_display(label, num_vectors):
    if num_vectors > 1:
        return '<a href="%s">%s</a>' % (url_for(".vector_listing", label=label), num_vectors)
    return str(num_vectors)


def braid_table_columns(signature, show_chars):
    """
    The braid classes of a topological class, each with its representative vector.
    """
    cols = [ProcessedCol("label", "curve.highergenus.aut.braid_equivalence", "Braid class", braid_label_display),
            MultiProcessedCol("num_vectors", "curve.highergenus.aut.generatingvector", "Vectors",
                              ["label", "num_vectors"], vector_count_display, align="center")]
    cols += [entry_col(i + 1, m) for i, m in enumerate(signature[1:])]
    if show_chars:
        cols.append(MultiProcessedCol("chars", "curve.highergenus.aut.characters", "Characters",
                                      ["group", "chars"],
                                      lambda group, chars: ", ".join(qchar_link(group, c) for c in chars)))
    return SearchColumns(cols)


def vector_table_columns(signature, show_braid):
    """
    A list of generating vectors (vector listing page).
    """
    cols = [ProcessedCol("label", "curve.highergenus.aut.generatingvector", "Vector", braid_label_display)]
    if show_braid:
        cols.append(ProcessedCol("braid", "curve.highergenus.aut.braid_equivalence", "Braid class",
                                 lambda lab: '<a href="%s">%s</a>' % (url_for(".by_label", label=lab), relative_label(lab))))
    cols += [entry_col(i + 1, m) for i, m in enumerate(signature[1:])]
    return SearchColumns(cols)


# The representative generating vector, one entry per row
rep_vector_columns = SearchColumns([
    MathCol("entry", "curve.highergenus.aut.generatingvector", "Entry"),
    ProcessedCol("perm", "group.permutation_representation", "Permutation", perm_display),
    MultiProcessedCol("cl", "group.conjugacy_class", "Conjugacy class", ["group", "cl"], cc_knowl),
    MathCol("order", "group.element_order", "Order")])

# Distinct tuples of conjugacy classes in a topological class
class_tuple_columns = SearchColumns([
    MultiProcessedCol("con", "curve.highergenus.aut.conjugacyclasses", "Class tuple", ["group", "con"], class_tuple_display),
    MathCol("num_vectors", "curve.highergenus.aut.generatingvector", "Vectors"),
    MathCol("num_braid", "curve.highergenus.aut.braid_equivalence", "Braid classes")])

# The group algebra decomposition of the Jacobian
jacobian_columns = SearchColumns([
    MultiProcessedCol("factor", "curve.highergenus.aut.groupalgebradecomp", "Factor", ["dim", "mult"],
                      lambda d, m: "$%s$" % jacobian_factor_tex(d, m)),
    MathCol("dim", "ag.dimension", "Dimension"),
    MathCol("mult", "curve.highergenus.aut.groupalgebradecomp", "Multiplicity"),
    MultiProcessedCol("char", "curve.highergenus.aut.characters", "Rational character", ["group", "char"],
                      qchar_link)])


# Downloads #################################################################

_curdir = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(_curdir, "code.yaml")) as _codefile:
    code_list = yaml.safe_load(_codefile)


class HGCWA_download(Downloader):
    table = db.hgcwa_topological
    title = "Topological classes of group actions on higher genus curves"
    short_name = "topological classes"
    var_name = "topological_classes"
    vector_languages = ["magma", "gap", "sage"]

    def postprocess(self, row, info, query):
        # The class tuple column of the family page needs the representative vector
        if row.get("rep_vector") and "rep_con" not in row:
            row["rep_con"] = db.hgcwa_vectors.lookup(row["rep_vector"], "con")
        return row

    def download_vectors(self, label, lang):
        """
        All generating vectors of a family, topological class or braid class,
        as permutations in the LMFDB permutation representation of the group.
        """
        if lang not in self.vector_languages:
            return abort(404, "Unknown download type %s" % lang)
        kind = label_kind(label)
        if kind not in ["family", "topological", "braid"]:
            return abort(404, "Invalid label %s" % label)
        family = split_label(label)[0]
        fam = db.hgcwa_families.lookup(family, ["genus", "group", "signature"])
        vecs = list(db.hgcwa_vectors.search({kind: label}, ["label", "con", "gen_vectors"]))
        if fam is None or not vecs:
            return abort(404, "No generating vectors for label %s" % label)
        vecs.sort(key=lambda rec: label_sort_key(rec["label"]))
        grp = db.gps_groups.lookup(fam["group"], ["representations"])
        d, gens = perm_rep(grp)
        if d is None:
            return abort(404, "The group %s has no permutation representation in the LMFDB" % fam["group"])
        language = self.languages[lang]
        s = language.comment(
            " Generating vectors (c_1, ..., c_r) of %s, up to simultaneous conjugation.\n"
            " The entries are permutations in the LMFDB permutation representation of the group %s,\n"
            " of degree %s.  Products are read from left to right (apply c_1 first),\n"
            " and c_1 c_2 ... c_r = 1.\n"
            " conjugacy_classes gives the LMFDB labels of the conjugacy classes of the entries.\n\n"
            % (label, fam["group"], d))
        s += language.assign("label", label)
        s += language.assign("group_label", fam["group"])
        s += language.assign("genus", fam["genus"])
        s += language.assign("signature", fam["signature"])
        s += language.assign("d", d)
        s += language.assign("gens", gens)
        s += language.assign("vector_labels", [rec["label"] for rec in vecs])
        s += language.assign("conjugacy_classes", [rec["con"] for rec in vecs])
        s += language.assign("vector_data", [[unrank_perm(c, d) for c in rec["gen_vectors"]] for rec in vecs])
        s += "\n" + language.comment(" The group and the vectors as permutations\n")
        for key in ["group", "vectors", "check"]:
            s += code_list[key][lang]
        return self._wrap(s, label + ".vectors", lang=lang,
                          title="Generating vectors of %s," % label)


class HGCWA_family_download(Downloader):
    table = db.hgcwa_families
    title = "Families of higher genus curves with automorphisms"
    short_name = "families"
    var_name = "families"


@higher_genus_w_automorphisms_page.route("/<label>/download/<download_type>")
def hgcwa_code_download(label, download_type):
    return HGCWA_download().download_vectors(label, download_type)


def vector_downloads(label):
    return [('Generating vectors to Magma', url_for(".hgcwa_code_download", label=label, download_type='magma')),
            ('Generating vectors to GAP', url_for(".hgcwa_code_download", label=label, download_type='gap')),
            ('Generating vectors to SageMath', url_for(".hgcwa_code_download", label=label, download_type='sage')),
            ('Underlying data', url_for(".hgcwa_data", label=label))]


# Browse, search, random, interesting, statistics ##########################

@higher_genus_w_automorphisms_page.route("/")
def index():
    info = to_dict(request.args, search_array=HGCWASearchArray(), stats=HGCWAstats())
    if request.args:
        search_type = info.get("search_type", info.get("hst", ""))
        if search_type in ["Families", "RandomFamily"]:
            return hgcwa_family_search(info)
        return higher_genus_w_automorphisms_search(info)
    small = db.hgcwa_families.distinct("group", {"group_order": {"$lte": 12}})
    cache = group_cache(small)
    small.sort(key=lambda lab: [int(x) if x.isdigit() else x for x in lab.split(".")])
    info["browse_groups"] = [(lab, group_tex(lab, cache)) for lab in small]
    info["browse_r"] = sorted(db.hgcwa_families.distinct("r"))
    info["browse_orders"] = ["2-9", "10-19", "20-49", "50-199", "200-999", "1000-"]
    return render_template("hgcwa-index.html",
                           title="Families of higher genus curves with automorphisms",
                           bread=get_bread(),
                           info=info,
                           learnmore=learnmore_list())


@higher_genus_w_automorphisms_page.route("/random")
@redirect_no_cache
def random_topological():
    label = db.hgcwa_topological.random()
    return url_for(".by_label", label=label)


@higher_genus_w_automorphisms_page.route("/interesting")
def interesting():
    return interesting_knowls(
        "curve.highergenus.aut",
        db.hgcwa_topological,
        url_for_label,
        label_col="label",
        regex=re.compile(r"^\d+\.\d+-[0-9a-z]+\.\d+\.\d+(-\d+)*(\.\d+)?$"),
        title=r"Some interesting higher genus families and topological classes",
        bread=get_bread("Interesting"),
        learnmore=learnmore_list()
    )


@higher_genus_w_automorphisms_page.route("/stats")
def statistics():
    title = 'Families of higher genus curves with automorphisms: statistics'
    bread = get_bread('Statistics')
    return render_template("hgcwa-stats.html", info=HGCWAstats(), title=title, learnmore=learnmore_list(), bread=bread)


def hgcwa_jump(info):
    label = clean_input(info['jump']).strip()
    if label_kind(label) is None:
        flash_error("%s is not a valid label for a family, topological class, braid class or generating vector.", label)
        return redirect(url_for(".index"))
    return redirect(url_for_label(label))


# Search parsing ############################################################

@search_parser(clean_info=True)
def parse_group(inp, query, qfield):
    labels = parse_group_labels(inp)
    query[qfield] = labels[0] if len(labels) == 1 else {"$in": labels}


@search_parser(clean_info=True)
def parse_signature_box(inp, query, qfield):
    query[qfield] = parse_signature(inp)


@search_parser(clean_info=True)
def parse_periods_box(inp, query, qfield):
    periods = parse_periods(inp)
    if isinstance(query.get(qfield), list):
        # an exact signature was also given
        if not set(periods).issubset(query[qfield][1:]):
            raise SearchParsingError("The signature does not contain all of the periods %s" % inp)
    else:
        query[qfield] = {"$contains": periods}


def _int_values(inp):
    """
    The values of an integer search input (an integer, a range or a list), as
    query values for a single column.
    """
    if not LIST_RE.match(inp):
        raise SearchParsingError("It needs to be an integer (such as 25), a range of integers (such as 2-10 or 2..10), or a comma-separated list of these (such as 4,9,16 or 4-25, 81-121).")
    parsed = parse_range2(inp, "x", int)
    if parsed[0] == "$or":
        return [D["x"] for D in parsed[1]]
    return [parsed[1]]


def _add_or(query, ors):
    # $and avoids mixing these alternatives with other $or conditions (from group_order, say)
    query.setdefault("$and", []).append({"$or": ors})


@search_parser(clean_info=True, prep_ranges=True)
def parse_aut_order(inp, query, qfield):
    """
    The order of the full automorphism group: the order of G for full actions,
    full_auto_order otherwise.
    """
    ors = []
    for v in _int_values(inp):
        ors.append({"is_full": True, "group_order": v})
        ors.append({"full_auto_order": v})
    _add_or(query, ors)


@search_parser(clean_info=True)
def parse_aut_group(inp, query, qfield):
    """
    The full automorphism group: G for full actions, full_auto otherwise.
    """
    labels = parse_group_labels(inp)
    v = labels[0] if len(labels) == 1 else {"$in": labels}
    _add_or(query, [{"is_full": True, "group": v}, {"full_auto": v}])


# Similar to parse_ints in lmfdb/utils, but also allows linear expressions in the genus g
LIST_G_RE = re.compile(r'^(\d+|(\d*-(\d+)?)|((\d*)\**(g((\+|\-)(\d*))*|\(g(\+|\-)(\d+)\))))(,(\d+|(\d*-(\d+)?)|((\d*)\**(g((\+|\-)(\d*))*|\(g(\+|\-)(\d+)\)))))*$')
GENUS_RE = re.compile(r'^(\d*)\**(g((\+|\-)(\d*))*|\(g(\+|\-)(\d+)\))$')
GROUP_ORDER_ERROR = ("It needs to be an integer (such as 25), a range of integers (such as 2-10 or 2..10), "
                     "a linear function of the variable g for genus (such as 84(g-1), 84g-84, 84g, or g-1), "
                     "or a comma-separated list of these (such as 4,9,16 or 4-25, 81-121).")


def parse_range2_extend(arg, key, parse_singleton=int, parse_endpoint=None, instance=1):
    if parse_endpoint is None:
        parse_endpoint = parse_singleton
    if isinstance(arg, str):
        arg = arg.replace(' ', '')
    if isinstance(arg, parse_singleton):
        return [key, arg]
    if ',' in arg:
        instance = len(arg.split(','))
        tmp = [parse_range2_extend(a, key, parse_singleton, parse_endpoint, instance) for a in arg.split(',')]
        ret = []
        for a in tmp:
            if a[0] == key:
                if len(a) == 3:
                    ret.append({a[0]: a[1], 'genus': a[2]})
                else:
                    ret.append({a[0]: a[1]})
            else:
                for ai in a:
                    ret.append({ai[0]: ai[1], 'genus': ai[2]})
        return ['$or', ret]
    elif 'g' in arg:  # linear function of variable g (ax+b)
        if GENUS_RE.match(arg):
            a = GENUS_RE.match(arg).groups()[0]
            genus_list = sorted(db.hgcwa_families.distinct('genus'))
            queries = []
            for g in range(genus_list[0], genus_list[-1] + 1):
                if '(' in arg:
                    b = int(GENUS_RE.match(arg).groups()[6])
                    if '+' in arg:  # a(g+b)
                        group_order = int(a) * (g + b)
                    else:  # a(g-b)
                        group_order = int(a) * (g - b)
                else:
                    if '+' in arg:
                        b = int(GENUS_RE.match(arg).groups()[4])
                        group_order = g + b if a == '' else int(a) * g + b
                    elif '-' in arg:
                        b = int(GENUS_RE.match(arg).groups()[4])
                        group_order = g - b if a == '' else int(a) * g - b
                    else:
                        group_order = g if a == '' else int(a) * g
                queries.append((group_order, g))
            if instance == 1:  # If there is only one linear function
                return ['$or', [{key: gp_ord, 'genus': g} for gp_ord, g in queries]]
            else:
                return [[key, gp_ord, g] for gp_ord, g in queries]  # Nested list
        else:
            raise ValueError(GROUP_ORDER_ERROR)
    elif '-' in arg:
        ix = arg.index('-', 1)
        start, end = arg[:ix], arg[ix + 1:]
        q = {}
        if start:
            q['$gte'] = parse_endpoint(start)
        if end:
            q['$lte'] = parse_endpoint(end)
        return [key, q]
    else:
        return [key, parse_singleton(arg)]


@search_parser(clean_info=True, prep_ranges=True)
def parse_group_order(inp, query, qfield, parse_singleton=int):
    if LIST_G_RE.match(inp):
        collapse_ors(parse_range2_extend(inp, qfield, parse_singleton), query)
    else:
        raise ValueError(GROUP_ORDER_ERROR)


def common_parse(info, query):
    """
    The inputs that are columns of both hgcwa_families and hgcwa_topological.
    """
    parse_ints(info, query, 'genus', name='Genus')
    parse_ints(info, query, 'dim', name='Dimension')
    parse_ints(info, query, 'r', name='Number of branch points')
    parse_group(info, query, 'group', name='Group')
    parse_group_order(info, query, 'group_order', name='Group order')
    parse_signature_box(info, query, 'signature', name='Signature')
    parse_periods_box(info, query, 'periods', qfield='signature', name='Periods')
    parse_ints(info, query, 'num_braid', name='Braid classes')
    parse_ints(info, query, 'num_vectors', name='Generating vectors')


# Inputs that only make sense for topological classes
TOPOLOGICAL_ONLY = [("aut_order", "Order of Aut"), ("aut_group", "Aut group"),
                    ("full_index", "[Aut:G]"), ("hyperelliptic", "Hyperelliptic"),
                    ("cyclic_trigonal", "Cyclic trigonal"), ("jac_dim", "Jacobian factor dimension")]

hgcwa_search_kwds = dict(
    per_page=50,
    url_for_label=url_for_label,
    postprocess=tex_cache_postprocess,
    bread=lambda: get_bread("Search results"),
    learnmore=learnmore_list)


@search_wrap(
    table=db.hgcwa_topological,
    title='Higher genus curve action search results',
    err_title='Higher genus curve action search input error',
    columns=hgcwa_columns,
    shortcuts={'jump': hgcwa_jump, 'download': HGCWA_download()},
    **hgcwa_search_kwds)
def higher_genus_w_automorphisms_search(info, query):
    common_parse(info, query)
    parse_bool(info, query, 'is_full', name='Full automorphism group')
    parse_ints(info, query, 'full_index', name='[Aut:G]')
    parse_aut_order(info, query, 'aut_order', name='Order of Aut')
    parse_aut_group(info, query, 'aut_group', name='Aut group')
    parse_bool(info, query, 'hyperelliptic', name='Hyperelliptic')
    parse_bool(info, query, 'cyclic_trigonal', name='Cyclic trigonal')
    parse_element_of(info, query, 'jac_dim', qfield='jacobian_factor_dims', split_interval=200,
                     name='Jacobian factor dimension')


@search_wrap(
    table=db.hgcwa_families,
    title='Higher genus curve family search results',
    err_title='Higher genus curve family search input error',
    columns=family_columns,
    shortcuts={'jump': hgcwa_jump, 'download': HGCWA_family_download()},
    **hgcwa_search_kwds)
def hgcwa_family_search(info, query):
    bad = [name for field, name in TOPOLOGICAL_ONLY if info.get(field)]
    if bad:
        msg = "%s can only be searched for topological classes; use the list of topological classes." % ", ".join(bad)
        flash_error(msg)
        raise SearchParsingError(msg)
    common_parse(info, query)
    # For a family: whether some (yes) or no (no) topological class has G as its full automorphism group
    parse_bool(info, query, 'is_full', qfield='num_topological_full', name='Full automorphism group',
               process=lambda b: {"$gt": 0} if b else 0)


class HGCWASearchArray(SearchArray):
    noun = "topological class"
    plural_noun = "topological classes"
    jump_example = "4.5-1.0.5-5-5-5.2"
    jump_egspan = "e.g. 2.10-2.0.2-5-10 (a family), 4.5-1.0.5-5-5-5.2 (a topological class) or 21.32-47.0.2-2-4-4-4.1.17 (a braid class)"
    jump_knowl = "curve.highergenus.aut.search_input"
    jump_prompt = "Label"
    label_knowl = "dq.curve.highergenus.aut.label"
    has_diagram = False
    null_column_explanations = {
        # null when G is not the full automorphism group: the data belongs to the full action
        'hyperelliptic': "hyperelliptic (recorded only when $G$ is the full automorphism group)",
        'cyclic_trigonal': "cyclic trigonal (recorded only when $G$ is the full automorphism group)",
        'full_auto': False,
        'full_auto_order': False,
    }
    sort_knowl = "curve.highergenus.aut.sort_order"
    _top_sorts = [
        ("", "genus", ['genus', 'group_order', 'group_counter', 'signature', 'top_index']),
        ("group_order", "group order", ['group_order', 'group_counter', 'genus', 'signature', 'top_index']),
        ("dim", "dimension", ['dim', 'genus', 'group_order', 'group_counter', 'signature', 'top_index']),
        ("full_index", "index [Aut:G]", ['full_index', 'genus', 'group_order', 'group_counter', 'signature', 'top_index']),
        ("num_braid", "braid classes", ['num_braid', 'genus', 'group_order', 'group_counter', 'signature', 'top_index']),
        ("num_vectors", "generating vectors", ['num_vectors', 'genus', 'group_order', 'group_counter', 'signature', 'top_index'])]
    _family_sorts = [
        ("", "genus", ['genus', 'group_order', 'group_counter', 'signature']),
        ("group_order", "group order", ['group_order', 'group_counter', 'genus', 'signature']),
        ("dim", "dimension", ['dim', 'genus', 'group_order', 'group_counter', 'signature']),
        ("num_topological", "topological classes", ['num_topological', 'genus', 'group_order', 'group_counter', 'signature']),
        ("num_vectors", "generating vectors", ['num_vectors', 'genus', 'group_order', 'group_counter', 'signature'])]
    sorts = {"": _top_sorts, "Random": _top_sorts,
             "Families": _family_sorts, "RandomFamily": _family_sorts}

    def __init__(self):
        genus = TextBox(
            name="genus",
            label="Genus",
            knowl="ag.curve.genus",
            example="3",
            example_span="4, or a range like 3..5")
        dim = TextBox(
            name="dim",
            label="Dimension",
            knowl="curve.highergenus.aut.dimension",
            example="1",
            example_span="1, or a range like 0..2")
        group = TextBox(
            name="group",
            label="Group",
            knowl="group.label",
            example="4.2",
            example_span="an LMFDB label, such as 4.2 or 32.47")
        group_order = TextBox(
            name="group_order",
            label="Group order",
            knowl="group.order",
            example="2..5",
            example_span="12, or a range like 10..20, or you may include the variable g for genus like 84(g-1)")
        signature = TextBox(
            name="signature",
            label="Signature",
            knowl="curve.highergenus.aut.signature",
            example="[0,2,3,7]",
            example_span="[0,2,3,7] or [0;2,3,7]")
        periods = TextBox(
            name="periods",
            label="Periods",
            knowl="curve.highergenus.aut.signature",
            example="5",
            example_span="5, or 2,4 for signatures containing all of them")
        r = TextBox(
            name="r",
            label="Branch points",
            knowl="curve.highergenus.aut.branchpoints",
            example="4",
            example_span="4, or a range like 4..6")
        is_full = YesNoBox(
            name="is_full",
            label="Full",
            knowl="curve.highergenus.aut.full")
        aut_order = TextBox(
            name="aut_order",
            label="Order of Aut",
            knowl="curve.highergenus.aut.full",
            example="48",
            example_span="48, or a range like 20..50")
        aut_group = TextBox(
            name="aut_group",
            label="Aut group",
            knowl="curve.highergenus.aut.full",
            example="48.29",
            example_span="an LMFDB label, such as 48.29")
        full_index = TextBox(
            name="full_index",
            label=r"Index $[\mathrm{Aut}:G]$",
            knowl="curve.highergenus.aut.full_index",
            example="2",
            example_span="1 for full actions, or a range like 2..")
        hyperelliptic = YesNoBox(
            name="hyperelliptic",
            label="Hyperelliptic",
            knowl="ag.hyperelliptic_curve")
        cyclic_trigonal = YesNoBox(
            name="cyclic_trigonal",
            label="Cyclic trigonal",
            knowl="ag.cyclic_trigonal")
        jac_dim = TextBox(
            name="jac_dim",
            label="Jacobian factor dim.",
            knowl="curve.highergenus.aut.groupalgebradecomp",
            example="1",
            example_span="1 (an elliptic factor), or a range like 1..2")
        num_braid = TextBox(
            name="num_braid",
            label="Braid classes",
            knowl="curve.highergenus.aut.braid_equivalence",
            example="2..",
            example_span="1, or a range like 2..10")
        num_vectors = TextBox(
            name="num_vectors",
            label="Generating vectors",
            knowl="curve.highergenus.aut.generatingvector",
            example="100..",
            example_span="1, or a range like 100..")
        count = CountBox()

        self.browse_array = [
            [genus, dim],
            [group, group_order],
            [signature, periods],
            [r, is_full],
            [aut_order, full_index],
            [aut_group, hyperelliptic],
            [jac_dim, cyclic_trigonal],
            [num_braid, num_vectors],
            [count]]

        self.refine_array = [
            [genus, dim, group, group_order, signature],
            [periods, r, num_braid, num_vectors, jac_dim],
            [is_full, aut_order, aut_group, full_index],
            [hyperelliptic, cyclic_trigonal]]

        self.family_refine_array = [
            [genus, dim, group, group_order],
            [signature, periods, r, is_full],
            [num_braid, num_vectors]]

    def main_array(self, info):
        if info is None:
            return self.browse_array
        if self._st(info) in ["Families", "RandomFamily"]:
            array = self.family_refine_array
        else:
            array = self.refine_array
        return array + [[SneakyTextBox(name="labels", label="Labels", knowl=self.label_knowl)]]

    def search_types(self, info):
        return self._search_again(info, [
            ('List', 'List of topological classes'),
            ('Families', 'List of families'),
            ('Random', 'Random topological class'),
            ('RandomFamily', 'Random family')])


# Object pages ##############################################################

@higher_genus_w_automorphisms_page.route("/<label>")
def by_label(label):
    clean_label = clean_input(label)
    if clean_label != label:
        return redirect(url_for(".by_label", label=clean_label), 301)
    kind = label_kind(label)
    if kind == "family":
        return render_family_page(label)
    if kind == "topological":
        return render_topological_page(label)
    if kind == "braid":
        # Braid classes are rows of the topological class page
        rec = db.hgcwa_braid.lookup(label, ["topological", "braid_index"])
        if rec is not None:
            return redirect(braid_row_url(rec["topological"], label, rec["braid_index"]))
        flash_error("No braid class with label %s was found in the database.", label)
        return redirect(url_for(".index"))
    if kind == "vector":
        rec = db.hgcwa_vectors.lookup(label, ["braid", "topological"])
        if rec is not None:
            braid = db.hgcwa_braid.lookup(rec["braid"], ["num_vectors", "braid_index"])
            if braid["num_vectors"] == 1:
                return redirect(braid_row_url(rec["topological"], rec["braid"], braid["braid_index"]))
            v = split_label(label)[1][-1]
            start = ((v - 1) // VECTORS_PER_PAGE) * VECTORS_PER_PAGE
            return redirect(url_for(".vector_listing", label=rec["braid"], start=start) + "#" + label)
        flash_error("No generating vector with label %s was found in the database.", label)
        return redirect(url_for(".index"))
    flash_error("No family with label %s was found in the database.", label)
    return redirect(url_for(".index"))


def braid_row_url(topological, braid, braid_index):
    if braid_index > BRAID_ROWS_SHOWN:
        return url_for(".by_label", label=topological, braids="all") + "#" + braid
    return url_for(".by_label", label=topological) + "#" + braid


class FamilyTopologicalArray(EmbeddedSearchArray):
    sorts = [("", "class number", ['top_index']),
             ("is_full", "full", [('is_full', -1), 'top_index']),
             ("full_index", "index [Aut:G]", ['full_index', 'top_index']),
             ("num_braid", "braid classes", ['num_braid', 'top_index']),
             ("num_vectors", "generating vectors", ['num_vectors', 'top_index'])]

    def __init__(self, show_filter=False):
        if show_filter:
            self.refine_array = [[YesNoBox(name="is_full", label="Full", knowl="curve.highergenus.aut.full")]]
        else:
            self.refine_array = []


def render_family_page(label):
    fam = db.hgcwa_families.lookup(label)
    if fam is None:
        flash_error("No family with label %s was found in the database.", label)
        return redirect(url_for(".index"))
    # The filter is useful for families with many topological classes, full and not full
    show_filter = fam["num_topological"] > 15 and 0 < fam["num_topological_full"] < fam["num_topological"]
    info = to_dict(request.args, search_array=FamilyTopologicalArray(show_filter))
    if info.get("download"):
        info["columns"] = family_top_columns
        return HGCWA_download()(info)
    cache = group_cache([fam["group"]])
    tex = group_tex(fam["group"], cache)
    gname = "$%s$" % tex if tex else fam["group"]
    sig = signature_display(fam["signature"])
    info["family"] = fam
    info["group_link"] = group_link(fam["group"], cache)
    info["signature_display"] = sig
    info["title"] = "Family %s of genus %s curves with automorphism group %s" % (label, fam["genus"], gname)
    info["titletag"] = "Family %s of genus %s curves with automorphism group %s" % (label, fam["genus"], fam["group"])
    info["bread"] = family_bread(label, tex)
    info["properties"] = [
        ('Label', label),
        ('Genus', r'\(%s\)' % fam["genus"]),
        ('Group', gname),
        ('Signature', r'\(%s\)' % sig),
        ('Dimension', r'\(%s\)' % fam["dim"]),
        ('Topological classes', r'\(%s\)' % comma(fam["num_topological"])),
        ('Braid classes', r'\(%s\)' % comma(fam["num_braid"])),
        ('Generating vectors', r'\(%s\)' % comma(fam["num_vectors"]))]
    friends = [("Group %s" % gname, url_for("abstract.by_label", label=fam["group"]))]
    if fam["genus"] == 2 and fam["group"] in G2_AUT_GROUPS:
        friends.append((r"Genus 2 curves over $\Q$", url_for("g2c.index_Q", geom_aut_grp_label=fam["group"])))
    info["friends"] = friends
    info["downloads"] = vector_downloads(label)
    info["KNOWL_ID"] = "curve.highergenus.aut.%s" % label
    return render_family(info)


@embed_wrap(
    table=db.hgcwa_topological,
    template="hgcwa-family.html",
    err_title="Higher genus curve family error",
    columns=family_top_columns,
    postprocess=family_top_postprocess,
    learnmore=learnmore_list,
    # Each of the following arguments is set here so that it is overridden when constructing template_kwds,
    # which prioritizes values found in info (which are set in render_family_page() before calling render_family)
    bread=lambda: None,
    properties=lambda: None,
    friends=lambda: None,
    downloads=lambda: None,
    KNOWL_ID=lambda: None,
    titletag=lambda: None,
)
def render_family(info, query):
    query["family"] = info["family"]["label"]
    parse_bool(info, query, "is_full", name="Full")


def render_topological_page(label):
    data = db.hgcwa_topological.lookup(label)
    if data is None:
        flash_error("No topological class with label %s was found in the database.", label)
        return redirect(url_for(".index"))
    family = data["family"]
    G = data["group"]
    sig = data["signature"]
    fam = db.hgcwa_families.lookup(family, ["num_topological"])
    num_top = fam["num_topological"] if fam else None
    cache = group_cache([G, data.get("full_auto")], extra=["representations"])
    d, gens = perm_rep(cache.get(G))
    tex = group_tex(G, cache)
    gname = "$%s$" % tex if tex else G

    # Braid classes and their representative vectors
    show_all = request.args.get("braids") == "all"
    limit = BRAID_ROWS_MAX if show_all else BRAID_ROWS_SHOWN
    braids = list(db.hgcwa_braid.search(
        {"topological": label}, ["label", "braid_index", "num_vectors", "rep_vector", "jacobian_decomp"],
        sort=["braid_index"], limit=limit))
    if data["num_braid"] <= len(braids):
        char_lists = []
        for rec in braids:
            if rec["jacobian_decomp"] not in char_lists:
                char_lists.append(rec["jacobian_decomp"])
    else:
        char_lists = db.hgcwa_braid.distinct("jacobian_decomp", {"topological": label})
    same_chars = len(char_lists) == 1
    rep_labels = [rec["rep_vector"] for rec in braids] + [data["rep_vector"]]
    vecs = {rec["label"]: rec for rec in db.hgcwa_vectors.search({"label": {"$in": rep_labels}}, ["label", "con", "gen_vectors"])}

    def entries(vec):
        return [(cl, entry_data(code, d)) for cl, code in zip(vec["con"], vec["gen_vectors"])]

    braid_rows = []
    for rec in braids:
        row = {"label": rec["label"], "num_vectors": rec["num_vectors"], "group": G,
               "chars": [c[2] for c in (rec["jacobian_decomp"] or [])]}
        vec = vecs.get(rec["rep_vector"])
        if vec is not None:
            for i, entry in enumerate(entries(vec)):
                row["c%s" % (i + 1)] = entry
        braid_rows.append(row)

    info = {"label": label, "data": data, "group": G, "degree": d, "show_all": show_all,
            "num_top": num_top, "same_chars": same_chars}
    info["group_link"] = group_link(G, cache)
    info["signature_display"] = signature_display(sig)
    info["perm_gens"] = ", ".join(perm_display(g) for g in gens)
    info["braid_table"] = {"columns": braid_table_columns(sig, not same_chars), "results": braid_rows}
    info["braid_rows_shown"] = BRAID_ROWS_SHOWN
    info["hyp_involution_display"] = cc_knowl(G, data.get("hyp_involution"))
    info["cinv_display"] = cc_knowl(G, data.get("cinv"))

    # The representative generating vector, one entry per row
    rep = vecs.get(data["rep_vector"])
    if rep is not None:
        info["rep_table"] = {"columns": rep_vector_columns,
                             "results": [{"entry": "c_{%s}" % (i + 1), "perm": perm, "cl": cl, "group": G, "order": m}
                                         for i, ((cl, perm), m) in enumerate(zip(entries(rep), sig[1:]))]}

    # The distinct tuples of conjugacy classes
    if data["num_vectors"] <= CLASS_TUPLE_MAX_VECTORS:
        tuples = class_tuple_summary(db.hgcwa_vectors.search({"topological": label}, ["braid", "con"]))
        for rec in tuples:
            rec["group"] = G
        info["tuple_table"] = {"columns": class_tuple_columns, "results": tuples}

    # The Jacobian
    if data.get("jacobian_dims"):
        info["jacobian"] = jacobian_tex(data["jacobian_dims"])
        if same_chars and char_lists[0]:
            info["jacobian_table"] = {"columns": jacobian_columns,
                                      "results": [{"dim": c[0], "mult": c[1], "char": c[2], "group": G}
                                                  for c in char_lists[0]]}

    # The full automorphism group
    if not data["is_full"]:
        info["full_auto_display"] = (group_link(data["full_auto"], cache) if data.get("full_auto")
                                     else "a group of order $%s$, not in the LMFDB" % data["full_auto_order"])
        if data.get("signH"):
            info["signH_display"] = signature_display(data["signH"])
        if data.get("full_gen_vector"):
            info["full_gen_vector"] = [cycle_string(p) for p in data["full_gen_vector"]]
            info["full_gen_degree"] = len(data["full_gen_vector"][0])

    # Previous and next topological class in the family
    T = data["top_index"]
    if T > 1:
        info["prev"] = "%s.%s" % (family, T - 1)
    if num_top is not None and T < num_top:
        info["next"] = "%s.%s" % (family, T + 1)

    if data["is_full"]:
        full_prop = "yes"
    elif data.get("full_auto"):
        ftex = group_tex(data["full_auto"], cache)
        full_prop = "no: $%s$" % ftex if ftex else "no: %s" % data["full_auto"]
    else:
        full_prop = r"no: order \(%s\)" % data["full_auto_order"]
    properties = [
        ('Label', label),
        ('Genus', r'\(%s\)' % data["genus"]),
        ('Group', gname),
        ('Signature', r'\(%s\)' % info["signature_display"]),
        ('Dimension', r'\(%s\)' % data["dim"]),
        ('Full', full_prop)]
    if data["is_full"] and data.get("hyperelliptic") is not None:
        properties.append(('Hyperelliptic', "yes" if data["hyperelliptic"] else "no"))
    properties += [
        ('Braid classes', r'\(%s\)' % comma(data["num_braid"])),
        ('Generating vectors', r'\(%s\)' % comma(data["num_vectors"]))]

    friends = [("Family %s" % family, url_for(".by_label", label=family)),
               ("Group %s" % gname, url_for("abstract.by_label", label=G)),
               ("Rational character table of %s" % gname, url_for("abstract.Qchar_table", label=G))]
    if not data["is_full"]:
        if data.get("full_label"):
            friends.append(("Family of the full action %s" % data["full_label"],
                            url_for(".by_label", label=data["full_label"])))
        if data.get("full_auto"):
            friends.append(("Full automorphism group %s" % data["full_auto"],
                            url_for("abstract.by_label", label=data["full_auto"])))
    elif data["genus"] == 2 and G in G2_AUT_GROUPS:
        friends.append((r"Genus 2 curves over $\Q$", url_for("g2c.index_Q", geom_aut_grp_label=G)))

    return render_template("hgcwa-topological.html",
                           title="Topological class %s" % label,
                           bread=family_bread(family, tex, [(str(T), " ")]),
                           info=info,
                           properties=properties,
                           friends=friends,
                           downloads=vector_downloads(label),
                           learnmore=learnmore_list(),
                           KNOWL_ID="curve.highergenus.aut.%s" % label)


@higher_genus_w_automorphisms_page.route("/<label>/vectors")
def vector_listing(label):
    """
    All generating vectors of a topological class or a braid class, one page at a time.
    Vectors are numbered within their braid class, so the labels of a page are
    determined by the braid classes, and they are sorted numerically here.
    """
    kind = label_kind(label)
    if kind not in ["topological", "braid"]:
        flash_error("%s is not the label of a topological class or a braid class.", label)
        return redirect(url_for(".index"))
    bquery = {"topological": label} if kind == "topological" else {"label": label}
    braids = list(db.hgcwa_braid.search(bquery, ["label", "num_vectors"], sort=["braid_index"]))
    if not braids:
        flash_error("No %s with label %s was found in the database.",
                    "braid class" if kind == "braid" else "topological class", label)
        return redirect(url_for(".index"))
    topological = label if kind == "topological" else parent_label(label)
    family = split_label(label)[0]
    data = db.hgcwa_topological.lookup(topological, ["group", "signature"])
    G = data["group"]
    cache = group_cache([G], extra=["representations"])
    d, _ = perm_rep(cache.get(G))
    labels = vector_labels(braids)
    try:
        start = max(0, int(request.args.get("start", 0)))
        count = min(1000, max(1, int(request.args.get("count", VECTORS_PER_PAGE))))
    except ValueError:
        start, count = 0, VECTORS_PER_PAGE
    page = labels[start:start + count]
    vecs = sorted(db.hgcwa_vectors.search({"label": {"$in": page}}, ["label", "braid", "con", "gen_vectors"]),
                  key=lambda rec: label_sort_key(rec["label"]))
    rows = []
    for vec in vecs:
        row = {"label": vec["label"], "braid": vec["braid"], "group": G}
        for i, (cl, code) in enumerate(zip(vec["con"], vec["gen_vectors"])):
            row["c%s" % (i + 1)] = (cl, entry_data(code, d))
        rows.append(row)
    info = {"label": label, "kind": kind, "topological": topological, "group": G, "degree": d,
            "number": len(labels), "start": start, "count": count,
            "upper": min(len(labels), start + count),
            "vector_table": {"columns": vector_table_columns(data["signature"], kind == "topological"), "results": rows}}
    if start > 0:
        info["prev_url"] = url_for(".vector_listing", label=label, start=max(0, start - count), count=count)
    if start + count < len(labels):
        info["next_url"] = url_for(".vector_listing", label=label, start=start + count, count=count)
    tex = group_tex(G, cache)
    tail = [(str(split_label(topological)[1][0]), url_for(".by_label", label=topological))]
    if kind == "braid":
        tail.append((str(split_label(label)[1][1]), braid_row_url(topological, label, split_label(label)[1][1])))
    tail.append(("Vectors", " "))
    noun = "braid class" if kind == "braid" else "topological class"
    return render_template("hgcwa-vectors.html",
                           title="Generating vectors in the %s %s" % (noun, label),
                           bread=family_bread(family, tex, tail),
                           info=info,
                           downloads=vector_downloads(label),
                           learnmore=learnmore_list())


@higher_genus_w_automorphisms_page.route("/data/<label>")
def hgcwa_data(label):
    kind = label_kind(label)
    if kind == "family":
        tables, labels, cols, sorts = ["hgcwa_families"], [label], ["label"], [None]
        title = f"Higher genus family - {label}"
    elif kind == "topological":
        tables = ["hgcwa_topological", "hgcwa_braid", "hgcwa_vectors"]
        labels, cols = [label] * 3, ["label", "topological", "topological"]
        sorts = [None, ["braid_index"], None]
        title = f"Higher genus topological class - {label}"
    elif kind == "braid":
        tables, labels, cols, sorts = ["hgcwa_braid", "hgcwa_vectors"], [label] * 2, ["label", "braid"], [None, None]
        title = f"Higher genus braid class - {label}"
    elif kind == "vector":
        tables, labels, cols, sorts = ["hgcwa_vectors"], [label], ["label"], [None]
        title = f"Higher genus generating vector - {label}"
    else:
        return abort(404, f"Invalid label {label}")
    bread = get_bread([(label, url_for_label(label)), ("Data", " ")])
    return datapage(labels, tables, title=title, bread=bread, label_cols=cols, sorts=sorts)


# Learn more ################################################################

@higher_genus_w_automorphisms_page.route("/Completeness")
def completeness_page():
    t = 'Completeness of higher genus curve with automorphisms data'
    bread = get_bread("Completeness")
    return render_template("single.html", kid='rcs.cande.curve.highergenus.aut',
                           title=t,
                           bread=bread,
                           learnmore=learnmore_list())


@higher_genus_w_automorphisms_page.route("/Labels")
def labels_page():
    t = 'Labels for higher genus curves with automorphisms'
    bread = get_bread("Labels")
    return render_template("single.html", kid='dq.curve.highergenus.aut.label',
                           learnmore=learnmore_list(),
                           title=t,
                           bread=bread)


@higher_genus_w_automorphisms_page.route("/Reliability")
def reliability_page():
    t = 'Reliability of higher genus curve with automorphisms data'
    bread = get_bread("Reliability")
    return render_template("single.html",
                           kid='rcs.rigor.curve.highergenus.aut',
                           title=t,
                           bread=bread,
                           learnmore=learnmore_list())


@higher_genus_w_automorphisms_page.route("/Source")
def how_computed_page():
    t = 'Source of higher genus curve with automorphisms data'
    bread = get_bread("Source")
    return render_template("multi.html",
                           kids=['rcs.source.curve.highergenus.aut',
                                 'rcs.ack.curve.highergenus.aut',
                                 'rcs.cite.curve.highergenus.aut'],
                           title=t,
                           bread=bread,
                           learnmore=learnmore_list())
