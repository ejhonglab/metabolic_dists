#!/usr/bin/env python3

# NOTE: pairwise requires python >= 3.10
from itertools import product, pairwise
from pprint import pprint, pformat
import time
from textwrap import dedent
import urllib
from typing import Optional

from neo4j import GraphDatabase
import pandas as pd
from tqdm import tqdm

from metabolike.db.neo4j import Neo4jClient
from metabolike.algo.compounds import CompoundMap
from metabolike.algo.routes import (get_reaction_route_between_compounds,
    get_high_degree_compound_nodes
)
# probably just run everything through pubchem compound ids for now...
from rdkit.Chem.inchi import InchiToInchiKey, MolFromInchi


# TODO ideally only instantiate db / id_table / cpds on first use (/ not on import)

# URI examples: "neo4j://localhost", "neo4j+s://xxx.databases.neo4j.io"
URI = 'neo4j://atlas'
# TODO don't hardcode
AUTH = ('neo4j', 'metabolike')

# TODO need to add port (7687) to end of URI?
db = Neo4jClient(URI, *AUTH)

cm = CompoundMap(db)

id_table = cm.id_table
cpds = cm.cpds

_meta2biocyc = id_table.dropna(subset='metaId').set_index('metaId',
    verify_integrity=True
)

def biocyc2metacyc_ids(biocyc):
    # TODO replace by setting index on biocyc id (w/ verify_integrity=True)
    # (once, globally, maybe on a copy of table)
    mdf = id_table.loc[
        id_table.biocyc == biocyc,
        [c for c in id_table.columns if c != 'biocyc']
    ]

    # TODO delete
    print(f'metaId(s) and compartments for biocyc ID {biocyc}')
    print(mdf.to_string(index=False))
    #

    return list(mdf.metaId)


def metacyc_id2compartment(metacyc_id):
    # will KeyError if metaId not in table. may want to handle some other way.
    return _meta2biocyc.at[metacyc_id, 'compartment']


def metacyc_id2biocyc(metacyc_id):
    # will KeyError if metaId not in table. may want to handle some other way.
    return _meta2biocyc.at[metacyc_id, 'biocyc']


# TODO TODO are rdf nodes and compound nodes 1:1 (how much can we rely on the various
# IDs?)

# TODO TODO see what fraction of chemicals, loading just the sbml xml, have
# inchikey/pubchem id defined (to pick one to use, maybe something else?)
# (may want to do this just by parsing the xml myself)
# TODO TODO and try w/ and w/o loading compounds.dat. is it removing any IDs (by
# accident?)?

# TODO count # of chemicals / reactions (and connected / unconnected) after each step of
# load process (or w/ loading increasing sets of files, starting w/ just sbml)
# TODO parse sbml another way and compare to # chemicals / reactions in there


def all_compound_metacyc_ids(tx):
    query = """
    MATCH (c:Compound)
    WHERE c.metaId IS NOT NULL
    RETURN DISTINCT(c.metaId) as meta_id
    """
    result = tx.run(query)
    records = list(result)  # a list of Record objects
    return records


def _query(tx, query, **kwargs):
    # TODO prefix 'InChI=' here?

    # TODO option to also check all but some stereochemistry layers for an inchi

    assert all([f'${k}' in query for k in kwargs.keys()])
    result = tx.run(query, **kwargs)
    records = list(result)  # a list of Record objects
    summary = result.consume()
    # TODO what is summary useful for anyway? anything?
    return records, summary


#def reactions_with_reactant(tx, inchi):
def reactions_with_reactant(tx, name):
    # TODO TODO TODO diff reaction id besides r.name?
    # canonicalId? metaId?
    #WHERE reactant.inchi = $inchi
    # was getting phenyl[ethanol/acetaldehyde]
    #WHERE reactant.name CONTAINS $name
    query = """
    MATCH (reactant:Compound)<-[x:hasLeft]-(r:Reaction)
    WHERE reactant.name = $name
    RETURN reactant.inchi as reactant_inchi, reactant.commonName AS reactant_name, r.name AS reaction_name
    """
    #return _query(tx, query, inchi=inchi)
    return _query(tx, query, name=name)


#def reactions_with_product(tx, inchi):
def reactions_with_product(tx, name):
    #WHERE product.inchi = $inchi
    #WHERE product.name CONTAINS $name
    query = """
    MATCH (r:Reaction)-[x:hasRight]->(product:Compound)
    WHERE product.name = $name
    RETURN product.inchi as product_inchi, product.commonName AS product_name, r.name AS reaction_name
    """
    #return _query(tx, query, inchi=inchi)
    return _query(tx, query, name=name)


# TODO TODO if lists of reactions share a reaction ID, set length to 1 and don't do next
# query (+ test this case)
# TODO option to also reverse inchi1 & inchi2, and return minimum (defined) distance
# TODO TODO investigate stuff where chemical shows up both as a reactant and as a
# product, and see if i need to handle these cases specially
def min_path_length_between_reactions(tx, start_reactions, end_reactions):
    # TODO more efficient way than looping over all pairs of reactions?
    # must be one...
    query = """
    UNWIND $start_reactions as starting
    UNWIND $end_reactions as ending
    MATCH (src :`Reaction` { name: starting}), (dest:`Reaction`{ name: ending}),
    p = allShortestPaths((src)-[:isPrecedingEvent*]->(dest))
    RETURN p as path, length(p) AS path_length
    """

    # TODO something wrong with this query / the database? why am i not getting any
    # paths?
    result = tx.run(query, start_reactions=start_reactions, end_reactions=end_reactions)
    records = list(result)

    # TODO delete?
    #summary = result.consume()

    min_path_length = min([x['path_length'] for x in records])
    min_length_paths = [x['path'] for x in records]
    return min_path_length, min_length_paths


def min_path_length_between_chemicals(session, name1: str, name2: str,
    either_direction=False) -> int:

    records, summary = session.execute_read(reactions_with_reactant, name1)
    print('reactants:')
    start_reactions = []
    for record in records:
        data = record.data()
        print(data)
        start_reactions.append(data['reaction_name'])

    records, summary = session.execute_read(reactions_with_product, name2)
    print('products:')
    end_reactions = []
    for record in records:
        data = record.data()
        print(data)
        end_reactions.append(data['reaction_name'])

    overlapping_reactions = set(start_reactions) & set(end_reactions)
    if len(overlapping_reactions) > 0:
        print(f'single reaction for reactant->product {overlapping_reactions}')
        return 1
    else:
        print('no single reaction from reactant->product')

    # TODO test w/ reactant from start and product from end

    # Should be length 3 (from looking at example data in neo4j, NOT from
    # groundtruth against metacyc interface)
    start_reactions = ['2PGADEHYDRAT-RXN']
    end_reactions = ['RXN3O-470']

    # Should have 2 paths of length 3.
    #start_reactions = ['PYRUFLAVREDUCT-RXN']
    #end_reactions = ['ACETOACETATE-DECARBOXYLASE-RXN']

    # TODO TODO test behavior if no path + have defined output
    # (that behavesa well w/ either_direction=True computation)

    min_path_length, min_length_paths = session.execute_read(
        min_path_length_between_reactions, start_reactions, end_reactions
    )

    #print('min length paths:')
    #for p in min_length_paths:
    #    pprint(p)

    print('min path length:', min_path_length)

    # + 1 since we want path length between chemicals, not between reactions
    # (which are interposed between chemicals)
    return min_path_length + 1


def read_db_records(db: Neo4jClient, query: str, **kwargs):

    def _read_records(tx, **kwargs):
        result = tx.run(query, **kwargs)
        records = list(result)
        return records

    return db.read_tx(_read_records, **kwargs)


def get_transport_reactions(db: Neo4jClient) -> set[str]:

    # NOTE: there seem to be 24 Reactions (having loaded tier1-tier2/metacyc) where
    # metaId is null. These seem to exactly account for difference between:
    # MATCH (r:Reaction) RETURN...
    # DISTINCT r.metaId vs DISTINCT r

    def print_reaction(reaction):
        print_metacyc_reaction_link(reaction)
        pprint(dict(reaction))

    # TODO delete
    def have_transport_types(reaction: dict) -> bool:
        if not 'types' in reaction:
            return False

        types = set(reaction['types'])
        if any([t.startswith('TR-') for t in types]):
            return True

        elif 'Transport-Reactions' in types:
            return True

        return False

    checks = False
    if checks:
        # TODO change query to also check that other IDs for chemicals (e.g. inchi)
        # IF AVAILABLE, are not the same (and hopefully this just excludes the ~4-6
        # things checks=True path code was turningup)
        #
        # If I didn't want to inspect c1, c2 (and if duplicate r bothered me), could do
        # `RETURN DISTINCT r`
        query_via_relationships = """
        MATCH (c1:Compound)<-[x:hasLeft]-(r:Reaction)-[y:hasRight]->(c2:Compound)
        WHERE c1.name = c2.name AND c1 <> c2
        RETURN r, c1, c2
        """
        records_via_rels = read_db_records(db, query_via_relationships)

        non_cross_compartment_reactions_ids = set()
        for record in records_via_rels:
            reaction = record['r']
            c1 = record['c1']
            c2 = record['c2']

            compartment1 = metacyc_id2compartment(c1['metaId'])
            compartment2 = metacyc_id2compartment(c2['metaId'])

            assert not (pd.isnull(compartment1) or pd.isnull(compartment2))

            # These seemed to all be things that had the same "name" but different (if
            # closely related) structures. Only a single digit number in
            # tier1-tier2/metacyc data.
            if compartment1 == compartment2:
                assert reaction['metaId'] not in non_cross_compartment_reactions_ids
                non_cross_compartment_reactions_ids.add(reaction['metaId'])

                print_reaction(reaction)

                print_metacyc_compound_link(c1)
                print_metacyc_compound_link(c2)
                print()
        print()

        lacking_transport_type = set()
        for record in records_via_rels:
            r = record['r']
            if not have_transport_types(r):
                assert r['metaId'] not in lacking_transport_type
                lacking_transport_type.add(r['metaId'])

                # TODO onnly actually print per canonical id (at least the URLs, since
                # that's all that is used in the URL?)
                print_reaction(r)
                print()

        print()
        assert non_cross_compartment_reactions_ids == lacking_transport_type

        print(f'{len(set(non_cross_compartment_reactions_ids))=}')
        pprint(non_cross_compartment_reactions_ids)
        print()

        via_relationships = [x['r'] for x in records_via_rels]

    # TODO need to deal w/ null types (doesn't seem like it?)? how?
    query_via_types = """
    MATCH (r:Reaction)
    WHERE ANY (
        t IN r.types WHERE (t STARTS WITH 'TR-') OR (t IN ['Transport-Reactions'])
    )
    RETURN r
    """
    via_types = [x['r'] for x in read_db_records(db, query_via_types)]

    # TODO are all things in transport but not transport2 things w/o .types defined?
    # (are there even any such things?)

    # TODO TODO compare via_types / via_relationships
    if checks:
        ids_via_rels =  {x['metaId'] for x in via_relationships}
        ids_via_types = {x['metaId'] for x in via_types}
        assert ids_via_rels - ids_via_types == non_cross_compartment_reactions_ids

        # TODO delete
        print(f'{len(ids_via_rels)=}')
        print(f'{len(ids_via_types)=}')

        # (same 4 as in non_cross... see assertion above)
        print(f'{len(ids_via_rels - ids_via_types)=}')

        print()

        # TODO also look at local graph for each of these (to see why other query not
        # getting them)?
        print('ids_via_types - ids_via_rels:')
        _canonical_id_set = set()
        from collections import defaultdict
        canonical_id2reactions = defaultdict(list)
        for meta_id in ids_via_types - ids_via_rels:
            #print_metacyc_reaction_link(meta_id)
            reaction = [x for x in via_types if x['metaId'] == meta_id][0]
            _canonical_id_set.add(reaction['canonicalId'])
            canonical_id2reactions[reaction['canonicalId']].append(reaction)

        # seems like if the website shows i/o names as the same (for a given
        # canonicalId, which can include a few particular reactions here), there are
        # usually slight mismatches in stereochem / similar across a few versions of a
        # reaction (and i/o doesn't match in each of the particular cases).
        # the differences often reflect even just ambiguity at one stereochemical locus
        # in one of the i/o, and probably treated just as if i/o were identical...
        for canonical_id in _canonical_id_set:
            print_metacyc_reaction_link(canonical_id)

            reactions = canonical_id2reactions[canonical_id]
            if len(reactions) <= 1:
                continue
            pprint([r['name'] for r in reactions])
            print()

        print()
        # This is just 660 / {15774 via_types, 15118 via_rels}
        print(f'{len(ids_via_types - ids_via_rels)=}')
        # This is just 102
        print(len(_canonical_id_set))

        # TODO maybe stuff that transports something BUT adds a phosphate/similar
        # shouldn't be treated any diff than stuff that doesn't change chemical at all?
        # still include the reaction to remove the phosphate / whatever when counting
        # path lengths that should have had transport collapsed?
        #
        # TODO why do these canonicalId(s) seem to have same name in i/o in website:
        # (see long comment a few lines above)
        # - TRANS-RXN-320
        # - RXN0-0
        # - TRANS-RXN-332
        # - ABC-2-RXN
        # - ABC-33-RXN
        # - TRANS-RXN-177
        # - TRANS-RXN66-1159 (tho only Na / phosphate i/o here...)
        # - TRANS-RXN-333
        # ...but they only turn up in via_types and not via_rels?

    # TODO delete
    # TODO maybe still check if any of these seem like true transport reactions, before
    # deleting?
    # ...how would i check though, just if the website somehow called it a transport
    # reaction?
    '''
    transport_str_query = """
    MATCH (r:Reaction)
    WHERE r.canonicalId STARTS WITH "TRANS-"
    RETURN r
    """
    transport_via_str = [x['r'] for x in read_db_records(db, transport_str_query)]

    for r in transport_via_str:
        if not have_transport_types(r):
            print_reaction(r)
            print()

    import ipdb; ipdb.set_trace()
    '''

    # TODO get set of all reaction types and see if it looks like i'm missing anything
    # (w/o filtering reactions first)

    return {x['metaId'] for x in via_types}


def get_blacklist_compound_metaids(db: Neo4jClient) -> list[str]:

    blacklist_compound_names = [
        'H',
        'H+',
        'ATP',
        'ADP',
        'GDP',
        'NADH',
        'NADPH',
        'NADP+',
        'NAD+',
        'phosphate',
        'diphosphate',
        'H2O',
        'CO2',
        'dioxygen',
    ]
    # fill in if need be
    blacklist_compound_metaids = set()

    blacklisted_name_metaids = db.read(
        # TODO TODO why did the distinct seem to matter here? match + uniqueness of
        # metaId (?) shouldn't already guarantee?
        """
        MATCH (c:Compound)<-[:hasLeft|hasRight]-(r:Reaction)
        WHERE c.name IN $blacklist_compound_names
        RETURN COLLECT(DISTINCT(c.metaId)) AS cpds
        """,
        blacklist_compound_names=blacklist_compound_names
    )
    assert len(blacklisted_name_metaids) == 1
    blacklisted_name_metaids = set(blacklisted_name_metaids[0]['cpds'])

    print('metaIds excluded because name in blacklist:')
    pprint(sorted(blacklisted_name_metaids))

    blacklist_compound_metaids.update(blacklisted_name_metaids)

    '''
    MATCH (c:Compound)<-[:hasLeft|hasRight]-(r:Reaction)
    WITH c as c, c.metaId as cpd, COUNT(*) AS cnt
    WHERE cnt >= $degree
    RETURN cnt, cpd, c
    ORDER BY cnt DESC
    '''
    # Looking through the sorted output, might not want to make a whitelist to not
    # blacklist: (count - metaId (synonym))
    # 450 - "AMMONIUM_c" (eh...)
    # 286 - "PYRUVATE_c"
    # 215 - "ACET_c" (acetate)
    # 190 - "FARNESYL__45__PP_c" (?)
    # 180 - "GLUTATHIONE_c" (?)
    # 171 - "GLY_c" (glycine)
    # 154 - "FORMALDEHYDE_c"
    # 147 - "HS_c" (hydrogen sulfide) (?)
    # 87  - "METOH_c" (methanol)
    # 79  - "ACETALD_c" (was there another one w/ higher degree?)
    # 72  - "GLYOX_c" (?)
    # 70  - "GERANYLGERANYL__45__PP_c" (?)
    # 69  - "GLYCEROL__45__3P_c"
    # TODO blacklist based on biocyc IDs? molecular weight? charge?
    # null/'H' chemicalFormula? stuff w/ enzyme in name?
    # all amino acids or something? probably not...
    whitelist_compound_metaids = {
        # TODO TODO probably whitelist other sugars in here too?
        # see what effect it actually has on computed paths
        'AMMONIUM_c',
        'PYRUVATE_c',
        'ACET_c',
        'FARNESYL__45__PP_c',
        'GLUTATHIONE_c',
        'GLY_c',
        'FORMALDEHYDE_c',
        'HS_c',
        'METOH_c',
        'ACETALD_c',
        'GLYOX_c',
        'GLYCEROL_c',
        'Glucopyranose_c',
        'GERANYL__45__PP_c',
        'GERANYLGERANYL__45__PP_c',
        'GLYCEROL__45__3P_c',
        'METHYLAMINE_c',
        'OXALACETIC_ACID_c',
        'ETHANOL__45__AMINE_c',
        'BUTANAL_c',
        'BUTANOL_c',
        'SUCROSE_c',
        'PHE_c',
    }
    degree = 45
    high_deg_metaids = db.read(
        """
        MATCH (c:Compound)<-[:hasLeft|hasRight]-(r:Reaction)
        WITH c.metaId as cpd, COUNT(*) AS cnt
        WHERE cnt >= $degree
        RETURN COLLECT(cpd) AS cpds
        """,
        degree=degree
    )
    assert len(high_deg_metaids) == 1
    high_deg_metaids = set(high_deg_metaids[0]['cpds'])

    high_deg_exclude_metaids = high_deg_metaids - whitelist_compound_metaids
    print(f'metaIds excluded because >={degree} reactions referencing '
        '(and not in whitelist):'
    )
    pprint(sorted(high_deg_exclude_metaids))

    blacklist_compound_metaids.update(high_deg_exclude_metaids)
    return sorted(blacklist_compound_metaids)


# TODO TODO in a more general db setup script, add 'compartment' property to each
# compound (to not need to figure out from suffix on metaId)?
def setup_canproduce_relationships(db: Neo4jClient) -> tuple[list[str], list[str]]:
    """
    For each Reaction node, create relationships between all Compounds in hasLeft and
    hasRight. Seems to facilitate faster/simpler finding of reaction paths between
    compounds.

    Returns compound metaId blacklist.
    """
    # TODO blacklist all / most deacetylase reactions?

    # TODO filter out all reactions w/ 'Electron-Transfer-Reactions' in types?
    # (wouldn't really expect them to be used anyway... could check)

    blacklist_compound_metaids = get_blacklist_compound_metaids(db)

    # TODO print summary of how many relationships were deleted / created
    # (prob need to not use db.write wrapper, so i can get summary object of result)

    # adding this to transaction below seemed to change behavior when i tried
    # interactively in web interface (was taking a long time)
    db.write('MATCH ()-[r:canProduce]-() DELETE r')

    db.write(
        """
        WITH $metaIdBlacklist as blacklist
        MATCH (a:Compound)<-[x:hasLeft]-(r:Reaction)-[y:hasRight]->(b:Compound)
        WHERE NOT(a.metaId IN blacklist) AND NOT(b.metaId in blacklist)
        CALL apoc.create.relationship(a, "canProduce", properties(r), b)
        YIELD rel
        RETURN null
        """,
        metaIdBlacklist=blacklist_compound_metaids
    )

    # TODO TODO could also define similar extra metadata just for transport reactions
    # TODO TODO define something in database to indicate blacklisted stuff, to aid in
    # manual queries using same blacklist?
    #
    # maybe some node that connects to blacklisted stuff (would change existing db
    # minimally and be easy to add/remove, but might make some queries more
    # complicated)?
    # https://stackoverflow.com/questions/44327169
    #
    # via an additional label to each node? (would have to slightly modify some of the
    # assertion on labels)
    # https://stackoverflow.com/questions/21625081
    #
    # an additional property?

    return blacklist_compound_metaids


def format_compound(compound: dict):
    # NOTE: I used the query:
    # MATCH (c:Compound)-[]-(r:Reaction)
    # WHERE c.name IS NULL
    # RETURN c
    # ...to verify that no compounds THAT ARE CONNECTED TO REACTIONS have missing .name

    assert compound is not None
    # This dict-like access works if compound is a Node type (from a Compound node)
    # (intending to mostly / always take dict input now though, converted from Node /
    # etc)
    if 'commonName' in compound:
        return compound['commonName']

    # TODO just always use compound['name'], for consistency? if more useful in lookups?
    elif 'name' in compound:
        return compound['name']

    # TODO nicer fallback here? metaId?
    else:
        # TODO delete
        pprint(compound)
        import ipdb; ipdb.set_trace()
        #
        return pformat(compound)


def print_metacyc_reaction_link(reaction: dict|str) -> None:

    if isinstance(reaction, str):
        canonical_id = reaction
    else:
        # This seems to correspond to "BioCyc ID" on the metacyc website pages for each
        # reaction.
        canonical_id = reaction['canonicalId']
        if canonical_id is None:
            return

    # TODO do different reactions w/ diff "name" in metabolike correspond to different
    # entities on metacyc website (for the same canonical ID)?
    url = ('https://metacyc.org/META/NEW-IMAGE?type=REACTION&object='
        f'{urllib.parse.quote(canonical_id)}'
    )
    print(url)


def print_biocyc_compound_link(biocyc_id: str) -> None:
    url = f'https://metacyc.org/compound?orgid=META&id={urllib.parse.quote(biocyc_id)}'
    print(url)


def print_metacyc_compound_link(compound: dict|str) -> None:
    if isinstance(compound, str):
        meta_id = compound
    else:
        meta_id = compound['metaId']
        assert meta_id is not None

    biocyc_id = metacyc_id2biocyc(meta_id)
    assert not pd.isnull(biocyc_id)

    print_biocyc_compound_link(biocyc_id)


def print_path(path: dict) -> None:
    compounds = path['compounds']
    reactions = path['reactions']
    assert len(reactions) + 1 == len(compounds)

    for i, reaction in enumerate(reactions):
        print(f'{i + 1}:')

        c1 = compounds[i]
        c2 = compounds[i + 1]

        if c1 is not None and c2 is not None:
            print(f'{format_compound(c1)} -> {format_compound(c2)}')

        elif c1 is not None:
            print(f'{format_compound(c1)} -> ?')
        elif c2 is not None:
            print(f'? -> {format_compound(c2)}')
        print_metacyc_reaction_link(reaction)
        # TODO only print if != canonicalId?
        print(f'name: {reaction["name"]}')

    print()


def print_paths(paths: list[dict], limit: Optional[int] = 5) -> None:
    for i, path in enumerate(paths):
        if limit is not None and i == limit:
            print(f'printed {limit=} paths ({len(paths) - limit} more)')
            return

        if len(paths) > 1:
            print(f'path {i}:')

        print_path(path)


# TODO TODO cache (across runs)
# TODO import correct Node type and hint that as return type
def compound_between_pathway_reactions(db: Neo4jClient, r1: dict, r2: dict,
    ignore_node_metaids=None):

    if ignore_node_metaids is None:
        ignore_node_metaids = []

    # TODO TODO if blacklist is ever generalized (in my usage) to ignore reactions too,
    # need to change query
    #
    # There should also be a isPrecedingEvent relationship from r1->r2, but that
    # probably shouldn't change the results.
    query = """
    MATCH (r1:Reaction)-[:hasRight]->(c:Compound)<-[:hasLeft]-(r2:Reaction)
    WHERE r1.metaId = $r1 AND r2.metaId = $r2 AND NOT(c.metaId IN $ignore_nodes)
    RETURN c
    """
    compounds = read_db_records(db, query, r1=r1['metaId'], r2=r2['metaId'],
        ignore_nodes=ignore_node_metaids
    )

    if len(compounds) == 0:
        # TODO if blacklist is ever generalized (in my usage) to ignore reactions too,
        # need to change query
        query = """
        MATCH (r1:Reaction)-[:hasRight|hasLeft]->(c:Compound)<-[:hasLeft|hasRight]-(r2:Reaction)
        WHERE r1.metaId = $r1 AND r2.metaId = $r2 AND NOT(c.metaId IN $ignore_nodes)
        RETURN c
        """
        compounds = read_db_records(db, query, r1=r1['metaId'], r2=r2['metaId'],
            ignore_nodes=ignore_node_metaids
        )

    # TODO TODO TODO fix. seems to fail in cases where isPrecedingEvent sequence
    # has an r1->r2 where >1 is reversible, and thus hasLeft/hasRight lose some meaning.
    # will always treating hasLeft|hasRight as interchangeable cause more problems?
    try:
        assert len(compounds) == 1
    except AssertionError:
        print()
        print(dedent(query.replace('$r1', r1['metaId']).replace('$r2', r2['metaId'])))

        print()
        print(f'{len(compounds)=}')
        print()

        print('r1: ', end='')
        print_metacyc_reaction_link(r1)
        pprint(r1)

        print('r2: ', end='')
        print_metacyc_reaction_link(r2)
        pprint(r2)

        import ipdb; ipdb.set_trace()

    compound = compounds[0]['c']

    # TODO delete
    print('r1: ', end='')
    print_metacyc_reaction_link(r1)

    print('compound:')
    pprint(dict(compound))

    print('r2: ', end='')
    print_metacyc_reaction_link(r2)
    print()
    #
    return compound


def min_pathway_reactions_between_chemicals(db: Neo4jClient, c1: str, c2: str,
    max_hops=12, ignore_node_metaids=None):

    if ignore_node_metaids is None:
        ignore_node_metaids = []

    # TODO probably replace w/ my code w/ pratyush, b/c this will not return
    # all shortest paths, if there are ever multiple of same length (and may
    # want to inspect each? not sure it really matters though)
    #
    # TODO what is the meaning of the `c2, ns, ...` in the second WITH
    # clause?
    #
    # Adapted from only_pathway_reactions=True path in
    # metabolike.algo.routes.get_reaction_route_between_compounds
    # I set bfs=True, so now this should find shortest paths.
    query = """
    MATCH (n)
    WHERE n.metaId IN $ignore_nodes
    WITH COLLECT(n) AS ns

    MATCH (r:Reaction)-[:hasRight]->(c2:Compound {metaId: $c2})
    WITH c2, ns, COLLECT(r) AS rxns
    MATCH (c1:Compound {metaId: $c1})
    CALL apoc.path.expandConfig(c1, {
        relationshipFilter: "<hasLeft|isPrecedingEvent>|isRelatedEvent>",
        labelFilter: "+Reaction",
        terminatorNodes: rxns,
        blacklistNodes: ns,
        bfs: true,
        limit: $num_routes,
        minLevel: 2,
        maxLevel: $max_hops
    })
    YIELD path
    RETURN c2, path, length(path) AS hops
    ORDER BY hops;
    """
    # In the Extended Data Figure 4 in:
    # "Metabolic activity organizes olfactory representations"
    # https://www.biorxiv.org/content/10.1101/2022.07.21.500995v3
    # their distribution of reaction distances tops out at 12 hops.
    routes = read_db_records(db, query, c1=c1, c2=c2, max_hops=12,
        # TODO TODO by default(?), ignore_nodes=blacklist_compound_metaids
        # (but maybe try w/o, as pathway stuff more curated anyway, and
        # certainly is a subset of the reactions)
        num_routes=1, ignore_nodes=[]
    )
    if len(routes) == 0:
        # Returning list of paths to be consistent w/ similar fns that actually can
        # return all shortest paths, or in case we want that behavior here in the
        # future.
        return float('inf'), []

    # with `limit: 1`, this should be the case
    assert len(routes) == 1
    route = routes[0]
    min_path_len = route['hops']
    assert min_path_len >= 1

    path = route['path']

    assert path.relationships[0].type == 'hasLeft'
    try:
        assert all([x.type == 'isPrecedingEvent' for x in path.relationships[1:]])
    except AssertionError:
        print([x.type == 'isPrecedingEvent' for x in path.relationships[1:]])


    reaction_nodes = path.nodes[1:]
    reactions = [dict(x) for x in path.nodes[1:]]
    assert len(reactions) == min_path_len
    assert set([x.labels for x in path.nodes[1:]]) == {frozenset({'Reaction'})}

    # some edge cases where this isn't working, and only needed it so i could use
    # compound sequences as a means of filtering out transport reactions. if i can
    # figure out another way to do that, i don't really need this
    #
    # TODO TODO possible to use set of reactions before/after ambiguous steps to
    # disambiguate which chemical to choose?
    '''
    # pairwise('ABCDEFG') --> AB BC ... FG
    intermediates = [
        compound_between_pathway_reactions(db, r1, r2,
            ignore_node_metaids=ignore_node_metaids
        )
        for r1, r2 in pairwise(reactions)
    ]
    '''
    intermediates = (min_path_len - 1) * [None]

    # TODO metabolike have something that can get intermediate compounds more easily?
    # TODO TODO possible to get intermediate compounds (without having to guess between
    # all hasLeft/hasRight for a given (pair of) reactions?). just want for
    # troubleshooting, so not necessarily that critical.
    # For now, just going to set all non-[first|last] elements in compounds to None
    # here, and won't print them downstream.
    #compound_nodes = [path.nodes[0]] + (min_path_len - 1) * [None] + [route['c2']]

    compound_nodes = [path.nodes[0]] + intermediates + [route['c2']]

    # TODO remove references to None if above intermediates calculation strategy works
    assert (
        set([x.labels for x in compound_nodes if x is not None]) ==
        {frozenset({'Compound'})}
    )
    compounds = [dict(x) if x is not None else x for x in compound_nodes]

    # TODO TODO possible to query for each compound based on which reactions it's
    # connected to (providing the two adjacent for each step in the pathway)?
    # would that always uniquely identify compounds?

    path = {
        'compounds': compounds,
        'reactions': reactions,
    }
    return min_path_len, [path]


def min_reactions_between_chemicals(db: Neo4jClient, c1: str, c2: str):
    # TODO limit needed (probably...)? any problems (shortest path finding not optimal
    # ig...)?
    query = """
    MATCH (c1:Compound {metaId: $c1}), (c2:Compound {metaId: $c2}),
        p = allShortestPaths((c1)-[:canProduce*]->(c2))

    RETURN p as path, length(p) as hops
    LIMIT 100
    """
    # doesn't seem like i need the relationships(p) call actually, now that i'm doing
    # the read differently?
    #RETURN relationships(p) as rs, p as path, length(p) as hops

    # The .data() call in db.read() was not giving me the
    # relationships(path) with the information I needed.
    routes = read_db_records(db, query, c1=c1, c2=c2)
    if len(routes) == 0:
        return float('inf'), []

    assert len(set([x['hops'] for x in routes])) == 1
    min_path_len = routes[0]['hops']

    paths = []
    for route in routes:
        path = route['path']

        assert all([
            x.labels == frozenset({'Compound'}) for x in path.nodes
        ])
        compounds = [dict(x) for x in path.nodes]

        assert set([x.type for x in path.relationships]) == {'canProduce'}
        reactions = [dict(x) for x in path.relationships]
        assert len(reactions) == min_path_len

        assert len(compounds) == len(reactions) + 1
        paths.append({
            'compounds': compounds,
            'reactions': reactions,
        })

    return min_path_len, paths


# TODO TODO TODO test in cases where first/last/middle reaction is a transport reaction
# TODO TODO replace w/ collapsing in db before searching / in query, in case there are
# weird cases where doing so changes which path is the shortest
def collapse_transport_rxns(transport_reaction_metaids: set[str], path: dict) -> dict:

    # TODO pathway reactions even have transport reactions in them?
    # TODO need to collapse pre/during shortest path search, or doesn't matter?

    assert len(path['compounds']) >= 2

    # TODO TODO TODO assert all reactions not determined to be transport reactions DONT
    # start with 'TRANS-' and all others DO (or same assertion for whatever other
    # criteria i might want to use)
    # TODO TODO do all transport reactions have 'TR-12' (or 'TR-<x>') in types?
    # TODO do all transport reactions have candonical ID / name starting with 'TRANS-'?
    # TODO and, if so, are all reactions with these properties transport reactions?
    # TODO all transport reactions have gibbs0 == 0.0? uniquely?
    compounds = []
    reactions = []

    # TODO delete
    _had_transport = False

    for (c1, c2), reaction in zip(pairwise(path['compounds']), path['reactions']):
        # TODO print what is collapsed?
        if reaction['metaId'] not in transport_reaction_metaids:
            compounds.append(c1)
            reactions.append(reaction)
        # TODO delete
        else:
            _had_transport = True

    # TODO delete
    # TODO warn if any of c1/c2 would be in blacklist (would be ignored w/ current
    # pathway-only implementation)
    '''
        reaction_id = reaction['canonicalId']

        c1_biocyc = metacyc_id2biocyc(c1['metaId'])
        c2_biocyc = metacyc_id2biocyc(c2['metaId'])

        if c1_biocyc == c2_biocyc:
            try:
                assert any([x.startswith('TR-') for x in reaction['types']])
            except AssertionError:
                print('reaction:')
                pprint(dict(reaction))
                import ipdb; ipdb.set_trace()
            #try:
            #    assert reaction_id.startswith('TRANS-'), reaction_id
            #except AssertionError as err:
            #    print(f'{reaction_id=}')
            #    pprint(reaction)
            #    import ipdb; ipdb.set_trace()

            # TODO assert only metaId key differs in dict of c1 vs c2?
        else:
            #assert not reaction_id.startswith('TRANS-'), reaction_id
            try:
                assert not any([x.startswith('TR-') for x in reaction['types']])

            except AssertionError:
                print(format_compound(c1))
                print(format_compound(c2))
                print_metacyc_reaction_link(reaction)
                print('reaction:')
                pprint(dict(reaction))
                import ipdb; ipdb.set_trace()

            compounds.append(c1)
            reactions.append(reaction)
    '''

    # TODO what if last reaction is a transport reaction? i think it still works...
    compounds.append(c2)

    assert len(compounds) == len(reactions) + 1
    new_path = {
        'compounds': compounds,
        'reactions': reactions,
    }

    # TODO delete / put behind checks flag
    if len(reactions) != len(path['reactions']):
        assert _had_transport

        # TODO want to also limit this somehow? would be easier to do that in wrapper,
        # but would be easier to compare to old path here...
        print('original path:')
        print_path(path)
        print('path with transport reactions collapsed:')
        print_path(new_path)

        print(f'old length: {len(path)}, after collapsing transport: {len(new_path)}')
        print()
    else:
        assert not _had_transport
        assert (
            [x['metaId'] for x in reactions] ==
            [x['metaId'] for x in path['reactions']]
        )
        assert (
            [x if x is None else x['metaId'] for x in compounds] ==
            [x if x is None else x['metaId'] for x in path['compounds']]
        )
    #

    return new_path


def collapse_transport_rxns_across_paths(transport_reaction_metaids: set[str],
    paths: list[dict]) -> tuple[int, list[dict]]:

    #  TODO refactor to share w/ below
    #old_min_len = min([len(p['reactions']) for p in paths])

    paths = [collapse_transport_rxns(transport_reaction_metaids, p) for p in paths]
    new_min_len = min([len(p['reactions']) for p in paths])

    #if old_min_len != new_min_len:

    paths = [
        p for p in paths if len(p['reactions']) == new_min_len
    ]
    return new_min_len, paths


def main():
    # TODO TODO use remy's spreadsheet that specifically includes metacyc IDs for the
    # hallem odors (partially manually looked up). and was her automated process more
    # succesful than mine?
    my_hdf = pd.read_csv('hallem_inchis.csv')
    hdf = pd.read_excel('hallem110_pubchempy.xlsx')

    # NOTE: column 1 starts w/ "PubChem:"for ones w/o successful lookup
    remy_metacyc_df = pd.read_excel('hallem110_from_pubchem__ALL__common_names.xlsx',
        sheet_name='hallem110_from_pubchem__ALL__common_names'
    )

    assert (
        hdf.biocyc_pubchem.apply(lambda x: x[len('PubChem:'):]).astype('int64') ==
        hdf.cid
    ).all()

    assert ('InChI=' + my_hdf.inchi == hdf.inchi).all()
    del my_hdf

    # TODO TODO modify code that computes cpds to also add inchi. currently have:
    #ipdb> pp [x for x in cpds.key_id.unique()]
    #['commonName',
    # 'smiles',
    # 'synonyms',
    # 'name',
    # 'drugbank',
    # 'inchikey',
    # 'chebi',
    # 'pubchemCompound',
    # 'keggGlycan',
    # 'lipidmaps',
    # 'metabolights',
    # 'hmdb',
    # 'cas',
    # 'keggCompound',
    # 'chemspider',
    # 'knapsack',
    # 'umbbdCompound']

    pubchem_df = cpds[cpds.key_id == 'pubchemCompound']
    assert all([str(int(x)) == x for x in pubchem_df.value])
    db_pubchem_ids = set(pubchem_df.value)

    # TODO TODO probably use 'AMMONIUM[_c]' for hallem ammonium hydroxide?

    cids_not_in_db = set()
    not_found_in_db = set()
    # TODO TODO refactor to module-level + add fn that calls this on an iterable of
    # input, and reports cids/not_found for all after (-> remove global set vars)
    def get_metacyc_ids(chemical) -> Optional[list[str]]:
        #print(chemical.name)

        cid = str(chemical.cid)
        if cid in db_pubchem_ids:
            # didn't want to risk matching against small integer CIDs, and having
            # partial / exact matches w/ the same integer used in a diff context
            # (could have probably just never used partial matches, but might still have
            # been some ambiguous cases...)
            sdf = pubchem_df[pubchem_df.value == cid]

            #print(sdf.to_string(index=False))
            #print()

            # would need to deal with otherwise...
            assert len(sdf) == 1
            biocyc = sdf.biocyc.iat[0]

            meta_ids = biocyc2metacyc_ids(biocyc)

            assert len(meta_ids) > 0
            return meta_ids

        cids_not_in_db.add(cid)

        inchikey = InchiToInchiKey(chemical.inchi)
        assert inchikey == chemical.inchikey

        # TODO smiles after inchikey? check if it's ever available when cid/inchikey
        # aren't
        additional_ids_to_try = ['inchikey', 'name', 'iupac_name']
        meta_ids = []
        for chem_id in additional_ids_to_try:
            # TODO edit type annotation of output of this fn to indicate it's optional
            hits = cm.compound_exact_match(getattr(chemical, chem_id))
            if hits is None:
                continue

            assert len(hits) > 0

            print(f'(get_metacyc_ids) {chem_id=}')
            print(f'(get_metacyc_ids) chemical.{chem_id}={getattr(chemical, chem_id)}')
            print(f'(get_metacyc_ids) {chem_id} matches: ', end='')
            pprint(hits)
            print()

            for biocyc in hits:
                hit_meta_ids = biocyc2metacyc_ids(biocyc)
                assert len(hit_meta_ids) > 0
                meta_ids.extend(hit_meta_ids)

            assert len(meta_ids) > 0
            return meta_ids

        # no way this could work as is. not in IDs metabolike computes.
        # InChI is in db (property under some Compounds, if i recall correctly), so
        # could change code to use.
        '''
        print('InChI result:')
        res = cm.search_compound_biocyc_id(chemical.inchi)
        pprint(res)
        print()
        '''

        not_found_in_db.add(chemical.name)

        print('(get_metacyc_ids) not found in db')
        print()
        return None

    transport_reaction_metaids = get_transport_reactions(db)

    # TODO call in separate script, maybe one that does other database initialization
    # stuff too (would need to trust that if we get the blacklist in the script
    # separately, for use in queries not involving canProduce relationships, that it
    # matches the blacklist used when setting up existing canProduce relationships...)
    #
    # TODO try to replace all relationship creation->blacklisting->deletion stuff with a
    # graph projection? would memory be an issue?
    blacklist_compound_metaids = setup_canproduce_relationships(db)

    # If False, blacklist will only apply to non-pathway-only searches.
    also_blacklist_in_pathway_search = True
    if also_blacklist_in_pathway_search:
        pathway_ignore_node_metaids = blacklist_compound_metaids
    else:
        pathway_ignore_node_metaids = []

    biocyc_id_set = set(id_table.biocyc)

    remy_ids_not_in_db = set()
    name2remy_biocyc_id = dict()

    # TODO TODO use this to print those w/ multiple biocyc ids
    # (to decide if i want to manually disambiguate any of them)
    name2biocyc_ids = dict()
    name2meta_ids = dict()
    for (_, remy_row), row in zip(remy_metacyc_df.iterrows(), hdf.itertuples()):
        name = row.name
        # +2 since spreadsheet rows start from 1, and then first row is header
        print(f'name: {name} (spreadsheet row {row.Index + 2})')
        meta_ids = get_metacyc_ids(row)
        name2meta_ids[name] = meta_ids

        remy_biocyc = remy_row['Object ID']
        if pd.notna(remy_biocyc):

            print("Remy's biocyc ID:")
            print_biocyc_compound_link(remy_biocyc)
            name2remy_biocyc_id[name] = remy_biocyc

            if remy_biocyc not in biocyc_id_set:
                remy_ids_not_in_db.add(remy_biocyc)
                remy_biocyc = None

            #assert remy_biocyc in biocyc_id_set

        # TODO TODO maybe ammonium hydroxide should really map to ammonium in the db
        # though? or both? remy has NH4OH in her spreadsheet
        if meta_ids is None:

            # refactor
            if pd.notna(remy_biocyc):
                meta_ids = biocyc2metacyc_ids(remy_biocyc)
                if meta_ids is not None:
                    name2meta_ids[name] = meta_ids

            print()
            print()
            continue

        biocyc_ids = {metacyc_id2biocyc(x) for x in meta_ids}
        name2biocyc_ids[name] = biocyc_ids
        for b in biocyc_ids:
            print('biocyc ID(s) from lookup:')
            print_biocyc_compound_link(b)

        if len(biocyc_ids) > 1:
            print(f'multiple biocyc IDs:\n{pformat(biocyc_ids)}\n...for metacyc IDs: '
                f'{pformat(meta_ids)}'
            )

        if pd.notna(remy_biocyc):
            #assert remy_biocyc in biocyc_ids, f'{remy_biocyc=}'
            if remy_biocyc not in biocyc_ids:
                print(f'{remy_biocyc=} not in output of lookup {biocyc_ids}!')
                import ipdb; ipdb.set_trace()

            # TODO hardcode using remy's here (if multiple biocyc_ids)?

            # TODO TODO separately categorize stuff where remy lists only one of the
            # biocyc IDs we have here?

            # TODO TODO summarize differences here

        '''
        try:
            assert len(biocyc_ids) == 1
        except:
            print(f'{len(biocyc_ids)=}')
            import ipdb; ipdb.set_trace()

        biocyc_id = biocyc_ids.pop()

        # gonna work? (no but at least initially b/c biocyc_ids sometimes > len 1...)
        assert remy_biocyc == biocyc_id
        '''
        print()
        print()

    print()
    print("IDs in Remy's spreadsheet, but not in db biocyc IDs:")
    # TODO confirm these not in db (via separate query? metaId STARTS WITH?
    # def not in _meta2biocyc right?
    pprint(remy_ids_not_in_db)

    print()
    print(f'{len(cids_not_in_db)} CIDs not in db')
    print(f'{len(not_found_in_db)} not found in db')
    print()
    import ipdb; ipdb.set_trace()

    names = list(name2meta_ids.keys())

    # TODO separate matrices pre-collapsing (esp until i'm confident in that...)?

    pathway_min_dists = pd.DataFrame(index=names, columns=names, data=float('inf'))
    pathway_min_dists.index.name = 'from'
    pathway_min_dists.columns.name = 'to'

    min_dists = pd.DataFrame(index=names, columns=names, data=float('inf'))
    min_dists.index.name = 'from'
    min_dists.columns.name = 'to'

    # TODO check these work
    pathway_min_paths = pd.DataFrame(index=names, columns=names, data=float('nan'),
        dtype='object'
    )
    pathway_min_paths.index.name = 'from'
    pathway_min_paths.columns.name = 'to'

    min_paths = pd.DataFrame(index=names, columns=names, data=float('nan'),
        dtype='object'
    )
    min_paths.index.name = 'from'
    min_paths.columns.name = 'to'

    # TODO need to manually specify len in tqdm? regardless, might want to so i can
    # advance in inner loop...
    for name1, meta_ids1 in tqdm(name2meta_ids.items()):

        if meta_ids1 is None:
            print(f'no metacyc IDs for {name1}')
            print()
            continue

        for name2, meta_ids2 in name2meta_ids.items():
            if name1 == name2:
                continue

            print(f'{name1=}')
            print(f'{name2=}')

            if meta_ids2 is None:
                print(f'no metacyc IDs for {name2}')
                print()
                continue

            # TODO profile
            for c1, c2 in product(meta_ids1, meta_ids2):
                print(f'getting reaction route between {c1} and {c2}')

                before = time.time()
                # TODO TODO check against my implementation w/ Pratyush (w/ the 3
                # separate queries)
                # TODO replace w/ version of pratyush's implementation of pathway route
                # searching, where multiple can be returned (or modify to start w/ low
                # maxHops and increase until we find stuff (and no limit on # paths
                # returned thre)
                pathway_len, pathway_paths = min_pathway_reactions_between_chemicals(
                    db, c1, c2, ignore_node_metaids=pathway_ignore_node_metaids
                )
                query_s = time.time() - before
                #print(f'pathway-only shortest-path query took {query_s:.2f}s')

                if len(pathway_paths) > 0:
                    # TODO flag to decide whether to collapse transport reactions?
                    # maybe they should be counted?

                    pathway_len, pathway_paths = collapse_transport_rxns_across_paths(
                        transport_reaction_metaids, pathway_paths
                    )
                    if pathway_len < pathway_min_dists.at[name1, name2]:
                        pathway_min_dists.at[name1, name2] = pathway_len
                        pathway_min_paths.at[name1, name2] = pathway_paths

                    # NOTE: current implementation will only ever return one path, even
                    # of there are two of same length
                    print('pathway routes:')
                    print_paths(pathway_paths)


                # TODO do all reactions crossing pathway boundaries at some point have a
                # transport reaction (only compartment -> compartment, w/o changing main
                # chemical)? are these kinds of reactions the only thing that cross
                # compartment boundaries?

                before = time.time()
                reaction_len, reaction_paths = min_reactions_between_chemicals(
                    db, c1, c2
                )
                query_s = time.time() - before
                #print(f'shortest-path query took {query_s:.2f}s')

                # TODO TODO TODO ensure any transport reactions that can move multiple
                # chemicals don't end up effectively converting one of the possible
                # transported chemicals into the other (whether collapsing or not)
                # (see that one putrescine/cadeverine one? or it was in path there?)
                # https://metacyc.org/META/NEW-IMAGE?type=REACTION&object=TRANS-RXN0-211

                # TODO query like that for via_relationships in get_transport..., but
                # with name1 <> name2?

                # TODO TODO TODO any other cases where reactions might incorrectly seem
                # like they can convert one chemical to another? any reaction ~class
                # (what was their term again?) in here, where it matters which
                # particular input you have as to which particular output you'll get?

                if len(reaction_paths) > 0:
                    reaction_len, reaction_paths = collapse_transport_rxns_across_paths(
                        transport_reaction_metaids, reaction_paths
                    )
                    if reaction_len < min_dists.at[name1, name2]:
                        min_dists.at[name1, name2] = reaction_len
                        min_paths.at[name1, name2] = reaction_paths

                    print(f'general routes ({len(reaction_paths)}):')
                    # TODO don't print any also in pathway routes (presumably those ones
                    # should also show up here? check that?)
                    print_paths(reaction_paths)

                # TODO TODO instead of collapsing transport reactions after, use more
                # flexible shortest path algorithm, and specify costs of transport
                # reactions as 0 (if they are not generic reactions w/ diff i/o)?
                # TODO TODO and set cost to inf instead of  0 for stuff where i/o
                # chemicals are not the same (where i=o and generic reaction has
                # multiple such pairs)?


                # took too long sometimes (at least on atlas w/ ~15Gb free memory and
                # w/ db populated using tier1-tier2/metacyc folder). unclear if it would
                # terminate. maybe blacklist would also help?
                #reaction_paths = get_reaction_route_between_compounds(db, c1, c2,
                #    only_pathway_reactions=False
                #)

                # TODO assert that if there are pathway paths, there are also general
                # paths (and they should all be <= min length of pathway paths) (?)


                # TODO delete
                '''
                # TODO convert this check into a query checking all (pathway) reactions
                # are only between compounds w/ metaIds referring to diff compartments
                #
                # To check if reactions are always compartment->compartment, b/c if so
                # we could simplify some of this code in the meantime, but we may also
                # want to modify the graph to ignore which compartment a reaction was
                # recorded in?
                compartment1 = metacyc_id2compartment(c1)
                compartment2 = metacyc_id2compartment(c2)
                across_compartments = compartment1 != compartment2
                if across_compartments:
                    assert len(pathway_paths) == 0

                # assertion not true
                #if across_compartments:
                #    assert len(reaction_paths) == 0
                '''

            # TODO TODO TODO also load these and resume from them

            print('saving pathway distances and paths...', flush=True, end='')
            pathway_min_dists.to_pickle('pathway_min_dists.p')
            pathway_min_paths.to_pickle('pathway_min_paths.p')
            print(' done', flush=True)

            print('saving distances and paths...', flush=True, end='')
            min_dists.to_pickle('min_dists.p')
            min_paths.to_pickle('min_paths.p')
            print(' done', flush=True)

    import ipdb; ipdb.set_trace()

    '''
    with GraphDatabase.driver(URI, auth=AUTH) as driver:
        driver.verify_connectivity()

        # TODO TODO TODO find some >=2-step tests
        # pyruvate to etoh? glucose (sterochem?) to etoh?

        with driver.session(database="neo4j") as session:
            meta_ids = session.execute_read(all_compound_metacyc_ids)

            # ethanol
            #'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
            name1 = 'ethanol'

            # acetaldehyde. should have 1 step to ethanol
            #'InChI=1S/C2H4O/c1-2-3/h2H,1H3'
            name2 = 'acetaldehyde'

            min_path_length = min_path_length_between_chemicals(session, name1, name2)

            import ipdb; ipdb.set_trace()
    '''


if __name__ == '__main__':
    main()

