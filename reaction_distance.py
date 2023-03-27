#!/usr/bin/env python3

from itertools import product
import time
from pprint import pprint, pformat

from neo4j import GraphDatabase
import pandas as pd

from metabolike.db.neo4j import Neo4jClient
from metabolike.algo.compounds import CompoundMap
from metabolike.algo.routes import (get_reaction_route_between_compounds,
    get_high_degree_compound_nodes
)
# probably just run everything through pubchem compound ids for now...
from rdkit.Chem.inchi import InchiToInchiKey, MolFromInchi

# TODO TODO TODO are rdf nodes and compound nodes 1:1 (how much can we rely on the
# various IDs?)

# TODO TODO debug queries for # of chemicals, reactions, and reactions with chemicals
# (& vice versa). also one to get unconnected chemicalis? to see what is adding those

# TODO TODO see what fraction of chemicals, loading just the sbml xml, have
# inchikey/pubchem id defined (to pick one to use, maybe something else?)
# (may want to do this just by parsing the xml myself)
# TODO TODO and try w/ and w/o loading compounds.dat. is it removing any IDs (by
# accident?)?
# TODO TODO some way i could modify database myself simply, to add relationships
# searchable like via isPrecedingEvent (but in a way that doesn't require loading
# reactions.dat?)

# TODO TODO compare number of isPrecedingEvent and hasLeft/hasRight defined? is the
# latter what i'd want to make my own "produces" relationships from (chem->chem), or
# something else? is reactions.dat what adds left/right? any directionality in sbml
# alone?
# TODO TODO count # of chemicals / reactions (and connected / unconnected) after each
# step of load process (or w/ loading increasing sets of files, starting w/ just sbml)
# TODO parse sbml another way and compare to # chemicals / reactions in there

# TODO TODO TODO change parser to use inchikey / whatever else might be more reliably in
# sbml


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
    #RETURN length(p) AS path_length

    # TODO something wrong with this query / the database? why am i not getting any
    # paths?
    result = tx.run(query, start_reactions=start_reactions, end_reactions=end_reactions)
    records = list(result)

    # TODO delete?
    #summary = result.consume()

    min_path_length = min([x['path_length'] for x in records])
    min_length_paths = [x['path'] for x in records]
    return min_path_length, min_length_paths


# TODO support both when reaction metadata is in a relationship (e.g. new canProduce)
# or in original reaction nodes
#def print_reaction_path()


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


def get_blacklist_compound_metaids(db):

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
    return blacklist_compound_metaids


def format_compound(compound):

    # This dict-like access works if compound is a Node type (from a Compound node)
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


def main():
    # URI examples: "neo4j://localhost", "neo4j+s://xxx.databases.neo4j.io"
    URI = 'neo4j://atlas'
    # TODO don't hardcode
    AUTH = ('neo4j', 'metabolike')

    # TODO need to add port (7687) to end of URI?
    db = Neo4jClient(URI, *AUTH)

    cm = CompoundMap(db)

    id_table = cm.id_table
    cpds = cm.cpds

    def biocyc2metacyc_ids(biocyc):
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
        meta_id_mask = id_table['metaId'] == metacyc_id
        assert meta_id_mask.sum() == 1
        row = id_table.loc[meta_id_mask.idxmax()]
        assert row['metaId'] == metacyc_id
        return row['compartment']

    my_hdf = pd.read_csv('hallem_inchis.csv')
    hdf = pd.read_excel('hallem110_pubchempy.xlsx')

    assert (
        hdf.biocyc_pubchem.apply(lambda x: x[len('PubChem:'):]).astype('int64') ==
        hdf.cid
    ).all()

    assert ('InChI=' + my_hdf.inchi == hdf.inchi).all()

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
    def get_metacyc_ids(chemical):
        print(chemical.name)

        cid = str(chemical.cid)
        if cid in db_pubchem_ids:
            # didn't want to risk matching against small integer CIDs, and having
            # partial / exact matches w/ the same integer used in a diff context
            # (could have probably just never used partial matches, but might still have
            # been some ambiguous cases...)
            sdf = pubchem_df[pubchem_df.value == cid]

            #print(sdf.to_string(index=False))

            # would need to deal with otherwise...
            assert len(sdf) == 1
            biocyc = sdf.biocyc.iat[0]

            meta_ids = biocyc2metacyc_ids(biocyc)

            # TODO check whether there are ever any distances actually defined across a
            # compartment boundary
            # TODO maybe augment database to allow searching across these boundaries
            # though? why do the chemicals have the compartment in their id...?
            print()

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

            print(f'{chem_id=}')
            print(f'chemical.{chem_id}={getattr(chemical, chem_id)}')
            print(f'{chem_id} matches:')
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

        print('not found in db')
        print()
        return None


    # TODO try to replace all relationship creation->blacklisting->deletion stuff with a
    # graph projection? would memory be an issue?
    # TODO TODO factor all blacklist creation / canProduce deletion/creation into
    # separate script, maybe one that does other database initialization stuff too

    # TODO blacklist all / most deacetylase reactions?

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
        # sorted(...) mainly just to convert to a list, as set not supported
        metaIdBlacklist=sorted(blacklist_compound_metaids)
    )


    name2meta_ids = dict()
    for row in hdf.itertuples():
        # TODO delete
        #if row.name not in ('E2-hexenal', '1-hexanol', 'ethanol', 'acetaldehyde'):
        #    continue

        # TODO delete
        #if row.name == 'ammonium hydroxide':
        #    continue

        meta_ids = get_metacyc_ids(row)
        name2meta_ids[row.name] = meta_ids

    print()
    print(f'{len(cids_not_in_db)} CIDs not in db')
    print(f'{len(not_found_in_db)} not found in db')
    print()

    for name1, meta_ids1 in name2meta_ids.items():
        for name2, meta_ids2 in name2meta_ids.items():
            if name1 == name2:
                continue

            print(f'{name1=}')
            print(f'{name2=}')

            if meta_ids1 is None or meta_ids2 is None:
                print('no metacyc IDs for one or more of these compounds')
                print()

                continue

            # TODO TODO only store if path is min length across all name1, name2 combos

            pathway_only_min_path_len = None
            pathway_only_min_path = None

            min_path_len = None
            min_path = None

            # TODO TODO TODO simplify to zipping over compartments, if there are no
            # cross-compartment reactions (tho may want to add my own cross-compartment
            # support later, if this is the case)
            for c1, c2 in product(meta_ids1, meta_ids2):
                print(f'getting reaction route between {c1} and {c2}')

                # To check if reactions are always compartment->compartment, b/c if so
                # we could simplify some of this code in the meantime, but we may also
                # want to modify the graph to ignore which compartment a reaction was
                # recorded in?
                compartment1 = metacyc_id2compartment(c1)
                compartment2 = metacyc_id2compartment(c2)
                across_compartments = compartment1 != compartment2

                before = time.time()

                # TODO TODO TODO replace w/ my own query (using ~same query str as
                # inside the only_pathway_reactions=True path here, maybe trying to only
                # use isPrecedingEvent in relationship filter?), and use read_tx
                # (to get output types same, for printing relationship info along path)
                #
                # TODO TODO by default, also pass blacklist of node ids into here
                # (but maybe try w/o, as pathway stuff more curated anyway, and
                # certainly is a subset of the reactions)
                pathway_routes = get_reaction_route_between_compounds(db, c1, c2,
                    only_pathway_reactions=True
                )
                query_s = time.time() - before
                #print(f'pathway-only shortest-path query took {query_s:.2f}s')

                # TODO delete
                if len(pathway_routes) > 0:
                    print('had pathway routes')
                    import ipdb; ipdb.set_trace()
                #

                # TODO convert this check into a query checking all (pathway) reactions
                # are only between compounds w/ metaIds referring to diff compartments
                if across_compartments:
                    assert len(pathway_routes) == 0

                before = time.time()

                # TODO factor to top-level. no need to define in loop.
                query = (
                    """
                    MATCH (c1:Compound {metaId: $c1}), (c2:Compound {metaId: $c2}),
                        p = allShortestPaths((c1)-[:canProduce*]->(c2))

                    RETURN relationships(p) as rs, p as path, length(p) as path_length
                    """
                )
                # The .data() call in db.read() was not giving me the
                # relationships(path) with the information I needed.
                def tx_fn(tx, **kwargs):
                    result = tx.run(query, **kwargs)
                    records = list(result)
                    return records

                routes = db.read_tx(tx_fn, c1=c1, c2=c2)

                query_s = time.time() - before
                #print(f'shortest-path query took {query_s:.2f}s')

                # assertion not true
                #if across_compartments:
                #    assert len(routes) == 0

                if len(routes) > 0:
                    assert len(set([x['path_length'] for x in routes])) == 1

                    print(f'{len(routes)=}')
                    for route in routes:
                        reactions = route['rs']
                        compounds = route['path'].nodes
                        assert len(reactions) + 1 == len(compounds)
                        # TODO would need to relax if input was from search where
                        # reactions are still represented as nodes...
                        # detect and only assert in case of my input? flag? broaden?
                        assert all([
                            x.labels == frozenset({'Compound'}) for x in compounds
                        ])

                        # TODO factor to a route-printing fn (and get routes from
                        # pathway only call into a suitable format for also passing to
                        # this)
                        for i, reaction in enumerate(reactions):
                            rxn_c1 = compounds[i]
                            rxn_c2 = compounds[i + 1]
                            print(f'{format_compound(rxn_c1)} -> '
                                f'{format_compound(rxn_c2)}'
                            )

                            # TODO also generate + print clickable links to metacyc
                            # (to check quickly)?
                            print(f'canonicalId: {reaction["canonicalId"]}')
                            print(f'name: {reaction["name"]}')

                        print()

                # took too long sometimes (at least on atlas w/ ~15Gb free memory and
                # w/ db populated using tier1-tier2/metacyc folder). unclear if it would
                # terminate.
                '''
                before = time.time()
                print('only_pathway_reactions=False')
                # TODO TODO does num_routes actually matter? bfs != shortest path(s)
                # here?
                # TODO max_hops (default=10)?
                # TODO TODO TODO pass blacklist as ignore_node_metaids
                routes = get_reaction_route_between_compounds(db, c1, c2,
                    only_pathway_reactions=False
                )
                query_s = time.time() - before
                print(f'took {query_s:.2f}s')
                print()
                '''

            # TODO TODO TODO actually reduce to single path / path-length output, across
            # calls

            # TODO compare only_pathway_reactions=True/False

    import ipdb; ipdb.set_trace()

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


if __name__ == '__main__':
    main()

