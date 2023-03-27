#!/usr/bin/env python3

# NOTE: didn't actually use this. was trying to adapt one of the metabolike scripts to
# load multiple input biocyc folders, but hadn't got it working / tested it

import logging
from pathlib import Path

import typer

from metabolike.config import load_config
from metabolike.db import MetacycClient
from metabolike.parser import MetacycParser

logging.basicConfig(
    format="[%(levelname)s] %(asctime)s - %(name)s:%(lineno)s:%(funcName)s - %(message)s",
    level=logging.INFO,
)
logger = logging.getLogger(__name__)


def main(
    config_file: Path = typer.Option(
        Path("metabolike-setup.yaml"), "--config", "-c", help="Path to the configuration file."
    ),
    biocyc_data_root: Path = typer.Option(
        Path("tier1-tier2"), "-d", help="Directory with biocyc data (containing *cyc/ directories)",
    ),
    create_db: bool = typer.Option(
        True, help="When database creation is not allowed, set this to False."
    ),
    #drop_if_exists: bool = typer.Option(
    #    False, "--drop-if-exists", "-f", help="Drop the database if it already exists."
    #),
):

    conf = load_config(config_file)

    logger.info("Connecting to neo4j database")
    db = MetacycClient(
        **conf.neo4j.dict(include={"uri", "database"}),
        neo4j_user=conf.neo4j.user,
        neo4j_password=conf.neo4j.password.get_secret_value(),
        create_db=create_db,

        drop_if_exists=False,
        # TODO test (+ fix if needed)

        # might cause problems on repeated calls...
        # (within one run of this script, across loop iterations)
        #drop_if_exists=drop_if_exists,
    )

    keys2filenames = {
        "sbml": ["metabolic-reactions.sbml", "metabolic-reactions.xml"],
        "compounds": "compounds.dat",

        # TODO need?
        "classes": "classes.dat",

        # TODO TODO ok to not use these?
        #
        # i think the code that merges nodes into "compound reaction" nodes has errors
        # (or at least isn't designed to be run more than once, but it could be
        # unrelated)
        #"pathways": "pathways.dat",

        # TODO need?
        "publications": "pubs.dat",

        # TODO TODO TODO need? test in yeast, that many don't lengths change!
        "reactions": "reactions.dat",

        # "atom_mapping" also supported in theory, but it led to some errors, and don"t
        # think it's needed
    }
    required_keys = [
        "sbml",
    ]
    # TODO define set of required (or optional) keys, if some input dirs are missing
    # some non-required

    # TODO particular order? if it does matter, probably indicates an error...
    for biocyc_dir in biocyc_data_root.glob("*cyc/"):

        logger.info(f"Adding BioCyc directory {biocyc_dir} to database")

        biocyc_files = dict()
        missing_required = False

        for key, filenames in keys2filenames.items():
            if isinstance(filenames, str):
                filenames = [filenames]

            found = False
            for filename in filenames:
                path = biocyc_dir / filename
                if path.exists():
                    biocyc_files[key] = path
                    found = True
                    break

            if not found:
                logger.warn(f"{biocyc_dir} missing {filenames}!")
                if key in required_keys:
                    missing_required = True
                    break

        if missing_required:
            logger.error(f"skipping all data in {biocyc_dir}!")
            continue

        meta = MetacycParser(**biocyc_files)

        logger.info(f"Adding to database using {type(meta)}")
        db.metacyc_to_graph(meta)

    db.close()


if __name__ == "__main__":
    typer.run(main)

