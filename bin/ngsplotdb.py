#!/usr/bin/env python
#
# ngsplotdb - PYTHON script to manipulate the ngs.plot database installation.
# Author: Li Shen, Mount Sinai School of Medicine
# Date created: Dec 28, 2012
# Last updated: Dec 28, 2012
#
# Functions implemented:
# 
# listlocal - List locally installed genomes.
# listremote - List available genomes on remote databases.
# check - Check the integrity of locally installed genomes.
# install - Install a genome from remote databases.
# remove - Remove a locally installed genome.
# update - Update all installed genomes with current releases.


def dir_to_idspec(dirn):
    """Convert an FTP directory name to species ID and Name.
       This works with both UCSC and Ensembl FTP.
    """
    parts = dirn.split('_')
    sid = parts[0][:3].lower() + parts[1][:3].lower()
    spec = str.title(parts[0] + ' ' + parts[1])
    return sid, spec

def listremote(database="refseq"):
    """List available genomes on a remote FTP site."""

    import urllib2

    url_ensembl = "ftp://ftp.ensembl.org/pub/current_gtf/"
    url_ucsc = "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/"

    if database == "refseq":
        ftp = urllib2.urlopen(url_ucsc)
    elif database == "ensembl":
        ftp = urllib2.urlopen(url_ensembl)
    else:
        pass

    # Read directory names from the remote FTP site.
    dir_names = ftp.read().splitlines()
    # Parse lines to obtain clean directory names only.
    if database == "refseq":
        dir_names = map(lambda x: x.split()[-3], dir_names)
    elif database == "ensembl":
        dir_names = map(lambda x: x.split()[-1], dir_names)
    else:
        pass
    # Filter directory names that do not correspond to a species.
    # On UCSC FTP, both "." and ".." will be retrieved from above and can 
    # cause troubles.

    id_spec_list = map(dir_to_idspec, dir_names)
    print "ID\tName"
    for r in id_spec_list:
        print r[0], "\t", r[1]


# Main procedure.
if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description="Manipulate ngs.plot's \
        annotation database")
    parser.add_argument("subcommand", help="Subcommand to execute.", type=str,
        choices=["listlocal", "listremote", "check", "install", "remove", 
                "update"])
    parser.add_argument("-db", "--database", type=str, choices=["refseq", 
        "ensembl"], default="refseq", help="Database name.")
    parser.add_argument("-gn", "--genome", type=str, help="Genome name \
                        concatenated by the first 3 letters of the two parts \
                        of the official name. For example, homo sapiens = \
                        homsap.")
    args = parser.parse_args()

    # Execute subcommand by calling functions.
    if args.subcommand == "listlocal":
        pass
    elif args.subcommand == "listremote":
        pass
    elif args.subcommand == "check":
        pass
    elif args.subcommand == "install":
        pass
    elif args.subcommand == "remove":
        pass
    elif args.subcommand == "update":
        pass
