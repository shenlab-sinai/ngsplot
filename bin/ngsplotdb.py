#!/usr/bin/env python
#
# ngsplotdb - PYTHON script to manipulate the ngs.plot database installation.
# Author: Li Shen, Mount Sinai School of Medicine
# Date created: Dec 28, 2012
# Last updated: May 28, 2013
#
# Functions implemented:
# 
# list - List locally installed genomes.
# check - Check the integrity of local database.
# install - Install a genome from ngsplotdb.XXX.tar.gz.
# remove - Remove a locally installed genome.


# def dir_to_idspec(dirn):
#     """Convert an FTP directory name to species ID and Name.
#        This works with both UCSC and Ensembl FTP.
#     """
#     parts = dirn.split('_')
#     sid = parts[0][:3].lower() + parts[1][:3].lower()
#     spec = str.title(parts[0] + ' ' + parts[1])
#     return sid, spec

# def listremote(database="refseq"):
#     """List available genomes on a remote FTP site."""

#     import urllib2

#     url_ensembl = "ftp://ftp.ensembl.org/pub/current_gtf/"
#     url_ucsc = "ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/"

#     if database == "refseq":
#         ftp = urllib2.urlopen(url_ucsc)
#     elif database == "ensembl":
#         ftp = urllib2.urlopen(url_ensembl)
#     else:
#         pass

#     # Read directory names from the remote FTP site.
#     dir_names = ftp.read().splitlines()
#     # Parse lines to obtain clean directory names only.
#     if database == "refseq":
#         dir_names = map(lambda x: x.split()[-3], dir_names)
#     elif database == "ensembl":
#         dir_names = map(lambda x: x.split()[-1], dir_names)
#     else:
#         pass
#     # Filter directory names that do not correspond to a species.
#     # On UCSC FTP, both "." and ".." will be retrieved from above and can 
#     # cause troubles.

#     id_spec_list = map(dir_to_idspec, dir_names)
#     print "ID\tName"
#     for r in id_spec_list:
#         print r[0], "\t", r[1]

def read_gnlist(root_path, mode):
    """Read installed genomes list.

       Args:
         root_path: ngs.plot installation directory.
         mode: row reading mode - "vector" or "hash", correspond to each record
               read as a vector of hash table.
       Returns: (header split vector, hash table of installed genomes, 
                 vector of column widths)
    """

    gn_list = root_path + "/database/gn_list.txt"
    try:
        db_f = open(gn_list, "r")
    except IOError:
        print "Open {0} error: your ngs.plot database may be corrupted.".\
              format(gn_list)
        sys.exit()

    header = db_f.readline()
    h_sp = header.rstrip().split("\t")
    v_cw = map(len, h_sp)  # column widths initialize to header widths.

    g_tbl = {}
    for rec in db_f:
        r_sp = rec.rstrip().split("\t")
        r_cw = map(len, r_sp)
        v_cw = map(max, v_cw, r_cw)

        r_tbl = dict(zip(h_sp, r_sp))
        if(mode == "vector"):
            g_tbl[r_tbl["ID"]] = r_sp
        elif(mode == "hash"):
            g_tbl[r_tbl["ID"]] = r_tbl
        else:
            pass

    db_f.close()

    return (h_sp, g_tbl, v_cw)


def update_gnlist(root_path, g_tbl, gn, assembly, species, ens_v, np_v):
    """Update/Add genomes in ngs.plot database.

       Args:
         root_path: ngs.plot installation directory.
         g_tbl: genome hash table. Each value is also a hash table of 
                genome info.
         gn: genome name.
         assembly: genome assembly name.
         species: species name.
         ens_v: new Ensembl version.
         np_v: new ngs.plot version.
       Returns: None
    """

    if gn in g_tbl:
        g_tbl[gn]["EnsVer"] = ens_v
        g_tbl[gn]["NPVer"] = np_v
    else:
        g_tbl[gn] = {"ID": gn, "Assembly": assembly, "Species": species, \
                     "EnsVer": ens_v, "NPVer": np_v}

    write_gnlist(root_path, g_tbl)


def write_gnlist(root_path, g_tbl):
    """Write genome hash table to database file: gn_list.txt.

       Args:
         root_path: ngs.plot installation directory.
         g_tbl: genome hash table. Each value is also a hash table of 
                genome info.
       Returns: None
    """

    output_order = ["ID", "Assembly", "Species", "EnsVer", "NPVer"]
    gn_list = root_path + "/database/gn_list.txt"

    try:
        db_f = open(gn_list, "w")
    except IOError:
        print "Open {0} error: your ngs.plot database may be corrupted.".\
              format(gn_list)
        sys.exit()

    db_f.write("\t".join(output_order) + "\n")

    for k in g_tbl.viewkeys():
        rec = []
        for col in output_order:
            rec.append(str(g_tbl[k][col]))
        db_f.write("\t".join(rec) + "\n")
    
    db_f.close()


def listgn(root_path, args):
    """List installed genomes in local database on screen.

       Args:
         root_path: ngs.plot installation directory.
         args: currently not used.
       Returns: None.
    """

    import math

    (h_sp, g_tbl, v_cw) = read_gnlist(root_path, "vector")

    # Format column widths to beautify output.
    tab_u = 4
    v_cw = map(lambda x: int(math.ceil(float(x) / tab_u) * tab_u + 1), v_cw)

    # List genomes to screen.
    print "".join(map(lambda x, y: x.ljust(y), h_sp, v_cw))  # header.
    for k in sorted(g_tbl.viewkeys(), key=str.lower):
        print "".join(map(lambda x, y: x.ljust(y), g_tbl[k], v_cw))


def install(root_path, args):
    """Interactive session for installing genome from package file.

       Args:
         root_path: ngs.plot installation directory.
         args.pkg: path of package file.
       Returns:
         None.
    """

    import tarfile
    import sys

    pkg_file = args.pkg
    # Package name: e.g. ngsplotdb_mm10_71_1.0.tar.gz.
    # Information about the genome is in mm10/.metainfo.
    try:
        pkg_f = tarfile.open(pkg_file, "r:gz")
    except tarfile.ReadError:
        print "Read package file {0} error.".format(pkg_file),
        print "The downloaded file may be corrupted."
        sys.exit()

    print "Extracting information from package...",
    sys.stdout.flush()
    pkg_files = pkg_f.getnames()
    g_folder = pkg_files[0]
    # Minus folder name and .metainfo file name.
    print "contains " + str(len(pkg_files) - 2) + " tables."
    # .metainfo file:
    # "ID"<TAB>mm10
    # "Assembly"<TAB>GRCm38
    # "Species"<TAB>mus_musculus
    # "EnsVer"<TAB>71
    # "NPVer"<TAB>1.0
    meta_f = pkg_f.extractfile(g_folder + "/.metainfo")
    meta_tbl = {}
    if meta_f:
        for rec in meta_f:
            rec_sp = rec.strip().split("\t")
            meta_tbl[rec_sp[0]] = rec_sp[1]  # format: Variable<TAB>Value
    else:
        print "No meta information found in the package file.",
        print "The downloaded file may be corrupted."
        sys.exit()
    meta_f.close()
    pkg_f.close()

    gn_inst = meta_tbl["ID"]
    assembly = meta_tbl["Assembly"]
    species = meta_tbl["Species"]
    ens_ver = float(meta_tbl["EnsVer"])
    np_ver = float(meta_tbl["NPVer"])

    # Database file.
    (h_sp, g_tbl, v_cw) = read_gnlist(root_path, "hash")

    if gn_inst in g_tbl:
        installed_ens = float(g_tbl[gn_inst]["EnsVer"])
        installed_np = float(g_tbl[gn_inst]["NPVer"])
        # Install a newer version.
        if ens_ver > installed_ens or \
           ens_ver == installed_ens and np_ver > installed_np:
            print "Will upgrade the genome {0}:".format(gn_inst)
            print "Ensembl: {0}==>{1}; ngs.plot: {2}==>{3}.".\
                  format(installed_ens, ens_ver, installed_np, np_ver),
            ans = raw_input("Continue?(y/n): ")
            while True:
                if(ans == "y" or ans == "Y" or ans == "n" or ans == "N"):
                    break
                else:
                    ans = raw_input("Answer must be y/Y or n/N: ")
            if(ans == "y" or ans == "Y"):
                install_pkg(root_path, pkg_file, gn_inst, gn_inst)
                update_gnlist(root_path, g_tbl, gn_inst, assembly, \
                              species, ens_ver, np_ver)
        # Install an older version.
        else:
            print "Will install the same or older version",
            print "of the genome {0}:".format(gn_inst)
            print "Ensembl: {0}==>{1}; ngs.plot: {2}==>{3}.".\
                  format(installed_ens, ens_ver, installed_np, np_ver),
            ans = raw_input("Continue?(y/n): ")
            while True:
                if(ans == "y" or ans == "Y" or ans == "n" or ans == "N"):
                    break
                else:
                    ans = raw_input("Answer must be y/Y or n/N: ")
            if(ans == "y" or ans == "Y"):
                install_pkg(root_path, pkg_file, gn_inst, gn_inst)
                update_gnlist(root_path, g_tbl, gn_inst, assembly, \
                              species, ens_ver, np_ver)
    else:
        print "Will install new genome {0}:".format(gn_inst),
        print "Ensembl=> v{0}; ngs.plot=> v{1}.".format(ens_ver, np_ver),
        ans = raw_input("Continue?(y/n): ")
                        
        while True:
            if(ans == "y" or ans == "Y" or ans == "n" or ans == "N"):
                break
            else:
                ans = raw_input("Answer must be y/Y or n/N:")
        if(ans == "y" or ans == "Y"):
            install_pkg(root_path, pkg_file, gn_inst)
            update_gnlist(root_path, g_tbl, gn_inst, assembly, \
                          species, ens_ver, np_ver)


def install_pkg(root_path, pkg_file, new_gn, old_gn=None):
    """Install genome from package file.

       Args:
         root_path: ngs.plot installation directory.
         pkg_file: package file name.
         old_gn: old genome name to be removed first.
       Returns:
         None.
    """

    import shutil
    import tarfile
    import sys

    if old_gn:
        print "Removing installed genome:{0}...".format(old_gn),
        sys.stdout.flush()
        shutil.rmtree(root_path + '/database/', old_gn)
        rm_dbtbl(root_path, old_gn)
        print "Done"

    print "Installing new genome...",
    sys.stdout.flush()
    try:
        tar_f = tarfile.open(pkg_file, "r:gz")
    except tarfile.ReadError:
        print "Read package file {0} error: ".format(pkg_file),
        print "downloaded file may be corrupted."
        sys.exit()

    try:
        tar_f.extractall(root_path + "/database")
    except tarfile.ExtractError:
        print "Extract files from package error.", 
        print "The downloaded file may be corrupted."

    add_dbtbl(root_path, new_gn)

    print "Done"


def add_dbtbl(root_path, gn):
    """Add a new genome into database meta tables.

       Args:
         root_path: ngs.plot installation directory.
         gn: genome name to add.
       Returns: None.
    """

    import os.path
    import os
    import glob
    import tempfile
    import subprocess

    (dblist_f, dblist_n) = tempfile.mkstemp(text=True)
    (gnlist_f, gnlist_n) = tempfile.mkstemp(text=True)

    # Obtain a list of .RData file names to extract info from.
    all_rdata = root_path + '/database/{0}/{1}*.RData'.format(gn, gn)
    all_set = set(map(os.path.basename, glob.glob(all_rdata)))
    exm_rdata = root_path + '/database/{0}/{1}*exonmodel*.RData'.format(gn, gn)
    exm_set = set(map(os.path.basename, glob.glob(exm_rdata)))
    nex_list = sorted(list(all_set - exm_set))

    # DB file list.
    gn_ens_list = {}  # Ensembl genome name unique list.
    desired_ncols = 6
    for fname in nex_list:
        tokens = fname.split('.')
        tokens = tokens[:6]
        if tokens[-1] == 'RData':
            tokens[-1] = 'NA'
        tokens.extend(['NA'] * (desired_ncols - len(tokens)))
        os.write(dblist_f, fname + "\t" + "\t".join(tokens) + "\n")
        if tokens[1] == 'ensembl':
            gn_ens_list[".".join(tokens[:3])] = 1
    os.close(dblist_f)

    # Ensembl genome name list.
    for gn in sorted(gn_ens_list.viewkeys()):
        os.write(gnlist_f, gn.replace(".", "\t") + "\n")
    os.close(gnlist_f)    

    # Call external program to create table default file.
    (default_f, default_n) = tempfile.mkstemp(text=True)
    subprocess.call(["setTableDefaults.py", gnlist_n, default_n])
    os.remove(gnlist_n)

    # Call external program to calcualte DB score and install the new entries 
    # into ngs.plot database.
    subprocess.call(["install.db.tables.r", default_n, dblist_n])
    os.remove(default_n)
    os.remove(dblist_n)


def remove(root_path, args):
    """Remove genome from ngs.plot database.

       Args:
         root_path: ngs.plot installation directory.
         args.gn: genome name.
       Returns:
         None.
    """

    import shutil
    import sys

    (h_sp, g_tbl, v_cw) = read_gnlist(root_path, "hash")
    gn = args.gn

    if gn in g_tbl:
        print "Will remove genome {0} from database.".format(gn),
        ans = raw_input("Continue?(y/n): ")
        while True:
            if ans == 'y' or ans == 'Y' or ans == 'n' or ans == 'N':
                break
            else:
                ans = raw_input("The answer must be y/Y or n/N: ")
        if ans == 'y' or ans == 'Y':
            folder_to_rm = root_path + "/database/" + gn
            print "Removing genome...",
            sys.stdout.flush()
            shutil.rmtree(folder_to_rm)
            del g_tbl[gn]
            write_gnlist(root_path, g_tbl)
            rm_dbtbl(root_path, gn)
            print "Done"
    else:
        print "Cannot find the genome in database. Nothing was done."


def rm_dbtbl(root_path, gn):
    """Add a new genome into database meta tables.

       Args:
         root_path: ngs.plot installation directory.
         gn: genome name to remove.
       Returns: None.
    """

    import subprocess

    subprocess.call(["remove.db.tables.r", gn])




###############################################
######## Main procedure. ######################
if __name__ == "__main__":

    import argparse
    import os
    import sys

    try:
        root_path = os.environ["NGSPLOT"]
    except KeyError:
        print "Environmental variable NGSPLOT is not set! Exit.\n"
        sys.exit()

    # Main parser.
    parser = argparse.ArgumentParser(description="Manipulate ngs.plot's \
                                                  annotation database",
                                     prog="ngsplotdb.py")
    subparsers = parser.add_subparsers(title="Subcommands",
                                       help="additional help")

    # ngsplotdb.py list parser.
    parser_list = subparsers.add_parser("list", help="List genomes installed \
                                                      in database")
    parser_list.set_defaults(func=listgn)

    # ngsplotdb.py install parser.
    parser_install = subparsers.add_parser("install", 
                                           help="Install genome from package \
                                                 file")
    parser_install.add_argument("pkg", help="Package file(.tar.gz) to install",
                                type=str)
    parser_install.set_defaults(func=install)

    # ngsplotdb.py remove parser.
    parser_remove = subparsers.add_parser("remove", 
                                          help="Remove genome from database")
    parser_remove.add_argument("gn", help="Name of genome to be \
                                           removed(e.g. hg18)", type=str)
    parser_remove.set_defaults(func=remove)


    args = parser.parse_args()
    args.func(root_path, args)

