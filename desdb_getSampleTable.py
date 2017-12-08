""" Simple script to get sample tables from the DESDB. It uses the user 
.desservices.ini file to get the connection information
"""

import os
import sys
import time
import socket
import argparse
import uuid
import logging
import easyaccess as ea
if False:
    import subprocess
    import shlex


class Query():
    def ea_query(self, db="dessci", tab=None, sample=None, 
                 cond1=None, cond2=None, outfile=None, ext=None):
        """ Easyaccess simple query to select columns (all or a sample), 
        applying some conditions plus the sampling parameter 
        Inputs
        - db: DB to connect to
        - tab: table name
        - sample: percentage to be used, in range [0.000001,100)
        - cond1: columns to be used, statement after SELECT ...
        - cond2: condition to the selection, statement after WHERE ...
        - outfile: 
        - ext: extension, '.csv', '.tab', '.h5', or '.fits'
        """
        # Naming
        if (outfile is None):
            auxid = str(uuid.uuid4())
            outfile = "{0}_{1}pcent_{2}{3}".format(tab, sample, auxid, ext)
        # review the below
        if (cond1 is None): 
            if (cond2 is None):
                q = "select * from {0} sample({1});".format(tab, sample)
                q += " > {0}".format(outfile)
            else:
                q = "select * from {0} where {1}".format(tab, cond2)
                q += " sample({0}); > {1}".format(sample, outfile)
        else:
            if (cond2 is None):
                q = "select {0} from {1}".format(cond1, tab)
                q += " sample({0}); > {1}".format(sample, outfile)
            else:
                q = "select {0} from {1} where {2}".format(cond1, tab, cond2)
                q += " sample({0}); > {1}".format(sample, outfile)
        connect =  ea.connect(db)
        cursor = connect.cursor()
        logging.info("Running query\n > {0}".format(q))
        cursor.execute(q)
        connect.close()
        logging.info("Output file: {0}".format(outfile))
        return True

if __name__ == "__main__":
    print("ERROR WHEN TRYING TO RUN THE CODE")
    print(" ORA-00933: SQL command not properly ended")

    descr = "Simple script to call easyaccess to get a sample of some DES"
    descr += " database table, able to add simples conditions"
    ap = argparse.ArgumentParser(description=descr)
    h0 = "Name of the table in which operate, as in \'FROM TABLE_NAME\'"
    ap.add_argument("table", help=h0, type=str)
    h1 = "Percentaje of the table to be sampled, in range [0.000001,100)"
    ap.add_argument("percent", help=h1, type=float)
    h2 = "Condition to be used to select columns, as in \'SELECT CONDITION\'."
    h2 += " If no value is provided, all columns are retrieved"
    ap.add_argument("--c1", help=h2, metavar="", type=str)
    h3 = "Condition to be used to constraint the selected columns, as in"
    h3 += " \'WHERE CONDITION\'. If no value is provided, no constraint is"
    h3 += " applied"
    ap.add_argument("--c2", help=h3, metavar="", type=str)
    h4 = "Filename for the output table. If no value is inputed, then"
    h4 += " {TABLE}_{SAMPLE}pcent_{RANDOM_STR}{.EXTENSION}"
    ap.add_argument("--out", help=h4, metavar="", type=str)
    h5 = "Extension to be used on the output file. Independent from the"
    h5 += " output filename. Default: \'.h5\'"
    ap.add_argument("--ext", help="h5", metavar="", type=str, default=".h5",
                    choices=[".csv", ".tab", ".h5", ".fits"])
    h6 = "Which database to use. Default is \'dessci\'"
    ap.add_argument("--db", help=h6, metavar="", type=str, default="dessci",
                    choices=["desoper", "destest", "dessci"])
    h7 = "In case an output file is needed, input the name to be used for it"
    ap.add_argument("--log", help=h7, metavar="", type=str)
    # Parse Args
    Nms = ap.parse_args()
    d = dict()
    d["db"] = Nms.db
    d["tab"] = Nms.table
    d["sample"] = Nms.percent
    d["cond1"] = Nms.c1 
    d["cond2"] = Nms.c2
    d["outfile"] = Nms.out
    d["ext"] = Nms.ext
    # Setup the logging
    if (Nms.log is None):
        logging.basicConfig( 
            format= "%(asctime)s - %(name)s - %(levelname)s - %(message)s", 
            level=logging.DEBUG
            )
    else:
        logging.basicConfig( 
            filename=Nms.log,
            format= "%(asctime)s - %(name)s - %(levelname)s - %(message)s", 
            level=logging.DEBUG
            )
    logging.info("Hostname: {0}".format(socket.gethostname()))
    # Run the query
    Q =  Query()
    Q.ea_query(**d)

