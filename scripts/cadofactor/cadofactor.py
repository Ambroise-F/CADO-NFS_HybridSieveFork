#!/usr/bin/env python3
import os
import sqlite3
import argparse
import logging
import cadologger
import wudb
import cadotask
import cadoparams

# FIXME: Is this a Task object? Probably not
# Should this be in cadotask or in cadofactor?
class CompleteFactorization(object):
    """ The complete factorization, aggregate of the individual tasks """
    def __init__ (self, db, *args, **kwargs):
        self.polysel = cadotask.PolyselTask(*args, db=db, **kwargs)
        self.fb = cadotask.FactorBaseTask(self.polysel, *args, db=db, **kwargs)
        self.freerel = cadotask.FreeRelTask(self.polysel, *args, db=db,
                                            **kwargs)
        self.sieving = cadotask.SievingTask(self.polysel, self.fb, *args, 
                                            db=db, **kwargs)
        self.dup1 = cadotask.Duplicates1Task(self.sieving, *args, db=db, 
                                             **kwargs)
        self.dup2 = cadotask.Duplicates2Task(self.dup1, *args, db=db, 
                                             **kwargs)
        self.purge = cadotask.PurgeTask(self.polysel, self.freerel, self.dup2,
                                       *args, db=db, **kwargs)
        self.merge = cadotask.MergeTask(self.purge, *args, db=db, **kwargs)
        self.linalg = cadotask.LinAlgTask(self.merge, *args, db=db, **kwargs)
        self.characters = cadotask.CharactersTask(self.polysel, self.purge, self.merge,
                                                  self.linalg, *args, db=db, **kwargs)
        self.sqrt = cadotask.SqrtTask(self.polysel, self.freerel, self.purge,
                                      self.merge, self.linalg, self.characters,
                                      *args, db=db, **kwargs)
    
    def run(self, *args, **kwargs):
        self.sqrt.run(*args, **kwargs)

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Integer Factorisation with the Number Field Sieve')
    parser.add_argument("--screenlog", help="Screen logging level, e.g., INFO/COMMAND/DEBUG", default="INFO", metavar="LEVEL")
    parser.add_argument("parameters", help="A file with the parameters to use")
    args = parser.parse_args()
    paramfile = args.parameters
    screenlvlname = args.screenlog
    
    parameters = cadoparams.Parameters()
    parameters._readfile(open(paramfile))
    tasksparams = parameters.myparams(("workdir", "name"), "tasks")
    
    screenlvl = cadologger.Logger.translate_level(screenlvlname)
    logger = cadologger.Logger()
    logger.addHandler(cadologger.ScreenHandler(lvl = screenlvl))
    cmdfilename = tasksparams["workdir"] + os.sep + tasksparams["name"] + ".cmd"
    logger.addHandler(cadologger.CmdFileHandler(cmdfilename))
    
    logger.debug("Root parameter dictionary:\n%s" % parameters)
    
     # logger.addHandler(cadologger.FileHandler(filename = "log", lvl = logging.DEBUG))
    wudb_file = tasksparams["workdir"] + os.sep + tasksparams["name"] + ".db"
    logger.info ('Opening database file "%s"', wudb_file)
    # dbconn = sqlite3.connect("wudb")
    logger.info ("Beginning factorization")
    factorjob = CompleteFactorization(wudb_file, parameters = parameters)
    factorjob.run()
