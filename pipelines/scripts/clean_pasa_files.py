import sys
import os

os.system("rm -r rnd1/__pasa_pasa.sqlite_SQLite_chkpts")
os.system("rm -r rnd1/__pasa_pasa.sqlite_SQLite_chkpts.cmds_log")
os.system("rm -r rnd1/pasa_run.log.dir")

os.system("rm -r rnd2/__pasa_pasa.sqlite_SQLite_chkpts")
os.system("rm -r rnd2/__pasa_pasa.sqlite_SQLite_chkpts.cmds_log")
os.system("rm -r rnd2/pasa_run.log.dir")

#os.system("rm -r rnd3/__pasa_pasa.sqlite_SQLite_chkpts")
#os.system("rm -r rnd3/__pasa_pasa.sqlite_SQLite_chkpts.cmds_log")
#os.system("rm -r rnd3/pasa_run.log.dir")


for file_path in os.listdir(sys.argv[1]):
    if file_path in ["alignAssembly.config", "annotationCompare.config", "rnd1", "rnd2", "rnd3", "run_alignAssembly.sh", "run_seqclean.sh"]:
        continue
    if "run_alignAssembly.sh" in file_path:
        continue
    if "run_seqclean.sh" in file_path:
        continue
    os.system("rm -r %s" % file_path)
