import ROOT as r
import os
import numpy as np
import uproot 

if __name__ == "__main__":

    dataDir = '/pnfs/icarus/scratch/users/micarrig/nueCCInclusiveSpine24/output/'

    completed_jobs = []
    completed_sys = []
    completed_var = []

    cutoff = 3749

    for filename in os.listdir(dataDir):
        if not filename.endswith('.root'):
            continue
        if 'systematics' in filename:
            id_ = filename.split('jobid')[-1].split('.root')[0]
            if int(id_) > cutoff:
                continue
            completed_sys.append(id_)
            continue
        elif 'varsys_jobid' in filename:
            id_ = filename.split('jobid')[-1].split('.root')[0]
            if int(id_) > cutoff:
                continue
            completed_var.append(id_)
            continue
        elif 'jobid' in filename:
            id_ = filename.split('jobid')[-1].split('.root')[0]
            if int(id_) > cutoff:
                continue
            completed_jobs.append(id_)
    
    completed_jobs.sort()
    completed_sys.sort()
    completed_var.sort()

    missing_sys = np.setdiff1d(completed_jobs, completed_sys)
    missing_var = np.setdiff1d(completed_jobs, completed_var)

    print("There are {} completed jobs".format(len(completed_jobs)))
    print("There are {} completed and {} missing systematics jobs".format(len(completed_sys), len(missing_sys)))
    print("There are {} completed and {} missing variations jobs".format(len(completed_var), len(missing_var)))

    # make list of var sys files to merge
    var_sys_files = [f'root://fndcadoor.fnal.gov:{dataDir.replace("/pnfs", "")}/output_varsys_jobid{id_}.root' for id_ in completed_sys]

    # make list of sys files to merge
    sys_files = [f'root://fndcadoor.fnal.gov:{dataDir.replace("/pnfs", "")}/output_systematics_jobid{id_}.root' for id_ in completed_sys]

    RED = "\033[91m"
    RESET = "\033[0m"

    mismatched_ids = set()
    total_var = 0
    total_sys = 0
    for id_, v, s in zip(completed_sys, var_sys_files, sys_files):
        try:
            sys_file = uproot.open(s)
            n_var = uproot.open(v)["events/NuMIFull/selected_variationTree"].num_entries
            n_sys = sys_file["events/NuMIFull/selected"].num_entries
            if n_var != n_sys:
                n_nonmatched = sys_file["events/NuMIFull/selected_nonmatched"].num_entries
                if n_sys + n_nonmatched != n_var:
                    print(f"{RED}{id_}: MISMATCH - var={n_var}, sys={n_sys}, nonmatched={n_nonmatched}{RESET}")
                    mismatched_ids.add(id_)
                    continue
                else:
                    print(f"{id_}: var={n_var}, sys={n_sys}, nonmatched={n_nonmatched} (sum matches)")
            else:
                print(f"{id_}: var={n_var}, sys={n_sys}")
            total_var += n_var
            total_sys += n_sys
        except Exception as e:
            print(f"{id_}: ERROR - {e}")
            mismatched_ids.add(id_)

    var_sys_files = [f for id_, f in zip(completed_sys, var_sys_files) if id_ not in mismatched_ids]
    sys_files     = [f for id_, f in zip(completed_sys, sys_files)     if id_ not in mismatched_ids]

    print(f"\nRemoved {len(mismatched_ids)} mismatched jobs. Writing {len(sys_files)} files.")

    with open('varsys_files.txt', 'w') as f:
        for file in var_sys_files:
            f.write(file + '\n')

    with open('sys_files.txt', 'w') as f:
        for file in sys_files:
            f.write(file + '\n')

    print(f"\nTotal entries:")
    print(f"  VarSys (selected_variationTree): {total_var}")
    print(f"  Sys    (selected):               {total_sys}")
