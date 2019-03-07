
def makebatchlist(command_list, basename, mem=16, time=50, nodes=1, queue="none", core=16, mode="sbatch", njobs_per_batch = 16, machine = "occigen", path2batch = "", path2run = ""):
    njobs = len(command_list)
    nbatch = njobs // njobs_per_batch;
    if njobs - nbatch*njobs_per_batch:
        nbatch = nbatch + 1

    for batch in range(nbatch):
        with open(path2batch + basename + "_{0}.sh".format(batch), 'w') as batchfile:
            if mode == "sbatch" :
                batchfile.write("#!/bin/bash\n")
                if queue != "none":
                    batchfile.write("#SBATCH --partition={0}\n".format(queue))
                batchfile.write("#SBATCH --time={0}:00:00\n".format(time))
                batchfile.write("#SBATCH --nodes={0}\n".format(nodes))
                batchfile.write("#SBATCH --ntasks-per-node={0}\n".format(core))
                batchfile.write("#SBATCH --mem={0}gb\n".format(mem))
                batchfile.write("#SBATCH -o {0}_{1}.out\n".format(basename, batch))
                batchfile.write("#SBATCH -e {0}_{1}.err\n".format(basename, batch))
                if machine == "p2chpd":
                    batchfile.write("module purge\n")
                    batchfile.write("module use /softs/modulefiles\n")
                    batchfile.write("module load gcc/4.7.3 mpi/openmpi/1.8.4/openmpi-gcc47\n")
                    # batchfile.write("export DIR={0}\n".format(path))
                if machine == "occigen":
                    batchfile.write("#SBATCH --constraint=HSW24\n")
                    batchfile.write("module purge\n")
                    batchfile.write("module load intel/18.1 gcc/6.2.0 openmpi/gnu/2.0.2\n")
                batchfile.write("\n")
            
            n = njobs_per_batch
            if batch == nbatch-1:
                n = njobs - (nbatch-1)*njobs_per_batch

            if n == 1:
                batchfile.write(path2run + command_list[batch*njobs_per_batch] + " \n")
            else:
                for i in range(n):
                    batchfile.write(path2run + command_list[batch*njobs_per_batch + i] + " &\n")
                batchfile.write("wait\n")

def makebatch(command, basename, mem=16, time=50, nodes=1, core=16, queue="none", mode="sbatch", machine = "occigen", path2batch = "", path2run = ""):
    with open(path2batch + basename + ".sh", 'w') as batchfile:
        if mode == "sbatch" :
            batchfile.write("#!/bin/bash\n")
            if queue != "none":
                batchfile.write("#SBATCH --partition={0}\n".format(queue))
            batchfile.write("#SBATCH --time={0}:00:00\n".format(time))
            batchfile.write("#SBATCH --nodes={0}\n".format(nodes))
            batchfile.write("#SBATCH --ntasks-per-node={0}\n".format(core))
            batchfile.write("#SBATCH --mem={0}gb\n".format(mem))
            batchfile.write("#SBATCH -o {0}.out\n".format(basename))
            batchfile.write("#SBATCH -e {0}.err\n".format(basename))
            if machine == "p2chpd":
                batchfile.write("module purge\n")
                batchfile.write("module use /softs/modulefiles\n")
                batchfile.write("module load gcc/4.7.3 mpi/openmpi/1.8.4/openmpi-gcc47\n")
                # batchfile.write("export DIR={0}\n".format(path))
            if machine == "occigen":
                batchfile.write("#SBATCH --constraint=HSW24\n")
                batchfile.write("module purge\n")
                batchfile.write("module load intel/18.1 gcc/6.2.0 openmpi/gnu/2.0.2\n")
            batchfile.write("\n")

        if machine == "p2chpd":
            batchfile.write("srun -n {0} -l --mpi=pmi2 --cpu_bind=Verbose {1}\n".format(nodes*core, path2run + command))
        if machine == "occigen":
            batchfile.write("srun --mpi=pmi2 -K1 --resv-ports -n {0} {1}\n".format(nodes*core, path2run + command))

if __name__ == "__main__":
    import sys
    makebatch(*sys.argv[1:])


