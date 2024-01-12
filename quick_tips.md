# Quick tips for interacting with the Alpine cluster

Full documentation [here](https://curc.readthedocs.io/en/latest/index.html).

## Login

Login at [this link](https://ondemand-rmacc.rc.colorado.edu). 

* Open a terminal on the head node by selecting *Clusters/Alpine*
* Start a VS Code interactive job *Interactive Apps/vS Code-Server (Presets)*

### ssh Login

Directions to set up ssh are [here](https://github.com/LRFreeborn/CURCDocumentation/blob/amc_login/docs/access/amc-access.md)

Command to log in using the created key

```bash
ssh -i ~/.ssh/id_rsa {user_name}@login-ci.rc.colorado.edu
```

## Helpful commands

Submit a script
```bash
sbatch job_script.sh
```

Cancel a job
```bash
scancel job_id
```

Start an interactive job
```bash
sinteractive --partition=amilan --time=00:10:00 --ntasks=1
```

To view all jobs

```bash
squeue
```

To view just your jobs

```bash
squeue --user {user_name}
```


## Transfer files
Set up [Globus](https://curc.readthedocs.io/en/latest/compute/data-transfer.html). Note, need to set up and use the `CU Boulder Research Computing ACCESS`.

[Link to transfer](https://app.globus.org/file-manager)

### Rsync


```bash
rsync -avz --progress path/to/local/files {user_name}@login-ci.rc.colorado.edu:/home/{user_name}/path/to/alpine/location
```

## Setting up Conda

To start conda

```bash
module load anaconda
conda init
```

By default conda will install in the home directory, but there is not enough room. Instead it needs to install into projects.

```bash
vim ~/.condarc
```

And add

```
pkgs_dirs:
  - /projects/$USER/.conda_pkgs
envs_dirs:
  - /projects/$USER/software/anaconda/envs
```

Now when you create a new environment, it will point to the correct location

Information for chainging the base directory of installing packages with conda is [here](https://curc.readthedocs.io/en/latest/software/python.html)

## Checking usage

```bash
module load slurmtools
```

See usage stats

```bash
suuser {user_name}
```

See usage stats for the last year

```bash
suuser {user_name} 365
```

This is important to check because
```
The cumulative computing allocations of a single research group across all projects may not exceed 5 M SU/year or 2.5 M SU over 6 months
```

Check on jobs that have run (including wait times)

```bash
jobstats {user_name}
```

Check priority (greater than 1 is higher than average, less than 1 is lower than average)

```bash
levelfs {user_name}
```

View resources requested by a job

```bash
scontrol show job {job_id}
```
