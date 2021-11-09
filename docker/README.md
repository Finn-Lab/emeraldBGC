# emeraldBGC container

#### Get InterProsScan data:
##### The size of download file is ~ 24G, the final directory is 12G. Be sure to have enough space
```bash
$ bash ./get_ips_slim.sh
```

#### Docker ready to use script:
##### Only works if "data/" and emeraldbgc_container.py are in the same directory
```bash
$ emeraldbgc_container.py --help
$ emeraldbgc_container.py [OPTIONS] ARGUMENTS
```

#### Docker image shell:
```bash
$ docker -it --entrypoint bash -v <path to emeraldBGC/docker>/data/:/opt/interproscan docker://santiagosanchezf/emeraldbgc:ips_nodata
$ emeraldbgc --help
$ emeraldbgc [OPTIONS] ARGUMENTS
```

#### Singularity image shell:
##### Get a copy of emeraldbgc container:
```bash
$ singularity build emeraldbgc.sif docker://santiagosanchezf/emeraldbgc:ips_nodata
```
##### Start container shell
```bash
$ singularity exec --bind $PWD:/home:rw --bind <path to emeraldBGC/docker>/data:/opt/interproscan/data emeraldbgc.sif 
```
##### Inside container shell
```bash
$ source ~/.bashrc
$ conda activate bgc
$ export PATH=/opt/interproscan/:$PATH
$ emeraldbgc --help
$ emeraldbgc [OPTIONS] ARGUMENTS
```
