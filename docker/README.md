# emeraldBGC container

#### Get InterProsScan data:
##### The size of download file is ~ 24G, the final directory is 16G. Be sure to have enough space
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
$ docker -it --entrypoint bash -v <path to emeraldBGC/docker>/data/:/opt/interproscan quay.io/repository/microbiome-informatics/emerald-bgc
$ emeraldbgc --help
$ emeraldbgc [OPTIONS] ARGUMENTS
```
