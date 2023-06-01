# fast5_exporter

## Description

Metrics exporter from set of FAST5 files. Exported metrics include: number of reads, average read length, number of reads per pore, etc.

## Getting Started

### Dependencies

* Golang
* Python3

### Compiling

To build the exporter you should clone the repository, enter to specified dir and compile the project:
```
git clone https://github.com/wkusmirek/fast5_exporter.git
cd fast5_exporter
make
```

### Executing program

To run the exporter on the default port and default path to FAST5 dir ('/tmp/fast5'), run compiled binary app:
```
cd fast5_exporter
./fast5_exporter
```
To check if exporter is running properly, please curl the port in the another terminal:
```
curl localhost:9307/metrics
```
### Executing program via docker

To run the exporter without building the binaries, you can use docker image:
```
docker pull wkusmirek/fast5-exporter:latest
```

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the Apache License 2.0 License - see the LICENSE file for details
