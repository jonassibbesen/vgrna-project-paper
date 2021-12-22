# Installation

## vg

Note that vg can take anywhere from 10 minutes to more than an hour to compile depending on your machine and the number of threads used. 

### Building on Linux

If you don't want to or can't use a pre-built release of vg, or if you want to become a vg developer, you can build it from source instead.

First, obtain the repo and its submodules:

    git clone https://github.com/vgteam/vg.git
    cd vg
    # the commit used in the paper
    git checkout ccf6e5d4420ab0e6c8a18aea18b03e9185d016f6
    git submodule update --init --recursive
    
Then, install VG's dependencies. You'll need the protobuf and jansson development libraries installed, and to run the tests you will need:
    * `jq`, `bc`, `rs`, and `parallel`
    * `hexdump` and `column` from `bsdmainutils`
    * [`npm` for testing documentation examples](https://github.com/anko/txm)).
On Ubuntu, you should be able to do:

    make get-deps
    
On other distros, you will need to perform the equivalent of:

    sudo apt-get install build-essential git cmake pkg-config libncurses-dev libbz2-dev  \
                         protobuf-compiler libprotoc-dev libprotobuf-dev libjansson-dev \
                         automake libtool jq bsdmainutils bc rs parallel npm curl unzip \
                         redland-utils librdf-dev bison flex gawk lzma-dev liblzma-dev \
                         liblz4-dev libffi-dev libcairo-dev libboost-all-dev
                         
Note that **Ubuntu 16.04** does not ship a sufficiently new Protobuf; vg requires **Protobuf 3** which will have to be manually installed.

At present, you will need GCC version 4.9 or greater, with support for C++14, to compile vg. (Check your version with `gcc --version`.)

Other libraries may be required. Please report any build difficulties.

Note that a 64-bit OS is required. Ubuntu 18.04 should work.

When you are ready, build with `. ./source_me.sh && make`, and run with `./bin/vg`.

You can also produce a static binary with `make static`, assuming you have static versions of all the dependencies installed on your system.

### Building on MacOS

#### Clone VG

The first step is to clone the vg repository:

    git clone https://github.com/vgteam/vg.git
    cd vg
    # the commit used in the paper
    git checkout ccf6e5d4420ab0e6c8a18aea18b03e9185d016f6
    git submodule update --init --recursive

#### Install Dependencies

VG depends on a number of packages being installed on the system where it is being built. Dependencies can be installed using either [MacPorts](https://www.macports.org/install.php) or [Homebrew](http://brew.sh/).

##### Using MacPorts

You can use MacPorts to install VG's dependencies:

    sudo port install libtool protobuf3-cpp jansson jq cmake pkgconfig autoconf automake libtool coreutils samtools redland bison gperftools md5sha1sum rasqal gmake autogen cairo libomp boost
    

##### Using Homebrew

Homebrew provides another package management solution for OSX, and may be preferable to some users over MacPorts. VG ships a `Brewfile` describing its Homebrew dependencies, so from the root vg directory, you can install dependencies, and expose them to vg, like this:

    # Install all the dependencies in the Brewfile
    brew bundle
    
    # Use GNU versions of coreutils over Apple versions
    export PATH="/usr/local/opt/coreutils/libexec/gnubin:/usr/local/bin:$PATH"

    # Force use of new version of bison
    brew link bison --force
    # NOTE! If brew says that it is refusing to link Bison, follow its suggested
    # instructions to put Bison on your PATH instead.

    # Use glibtool/ize
    export LIBTOOL=glibtool
    export LIBTOOLIZE=glibtoolize

    # Use installed libraries
    export LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH;
    export LIBRARY_PATH=$LD_LIBRARY_PATH;

#### Build

With dependencies installed, VG can now be built:

    . ./source_me.sh && make
    
**Note that static binaries cannot yet be built for Mac.**

Our team has successfully built vg on Mac with GCC versions 4.9, 5.3, 6, 7, and 7.3, as well as Clang 9.0.

## rpvg

### Compilation

*rpvg* requires that [CMake](https://cmake.org) (3.10 or higher), [protobuf](https://github.com/protocolbuffers/protobuf) (version 3) and OpenMP are installed before compilation. 

1. `git clone https://github.com/jonassibbesen/rpvg.git`
2. `cd rpvg`
3. `git checkout ad52c1e22934c646c6614fc1fe2c4e239f28f3c8`
4. `git submodule update --init --recursive`
5. `mkdir build && cd build`
6. `cmake ..`
7. `make -j <threads>` 

Compiling *rpvg* should take 5-10 minutes using 4 threads (`-j`). *rpvg* has been successfully built and tested on Linux (CentOS Linux 7 with GCC 8.1.0 and Ubuntu 18.04 with GCC 7.5.0) and Mac (macOS 10.14.6 with Clang 10.0.1). 

### Docker container

A docker container of the latest commit to master is available on [Docker Hub](https://hub.docker.com/repository/docker/jsibbesen/rpvg) and [Quay](https://quay.io/repository/jsibbesen/rpvg).

# Demo

`vg` has a number of tools to support transcriptomic analyses with spliced graphs (i.e. graphs that have annotated splice junctions added as edges into the graph). These edges can be added into an existing graph using `vg rna`. We can then perform splice-aware mapping to these graphs using `vg mpmap`. The mapped reads can be used for haplotype-aware transcript quantification based on these tools in `rpvg`. 

### Construct and index spliced pangenome graph

The easiest way to start this pipeline is to use the `vg autoindex` subcommand to make indexes for `vg mpmap`. `vg autoindex` creates indexes for mapping from common interchange formats like FASTA, VCF, and GTF. It effectively combines the `vg rna` step and the indexing for `vg mpmap`.

More information is available in the [wiki page on transcriptomics](https://github.com/vgteam/vg/wiki/Transcriptomic-analyses).

Working from this directory, the following example shows how to create a spliced pangenome graph and indexes using `vg autoindex` with 4 threads:

```
# Create spliced pangenome graph and indexes for vg mpmap
vg autoindex --workflow mpmap -t 4 --prefix vg_rna --ref-fasta example_data/x.fa --vcf example_data/x.vcf.gz --tx-gff example_data/x.gtf
```

This will create several files with the prefix `vg_rna`, which can be used in `rpvg` and `vg mpmap`.

### Map reads to the spliced pangenome graph

RNA-seq reads can be mapped to the spliced pangenome graph using `vg mpmap` with 4 threads:

```sh
# Map simulated RNA-seq reads using vg mpmap
vg mpmap -t 4 -x vg_rna.spliced.xg -g vg_rna.spliced.gcsa -d vg_rna.spliced.dist -f example_data/x_rna_1.fq -f example_data/x_rna_2.fq > mpmap.gamp
```

This will create a multipath alignment file called `mpmap.gamp`.

### Perform haplotype-specific expression inferrence

rpvg* requires the following five arguments:

```
rpvg -g graph.xg -p paths.gbwt -a alignments.gamp -o rpvg_results -i <inference-model>
```

The prefix used for all output files are given using `-o`. The number of threads can be given using `-t`. 

A small example dataset containing 36,120 haplotype-specific transcripts and 100,000 read pairs is available as example data.

```
rpvg -t 4 -g example_data/graph.xg -p example_data/pantranscriptome.gbwt -f <(zcat example_data/pantranscriptome.txt.gz) -a example_data/mpmap_align.gamp -o rpvg --inference-model haplotype-transcripts
```

This should take less than a minute to run and will create two files: 

* *rpvg.txt*: Contains the estimated haplotype probability and transcript expression for each haplotype-specific transcript in the pantranscriptome.
* *rpvg_haps.txt*: Contains the estimated probability of each haplotype combination (e.g. diplotype) for all transcripts (only combinations with a probability above zero are shown).



