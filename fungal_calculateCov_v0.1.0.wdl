
workflow analysis {
    File faidx
    File dat
    String prefix
    Array[String] sample_names
    Array[File] input_bams

    Int mem_size
    Int disk_size

    scatter(i in range(length(sample_names))) {
        String sample_name = sample_names[i]
        File input_bam = input_bams[i]

        call run_ploidy {
            input:
            in_bam = input_bam,
            in_faidx = faidx,
            sample_name = sample_name,
            mem_size = mem_size,
            disk_size = disk_size,
        }
    }

    call coverage_analysis {
        input:
        dat = dat,
        list = run_ploidy.out1,
        sam = sample_names,
        prefix = prefix,
        mem_size = mem_size,
        disk_size = disk_size,
    }

}


## TASK DEFINITIONS
task run_ploidy {
    File in_bam
    File in_faidx
    String sample_name

    Int mem_size
    Int disk_size

    command {
        mkdir ${sample_name}
        run_ploidy.py -i ${in_bam} --faidx ${in_faidx} -n -o $PWD/${sample_name}/
    }

    output {
        String done = "Done"
        File out1 = "${sample_name}/${sample_name}.depth.tsv"
    }

    runtime {
        docker: "lizh666/funpipe:0.1.3"
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}

task coverage_analysis {
    File dat
    Array[String] sam
    Array[File] list
    String prefix

    Int mem_size
    Int disk_size

    command {

        paste ${write_lines(sam)} ${write_lines(list)} > file.list
        coverage_analysis.py -c ${dat} -i file.list -p ${prefix} --no_sub
    }

    output {
        File out = "${prefix}_cov.tsv"
        File png1 = "${prefix}_clustermap.png"
    }
    runtime {
        docker: "lizh666/funpipe:0.1.3"
        memory: mem_size + " GB"
        disks: "local-disk " + disk_size + " HDD"
    }
}
