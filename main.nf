#!/usr/bin/env nextflow
// -*- mode:groovy -*-

params.input="$baseDir/input" //


//chunk size test




process finemap {
    scratch true
    label 'finemap'
    input:
    output:
    file(finemap_output) into finemap_output
    script:

    """

    """ 
}


finemap_output.view()

workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following directory for the outputs --> $params.output\n" : "Oops .. something went wrong" )
}