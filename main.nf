#!/usr/bin/env nextflow

params.expected    = "${baseDir}/tutorial/output/ex00_ngsp_expected.txt"
params.finalDir    = "${baseDir}/tutorial"
params.experiment  = "ex00"
params.resolution  = 2

observed = Channel.fromPath("${params.finalDir}/*.txt").map { path -> tuple(sample(path), path )}
observedFiles = Channel.fromPath("${params.finalDir}/*.txt")


process validateInterpretation {
	tag{ s }

  input:
    set s, file(observed_file) from observed

  output:
    file {"${s}_ngsp_validate.txt"} into validatedFile

  """
  	ngs-validate-interpretation  -e ${params.expected} -b ${observed_file} -r ${params.resolution} > "${s}_ngsp_validate.txt"
  """

}

validated_file = validatedFile.collectFile(name:"${params.finalDir}/output/ex00_ngsp_validated.txt")
observed_file  = observedFiles.collectFile(name:"${params.finalDir}/output/ex00_ngsp_observed.txt")

process validationReport{
	
	input:
		file observed_file
		file validated_file

	"""
		ngs-validation-report -e ${params.expected} -o ${observed_file} -l ${validated_file} -p ${params.finalDir} -f 1 -v 1
    """
}


def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf("."))
}

def locus(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf(".")+1)
}