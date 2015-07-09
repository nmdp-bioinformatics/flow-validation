#!/usr/bin/env nextflow
/*

    flow-validation  Validation workflow.
    Copyright (c) 2014-2015 National Marrow Donor Program (NMDP)

    This library is free software; you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published
    by the Free Software Foundation; either version 3 of the License, or (at
    your option) any later version.

    This library is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this library;  if not, write to the Free Software Foundation,
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.

    > http://www.gnu.org/licenses/lgpl.html

*/
params.path        = "${baseDir}/tutorial"
params.experiment  = "ex00"
params.resolution  = 2

reportDir=file("${params.path}/${params.experiment}/reports")
observedDir=file("${params.path}/${params.experiment}/final")
expected=file("${params.path}/${params.experiment}/hml/${params.experiment}.hml")
hmlDir=file("${params.path}/${params.experiment}/hml")
reportDir.mkdirs()

observed      = Channel.fromPath("${observedDir}/*.observed.gz").map { path -> tuple(sample(path), locus(path), path )}
observedFiles = Channel.fromPath("${observedDir}/*.observed.gz")
observed_file = observedFiles.collectFile(name:"${observedDir}/${params.experiment}_ngsp_observed.gz")

/*
  Validating the typing for each subject at each locus
  - This skips glstrings that were not interpreted 
*/
process validateInterpretation {
	tag{ "${s} ${locus}" }

  input:
    set s, locus, file(observed_gz) from observed

  output:
    file {"${s}_${locus}_validate.txt"} into validatedFile

  """
    gzcat ${observed_gz} | perl -ne 'print \$_ if \$_ =~ /HLA-\\D+\\d{0,1}\\*/;' > new_observed.out
    ngs-extract-expected-haploids -i ${expected} | ngs-validate-interpretation -r ${params.resolution} -l ${locus} -b new_observed.out > "${s}_${locus}_validate.txt"
  """

}
validated_file = validatedFile.collectFile(name:"${observedDir}/${params.experiment}_ngsp_validated.txt")

/*
  Generating the validation report from the expected, validated and observed files
*/
process validationReport{
	tag{ validated_file }

  maxForks 1

	input:
		file observed_file
		file validated_file

  output:
    file {"failedSubjects.csv"} into failed
    file {"passedSubjects.csv"} into passed
    file {"validationReport.tar.gz"} into compressedReport

	"""
		ngs-validation-report -e ${expected} -o ${observed_file} -l ${validated_file} -f -r -v 1
    tar -zcvf validationReport.tar.gz report
  """
}


failedSubjects = failed.splitCsv()
.collectFile() { row ->
       [ "${row[0]}.txt", row[1] + '\n' ]
   }
.map{path -> 
  tuple(sample(path), path, file(expected) )
}

passedSubjects = passed.splitCsv()
.collectFile() { row ->
       [ "${row[0]}.txt", row[1] + '\n' ]
   }
.map{path -> 
  tuple(sample(path), path, file(expected)  )
}

/*
  Filtering out the failed subjects
*/
process filterFailedHml{
  tag{ exp }

  input:
    set exp, file(failedFile), file(expectedFile)  from failedSubjects

  output:
    file{"${exp}_failed.hml"} into failedHmlFile

  """
    ngs-filter-samples -i ${expectedFile} -s ${failedFile} > ${exp}_failed.hml
  """
}

/*
  Filtering out the passed subjects
*/
process filterPassedHml{
  tag{ exp }

  input:
    set exp, file(failedFile), file(expectedFile)  from passedSubjects

  output:
    file{"${exp}_passed.hml"} into passedHmlFile

  """
    ngs-filter-samples -i ${expectedFile} -s ${failedFile} > ${exp}_passed.hml
  """
}

//Copying the final hml files and the report to the final directory
passedHmlFile.subscribe {    file -> copyToFinalDir(file) }
failedHmlFile.subscribe {    file -> copyToFinalDir(file) }
compressedReport.subscribe { file -> copyToReportDir(file) }

def copyToFinalDir (file) { 
  log.info "Copying ${file.name} into: $hmlDir"
  file.copyTo(hmlDir)
}

def copyToReportDir (file) { 
  log.info "Copying ${file.name} into: $reportDir"
  file.copyTo(reportDir)
}

def sample(Path path) {
  def name = path.getFileName().toString()
  int start = Math.max(0, name.lastIndexOf('/'))
  return name.substring(start, name.indexOf("."))
}

def locus(Path path) {
  def name = path.getFileName().toString()
  loc = name =~ /\.(HLA-\D{1,3}\d{0,1}).observed.gz$/
  return loc[0][1]
}





