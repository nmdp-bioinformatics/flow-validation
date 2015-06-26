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
params.expected    = "${baseDir}/tutorial/output/ex00_ngsp_expected.xml"
params.finalDir    = "${baseDir}/tutorial"
params.outputDir   = "${baseDir}"
params.experiment  = "ex00"
params.resolution  = 2

observed      = Channel.fromPath("${params.finalDir}/*.gz").map { path -> tuple(sample(path), locus(path), path )}
observedFiles = Channel.fromPath("${params.finalDir}/*.gz")
observed_file = observedFiles.collectFile(name:"${params.finalDir}/output/${params.experiment}_ngsp_observed.gz")

process validateInterpretation {
	tag{ s }

  input:
    set s, locus, file(observed_gz) from observed

  output:
    file {"${s}_${locus}_validate.txt"} into validatedFile

  """
    ngs-extract-expected-haploids -i ${params.expected} | ngs-validate-interpretation -l ${locus} -b ${observed_gz} > "${s}_${locus}_validate.txt"
  """

}

validated_file = validatedFile.collectFile(name:"${params.finalDir}/output/${params.experiment}_ngsp_validated.txt")


process validationReport{
	tag{ validated_file }

	input:
		file observed_file
		file validated_file

	"""
		ngs-validation-report -e ${params.expected} -o ${observed_file} -l ${validated_file} -p ${params.outputDir} -f 1 -v 1 
  """
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





