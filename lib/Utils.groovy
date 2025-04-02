import nextflow.Nextflow

/**
 * Utility functions for ews-nf
 */
class Utils {

  // grabs all processed SCE files present in input directory
  static def getLibraryFiles(input_dir){
    input_dir = Nextflow.file(input_dir, type: 'dir')
    files = Nextflow.files(input_dir / "**_processed.rds").findAll{it.size() > 0}
    return files.collect{new Tuple(it.name, it)}
 }

}