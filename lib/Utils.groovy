import nextflow.Nextflow

/**
 * Utility functions for ews-nf
 */
class Utils {

  // grabs all processed SCE files present in input directory
  static def getLibraryFiles(input_dir){
    def sample_dirs = Nextflow.file(input_dir, type: 'dir')
    def files = Nextflow.files(sample_dirs / "**_processed.rds")
    files = files.findAll{it.size() > 0}
    return files.collect{new Tuple(it.name, it)}
 }

}