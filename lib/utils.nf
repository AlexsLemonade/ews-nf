// Utility functions for ews-nf

/**
  * Parses the processed file path to extract unique id, sample id, and number of processed cells from metadata json file
  * 
  * @param processed_file path to processed file (either .rds or .h5ad)
  * @param extension file extension to remove to get unique id (default is "_processed.rds")
  * @return tuple of unique id, sample id, path to processed file, and number of processed cells
**/


def parseProcessedFile(processed_file, extension = "_processed.rds") {

  def file_id = processed_file.name.replace(extension, "") // extract unique id which will be either library_id or <library_id>-<sample_id>
  def sample_id = processed_file.parent.name // extract sample id from folder name
  def unique_id = file_id.contains("-") ? file_id : "${file_id}-${sample_id}" // if unique id doesn't already contain sample id, append it

  def metadata_json = new File(processed_file.toString().replace(extension, "_metadata.json")) // construct path to metadata json file based on processed file path
  def meta = new groovy.json.JsonSlurper().parse(metadata_json) // read in metadata json file as groovy object
  def num_processed_cells = meta["processed_cells"] // pull out number of processed cells from metadata

  return [unique_id, sample_id, processed_file, num_processed_cells]
}
