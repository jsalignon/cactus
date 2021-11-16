import static nextflow.Nextflow.file
import nextflow.Channel

class atacUtils {

  //////////////////////////////////////////////////////////////////////////////
  //// extracting keys from the grouped key

  static def getFDR(file_obj) {
      def result = (file_obj.name =~ /.*(FDR_.*?)__.*/)
      return(result [0][1])
  }

  static def getOther(file_obj, keytype) {
      def result = (file_obj.name =~ /(${keytype}_\w[^_]*)/)
      return(result [0][1])
  }

  static def getFC(file_obj) {
      return(getOther(file_obj, 'FC'))
  }

  static def getPA(file_obj) {
      return(getOther(file_obj, 'PA'))
  }

  static def getET(file_obj) {
      return(getOther(file_obj, 'ET'))
  }

  static def getIdComp(file_obj) {
    def result = (file_obj.name =~ /([^_]*_vs_\w[^_]*)/)
    return(result [0][1])
  }

}


// trying these commands to insert into the script but it failed
evaluate(new File("${atacdir}/src/nextflow/lib/atac_utils.groovy"))

def script = new GroovyScriptEngine( '.' ).with {
  loadScriptByName( "${atacdir}/src/nextflow/lib/atac_utils.groovy" )
}

// examples
// Channel.from([ file('/home/salignon//peak_set__g4d3_vs_g4d1__FDR_1.3__FC_up__PA_genProm__ET_both.bed'), file('/home/salignon/peak_set__g4d3_vs_g4d1__FDR_1000__FC_up__PA_all__ET_both.bed')])
//     .map{ [ getIdComp(it), getFDR(it), getFC(it), getPA(it), getET(it) ] }
//     .view()




////////////////////////////////////////////////////////////////////////////
// creating a function that puts a file in a folder with the name of its extension

// commands to get file extension
// file.name.lastIndexOf('.').with {it != -1 ? file.name[0..<it] : file.name}
//
// file.name.replaceFirst(~/\.[^\.]+$/, '')
// String fileWithoutExt = file.name.take(file.name.lastIndexOf('.'))

// static def getExtensionAndFile(def filename) {
//   def extension = filename.drop(filename.lastIndexOf('.'))
//   return( extension/filename )
// }

// these 3 functions works. The only problem I have is to create variables
static def getExtensionAndFile(filename) {
  return( [ filename.drop(filename.lastIndexOf('.') + 1), filename].join('/') )
}

static def getExtensionAndFile(filename) {
  return( filename.drop(filename.lastIndexOf('.') + 1) + '/' + filename )
}

static def getExtensionAndFile(filename) {
  return( "${filename.drop(filename.lastIndexOf('.') + 1)}/${filename}" )
}

Channel.from([ file('test1.bed'), file('test2.bam'), file('test3.bai')])
    .map{ it.name }
    .map{ getExtensionAndFile(it) }
    .view()


// class AtacUtils {
//     static def getExtensionAndFile(def filename) {
//       def extension = filename.drop(filename.lastIndexOf('.'))
//       def path = extension + "/" + filename
//       return( path )
//     }
//   }
//
// Channel.from([ file('test1.bed'), file('test2.bam'), file('test3.bai')])
//     .map{ AtacUtils.getExtensionAndFile(it) }
//     .view()


// static def getExtension(it) {
//   def file = it.name
//   return( file.drop(file.lastIndexOf('.')) )
// }
//
// Channel.from([ file('/home/salignon/test1.bed'), file('/home/salignon/test2.bam'), file('/home/salignon/test3.bai')])
//     .map{ getExtension(it) }
//     .view()
//



