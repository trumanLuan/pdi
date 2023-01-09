tomtom<-function(pwm.list1,pwm.list2,output.dir="./tmp_tomtom",tomtom.option=""){

tomtom.bin="/usr/local/meme/bin/tomtom"

m.meme1=matrixToMemeText(pwm.list1)
m.meme2=matrixToMemeText(pwm.list2)
#
dir.create(output.dir,recursive=TRUE)
write.table(m.meme1,file.path(output.dir,"m.meme1"),quote=F,row.names=F,col.names=F)
write.table(m.meme2,file.path(output.dir,"m.meme2"),quote=F,row.names=F,col.names=F)

system(paste(tomtom.bin,tomtom.option, " -oc", output.dir,file.path(output.dir,"m.meme1"), file.path(output.dir,"m.meme2")))

tomtom.output=read.csv(file.path(output.dir,"tomtom.txt"),header=T,as.is=T,sep="\t")
return(tomtom.output)
}# end of function


#-------------------------------------------------------------------------------
# http://stuff.mit.edu/afs/athena/software/meme_v3.5.4/etc/meme-explanation.html
# The motif itself is a position-specific probability matrix giving, for each
# position in the pattern, the observed  frequency ('probability') of each
# possible letter. The probability matrix is printed 'sideways'--columns
# correspond  to the letters in the alphabet (in the same order as shown in
# the simplified motif) and rows corresponding to the  positions of the motif,
# position one first. The motif is preceded by a line starting with
# 'letter-probability matrix:' and containing the length of the alphabet,
# width of the motif, number of occurrences of the motif, and the E-value
# of the motif.
matrixToMemeText = function (matrices)
{
  matrix.count = length (matrices)

    # incoming matrices have nucleotide rows, position columns.  meme
    # format, however, requires position-rows, and nucleotide-columns
    # calculate the number of lines of text by counting columns
  total.transposed.matrix.rows = sum (as.integer (sapply (matrices, ncol)))
  predicted.line.count = 12 + (3 * length (matrices)) +
                           total.transposed.matrix.rows
  #s = vector ('character', predicted.line.count)
  s = character (predicted.line.count)

  s [1] = 'MEME version 4'
  s [2] = ''
  s [3] = 'ALPHABET= ACGT'
  s [4] = ''
  s [5] = 'strands: + -'
  s [6] = ''
  s [7] = 'Background letter frequencies'
  s [8] = 'A 0.250 C 0.250 G 0.250 T 0.250 '
  s [9] = '' 

  index = 10
  for (name in names (matrices)) {
       # transpose the frequency matrix version of the incoming matrix,
       # hence 'tfMat'
    tfMat = transformMatrixToMemeRepresentation (matrices [[name]])
       # meme output may be used by tomtom, which uses matrix names as
       # part of image filenames. removed all file-system-unfriendly
       # characters here
    fixed.name = gsub ('\\/', '_', name)
    s [index] = sprintf ('MOTIF %s', fixed.name)
    index = index + 1
    new.line =
       sprintf ('letter-probability matrix: alength= 4 w= %d nsites= %d E=8.1e-020',
          nrow (tfMat), 45, 8.1e-020)
    s [index] =  new.line
    index = index + 1
    for (r in 1:nrow (tfMat)) {
      new.row = sprintf (' %12.10f  %12.10f  %12.10f  %12.10f', tfMat [r,1],
                          tfMat [r,2], tfMat [r,3], tfMat [r,4])
      s [index] = new.row
      index = index + 1
      }
    s [index] = ''
    index = index + 1
    } # for name

  return(s)

} # matrixToMemeText

transformMatrixToMemeRepresentation = function (m)
{
  return (t (m))

} # transformMatrixToMemeRepresentation

