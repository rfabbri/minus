# generate 3 random frames and 3 random point ids
#
# each frame is 0-99 (number of frames in synthdata is 100)
# each line is 0-5116 (number of points in synthdata is 5117)
#
# not used

f=/tmp/minus-100-configs.$$;

seq 0 99 | gshuf > $f
seq 0 99 | paste $f -
