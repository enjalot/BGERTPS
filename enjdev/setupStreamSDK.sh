
# Select all display GPUs as candidates for computation (:0.0 would specify only the first GPU+CPU)

# Cut any .0 or .1 that is appended to an existing DISPLAY var. 
if [ -z $DISPLAY ]
then 
#	echo "DISPLAY UNSET"
	export DISPLAY=:0
else 
#	echo "DISPLAY SET"
	export DISPLAY=${DISPLAY%.*}
fi

# Where this script should reside
export ATIOPENCLHOME=~efb06/research/ati_opencl

# Variables required by the SDK and many examples
export ATISTREAMSDKROOT=$ATIOPENCLHOME/ati-stream-sdk-v2.2-lnx64
export ATISTREAMSDKSAMPLESROOT=$ATIOPENCLHOME/ati-stream-sdk-v2.2-lnx64/samples

# Update our lib search dirs so we find libOpenCL, etc
export LD_LIBRARY_PATH=$ATISTREAMSDKROOT/lib/x86_64:$ATISTREAMSDKROOT/lib/x86:/usr/lib/fglrx:/usr/lib32/fglrx:$LD_LIBRARY_PATH
