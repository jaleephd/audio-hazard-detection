# Audio Hazard Detection

_Algorithm to detect vehicle in vicinity of pedestrian through audio processing alone._

This algorithm addresses the problem of detection of vehicles approaching a pedestrian
by a novel nonresource intensive acoustic method. The method uses a set of
statistical signal processing tools to identify signal features. The audio features
are adaptively thresholded for relevance and classified with a three-component heuristic.
The resulting acoustic hazard detection system has a very low false-positive detection rate.

This code was used in the following publication:

"Acoustic hazard detection for pedestrians with obscured hearing"
Justin A. Lee & Andry Rakotonirainy (2011),
IEEE Transactions on Intelligent Transportation Systems, 12(4), pp. 1640-1649.
Available at <http://eprints.qut.edu.au/48109/>

## Development Environment and Disclaimer

This code is provided for illustrating how vehicles can be accurately identified acoustically,
as per the referenced paper above. It is research code, not production level code.

The code was written in Matlab in 2009-2010, probably with Matlab 7, on Windows XP Professional
(Cygwin was also utilised for Unix-like command-line functionality, but shouldn't be needed for
running the Matlab code.) It may or may not run "out of the box" on current versions of Matlab
(I haven't run it since 2010).


## Running the code

The main Matlab script to run the algorithm is `AHD.m`, and the file to be processed is defined
near the top of the script with `fileName=`. Other parameters can be modified as required.
The called functions are identified by their filename.

The code will output signal feature and detection score graphs (optionally these can be saved),
and a textual output of detection events.
