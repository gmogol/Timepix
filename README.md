# Timepix
This is a repository to preprocess TpX3cam data and analyse it. The camera has 3 outputs; camera pixel outputs and 2 time to digitial converter (TDC) outputs.
TDC1 is triggered by the laser pulse and has 1.5ns precision. TDC2 is the high precision TDC with 0.25ns and is triggered by the phospor plate and is used to
get time of flight of electrons. Outputs of these 3 devices are stored together in one file and have the form:

* Pixel Data:  ```TOA: 0.0007992671875, ToT: 150, x: 157, y: 32```
* TDC Data: ```TDC2: 0.00103862551879883```

Here the pixel data has time of arrival (TOA) and time over threshold (ToT) information. The former is used to associate the pixel data with corresponding TDC1
and TDC2 triggers. ToT is a measure of intensity of the pixel and is used center of mass centroiding.

First we create a seperate dataset for TDC data alone so that we can iterate over it. Pixel data is associated to TDC data by thresholding the 
difference between TOA of a pixel and TDCs.
We call the set of pixels that are associated to a single TDC2 trigger a batch and each batch is thresholded in x and y coordinates. This removes almost all of the 
many electron events.
