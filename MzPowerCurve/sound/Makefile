
all: sound analysis


allb: soundb analysisb

clean:

	rm -i sinetest-raw.txt sinetest.wav

sound:
	/var/www/websites/sig/bin/osc -d 5 -f 440 -a "0 0 1 1 2 0" sinetest.wav

soundb:
	/var/www/websites/sig/bin/osc -d 5 -f 440 -a "0 0 1 0.5 2 0" sinetestb.wav

analysis:

	../../vamp/sdk/host/simplehost  \
	../../vamp/sdk/examples/plugins.so:mzpowercurve \
	sinetest.wav 0 | sed 's/000000: /\t/'  > sinetest-raw.txt

analysisx:

	../../vamp/sdk/host/simplehost  \
	../../vamp/sdk/examples/plugins.so:mzpowercurve \
	sinetest.wav 0  > sinetest-rawx.txt


analysis2:
	/var/www/websites/soundfile/bin/smoothpower -SR sinetest.wav > sinetest-raw2.txt


analysisb:

	../../vamp/sdk/host/simplehost  \
	../../vamp/sdk/examples/plugins.so:mzpowercurve \
	sinetest.wav 0 | sed 's/000000: /\t/'  > sinetest-rawb.txt


smoothanalysis:

	../../vamp/sdk/host/simplehost  \
	../../vamp/sdk/examples/plugins.so:mzpowercurve \
	sinetest.wav 1 | sed 's/000000: /\t/' > sinetest-smooth.txt

smoothanalysisb:

	/var/www/websites/soundfile/bin/smoothpower -R sinetest.wav > sinetest-smoothb.txt


slopeanalysis:

	../../vamp/sdk/host/simplehost  \
	../../vamp/sdk/examples/plugins.so:mzpowercurve \
	sinetest.wav 2 | sed 's/000000: /\t/'  > slope.txt

productanalysis:

	../../vamp/sdk/host/simplehost  \
	../../vamp/sdk/examples/plugins.so:mzpowercurve \
	sinetest.wav 3 | sed 's/000000: /\t/'  > product.txt





