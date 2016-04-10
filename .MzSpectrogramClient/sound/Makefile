
all: sound analysis

sound:

	faad example1.mp4

analysis:

	../../sdk/vamp-plugin-sdk/host/simplehost  \
	../../sdk/vamp-plugin-sdk/examples/plugins.so:mzspectrogramclient \
	example1.wav 0 | sed 's/.*: //' > example1.txt

analysiscompare:

	time ../../sdk/vamp-plugin-sdk/host/simplehost  \
	../../sdk/vamp-plugin-sdk/examples/plugins.so:mzspectrogramhost \
	example1.wav 0 | sed 's/.*: //' > example1a.txt

	time ../../sdk/vamp-plugin-sdk/host/simplehost  \
	../../sdk/vamp-plugin-sdk/examples/plugins.so:mzspectrogramclient \
	example1.wav 0 | sed 's/.*: //' > example1b.txt


clean:

	rm -i example1.wav example1.txt


sinecycle:

	/var/www/websites/sig/bin/osc -s 1024 -f 43.06640625 -a 0.5 sinecycle.wav

sinecycleanalysis:

	time ../../sdk/vamp-plugin-sdk/host/simplehost  \
	../../sdk/vamp-plugin-sdk/examples/plugins.so:mzspectrogramclient \
	sinecycle.wav 0 | sed 's/.*: //' > sinecycleclient.txt

	time ../../sdk/vamp-plugin-sdk/host/simplehost  \
	../../sdk/vamp-plugin-sdk/examples/plugins.so:mzspectrogramhost \
	sinecycle.wav 0 | sed 's/.*: //' > sinecyclehost.txt
