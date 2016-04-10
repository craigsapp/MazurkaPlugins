
all: sound analysis

sound:

	faad example1.mp4

analysis:

	../../sdk/vamp-plugin-sdk/host/simplehost  \
	../../sdk/vamp-plugin-sdk/examples/plugins.so:mzspectrogramhost \
	example1.wav 0 | sed 's/.*: //' > example1.txt


clean:
	rm -i example1.wav example1.txt

